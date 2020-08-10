import numpy as np
import collections
from itertools import combinations
import yaml
from astropy.time import Time
import astropy.units as u
import astropy.constants as constants
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation
from mpltools import color as mpl_color
import matplotlib.pyplot as plt
import ipyvolume as ipv
import h5py

from .effective_area import EffectiveArea
from .geometry import Pointing, DetectorLocation, Location
from .lightcurve import LightCurve
from .grb import GRB
from .detector import Detector

from .io.plotting.projection import *
from .io.plotting.projection import create_skw_dict
from .io.plotting.spherical_circle import SphericalCircle, get_3d_circle
from .utils.hdf5_utils import (
    recursively_load_dict_contents_from_group,
    recursively_save_dict_contents_to_group,
)

from .utils.timing import (
    compute_annulus_from_time_delay,
    calculate_distance_and_norm,
    theta_from_time_delay,
)


class Universe(object):
    def __init__(self, grb, yaml_dict=None, locked=False):
        """FIXME! briefly describe function

        :param grb: 
        :returns: 
        :rtype: 

        """

        self._detectors = collections.OrderedDict()
        self._grb = grb

        self._grb_radius = 1e6

        self._time_differences = (
            None  # array of time differences ordered like _detectors
        )
        self._T0 = None  # array of times at which detectors get hit by GRB ordered like _detectors
        self._light_curves = None
        self._n_detectors = 0

        self._locked = locked
        self._yaml_dict = yaml_dict

    def register_detector(self, detector):
        """FIXME! briefly describe function

        :param detector: 
        :returns: 
        :rtype: 

        """

        self._detectors[detector.name] = detector

        self._n_detectors += 1

        assert self._n_detectors == len(self._detectors.keys())

    @property
    def grb(self):
        return self._grb

    @property
    def T0(self):
        return self._T0

    @property
    def detectors(self):
        return self._detectors

    @property
    def light_curves(self):
        return self._light_curves

    @property
    def grb_radius(self):
        return self._grb_radius

    def explode_grb(self, tstart, tstop, verbose=True):
        """FIXME! briefly describe function

        :param verbose: 
        :returns: 
        :rtype: 

        """

        if not self._locked:

            self._compute_time_differences()

            self._create_light_curves(tstart, tstop)
        else:

            print("This universe is locked")

    def _compute_time_differences(self):
        """FIXME! briefly describe function

        :returns: 
        :rtype: 

        """
        # compute which detector sees the GRB first
        ltd = []

        norm_grb_vec = self._grb.location.get_norm_vec(
            u.km
        )  # normalized vector towards GRB

        for name, detector in self._detectors.items():

            # calculate closest distanistancece to wavefront when the GRB reaches the detector
            # (negative sign for right order)
            ltd.append(
                -norm_grb_vec.dot(detector.location.get_cartesian_coord().xyz)
                .to("km")
                .value
            )

        # rank the distances in ascending order

        self._ltd_rank = np.argsort(ltd)
        unsort = self._ltd_rank.argsort()

        # for now compute considering all detectors are static
        # the TOA difference of each detector
        ltd = np.array(ltd)[self._ltd_rank]

        self._time_differences = [0.0]
        self._T0 = [0.0]
        T0 = 0.0
        for i in range(len(ltd) - 1):

            dt = (
                ((ltd[i + 1] - ltd[i]) * u.km / constants.c).decompose().to("s").value
            )  # time in seconds
            assert (
                dt >= 0
            ), "The time diferences should be positive if the ranking worked!"

            T0 += dt
            self._T0.append(T0)
            self._time_differences.append(dt)

        self._T0 = np.array(self._T0)
        self._time_differences = np.array(self._time_differences)  # time in s

        self._T0 = self._T0[unsort]
        self._time_differences[unsort]

    def _create_light_curves(self, tstart, tstop):
        """FIXME! briefly describe function

        :returns: 
        :rtype: 

        """

        self._light_curves = collections.OrderedDict()

        for t0, (name, detector) in zip(self._T0, self._detectors.items()):

            self._light_curves[name] = detector.build_light_curve(
                self._grb, t0, tstart, tstop
            )

    def write_to(self, file_name):

        src_lc = []
        bkg_lc = []

        for k, v in self._light_curves.items():

            src_lc.append(v.source_arrival_times)
            bkg_lc.append(v.bkg_arrival_times)

        uni_save = UniverseSave(self._yaml_dict, src_lc, bkg_lc)

        uni_save.write_to(file_name)

    @classmethod
    def from_dict(cls, d, locked=False):

        grb_params = d["grb"]

        if "t_start" in grb_params:

            t_start = grb_params["t_start"]
        else:

            t_start = None

        if isinstance(grb_params["K"], list):

            grb_params["K"] = np.array(grb_params["K"])
            grb_params["t_rise"] = np.array(grb_params["t_rise"])
            grb_params["t_decay"] = np.array(grb_params["t_decay"])
            grb_params["t_start"] = np.array(grb_params["t_start"])

        grb = GRB(
            grb_params["ra"],
            grb_params["dec"],
            grb_params["distance"] * u.Mpc,
            grb_params["K"],
            grb_params["t_rise"],
            grb_params["t_decay"],
            t_start,
        )

        universe = cls(grb, yaml_dict=d, locked=locked)

        for name, value in d["detectors"].items():

            eff_area = EffectiveArea(value["effective_area"])

            time = Time(value["time"])

            location = DetectorLocation(
                value["ra"], value["dec"], value["altitude"] * u.km, time
            )

            pointing = Pointing(value["pointing"]["ra"], value["pointing"]["dec"])

            det = Detector(location, pointing, eff_area, name)

            universe.register_detector(det)

        return universe

    @classmethod
    def from_yaml(cls, yaml_file):
        """
        Create a universe from a yaml file

        :param cls: 
        :param yaml_file: 
        :returns: 
        :rtype: 

        """

        with open(yaml_file, "r") as f:

            setup = yaml.load(f, Loader=yaml.SafeLoader)

            universe = cls.from_dict(setup)

        return universe

    @classmethod
    def from_save_file(cls, file_name):

        uni_save = UniverseSave.from_file(file_name)

        universe = cls.from_dict(uni_save.yaml_dict, locked=True)

        universe._light_curves = collections.OrderedDict()

        for i, (k, v) in enumerate(universe.detectors.items()):

            lc = LightCurve(uni_save.source_lightcurves[i], uni_save.bkg_lightcurves[i])

            universe._light_curves[k] = lc

        universe._compute_time_differences()

        return universe

    def calculate_annulus(self, detector1, detector2):
        """FIXME! briefly describe function

        :param detector1: 
        :param detector2: 
        :returns: 
        :rtype: 

        """

        d1, d2 = self._detectors[detector1], self._detectors[detector2]

        distance, norm_d, ra, dec = calculate_distance_and_norm(d1, d2)

        dt = (
            self._T0[list(self._detectors.keys()).index(detector1)]
            - self._T0[list(self._detectors.keys()).index(detector2)]
        ) * u.second  # seconds
        # rounding to 15th decimal because small numerical errors cause issues with numbers slightly over 1

        theta = theta_from_time_delay(dt, distance)

        return (norm_d, np.array([ra.value, dec.value]) * ra.unit, theta * u.rad)

    def plot_annulus(
        self,
        detector1,
        detector2,
        projection="astro degrees mollweide",
        ax=None,
        radius=None,
        center=None,
        threeD=True,
        **kwargs,
    ):

        if not threeD:
            if ax is None:

                skw_dict = create_skw_dict(projection, center, radius)

                fig, ax = plt.subplots(subplot_kw=skw_dict)

            else:

                fig = ax.get_figure()

        # compute the annulus for this set of detectors
        cart_vec, spherical_vec, theta = self.calculate_annulus(detector1, detector2)

        if not threeD:
            circle = SphericalCircle(
                spherical_vec,
                theta,
                vertex_unit=u.deg,
                resolution=5000,
                #            edgecolor=color,
                fc="none",
                transform=ax.get_transform("icrs"),
                **kwargs,
            )

            ax.add_patch(circle)

            return fig

        else:

            # get all the threeD point

            xyz = get_3d_circle(
                spherical_vec, theta, radius=self._grb_radius, resolution=1000
            )

            ipv.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], **kwargs)

    def to_stan_data(self, tstart, tstop, dt=0.2, k=50, n_cores=1, factor=0.5):

        n_dets = len(self._detectors)

        counts = []
        times = []
        exposures = []
        sc_pos = np.empty((n_dets, 3))

        n_time_bins = []

        tstart = np.atleast_1d(tstart)
        tstop = np.atleast_1d(tstop)
        dt = np.atleast_1d(dt)
        
        for n, (det_nam, v) in enumerate(self._detectors.items()):

            lc = self._light_curves[det_nam]

            _, t, c = lc.get_binned_light_curve(tstart[n], tstop[n], dt[n])

            mid = np.mean([t[:-1], t[1:]], axis=0)
            e = t[1:] - t[:-1]

            counts.append(c)
            times.append(mid)
            exposures.append(e)
            n_time_bins.append(len(c))

            xyz = v.location.get_cartesian_coord().xyz.value
            sc_pos[n] = xyz

        max_n_time_bins = max(n_time_bins)

        counts_stan = np.zeros((n_dets, max_n_time_bins), dtype=int)
        times_stan = np.zeros((n_dets, max_n_time_bins))
        exposure_stan = np.zeros((n_dets, max_n_time_bins))

        for n in range(n_dets):

            counts_stan[n, : n_time_bins[n]] = counts[n]
            times_stan[n, : n_time_bins[n]] = times[n]
            exposure_stan[n, : n_time_bins[n]] = exposures[n]

        #     data = dict(N_detectors=n_dets,
        #                 N_time_bins = n_time_bins[::-1],
        #                 max_N_time_bins = max_n_time_bins,
        #                 counts = counts_stan[::-1,:],
        #                 time = times_stan[::-1,:],
        #                 exposure = exposure_stan[::-1,:],
        #                 sc_pos = sc_pos[::-1,:],
        #                 k=k,
        #                 grainsize=1,
        #                 bw=1. )

        grainsize = []
        for n in n_time_bins:

            grainsize.append(int(np.round(n / n_cores) * factor))
            # grainsize.append(1)

        data = dict(
            N_detectors=n_dets,
            N_time_bins=n_time_bins,
            max_N_time_bins=max_n_time_bins,
            counts=counts_stan,
            time=times_stan,
            exposure=exposure_stan,
            sc_pos=sc_pos,
            k=k,
            grainsize=grainsize,
            bw=1.0,
        )

        return data

    def plot_all_annuli(
        self,
        projection="astro degrees mollweide",
        radius=None,
        center=None,
        cmap="Set1",
        threeD=True,
        **kwargs,
    ):

        if not threeD:

            assert projection in [
                "astro degrees aitoff",
                "astro degrees mollweide",
                "astro hours aitoff",
                "astro hours mollweide",
                "astro globe",
                "astro zoom",
            ]

            skw_dict = dict(projection=projection)

            if projection in ["astro globe", "astro zoom"]:

                assert center is not None, "you must specify a center"

                skw_dict = dict(projection=projection, center=center)

            if projection == "astro zoom":

                assert radius is not None, "you must specify a radius"

                skw_dict = dict(projection=projection, center=center, radius=radius)

            fig, ax = plt.subplots(subplot_kw=skw_dict)

        else:

            fig = ipv.figure()
            ipv.pylab.style.box_off()
            ipv.pylab.style.axes_off()
            ax = None

        # get the colors to use

        n_verts = self._n_detectors * (self._n_detectors - 1) / 2

        colors = mpl_color.colors_from_cmap(int(n_verts), cmap=cmap)

        for i, (d1, d2) in enumerate(combinations(self._detectors.keys(), 2)):

            _ = self.plot_annulus(
                d1,
                d2,
                projection=projection,
                center=center,
                radius=radius,
                ax=ax,
                edgecolor=colors[i],
                threeD=threeD,
                color=colors[i],
                **kwargs,
            )

            if threeD:

                loc1 = self._detectors[d1].location.get_cartesian_coord().xyz.value
                loc2 = self._detectors[d2].location.get_cartesian_coord().xyz.value

                ipv.plot(
                    np.array([loc1[0], loc2[0]]),
                    np.array([loc1[1], loc2[1]]),
                    np.array([loc1[2], loc2[2]]),
                    color=colors[i],
                )

        if threeD:

            ipv.scatter(
                *(
                    self._grb_radius
                    * self._grb.location.get_cartesian_coord().xyz.value
                    / np.linalg.norm(self._grb.location.get_cartesian_coord().xyz.value)
                )[np.newaxis].T,
                marker="sphere",
                color="green",
            )

            ipv.show()

        return fig

    def localize_GRB(self):
        M = []
        b = []

        """
        build matrix M consisting of connection vectors between two satellites
        and vector b containing corresponding cos of annulus angles
        """
        for (d0, d1) in combinations(self._detectors.keys(), 2):
            (cart_vec, spherical_vec, theta) = self.calculate_annulus(d0, d1)
            M.append(cart_vec.value)
            b.append(np.array([np.cos(theta.value)]))

        M = np.array(M)
        b = np.array(b)

        g = np.linalg.lstsq(M, b, rcond=None)
        grb_loc = Location(
            SkyCoord(
                x=g[0][0][0],
                y=g[0][1][0],
                z=g[0][2][0],
                representation_type="cartesian",
                unit="km",
            )
        )
        norm_grb_loc = grb_loc.get_norm_vec(u.km)
        return grb_loc


class UniverseSave(object):
    def __init__(self, yaml_dict, source_lightcurves, bkg_lightcurves):

        self._yaml_dict = yaml_dict
        self._source_lightcurves = source_lightcurves
        self._bkg_lightcurves = bkg_lightcurves

        assert len(source_lightcurves) == len(bkg_lightcurves)

        self._n_light_curves = len(source_lightcurves)

    def write_to(self, file_name):

        with h5py.File(file_name, "w") as f:

            recursively_save_dict_contents_to_group(f, "yaml", self._yaml_dict)

            f.attrs["n_dets"] = self._n_light_curves

            for i in range(self._n_light_curves):

                f.create_dataset(
                    f"src_lc{i}", data=self._source_lightcurves[i], compression="gzip"
                )
                f.create_dataset(
                    f"bkg_lc{i}", data=self._bkg_lightcurves[i], compression="gzip"
                )

    @classmethod
    def from_file(cls, file_name):

        with h5py.File(file_name, "r") as f:

            yaml_dict = recursively_load_dict_contents_from_group(f, "yaml")

            src_lcs = []
            bkg_lcs = []

            n_light_curves = f.attrs["n_dets"]

            for i in range(n_light_curves):

                src_lcs.append(f[f"src_lc{i}"][()])
                bkg_lcs.append(f[f"bkg_lc{i}"][()])

        return cls(yaml_dict, src_lcs, bkg_lcs)

    @property
    def yaml_dict(self):

        return self._yaml_dict

    @property
    def source_lightcurves(self):

        return self._source_lightcurves

    @property
    def bkg_lightcurves(self):

        return self._bkg_lightcurves
