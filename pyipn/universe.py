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

from .effective_area import EffectiveArea
from .geometry import Pointing, DetectorLocation, Location

from .grb import GRB
from .detector import Detector

from .io.plotting.projection import *
from .io.plotting.spherical_circle import SphericalCircle


class Universe(object):
    def __init__(self, grb):
        """FIXME! briefly describe function

        :param grb: 
        :returns: 
        :rtype: 

        """

        self._detectors = collections.OrderedDict()
        self._grb = grb

        self._time_differences = (
            None
        )  # array of time differences ordered like _detectors
        self._T0 = (
            None
        )  # array of times at which detectors get hit by GRB ordered like _detectors
        self._light_curves = None
        self._n_detectors = 0

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
    def detectors(self):
        return self._detectors

    @property
    def light_curves(self):
        return self._light_curves

    def explode_grb(self, tstart, tstop, verbose=True):
        """FIXME! briefly describe function

        :param verbose: 
        :returns: 
        :rtype: 

        """

        self._compute_time_differences()

        self._create_light_curves(tstart, tstop)

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

            # calculate closest distance to wavefront when the GRB reaches the detector
            # (negative sign for right order)
            ltd.append(
                -norm_grb_vec.dot(detector.location.get_cartesian_coord().xyz).value
            )

        # rank the distances in ascending order

        self._distance_rank = np.argsort(ltd)
        unsort = self._distance_rank.argsort()

        # for now compute considering all detectors are static
        # the TOA difference of each detector
        ltd = np.array(ltd)[self._distance_rank]

        self._time_differences = [0.0]
        self._T0 = [0.0]
        T0 = 0.0
        for i in range(len(ltd) - 1):

            dt = (
                ((ltd[i + 1] - ltd[i]) * u.km / constants.c).decompose().value
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

            grb_params = setup["grb"]

            grb = GRB(
                grb_params["ra"],
                grb_params["dec"],
                grb_params["distance"] * u.Mpc,
                grb_params["K"],
                grb_params["t_rise"],
                grb_params["t_decay"],
            )

            universe = cls(grb)

            for name, value in setup["detectors"].items():

                eff_area = EffectiveArea(value["effective_area"])

                time = Time(value["time"])

                location = DetectorLocation(
                    value["ra"], value["dec"], value["altitude"] * u.km, time
                )

                pointing = Pointing(value["pointing"]["ra"], value["pointing"]["dec"])

                det = Detector(location, pointing, eff_area, name)

                universe.register_detector(det)

            return universe

    def calculate_annulus(self, detector1, detector2):
        """FIXME! briefly describe function

        :param detector1: 
        :param detector2: 
        :returns: 
        :rtype: 

        """

        d1, d2 = self._detectors[detector1], self._detectors[detector2]
        dxyz = (
            d2.location.get_cartesian_coord().xyz
            - d1.location.get_cartesian_coord().xyz
        )

        # calculate ra and dec of vector d  pointing from detector1 to detector2
        dcart = Location(
            SkyCoord(
                x=dxyz[0],
                y=dxyz[1],
                z=dxyz[2],
                representation_type="cartesian",
                unit="km",
            )
        )

        norm_d = dcart.get_norm_vec(u.km)
        ra = dcart.coord.represent_as(UnitSphericalRepresentation).lon
        dec = dcart.coord.represent_as(UnitSphericalRepresentation).lat

        # calculate angle theta between center point d and annulus
        distance = np.linalg.norm(dxyz) * u.km
        dt = (
            self._T0[list(self._detectors.keys()).index(detector1)]
            - self._T0[list(self._detectors.keys()).index(detector2)]
        ) * u.s
        # rounding to 15th decimal because small numerical errors cause issues with numbers slightly over 1
        theta = np.arccos(
            np.around((constants.c * dt / distance).decompose().value, 15)
        )

        return (norm_d, np.array([ra.value, dec.value]) * ra.unit, theta * u.rad)

    def plot_annulus(
        self,
        detector1,
        detector2,
        projection="astro degrees mollweide",
        ax=None,
        radius=None,
        center=None,
        **kwargs
    ):

        if ax is None:

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

            fig = ax.get_figure()

        # compute the annulus for this set of detectors
        cart_vec, spherical_vec, theta = self.calculate_annulus(detector1, detector2)

        circle = SphericalCircle(
            spherical_vec,
            theta,
            vertex_unit=u.deg,
            resolution=5000,
            #            edgecolor=color,
            facecolor="none",
            transform=ax.get_transform("icrs"),
            **kwargs
        )

        ax.add_patch(circle)

        return fig

    def plot_all_annuli(
        self,
        projection="astro degrees mollweide",
        radius=None,
        center=None,
        cmap="Set1",
        **kwargs
    ):

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

        
        # get the colors to use

        n_verts = self._n_detectors * (self._n_detectors - 1)/2
        
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
                **kwargs
            )

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
