#
# pieces of this code were taken from and inspired by
#

import collections
import urllib
from math import pi

import astropy.time as atime
import ipyvolume as ipv
import numpy as np
import requests
import yaml
from astropy.coordinates import SkyCoord
from bs4 import BeautifulSoup
from orbital import KeplerianElements, earth
from tletools import TLE

from vangogh import Earth

heavenly_body_radius = {
    "earth": 6371,
    "luna": 1737,
    "mars": 3390,
    "venus": 6052,
    "mercury": 2440,
    "sol": 695700,
    "jupiter": 69911,
    "saturn": 58232,
    "uranus": 25362,
    "neptune": 24622,
    "pluto": 1188,
}


class SatelliteCollection(object):
    def __init__(self, *satellites, normal_pointing=False):
        """
        A collection of satellites either in a constellation
        or from TLEs. Can be used to construct YAML files for 
        simulations

        :param normal_pointing: Give the satellites pointing normal to earth
        :returns: 
        :rtype: 

        """

        self._n_satellties = 0
        self._satellites = collections.OrderedDict()

        # chech that everything is a satellite and
        # build up the dictionary

        for sat in satellites:

            self.add_satellite(sat)

        self._n_satellties = len(self._satellites)

        self._normal_pointing = normal_pointing

    @property
    def satellites(self):

        return self._satellites

    @property
    def n_satellites(self):

        return self._n_satellties

    @classmethod
    def from_tle_file(cls, tle_file, normal_pointing=False):
        """
        Construct a collection from a TLE file

        :param cls: 
        :param tle_file: 
        :param normal_pointing: 
        :returns: 
        :rtype: 

        """

        tles = TLE.load(tle_file)

        sats = []

        for tle in tles:

            sat = Satellite.from_TLE(tle)
            sats.append(sat)

        return cls(*sats)

    @classmethod
    def from_celestrack(cls, satellite_group):

        try:

            req = urllib.request.urlopen(
                f"https://www.celestrak.com/NORAD/elements/{satellite_group}.txt")

            x = req.read().decode()

            sats = []

            tles = TLE.loads(x)

            for tle in tles:

                sat = Satellite.from_TLE(tle)
                sats.append(sat)

            return cls(*sats)

        except:

            url = 'https://www.celestrak.com/NORAD/elements'

            page = requests.get(url).text
            soup = BeautifulSoup(page, 'html.parser')
            xx = [node.get('href') for node in soup.find_all('a')]

            xxx = []

            print("satellite_group must be one of the following:")

            for x in xx:

                if x is not None and (".txt" in x) and ("https" not in x):

                    print(x[:-4])

    def add_satellite(self, satellite):
        """
        Add a new satellite

        :param satellite: 
        :returns: 
        :rtype: 

        """

        assert isinstance(satellite, Satellite)
        self._satellites[satellite.name] = satellite
        self._n_satellties += 1

    def add_satellite_from_tle(self, tle):
        """
        Add a satellite from a TLE

        :param tle: 
        :returns: 
        :rtype: 

        """

        assert isinstance(tle, TLE)

        self.add_satellite(Satellite.from_TLE(tle))

    @classmethod
    def from_constellation(cls, num_sats, num_planes, phasing, inclination, altitude, eccentricity, name="det", normal_pointing=False):
        """
        Build a satellite collection from a constellation of satellites

        :param cls: 
        :param num_sats: 
        :param num_planes: 
        :param phasing: 
        :param inclination: 
        :param altitude: 
        :param name: 
        :param normal_pointing: 
        :returns: 
        :rtype: 

        """

        constellation = Constellation(
            num_sats=num_sats, num_planes=num_planes, phasing=phasing, inclination=inclination, altitude=altitude, name=name, eccentricity=eccentricity)

        return cls(*constellation.satellites, normal_pointing=normal_pointing)

    def as_dict(self, names=None):
        """
        Build a dict of the satellites for
        a Universe object

        :returns: 
        :rtype: 

        """

        group_dict = {}
        group_dict["seed"] = 1234

        grb_dict = {}

        Grb_dict["ra"] = 80
        grb_dict["dec"] = -30
        grb_dict["distance"] = 500
        grb_dict["K"] = 500
        grb_dict["t_rise"] = 1
        grb_dict["t_decay"] = 7

        group_dict["grb"] = grb_dict

        det_dict = {}

        if names is None:

            names = list(self._satellites.keys())

        for name, sat in self._satellites.items():

            if name in names:

                sat_dict = {}
                sat_dict["ra"] = sat.ra
                sat_dict["dec"] = sat.dec
                sat_dict["altitude"] = sat.altitude
                sat_dict["time"] = '2010-01-01T00:00:00'

                if self._normal_pointing:

                    sat_dict["pointing"] = dict(ra=sat.ra, dec=sat.dec)

                else:

                    sat_dict["pointing"] = dict(ra=80, dec=-30)

                sat_dict["effective_area"] = 1.

                det_dict[name] = sat_dict

        group_dict["detectors"] = det_dict

        return group_dict

    def __add__(self, other):

        sats = list(other.satellites.values())
        sats.extend(list(self._satellites.values()))

        return SatelliteCollection(*sats)

    def write_to(self, file_name, names=None):
        """
        write out a simulation file with all satellites
        and a proxy GRB to YAML file

        :param file_name: 
        :returns: 
        :rtype: 

        """

        with open(file_name, "w") as f:

            yaml.dump(stream=f, data=self.as_dict(
                names), Dumper=yaml.SafeDumper)

    def display(self, earth_time="day", obs_time='2010-01-01T00:00:00', names=None, size=10, color="yellow"):

        fig = ipv.figure()
        ipv.pylab.style.box_off()
        ipv.pylab.style.axes_off()
        ipv.pylab.style.set_style_dark()
        ipv.pylab.style.background_color("black")

        tt = atime.Time(obs_time)

        earth = Earth(
            earth_time=earth_time, realistic=True, astro_time=tt,
        )

        earth.plot()

        x = []
        y = []
        z = []

        distances = []

        for name, sat in self._satellites.items():

            add_sat = False

            if names is not None:

                if name in names:

                    add_sat = True
            else:

                add_sat = True

            if add_sat:

                x.append(sat.xyz[0])
                y.append(sat.xyz[1])
                z.append(sat.xyz[2])

                distances.append(sat.true_alt)

        ipv.pylab.scatter(np.array(x), np.array(
            y), np.array(z), marker='sphere', color=color, size=size)

        ipv.xyzlim(max(distances))

        ipv.show()

        return fig


class Constellation(object):
    """
    Class for describing and holding a constellation of satellites
    """

    def __init__(self, num_sats, num_planes, phasing, inclination, altitude,
                 eccentricity, name="det"):

        focus = "earth"
        starting_number = 0

        self._num_sats = num_sats
        self._num_planes = num_planes
        self._phasing = phasing
        self._inclination = inclination
        self._altitude = altitude
        self._e = eccentricity

        self._focus = focus
        self._start_num = starting_number
        self._constellation_name = name

        self._sats_per_plane, self._correct_phasing = self._corrected_planes()
        self._perigee_positions = self._perigee_positions()
        self._raan = self._calculate_raan()
        self._ta = self._calculate_ta()
        self._satellites = self._build_satellites()

    def _corrected_planes(self):
        sats_per_plane = int(self._num_sats / self._num_planes)
        corrected_phasing = 360 * self._phasing / self._num_sats
        return sats_per_plane, corrected_phasing

    def _perigee_positions(self):
        perigees = list(range(0, 360, int(360/self._sats_per_plane)))
        all_perigees = []
        for i in range(self._num_sats):
            all_perigees.extend(perigees)
        return all_perigees

    def _calculate_raan(self):
        raan = [0] * self._num_sats
        for i in range(self._sats_per_plane, self._num_sats):
            raan[i] = raan[i - self._sats_per_plane] + 360 / self._num_planes
        return raan

    def _calculate_ta(self):
        ta = [0] * self._num_sats
        for i in range(self._sats_per_plane, self._num_sats):
            ta[i] = ta[i - self._sats_per_plane] + self._correct_phasing
        return ta

    def _build_satellites(self):
        satellites = []
        for i in range(self._num_sats):
            sat_num = i + self._start_num + 1
            sat_name = f"{self._constellation_name}{sat_num}"
            satellites.append(Satellite(sat_name, self._altitude, self._e, self._inclination, self._raan[i],
                                        self._perigee_positions[i], self._ta[i], focus=self._focus))
        return satellites

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}, name={7}, starting_number={8}".format(self._num_sats, self._num_planes,
                                                                                         self._phasing, self._inclination,
                                                                                         self._altitude, self._e,
                                                                                         self._constellation_name,
                                                                                         self._start_num)

    @property
    def satellites(self):

        return self._satellites

    def __str__(self):
        sat_string = ""
        for sat in self.satellites:
            sat_string += sat.__str__() + '\n'

        return sat_string.rstrip()


class Satellite(object):

    def __init__(self, name, altitude, eccentricity, inclination, right_ascension, perigee, ta,
                 focus="earth", rads=True):
        """FIXME! briefly describe function

        :param name: 
        :param altitude: 
        :param eccentricity: 
        :param inclination: 
        :param right_ascension: 
        :param perigee: 
        :param ta: 

        :param focus: 
        :param rads: 
        :returns: 
        :rtype: 

        """

        self._name = name
        self._altitude = altitude
        self._focus = focus
        self._true_alt = self.altitude + self._get_radius()
        self._eccentricity = eccentricity

        if not rads:
            self._inclination = inclination
            self._right_ascension = right_ascension
            self._perigee = perigee
            self._ta = ta
            self._inclination_r, self._right_ascension_r, self._perigee_r, self._ta_r = self._convert_to_rads()
        else:
            self._inclination_r = inclination
            self._right_ascension_r = right_ascension
            self._perigee_r = perigee
            self._ta_r = ta
            self._inclination, self._right_ascension, self._perigee, self._ta = self._convert_to_degs()

        self._compute_position()

    def _compute_position(self):

        # add this to orbital

        ke = KeplerianElements.with_altitude(self._altitude*1000.,
                                             e=self._eccentricity,
                                             i=self._inclination_r,
                                             arg_pe=self._perigee_r,
                                             raan=self._right_ascension_r,
                                             body=earth
                                             )

        # natural output is in meters so convert

        coord = SkyCoord(x=ke.r.x/1000., y=ke.r.y/1000., z=ke.r.z /
                         1000., unit='km', frame='gcrs', representation_type='cartesian')

        self._ra, self._dec = float(
            coord.spherical.lon.deg), float(coord.spherical.lat.deg)

        self._xyz = np.array([ke.r.x, ke.r.y, ke.r.z])/1000.

    @classmethod
    def from_TLE(cls, tle):

        name = clean_name(tle.name)

        orbit = tle.to_orbit()

        altitude = np.sqrt((orbit.r**2).sum()).to('km').value - 6371

        return cls(name=name,
                   altitude=altitude,  # this is in km
                   eccentricity=tle.ecc,
                   inclination=tle.inc,
                   right_ascension=tle.raan,
                   perigee=tle.argp,
                   ta=0,
                   rads=False


                   )

    @property
    def name(self):
        return str(self._name)

    @property
    def altitude(self):
        return float(self._altitude)

    @property
    def ra(self):
        return float(self._ra)

    @property
    def dec(self):
        return float(self._dec)

    @property
    def xyz(self):
        """
        cartesian coordinates in km

        :returns: 
        :rtype: 

        """

        return self._xyz

    @property
    def true_alt(self):
        return float(self._true_alt)

    @property
    def eccentricity(self):
        return float(self._eccentricity)

    def _convert_to_rads(self, value=None):
        to_rad = pi / 180
        if value:
            return value * to_rad
        else:
            return self._inclination * to_rad, self._right_ascension * to_rad, self._perigee * to_rad, self._ta * to_rad

    def _convert_to_degs(self, value=None):
        to_deg = 180 / pi
        if value:
            return value * to_deg
        else:
            return self._inclination_r * to_deg, self._right_ascension_r * to_deg, self._perigee_r * to_deg, \
                self._ta_r * to_deg

    def _get_radius(self):
        return heavenly_body_radius[self._focus.lower()]

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}".format(self._name, self._altitude, self._eccentricity,
                                                          self._inclination, self._right_ascension, self._perigee, self._ta)

    def __str__(self):
        return "Satellite Name: {0}, Alt: {1}, e: {2}, " \
               "Inclination: {3}, RA: {4}, Periapsis: {5}, Anomaly: {6}".format(self._name, self._altitude,
                                                                                self._eccentricity, self._inclination,
                                                                                self._right_ascension, self._perigee,
                                                                                self._ta)


def clean_name(name):
    return name.replace(' ', '_').replace('(', '').replace(')', '')
