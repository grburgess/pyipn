import collections
from math import pi

import numpy as np
from astropy.coordinates import SkyCoord
from orbital import KeplerianElements, earth

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


class Constellation(object):
    """
    Class for describing and holding a constellation of satellites
    """

    def __init__(self, num_sats, num_planes, phasing, inclination, altitude,
                 eccentricity, name="det", focus="earth", starting_number=0, normal_pointing=False):
        self._num_sats = num_sats
        self._num_planes = num_planes
        self._phasing = phasing
        self._inclination = inclination
        self._altitude = altitude
        self._e = eccentricity
        self._normal_pointing = normal_pointing

        self._focus=focus
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

    def __str__(self):
        sat_string = ""
        for sat in self.satellites:
            sat_string += sat.__str__() + '\n'

        return sat_string.rstrip()

    def as_dict(self):
        constellation = {}
        constellation["seed"] = 1234

        grb_dict = {}

        grb_dict["ra"] = 80
        grb_dict["dec"] = -30
        grb_dict["distance"] = 500
        grb_dict["K"] = 500
        grb_dict["t_rise"] = 1
        grb_dict["t_decay"] = 7

        constellation["grb"] = grb_dict

        det_dict = {}

        for sat in self._satellites:
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

            det_dict[sat.name] = sat_dict

        constellation["detectors"] = det_dict
            
        return constellation


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


        
        self._ra, self._dec = float(coord.spherical.lon.deg), float(coord.spherical.lat.deg)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def altitude(self):
        return self._altitude


    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    
    @property
    def true_alt(self):
        return self._true_alt

    @property
    def eccentricity(self):
        return self._eccentricity

    @eccentricity.setter
    def eccentricity(self, new_e):
        if new_e < 0:
            return ValueError("Eccentricity can't be set below a perfect circle.")
        else:
            self._eccentricity = new_e

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
        return "{0}, {1}, {2}, {3}, {4}, {5}, {6}".format(self.name, self.altitude, self.eccentricity,
                                                          self.inclination, self.right_ascension, self.perigee, self.ta)

    def __str__(self):
        return "Satellite Name: {0}, Alt: {1}, e: {2}, " \
               "Inclination: {3}, RA: {4}, Periapsis: {5}, Anomaly: {6}".format(self.name, self.altitude,
                                                                                self.eccentricity, self.inclination,
                                                                                self.right_ascension, self.perigee,
                                                                                self.ta)
