import numpy as np
import astropy.constants as constants


from .lightcurve import LightCurve
from .possion_gen import source_poisson_generator, background_poisson_generator


class Detector(object):
    def __init__(self, location, pointing, effective_area, name):
        """
        A detector at a certain location in space,
        with a certain spatial pointing, 
        and a given effective area.

        :param location: a detector location
        :param pointing: a detector pointing
        :param effective_area: a detector effective area
        :param name: the detector's name
        :returns: 
        :rtype: 

        """

        self._effective_area = effective_area
        self._pointing = pointing
        self._location = location
        self._background_slope = 0.0
        self._background_norm = 500.0

        self._name = name

    @property
    def name(self):
        """
        The detector's name
        :returns: 
        :rtype: 

        """

        return self._name

    @property
    def effective_area(self):
        """
        the detector's effective area
        :returns: 
        :rtype: 

        """

        return self._effective_area

    @property
    def pointing(self):
        """
        the pointing of the detector

        :returns: 
        :rtype: 

        """

        return self._pointing

    @property
    def location(self):
        """
        detector's location 

        :returns: 
        :rtype: 

        """

        return self._location

    def light_travel_time(self, grb):
        """
        The light travel time to the GRB

        :param grb: A GRB
        :returns: 
        :rtype: 

        """

        # get the 3D seperation

        distance = self._location.coord.separation_3d(grb.location.coord)

        ltt = distance / constants.c

        return ltt.to("second")

    def angular_separation(self, grb):

        return self._location.coord.separation(grb.location.coord)

    def build_light_curve(self, grb, T0, tstart, tstop):
        """
        Build the light curve observed from the GRB

        :param grb: 
        :param T0: 
        :param tstart: 
        :param tstop: 
        :returns: 
        :rtype: 

        """

        # first get the seperation angle
        # from the detector pointing and
        # the grb

        seperation_angle = self._pointing.get_seperation_angle(grb)

        self._effective_area.set_seperation_angle(seperation_angle)

        # now get the grb's pulse parameters

        K, t_rise, t_decay = grb.pulse_parameters

        # scale the GRB by the effective area

        observed_intensity = K * self._effective_area.effective_area

        # compute the arrival times
        
        source_arrival_times = source_poisson_generator(
            tstart, tstop, observed_intensity, T0, t_rise, t_decay
        )

        bkg_arrival_times = background_poisson_generator(
            tstart, tstop, self._background_slope, self._background_norm
        )

        # return a light curve

        return LightCurve(source_arrival_times, bkg_arrival_times)
