import numpy as np
import astropy.constants as constants


from .lightcurve import LightCurve
from .possion_gen import source_poisson_generator, background_poisson_generator


class Detector(object):

    def __init__(self, location, pointing, effective_area, name):


        self._effective_area = effective_area
        self._pointing = pointing
        self._location = location
        self._background_slope = 0.
        self._background_norm = 50.
        
        self._name = name


    @property
    def name(self):

        return self._name
        
    @property
    def effective_area(self):

        return self._effective_area

    @property
    def pointing(self):

        return self._pointing

    @property
    def location(self):

        return self._pointing


    def light_travel_time(self, grb):


        # get the 3D seperation

        distance = self._location.coord.separation_3d(grb.location.coord)

        ltt = distance/constants.c

        return ltt
        
    
    def build_light_curve(self, grb, T0, tstart, tstop):


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
        
        source_arrival_times = source_poisson_generator(tstart, tstop, observed_intensity, T0, t_rise, t_decay)
        

        bkg_arrival_times = background_poisson_generator(tstart, tstop, self._background_slope, self._background_norm)


        # return a light curve

        return LightCurve(source_arrival_times, bkg_arrival_times)
        
        
