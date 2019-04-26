import numpy as np

from .geometry import Location

class GRB(object):

    def __init__(self, location, K, t_rise, t_decay):

        self._location = location

        self._K = K
        self._t_rise  = t_rise
        self._t_decay = t_decay


    @property
    def pulse_parameters(self):

        return self._K, self._t_rise, self._t_decay
        


    @property
    def location(self):
        return self._location
