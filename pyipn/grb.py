import numpy as np
try:
    from .geometry import GRBLocation
except Exception: #ImportError
    from geometry import GRBLocation


class GRB(object):
    def __init__(self, ra, dec, distance, K, t_rise, t_decay):
        """

        A GRB that emits a spectrum as a given location

        :param ra: RA of the GRB
        :param dec: DEC of the GRB
        :param distance: distance to the GRB
        :param K: normalization of the flux
        :param t_rise: rise time of the flux
        :param t_decay: decay time of the flux
        :returns: 
        :rtype: 

        """

        # create a GRB location

        self._location = GRBLocation(ra, dec, distance)

        self._K = K
        self._t_rise = t_rise
        self._t_decay = t_decay

    @property
    def pulse_parameters(self):
        """
        The temporal flux parameters
        :returns: (K, t_rise, t_decay) 
        :rtype: tuple

        """

        return self._K, self._t_rise, self._t_decay

    @property
    def location(self):
        return self._location
