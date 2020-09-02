import numpy as np
import collections
import pandas as pd

from .geometry import GRBLocation


class GRB(object):
    def __init__(self, ra, dec, distance, K, t_rise, t_decay, t_start=None):
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
        self._t_start = t_start

    @property
    def pulse_parameters(self):
        """
        The temporal flux parameters
        :returns: (K, t_rise, t_decay) 
        :rtype: tuple

        """

        return self._K, self._t_rise, self._t_decay, self._t_start

    @property
    def location(self):
        return self._location

    @property
    def table(self):

        output = collections.OrderedDict()

        position = f"{self._location.coord.ra.deg},{self._location.coord.dec.deg}"

        
        
        output["location"] = position
        output["K"] = ",".join([f"{x}" for x in  np.atleast_1d(self._K)])
        if self._t_start is not None:
            output[r"$t_s$"] = ",".join([f"{x}" for x in  np.atleast_1d(self._t_start)])
        
        output[r"$t_r$"] = ",".join([f"{x}" for x in  np.atleast_1d(self._t_rise)])
        output[r"$t_d$"] = ",".join([f"{x}" for x in  np.atleast_1d(self._t_decay)])

        return pd.Series(output)
