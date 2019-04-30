import numpy as np
import matplotlib.pyplot as plt


class LightCurve(object):
    def __init__(self, source_arrival_times, bkg_arrival_times):

        self._arrival_times = np.concatenate([source_arrival_times, bkg_arrival_times])

        self._source_arrival_time = source_arrival_times

        self._bkg_arrival_times = bkg_arrival_times

    def get_binned_light_curve(self, tstart, tstop, dt):
        """FIXME! briefly describe function

        :param tstart: 
        :param tstop: 
        :param dt: 
        :returns: 
        :rtype: 

        """

        # create the bins

        bins = np.arange(tstart, tstop, dt)

        # histogram the counts

        hist, edges = np.histogram(self._arrival_times, bins=bins)

        # compute the rate

        rate = hist / dt

        return rate, edges

    def display(self, tstart, tstop, dt, **kwargs):
        """FIXME! briefly describe function

        :param tstart: 
        :param tstop: 
        :param dt: 
        :returns: 
        :rtype: 

        """

        rate, edges = self.get_binned_light_curve(tstart, tstop, dt)

        err = np.sqrt(rate)

        mid_points = np.mean([edges[1:], edges[:-1]], axis=0)

        fig, ax = plt.subplots()

        ax.hlines(rate, edges[:-1], edges[1:], **kwargs)

        ax.vlines(mid_points, rate - err, rate + err, **kwargs)

        ax.set_ylabel("rate")
        ax.set_xlabel("time")
