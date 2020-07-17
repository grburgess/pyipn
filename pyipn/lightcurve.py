import numpy as np
import matplotlib.pyplot as plt


class LightCurve(object):
    def __init__(self, source_arrival_times, bkg_arrival_times):
        """
        A light curve containing events

        :param source_arrival_times: the arrival times 
        :param bkg_arrival_times: the arrival times from the background
        :returns: 
        :rtype: 

        """

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

        counts, edges = np.histogram(self._arrival_times, bins=bins)

        # compute the rate

        rate = counts / dt

        return rate, edges, counts

    def display(self, tstart, tstop, dt, ax=None,**kwargs):
        """FIXME! briefly describe function

        :param tstart: 
        :param tstop: 
        :param dt: 
        :returns: 
        :rtype: 

        """

        rate, edges, _ = self.get_binned_light_curve(tstart, tstop, dt)

        err = np.sqrt(rate)

        mid_points = np.mean([edges[1:], edges[:-1]], axis=0)

        if ax is None:
        
            fig, ax = plt.subplots()

        else:

            fig = ax.get_figure()
            
        ax.hlines(rate, edges[:-1], edges[1:], **kwargs)

        ax.vlines(mid_points, rate - err, rate + err, **kwargs)

        ax.set_ylabel("rate")
        ax.set_xlabel("time")

        return fig
