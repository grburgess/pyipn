import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u

class LightCurve(object):
    def __init__(self, source_arrival_times, bkg_arrival_times):
        """
        A light curve containing events

        :param source_arrival_times: the arrival times 
        :param bkg_arrival_times: the arrival times from the background
        :returns: 
        :rtype: 

        """

        self._arrival_times = np.concatenate(
            [source_arrival_times, bkg_arrival_times])

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

        rate = counts / float(dt)

        return rate, edges, counts

    def display(self, tstart, tstop, dt, ax=None, step=False, **kwargs):
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

        # ax.hlines(rate, edges[:-1], edges[1:], **kwargs)

        # ax.vlines(mid_points, rate - err, rate + err, **kwargs)

        if step:
            ax.step(mid_points, rate, **kwargs)
        else:
            ax.scatter(mid_points, rate, **kwargs)

        ax.set_ylabel("rate")
        ax.set_xlabel("time")

        return fig

    @property
    def source_arrival_times(self):

        return self._source_arrival_time

    @property
    def bkg_arrival_times(self):

        return self._bkg_arrival_times


class BinnedLightCurve(object):
    def __init__(self, counts, time_bins, tstart, tstop, dt):

        assert len(counts) == len(time_bins) - 1
        assert dt > 0

        self._counts = counts
        self._time_bins = time_bins
        self._dt = dt
        self._dt_ms = (dt*u.s).to("millisecond").value
        
        self._tstart = tstart
        self._tstop = tstop
        self._n_bins = len(counts)

    @property
    def counts(self):
        return self._counts

    @property
    def time_bins(self):
        return self._time_bins

    @property
    def dt(self):
        return self._dt

    @property
    def res_ms(self):
        return self._dt_ms
    
    @property
    def tstart(self):
        return self._tstart

    @property
    def tstop(self):
        return self._tstop

    @property
    def n_bins(self):
        return int(self._n_bins)

    def time2idx(self, t):
        return np.searchsorted(self._time_bins,t)


    def get_src_counts(self, bkg_rate=500.):

        return self._counts - bkg_rate * self._dt
    
    def get_max_sn(self, dt_ms, bkg_rate=500):
        """
        From D. Svinkn

        :param dt_ms: 
        :returns: 
        :rtype: 

        """

        len_int =  int(dt_ms//self._dt_ms)
        #print(dt_ms, len_int)
        counts_src = self.get_src_counts(bkg_rate)
        cs = np.cumsum(counts_src)
  
        fmax = 0.0
        i1 = 0
        i2 = cs.size
        for i in range(cs.size-len_int):
            dif = cs[i+len_int] - cs[i]
            if dif > fmax:
                fmax = dif
                i1 = i
                i2 = i+len_int

        #print(fmax)

        return i1, i2
        
    
    @classmethod
    def from_lightcurve(cls, lightcurve, tstart, tstop, dt):

        _, t, c = lightcurve.get_binned_light_curve(tstart, tstop, dt)

        return cls(c, t, tstart, tstop, dt)
