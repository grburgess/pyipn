import arviz as av
import numpy as np
import numba as nb
from pyipn.rff import RFF


class Fit(object):
    def __init__(self, inference_data, universe_save=None):

        self._beta1 = inference_data.posterior.beta1.stack(
            sample=("chain", "draw")
        ).values
        self._beta2 = inference_data.posterior.beta2.stack(
            sample=("chain", "draw")
        ).values

        self._omega1 = inference_data.posterior.omega.stack(
            sample=("chain", "draw")
        ).values[0]

        self._omega2 = inference_data.posterior.omega.stack(
            sample=("chain", "draw")
        ).values[0]

        self._amplitude = inference_data.posterior.amplitude.stack(
            sample=("chain", "draw")
        ).values
        self._background = inference_data.posterior.background.stack(
            sample=("chain", "draw")
        ).values

        try:

            self._dt = inference_data.posterior.dt.stack(
                sample=("chain", "draw")
            ).values

            self._grb_theta = inference_data.posterior.grb_theta.stack(
                sample=("chain", "draw")
            ).values
            self._grb_phi = inference_data.posterior.grb_phi.stack(
                sample=("chain", "draw")
            ).values

            self._is_dt_fit = True

        except:

            self._is_dt_fit = False

            self._dt = None
            self._grb_theta = None
            self._grb_phi = None

        self._n_samples = self._beta1.shape[-1]

    def _detector_check(self, det_number):

        assert det_number in range(self._n_dets)

    @classmethod
    def from_cmdstanpy(cls, fit):

        inference_data = av.from_cmdstanpy(fit)

        return cls(inference_data)

    @classmethod
    def from_netcdf(cls, file_name):

        inference_data = av.from_netcdf(file_name)

        return cls(inference_data)

    def set_universe(self, universe):

        self._universe = universe

    def expected_rate(self, time, detector):

        self._detector_check(detector)

        if self._is_dt_fit:

            if detector > 0:

                dt = self._dt[detector]

        else:

            dt = np.zeros(self._n_samples)

        out = _expeced_rate(
            time,
            self._omega1,
            self._omega2,
            self._beta1,
            self._beta2,
            self._bw,
            self._scale,
            self._amplitude,
            dt,
            self._n_samples,
        )

        return out

    def plot_location_fit(self):
        pass

    def plot_light_curve_fit(self, detector,tstart, tstop, dt=.2):

        self._detector_check(detector)
        assert self._has_universe

        lc = self._universe.lightcurves[list(self._lightcurves.keys())[detector]]

        rate, edges, counts = lc.get_binned_light_curve(tstart, tstop, dt)

        exposure = np.diff(edges)

        mid_points = 0.5 * (edges[:-1] + edges[1:])

        fig, ax = plt.subplots()


        ax.scatter(mid_points, rate,fc="none", ec="k")


        pred_rate = self.expected_rate(mid_points, time)

        for i in range(self._n_samples):

            ax.plot(mid_points, pred_rate[i] + self._background[i])
        

        return fig

    @property
    def beta1(self):

        return self._beta1

    @property
    def beta2(self):

        return self._beta1

    @property
    def omega1(self):

        return self._omega1

    @property
    def omega2(self):

        return self._omega1

    @property
    def dt(self):

        return self._dt

    @property
    def grb_theta(self):

        return self._grb_theta

    @property
    def grb_phi(self):

        return self._grb_phi


@nb.njit()
def _expeced_rate(time, omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N):

    out = np.empty((N, len(time)))

    for n in range(N):

        out[n] = amplitude[n] * RFF(
            time - dt[n], omega1[n], omega2[n], beta1[n], beta2[n], bw[n], scale[n]
        )

    return out
