import arviz as av
import numpy as np
import numba as nb
import matplotlib.pyplot as plt
import astropy.units as u

from pyipn.rff import RFF, RFF_multiscale
from .utils.timing import compute_annulus_from_time_delay
from .io.plotting.projection import *
from .io.plotting.projection import create_skw_dict
from .io.plotting.markers import reticle


class Fit(object):
    def __init__(self, inference_data, universe_save=None):

        self._posterior = inference_data

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
        ).values[1]

        self._amplitude = inference_data.posterior.amplitude.stack(
            sample=("chain", "draw")
        ).values

        self._background = inference_data.posterior.bkg.stack(
            sample=("chain", "draw")
        ).values

        self._scale = inference_data.posterior.scale.stack(
            sample=("chain", "draw")
        ).values

        if self._scale.shape[0] == 2:

            self._multi_scale = True

        else:

            self._multi_scale = False

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

            self._n_dets = self._background.shape[0]

        except:

            self._is_dt_fit = False

            self._dt = None
            self._grb_theta = None
            self._grb_phi = None

            self._n_dets = 1

        self._n_samples = self._beta1.shape[-1]

        self._use_bw = False

        try:

            self._bw1 = inference_data.posterior.bw1.stack(
                sample=("chain", "draw")
            ).values

            self._bw2 = inference_data.posterior.bw2.stack(
                sample=("chain", "draw")
            ).values

            self._multi_bw = True

        except:

            try:

                self._bw = inference_data.posterior.bw.stack(
                    sample=("chain", "draw")
                ).values

                if self._bw.shape[0] == 2:

                    self._multi_bw = True

                else:

                    self._multi_bw = False

            except:

                self._bw = inference_data.posterior.bw_out.stack(
                    sample=("chain", "draw")
                ).values

                self._use_bw = False
                self._multi_bw = True

        self.grb_color = "#06DC94"
        self._grb_style = "lrtb"
        self._has_universe = False

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
        self._has_universe = True

    def expected_rate(self, time, detector):

        self._detector_check(detector)

        if self._is_dt_fit and (detector > 0):

            dt = self._dt[detector - 1]

        else:

            dt = np.zeros(self._n_samples)

        if not self._is_dt_fit:

            amp = self._amplitude

        else:

            amp = self._amplitude[detector]

        if self._use_bw:

            bw = self._bw

        else:

            bw = np.ones(self._n_samples)

        if not self._multi_scale:
            out = _expeced_rate(
                time,
                self._omega1,
                self._omega2,
                self._beta1,
                self._beta2,
                bw,
                self._scale,
                self._amplitude[detector],
                dt,
                self._n_samples,
            )
        else:

            out = _expeced_rate_multiscale(
                time,
                self._omega1,
                self._omega2,
                self._beta1,
                self._beta2,
                bw,
                self._scale,
                self._amplitude[detector],
                dt,
                self._n_samples,
            )

        return out

    def _contour_two_detectors(
        self, levels=[0.68], colors=["green"], ax=None, **kwargs
    ):

        dt = self._dt[0]

        assert len(levels) == len(colors)

        dkey = list(self._universe.detectors.keys())

        d1 = self._universe.detectors[dkey[0]]
        d2 = self._universe.detectors[dkey[1]]

        for i, level in enumerate(levels):

            dt1, dt2 = av.hdi(dt, hdi_prob=level)

            compute_annulus_from_time_delay(
                dt1*u.s, dt2*u.s, d1, d2, color=colors[i], ax=ax, **kwargs
            )

    def _show_grb(self, ax):

        ax.plot(
            self._universe.grb.location.coord.ra.deg,
            self._universe.grb.location.coord.dec.deg,
            transform=ax.get_transform("icrs"),
            marker=reticle(inner=0.4, which=self._grb_style),
            markersize=30,
            markeredgewidth=2,
            color=self.grb_color,
        )

    def location_contour(
        self,
        levels=[0.68],
        colors=["green"],
        ax=None,
        projection="astro degrees mollweide",
        center=None,
        radius=None,
        show_grb=True,
        **kwargs
    ):

        if ax is None:

            skw_dict = create_skw_dict(projection, center, radius)

            fig, ax = plt.subplots(subplot_kw=skw_dict)

        else:

            fig = ax.get_figure()

        if self._n_dets == 2:

            self._contour_two_detectors(levels, colors, ax=ax, **kwargs)

        if show_grb and self._has_universe:

            self._show_grb(ax)

        return fig

    def location_scatter(
        self,
        color="green",
        projection="astro degrees mollweide",
        center=None,
        radius=None,
        show_grb=True,
        **kwargs
    ):

        skw_dict = create_skw_dict(projection, center, radius)

        fig, ax = plt.subplots(subplot_kw=skw_dict)

        theta = np.rad2deg(self._grb_theta)
        phi = np.rad2deg(self._grb_phi)

        idx = phi <= 0

        phi[idx] += 360

        ax.scatter(
            phi, theta, transform=ax.get_transform("icrs"), color=color, **kwargs
        )

        if self._has_universe and show_grb:
            self._show_grb(ax)

        return fig

    def plot_light_curve_fit(self, detector, tstart, tstop, dt=0.2, thin=1):

        self._detector_check(detector)
        assert self._has_universe

        lc = self._universe.light_curves[
            list(self._universe.light_curves.keys())[detector]
        ]

        rate, edges, counts = lc.get_binned_light_curve(tstart, tstop, dt)

        exposure = np.diff(edges)

        mid_points = 0.5 * (edges[:-1] + edges[1:])

        fig, ax = plt.subplots()

        pred_rate = self.expected_rate(mid_points, detector)

        if self._is_dt_fit:

            bkg = self._background[detector]

        else:

            bkg = self._background

        for i in range(self._n_samples)[::thin]:

            ax.plot(mid_points, pred_rate[i] + bkg[i], color="r", alpha=0.05)

        ax.scatter(mid_points, rate, fc="none", ec="k")

        ax.set(xlabel="time (s)", ylabel="rate (cnts/s)")

        return fig

    @property
    def beta1(self):

        return self._beta1

    @property
    def beta2(self):

        return self._beta2

    @property
    def omega1(self):

        return self._omega1

    @property
    def omega2(self):

        return self._omega2

    @property
    def bw(self):

        return self._bw

    @property
    def bw1(self):

        return self._bw1

    @property
    def bw2(self):

        return self._bw2

    @property
    def scale(self):
        return self._scale

    @property
    def dt(self):

        return self._dt

    @property
    def grb_theta(self):

        return self._grb_theta

    @property
    def grb_phi(self):

        return self._grb_phi

    @property
    def amplitude(self):

        return self._amplitude

    @property
    def background(self):

        return self._background

    @property
    def posterior(self):

        return self._posterior


@nb.njit()
def _expeced_rate(time, omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N):

    out = np.empty((N, len(time)))

    for n in range(N):

        out[n] = amplitude[n] * RFF(
            time - dt[n],
            omega1[:, n],
            omega2[:, n],
            beta1[:, n],
            beta2[:, n],
            bw[n],
            scale[n],
        )

    return out


@nb.njit()
def _expeced_rate_multiscale(
    time, omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N
):

    out = np.empty((N, len(time)))

    for n in range(N):

        out[n] = amplitude[n] * RFF_multiscale(
            time - dt[n],
            omega1[:, n],
            omega2[:, n],
            beta1[:, n],
            beta2[:, n],
            bw[n],
            scale[0, n],
            scale[1, n],
        )

    return out
