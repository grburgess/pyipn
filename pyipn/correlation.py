import matplotlib.pyplot as plt
import numba as nb
import numpy as np
from scipy.stats import chi2, norm, poisson

from .lightcurve import BinnedLightCurve


class Correlator(object):
    """
    Cross-correlation metodology is descibed in Pal'shin et al. ApJS 207, 38, 2013
    Original C++ code was developed by Valentin Pal'shin.
    Python translation by D. Svinkin.
    """

    def __init__(
        self,
        lc_1,
        lc_2,
        idx_lc_1,
        idx_beg_lc2,
        idx_end_lc2,
        cl_sigma=[1, 2, 3],
        bkg_rate=500,
    ):
        """FIXME! briefly describe function

        :param lc_1: 
        :param lc_2: 
        :param idx_lc_1: 
        :param idx_beg_lc2: 
        :param idx_end_lc2: 
        :param cl_sigma: 
        :param bkg_rate: 
        :returns: 
        :rtype: 

        """

        assert isinstance(lc_1, BinnedLightCurve)
        assert isinstance(lc_2, BinnedLightCurve)

        self._idx_beg_lc2, self._idx_end_lc2 = idx_beg_lc2, idx_end_lc2
        self.lc_1, self.lc_2 = lc_1, lc_2
        self._idx_lc_1 = idx_lc_1

        dt_grb_ms = (idx_end_lc2 - idx_beg_lc2) * self.lc_2.res_ms
        self._idx_beg_lc1, self._idx_end_lc1 = self.lc_1.get_max_sn(dt_grb_ms, bkg_rate)

        self._fscale = np.sum(
            lc_1.get_src_counts()[self._idx_beg_lc1 : self._idx_end_lc1 + 1]
        ) / np.sum(lc_2.get_src_counts()[self._idx_beg_lc2 : self._idx_end_lc2 + 1])

        self._n_max = int(
            lc_1.n_bins
            - idx_lc_1
            - (idx_end_lc2 - idx_beg_lc2 + 1) * lc_2.res_ms // lc_1.res_ms
        )

        self._cl_sigma = cl_sigma

        (
            self._arr_dt,
            self._arr_chi,
            self._nDOF,
            self._fRijMin,
            self._dTmin,
            self._iMin,
            self._nMin,
        ) = correlate(
            lc_1.time_bins,
            lc_1.get_src_counts(bkg_rate),
            lc_1.counts,
            lc_2.time_bins,
            lc_2.get_src_counts(bkg_rate),
            lc_2.counts,
            idx_beg_lc2,
            idx_end_lc2,
            idx_lc_1,
            self._n_max,
            self._fscale,
            lc_1.res_ms,
            lc_2.res_ms,
        )

        self._dTlower = []
        self._dTupper = []
        self._fSigma = []

        for sigma in self._cl_sigma:

            dTlower, dTupper, fSigma = self._get_dTcc(nSigma=sigma)

            self._dTlower.append(dTlower)
            self._dTupper.append(dTupper)
            self._fSigma.append(fSigma)

    def info(self):

        dT2 = (self._idx_end_lc2 - self._idx_beg_lc2) * self.lc_2.res_ms
        dT1 = (self._idx_end_lc1 - self._idx_beg_lc1) * self.lc_1.res_ms

        str_out = "Cross-correlation info:\n"
        str_out += "CC interval for lc2 (i_beg, i_end, dT): {:4d} {:4d} {:4d} ms\n".format(
            self._idx_beg_lc2, self._idx_end_lc2, dT2
        )
        str_out += "CC interval for lc1 (i_beg, i_end, dT): {:4d} {:4d} {:4d} ms\n".format(
            self._idx_beg_lc1, self._idx_end_lc1, dT1
        )
        str_out += "Lc scale factor: {:.3f}\n".format(self._fscale)
        str_out += "Lc 1 i_start: {:d}\n".format(self._idx_lc_1)
        str_out += "Lc 1 n_max: {:d}\n".format(self._n_max)

        # str_out +="dTmin chiMin: {:8.3f} {:8.3f}\n".format(self._dTmin, self._fRijMin)
        # str_out +="dTlower dTupper fSigma: {:8.3f} {:8.3f} {:8.3f}\n".format(self._dTlower, self._dTupper, self._fSigma)

        for i, sigma in enumerate(self._cl_sigma):
            str_out += "Calculated time delay and its {:.1f} sigma CI:\n".format(sigma)
            str_out += "dTcc dTcc- dTcc+: {:.3f} {:+.3f} {:+.3f}\n".format(
                self._dTmin,
                self._dTlower[i] - self._dTmin,
                self._dTupper[i] - self._dTmin,
            )

        print(str_out)

    def _get_dTcc(self, nSigma=3):
        """Calculate cross-correlation lag confidence interval

        Parameters
        ----------
        arr_dt 
        self._arr_chi 
        self_nDOF 
        self._fRijMin 
        self._nMin 
        nSigma

        Returns
        -------
        dTlower, dTupper - lower and upper boundary of the interval
        fSigma - chi2/dof level for nSigma
        """

        P0 = norm.cdf(nSigma) - norm.cdf(-nSigma)

        fSigma = chi2.ppf(P0, self._nDOF - 1) / (self._nDOF - 1) - 1.0 + self._fRijMin
        fSigmaSimple = self._fRijMin + nSigma ** 2 / self._nDOF

        n = self._arr_chi.size

        # search upper dTcc 3 sigma
        i = self._nMin
        fRij = self._arr_chi[i]

        while (i > 1) and (fRij < fSigma):
            i = i - 1
            fRij = self._arr_chi[i]

        dTupper = self._arr_dt[i]

        # search lower dTcc 3 sigma
        i = self._nMin
        fRij = self._arr_chi[i]
        while (i < n - 1) and (fRij < fSigma):
            i = i + 1
            fRij = self._arr_chi[i]

        dTlower = self._arr_dt[i]

        return dTlower, dTupper, fSigma

    @property
    def dt_min(self):
        return self._dTmin

    @property
    def dt_lower(self):
        return self._dTlower

    @property
    def dt_upper(self):
        return self._dTupper

    @property
    def f_sigma(self):

        return self._fSigma

    def plot(self, ax=None):

        if ax is None:

            fig, ax = plt.subplots()

        else:

            fig = ax.get_figure()

        ax.set_ylabel("$\chi^2_r$")
        ax.set_xlabel(r"$\delta T$ (s)")

        ax.plot(
            self._arr_dt,
            self._arr_chi,
            marker="s",
            color="r",
            linewidth=0.5,
            markersize=5,
        )

        for i in range(len(self._cl_sigma)):

            ax.axhline(self._fSigma[i], color="r", linewidth=0.5)
            ax.vlines(
                [self._dTlower[i], self._dTupper[i]],
                [0, 0],
                [self._fSigma[i] + 3, self._fSigma[i] + 3],
                linestyles="dashed",
                color="k",
                linewidth=0.5,
            )

        ax.axvline(x=self._dTmin, linestyle="dashed", color="b", linewidth=0.5)

        idx = np.argmax(self._cl_sigma)

        dT_err = (self._dTupper[idx] - self._dTlower[idx]) / 2.0
        ax.set_xlim(self._dTlower[idx] - dT_err, self._dTupper[idx] + dT_err)

        # ax.set_yticks(np.arange(0, self.y_max+1, 1))

        # minorLocator_x = MultipleLocator(0.1)
        # minorLocator_y = MultipleLocator(0.2)
        # ax.xaxis.set_minor_locator(minorLocator_x)   # x minor ticks
        # ax.yaxis.set_minor_locator(minorLocator_y)  # y minor ticks

        # ax.set_ylim(0, self.y_max)

        # if not self.caption is None:
        #     ax.set_title(caption, fontsize=18)

        # if not self.fig_file_name is None:
        #      plt.savefig(fig_f


@nb.njit(fastmath=True)
def correlate(
    arr_t_1,
    arr_cnts_1,
    arr_cnts_err_1,
    arr_t_2,
    arr_cnts_2,
    arr_cnts_err_2,
    i_beg_2,
    i_end_2,
    i_beg_1,
    n_max_1,
    fscale,
    dt_1,
    dt_2,
    n_sum_2=1,
):
    """ Calculate chi2/dof(dTcc)

    Parameters
    ----------
    arr_t_1         array with TH 1 bin starts
    arr_cnts_1      array with TH 1 bg sub counts
    arr_cnts_err_1  array with TH 1 bg sub counts err 
    arr_t_2         array with TH 2 bin starts
    arr_cnts_2      array with TH 2 bg sub counts
    arr_cnts_err_2  array with TH 2 bg sub counts err 
    i_beg_2         start bin of TH 2 to cross-correlate
    i_end_2         end bin of TH 2 to cross-correlate
    i_beg_1         start bin of TH 1
    n_max_1         maximum shift of start of TH 1
    fscale          scale TH 2 to TH 1
    dt_1            time resolution of sc1, ms!
    dt_2            time resolution of sc2, ms!
    n_sum_2=1       number of bins of TH 2 to sum, default = 1

    Returns
    -------
    arr_dt, arr_chi, self._nDOF, self._fRijMin, dTmin, iMin, nMin

    Notes
    -----
    We shift the time history 1 !!!!
    the time history 2 is in rest
    """

    arr_dt = np.zeros(n_max_1)
    arr_chi = np.zeros(n_max_1)
    arr_sigma3 = np.zeros(n_max_1)
    arr_sigma2 = np.zeros(n_max_1)
    arr_sigma3simple = np.zeros(n_max_1)
    arr_sigma2simple = np.zeros(n_max_1)

    n = 0
    fRijMin = 1e4
    i_min = 0
    n_min = 0

    nDOF = (i_end_2 - i_beg_2) / n_sum_2  # number degrees of freedom
    fRij = 0

    for i in range(i_beg_1, i_beg_1 + n_max_1):

        fRij = 0

        fT1 = dt_1
        fT2 = 0

        l = i

        fdCnt1 = 0
        fdErr1 = 0

        for k in range(i_beg_2, i_end_2 + 1, n_sum_2):

            fCnt2 = 0
            fErr2 = 0

            for iSum2 in range(k, k + n_sum_2):
                fCnt2 += arr_cnts_2[iSum2]
                fErr2 += arr_cnts_err_2[iSum2]

            fCnt2 *= fscale
            fErr2 *= fscale * fscale

            fCnt1 = fdCnt1
            fErr1 = fdErr1

            while fT1 <= (fT2 + dt_2 * n_sum_2):
                fCnt1 += arr_cnts_1[l]
                fErr1 += arr_cnts_err_1[l]
                fT1 += dt_1
                l = l + 1

            fdT = (fT1 - (fT2 + dt_2 * n_sum_2)) / dt_1
            fdCnt1 = fdT * arr_cnts_1[l]
            fCnt1 += (1.0 - fdT) * arr_cnts_1[l]

            fdErr1 = fdT * arr_cnts_err_1[l]
            fErr1 += (1.0 - fdT) * arr_cnts_err_1[l]

            l = l + 1
            fT1 += dt_1
            fT2 += dt_2 * n_sum_2

            fDif = fCnt1 - fCnt2
            fDif *= fDif

            fDif /= fErr1 + fErr2
            fRij += fDif

        fRij /= nDOF - 1

        dT = arr_t_2[i_beg_2] - arr_t_1[i]
        arr_dt[n] = dT
        arr_chi[n] = fRij

        if fRij < fRijMin:
            fRijMin = fRij
            iMin = i
            nMin = n
            dTmin = dT

        n = n + 1

    return arr_dt, arr_chi, nDOF, fRijMin, dTmin, iMin, nMin
