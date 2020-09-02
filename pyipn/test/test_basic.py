import numpy.testing as nt
import pytest
from astropy.units import Quantity

from pyipn import BinnedLightCurve, Correlator, Universe, copy_template
from pyipn.utils.timing import calculate_distance_and_norm


def test_calculate_annulus(universe):

    d1 = universe.detectors["det1"]
    d2 = universe.detectors["det2"]

    distance, norm_d, ra, dec = calculate_distance_and_norm(d1, d2)

    assert isinstance(distance, Quantity)

    assert isinstance(norm_d, Quantity)

    print(distance)


def test_simple(universe):

    for det, lc in universe.light_curves.items():

        lc.display(-10, 50, 1.0)
    print(universe.calculate_annulus("det1", "det2"))

    print(universe.T0)

    print(universe.table)

    universe.plot_all_annuli()


def test_save_universe(universe):

    universe.write_to("test_universe_save.h5")

    universe_reloaded = Universe.from_save_file("test_universe_save.h5")

    for det, lc in universe_reloaded.light_curves.items():

        lc.display(-10, 50, 1.0)

        rate, edges, counts = lc.get_binned_light_curve(-10, 50, 1)

        lc_original = universe.light_curves[det]

        rate2, edges2, counts2 = lc_original.get_binned_light_curve(-10, 50, 1)

        nt.assert_array_equal(rate, rate2)

        nt.assert_array_equal(counts, counts2)


def test_read_fit(fit, universe):

    assert not fit._has_universe

    fit.set_universe(universe)

    assert fit._has_universe

    fit.location_contour()

    fit.location_scatter()

    fit.plot_light_curve_ppcs(0, 0, 20, dt=1.)

    fit.plot_light_curve_fit(0, 0, 20, dt=1)


def test_correlation(universe):

    lc1 = BinnedLightCurve.from_lightcurve(
        universe.light_curves["det1"], -10, 50, .1)
    lc2 = BinnedLightCurve.from_lightcurve(
        universe.light_curves["det2"], -10, 50, .1)

    t_cc_beg_1 = -5  # for lc 1
    t_cc_beg_2 = 0.  # for lc 2
    t_cc_end_2 = 20.  # for lc 2

    idx_lc_1 = lc1.time2idx(t_cc_beg_1)
    idx_beg_lc2 = lc2.time2idx(t_cc_beg_2)
    idx_end_lc2 = lc2.time2idx(t_cc_end_2)

    cc = Correlator(lc1, lc2, idx_lc_1, idx_beg_lc2,
                    idx_end_lc2, cl_sigma=[1, 2, 3])

    cc.dt_min

    assert len(cc.dt_lower) == 3
    assert len(cc.dt_upper) == 3
