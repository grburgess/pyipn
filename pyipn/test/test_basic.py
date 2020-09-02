from pyipn import Universe, copy_template
import pytest


def test_simple(universe):


    for det, lc in universe.light_curves.items():

        lc.display(-10, 50, 1.0)
    print(universe.calculate_annulus("det1", "det2"))

    print(universe.T0)

    print(universe.table)

    universe.plot_all_annuli()
