import pytest
import os
from pyipn import Universe, copy_template, Fit
from pyipn.io.package_utils import get_path_of_data_file


@pytest.fixture(scope="session")
def universe():

    copy_template()
    uni = Universe.from_yaml("template_config.yaml")
    uni.explode_grb(tstart=-50, tstop=50)

    yield uni

    os.remove("template_config.yaml")

    os.remove("test_universe_save.h5")


@pytest.fixture(scope="session")
def fit():

    fit = Fit.from_netcdf(get_path_of_data_file("test_fit.h5"))

    yield fit
