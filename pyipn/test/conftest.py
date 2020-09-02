import pytest

from pyipn import Universe, copy_template




@pytest.fixture(scope="session")
def universe():

    copy_template()
    uni = Universe.from_yaml("template_config.yaml")
    uni.explode_grb(tstart=-50, tstop=50)

    yield uni
