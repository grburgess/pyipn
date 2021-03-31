from pyipn.constellation import SatelliteCollection


def test_build_celestar():

    sc1 = SatelliteCollection.from_celestrack("galileo")
    sc1.display()

    assert sc1.n_satellites == 26

    sc2 = SatelliteCollection.from_celestrack("cubesat")

    sc3 = sc1 + sc2

    sc3.display()

    sc3.satellites


def test_constellation():

    sc = SatelliteCollection.from_constellation(
        num_sats=10,
        num_planes=1,
        phasing=1,
        inclination=50,
        altitude=500,
        eccentricity=1,
    )

    assert sc.n_satellites == 10

    sc.as_dict()
