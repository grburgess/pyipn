from pyipn.io.orbits.tle import position_pyorbital, position_skyfield
from pyipn.geometry import Location

from pathlib import Path
import datetime as dt

def test_position_pyorbital():
	tle_path = Path(__file__).resolve().parent.parent.as_posix() +'/data/GLAST2018-01-01 00:00:00--2018-12-31 00:00:00_tle.txt'
	date = dt.datetime(2018, 6, 13)

	pos  = position_pyorbital(date, tle_path, 'GLAST')

def test_position_skyfield():
	tle_path_G = Path(__file__).resolve().parent.parent.as_posix() +'/data/GLAST2018-01-01 00:00:00--2018-12-31 00:00:00_tle.txt'
	tle_path_I = Path(__file__).resolve().parent.parent.as_posix() +'/data/INTEGRAL2018-01-01 00:00:00--2018-12-31 00:00:00_tle.txt'
	date = dt.datetime(2018, 6, 13)

	pos_G = position_skyfield(date, tle_path_G)
	pos_I = position_skyfield(date, tle_path_I)
	loc_G = Location.from_GCRS(pos_G)
	loc_I = Location.from_GCRS(pos_I)


