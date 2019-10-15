from pathlib import Path

from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

'''
Send query request to space-track.org for a specified timerange (drange) 
and satellite id and write the corresponding two-line-elements into a data file.
'''

def write_tle(st, sat_id, drange, name):
	lines = st.tle(iter_lines=True, norad_cat_id=sat_id, epoch=drange, format='tle')
	data_path = Path().resolve().parent.parent.as_posix() +'/data/'+'fermi('+drange+')_tle.txt'

	with open(data_path,'w') as fp:
		for line in lines:
			fp.write(line + "\n")

	return(Path().resolve().parent.parent.as_posix() +'/data/'+'fermi('+drange+')_tle.txt')


def convert_to_decimal_days(dt):
	dd = dt.timetuple().tm_yday + dt.hour/24 + dt.minute/(24*60) + dt.second/(24*3600) + dt.microsecond/(24*3600*1000000)
	return dd

'''
find closest epoch TLE to correstponding datetime from specified file and return the two lines
'''
def find_closest_epoch(dt, tle_path):
	day = convert_to_decimal_days(dt)
	diff = []

	with open(tle_path, 'r') as fp:
		lines = fp.readlines()
		for line in lines:
			elem = line.split()
			if elem[0] == '1':
				epoch_day = float(elem[3][2:])
				diff.append(abs(epoch_day - day))

		val, idx = min((val, idx) for (idx, val) in enumerate(diff))
		line1 = lines[idx*2]
		line2 = lines[idx*2+1]

	return(line1, line2)

'''
return position of satellite at time dt based on TLE file
'''
def position(dt, tle_path):
	line1, line2 = find_closest_epoch(dt, tle_path)

	satellite = twoline2rv(line1, line2, wgs72)
	position, velocity = satellite.propagate(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)

	return position

