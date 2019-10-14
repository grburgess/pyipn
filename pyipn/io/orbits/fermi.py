#package to calculate satellite positions based on TLE  
from skyfield.api import Topos, load

#package to send query requests to space-track.org
from spacetrack import SpaceTrackClient
st = SpaceTrackClient('singhartinger.moritz@gmail.com', 'WhydoIneedsuchalongpassword')

import spacetrack.operators as op
import datetime as dt
drange = op.inclusive_range(dt.datetime(2018, 1, 1),
                            dt.datetime(2018, 1, 31))

print(st.tle(norad_cat_id=33053, epoch=drange, format='tle'))

