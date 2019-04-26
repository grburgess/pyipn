import astropy.units as u
from astropy.coordinates import SkyCoord


class Pointing(object):

    def __init__(self, ra, dec):
        """
        The directional pointing of a detector.
        Currently, this is a static property, but could
        be extended to include time in the future

        :param ra: the RA of the pointing
        :param dec: the DEC of the pointing
        :returns: 
        :rtype: 

        """

        # build a sky_coord

        self._skycoord = SkyCoord(ra, dec, unit='deg', frame='icrs')


    def get_seperation_angle(self, grb):
        """
        get the angular seperation of this 
        pointing to the GRB

        :param grb: a GRB object
        :returns: 
        :rtype: 

        """
        

        return self._skycoord.separation(grb.location.coord).rad
        
        
        


class Location(object):
    def __init__(self, sky_coord):
        
        self._skycoord

    def get_light_travel_time(self, other_location):

        pass


    @property
    def coord(self):

        return self._skycoord

class GRBLocation(Location):
    def __init__(self, ra, dec, distance):

        sky_coord = SkyCoord(ra, dec, unit="deg", distance=distance, frame="icrs")

        super(DetectorLocation, self).__init__(sky_coord)


class DetectorLocation(Location):
    _EARTH_RADIUS = 6700 * u.km

    def __init__(self, lat, lon, altitude, obs_time):
        """FIXME! briefly describe function

        :param lat: 
        :param lon: 
        :param altitude: 
        :param obs_time: 
        :returns: 
        :rtype: 

        """

        assert isinstance(altitude, u.Quantity), "Altitude must be a qauntity"

        # compute the distance of the detector from the center of
        # the earth via it's altitude

        distance = DetectorLocation._EARTH_RADIUS + altitude

        # create a sky coordinate for the detector

        sky_coord = SkyCoord(
            lat,
            long,
            unit="deg",
            distance=distance,
            equinox=obs_time,
            frame="geocentric",
        )

        super(DetectorLocation, self).__init__(sky_coord)
