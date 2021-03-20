import astropy.units as u
from astropy.coordinates import SkyCoord, CartesianRepresentation
import numpy as np


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

        self._skycoord = SkyCoord(ra, dec, unit="deg", frame="icrs")

    def get_seperation_angle(self, grb):
        """
        get the angular seperation of this 
        pointing to the GRB

        :param grb: a GRB object
        :returns: 
        :rtype: 

        """

        return self._skycoord.separation(grb.location.coord).rad

    @property
    def cartesian(self):
        return self._skycoord.cartesian.xyz.value

    @property
    def coord(self):
        return self._skycoord


# why pointing and location class different
class Location(object):
    def __init__(self, sky_coord, base_frame="icrs"):

        self._skycoord = sky_coord
        self._base_frame = base_frame

    def get_light_travel_time(self, other_location):

        pass

    @property
    def coord(self):

        return self._skycoord

    def get_cartesian_coord(self, frame=None):

        if frame is None:

            frame = self._base_frame

        return self._skycoord.transform_to(frame).represent_as(CartesianRepresentation)

    def get_norm_vec(self, unit, frame=None):
        assert isinstance(
            unit, u.Unit
        ), "no (astropy) unit provided to get_morm_vec function!"
        vec = self.get_cartesian_coord(frame).xyz.to(unit)
        norm_vec = vec / np.linalg.norm(vec)
        return norm_vec


class GRBLocation(Location):
    def __init__(self, ra, dec, distance):

        sky_coord = SkyCoord(ra * u.deg, dec * u.deg, distance=distance, frame="icrs")

        super(GRBLocation, self).__init__(sky_coord, base_frame="icrs")


class DetectorLocation(Location):
    _EARTH_RADIUS = 6371 * u.km

    def __init__(self, ra, dec, altitude, obs_time):
        """FIXME! briefly describe function

        :param ra: 
        :param dec: 
        :param altitude: 
        :param obs_time: 
        :returns: 
        :rtype: 

        """

        assert isinstance(altitude, u.Quantity), "Altitude must be a qauntity"

        # compute the distance of the detector from the center of
        # the earth via it's altitude

        distance = DetectorLocation._EARTH_RADIUS + altitude

        self._altitude = altitude

        # create a sky coordinate for the detector

        sky_coord = SkyCoord(
            ra * u.deg, dec * u.deg, distance=distance, equinox=obs_time, frame="gcrs"
        )

        super(DetectorLocation, self).__init__(sky_coord, base_frame="gcrs")

    @property
    def altitude(self):
        return self._altitude
