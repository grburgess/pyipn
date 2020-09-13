from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.time import Time
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes.formatter_locator import AngleFormatterLocator
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from astropy import units as u
from matplotlib import rcParams
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, PathPatch
from matplotlib.projections import projection_registry
import numpy as np
from reproject import reproject_from_healpix
from scipy.ndimage import gaussian_filter
import scipy.optimize
from .angle import reference_angle_deg

__all__ = (
    "AstroDegreesAitoffAllSkyAxes",
    "AstroDegreesMollweideAllSkyAxes",
    "AstroHoursAitoffAllSkyAxes",
    "AstroHoursMollweideAllSkyAxes",
    "AutoScaledWCSAxes",
    "GeoDegreesAitoffAllSkyAxes",
    "GeoDegreesMollweideAllSkyAxes",
    "GeoHoursAitoffAllSkyAxes",
    "GeoHoursMollweideAllSkyAxes",
    "GlobeAxes",
    "ScaleBar",
    "ZoomSkyAxes",
)


def create_skw_dict(projection, center=None, radius=None, rotate=None):
    assert projection in [
        "astro degrees aitoff",
        "astro degrees mollweide",
        "astro hours aitoff",
        "astro hours mollweide",
        "astro globe",
        "geo globe",
        "astro zoom",
    ]

    skw_dict = dict(projection=projection)

    if projection in ["astro globe", "astro zoom", "geo globe"]:

        assert center is not None, "you must specify a center"

        skw_dict = dict(projection=projection, center=center)

    if projection == "astro zoom":

        assert radius is not None, "you must specify a radius"

        skw_dict = dict(
            projection=projection, center=center, radius=radius, rotate=rotate
        )

    elif center is not None:

        skw_dict["center"] = center

    return skw_dict
