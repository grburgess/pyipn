import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const
import astropy.units as u
from astropy.coordinates import SkyCoord, UnitSphericalRepresentation
from pyipn.io.plotting.projection import create_skw_dict

from pyipn.io.plotting.spherical_circle import SphericalCircle


from pyipn.geometry import Location


def theta_from_time_delay(dt, distance):

    return np.arccos(
        np.around(
            ((const.c * dt).to("km") / distance.to("km"))
            .to(u.dimensionless_unscaled)
            .value,
            15,
        )
    )


def calculate_distance_and_norm(d1, d2):

    dxyz = d2.location.get_cartesian_coord().xyz - d1.location.get_cartesian_coord().xyz

    # calculate ra and dec of vector d  pointing from detector1 to detector2
    dcart = Location(
        SkyCoord(
            x=dxyz[0],
            y=dxyz[1],
            z=dxyz[2],
            representation_type="cartesian",
            unit="km",
            frame="icrs",
        )
    )

    norm_d = dcart.get_norm_vec(u.km)
    ra = dcart.coord.represent_as(UnitSphericalRepresentation).lon
    dec = dcart.coord.represent_as(UnitSphericalRepresentation).lat
    distance = np.linalg.norm(dxyz).to("km")

    return distance, norm_d, ra, dec


def compute_annulus_from_time_delay(
    dt1,
    dt2,
    detector1,
    detector2,
    ax=None,
    color="green",
    projection=None,
    center=None,
    radius=None,
    **kwargs
):

    distance, norm_d, ra, dec = calculate_distance_and_norm(detector1, detector2)

    if ax is None:

        skw_dict = create_skw_dict(projection, center, radius)

        fig, ax = plt.subplots(subplot_kw=skw_dict)

    thetas = []

    time_delays = [dt1, dt2]

    for dt in time_delays:

        theta = theta_from_time_delay(-dt, distance * u.km )

        thetas.append(theta)

    thetas = np.array(thetas)

    circle1 = SphericalCircle(
        np.array([ra.value, dec.value]) * ra.unit,
        thetas[0] * u.rad,
        vertex_unit=u.deg,
        resolution=5000,
        edgecolor=color,
        facecolor="none",
        transform=ax.get_transform("icrs"),
        **kwargs
    )

    #     circle1 = Circle(
    #         np.array([ra.deg, dec.deg]), #* ra.unit,
    #         np.rad2deg(thetas[idx[-1]]),# * u.rad,
    #  #       vertex_unit=u.deg,
    #   #      resolution=5000,
    #         #edgecolor=color,
    #         facecolor=color,
    #     transform=ax.get_transform("icrs"),
    #         **kwargs
    #     )

    ax.add_patch(circle1)

    circle2 = SphericalCircle(
        np.array([ra.value, dec.value]) * ra.unit,
        thetas[1] * u.rad,
        vertex_unit=u.deg,
        resolution=5000,
        edgecolor=color,
        facecolor="none",
        # visible=True,
        transform=ax.get_transform("icrs"),
        **kwargs
    )

    ax.add_patch(circle2)
