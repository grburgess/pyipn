import numpy as np
import astropy.units as u
import astropy.constants as constants


from .lightcurve import LightCurve
from .possion_gen import source_poisson_generator, background_poisson_generator


class Detector(object):
    def __init__(self, location, pointing, effective_area, name):
        """
        A detector at a certain location in space,
        with a certain spatial pointing, 
        and a given effective area.

        :param location: a detector location
        :param pointing: a detector pointing
        :param effective_area: a detector effective area
        :param name: the detector's name
        :returns: 
        :rtype: 

        """

        self._effective_area = effective_area
        self._pointing = pointing
        self._location = location
        self._background_slope = 0.0
        self._background_norm = 500.0

        self._name = name

    @property
    def name(self):
        """
        The detector's name
        :returns: 
        :rtype: 

        """

        return self._name

    @property
    def effective_area(self):
        """
        the detector's effective area
        :returns: 
        :rtype: 

        """

        return self._effective_area

    @property
    def pointing(self):
        """
        the pointing of the detector

        :returns: 
        :rtype: 

        """

        return self._pointing

    @property
    def location(self):
        """
        detector's location 

        :returns: 
        :rtype: 

        """

        return self._location

    def light_travel_time(self, grb):
        """
        The light travel time to the GRB

        :param grb: A GRB
        :returns: 
        :rtype: 

        """

        # get the 3D seperation

        distance = self._location.coord.separation_3d(grb.location.coord)

        ltt = distance / constants.c

        return ltt.to("second")

    def angular_separation(self, grb):

        return self._location.coord.separation(grb.location.coord)

    def _check_earth_blockage(self, grb):

        earth_radius = 6371.0 * u.km

        grb_xyz = grb.location.get_cartesian_coord("gcrs").xyz / np.linalg.norm(
            grb.location.get_cartesian_coord("gcrs").xyz
        )

        horizon_angle = (
            0.5 * np.pi
            - np.arccos(earth_radius / self._location.altitude).to(u.rad).value
        )

        sc_pos = self._location.get_cartesian_coord("gcrs").xyz / np.linalg.norm(
            self._location.get_cartesian_coord("gcrs").xyz
        )

        angle = np.arccos(grb_xyz.value.dot(-sc_pos.value))

        if angle < horizon_angle:

            return True

        else:

            return False

    def build_light_curve(
        self, grb, T0, tstart, tstop, earth_blockage=False, seed=1234
    ):
        """
        Build the light curve observed from the GRB

        :param grb: 
        :param T0: 
        :param tstart: 
        :param tstop: 
        :returns: 
        :rtype: 

        """

        # first get the seperation angle
        # from the detector pointing and
        # the grb

        seperation_angle = self._pointing.get_seperation_angle(grb)

        self._effective_area.set_seperation_angle(seperation_angle)

        earth_occulted = False

        if earth_blockage:

            earth_occulted = self._check_earth_blockage(grb)

        # now get the grb's pulse parameters

        K, t_rise, t_decay, t_start = grb.pulse_parameters

        try:

            len(K)

            assert len(K) == len(t_rise)
            assert len(K) == len(t_decay)
            assert len(K) == len(t_start)

            if earth_occulted:
                K = np.array([0.0 for x in K])

                print(f"{self._name} is earth occulted")

            is_multi_pulse = True

        except:

            if earth_occulted:

                print(f"{self._name} is earth occulted")

                K = 0.0

            is_multi_pulse = False

        if not is_multi_pulse:

            # scale the GRB by the effective area

            observed_intensity = K * self._effective_area.effective_area

            # compute the arrival times

            source_arrival_times = source_poisson_generator(
                tstart, tstop, observed_intensity, T0, t_rise, t_decay, seed
            )

        else:

            K = np.array(K)
            t_rise = np.array(t_rise)
            t_decay = np.array(t_decay)
            t_start = np.array(t_start)

            observed_intensity = K * self._effective_area.effective_area

            # scale all the pulses

            t_start += T0

            source_arrival_times = source_poisson_generator(
                tstart, tstop, observed_intensity, t_start, t_rise, t_decay, seed
            )

        bkg_arrival_times = background_poisson_generator(
            tstart, tstop, self._background_slope, self._background_norm, seed + 1
        )

        # return a light curve

        return LightCurve(source_arrival_times, bkg_arrival_times)
