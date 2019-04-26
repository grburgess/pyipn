class Universe(object):
    def __init__(self, grb):
        """FIXME! briefly describe function

        :param grb: 
        :returns: 
        :rtype: 

        """
        
        self._detectors = collections.OrderedDict()
        self._grb = grb

        self._time_differences = None
        self._T0 = None

    def register_detector(self, detector):
        """FIXME! briefly describe function

        :param detector: 
        :returns: 
        :rtype: 

        """

        self._detectors[detector.name] = detector

    @property
    def detectors(self):
        return self._detectors

    def explode_grb(self, verbose=True):
        """FIXME! briefly describe function

        :param verbose: 
        :returns: 
        :rtype: 

        """

        self._compute_time_differences()

        self._create_light_curves()

    def _compute_time_differences(self):
        """FIXME! briefly describe function

        :returns: 
        :rtype: 

        """

        # compute which detector sees the GRB first
        ltt = []
        for name, detector in self._detectors.items():

            ltt.append(detector.light_travel_time(self._grb))

        # rank the distances in ascending order
        self._distance_rank = np.argsort(distances)
        unsort = self._distance_rank.argsort()

        # for now compute considering all detectors are static
        # the TOA difference of each detector
        ltt = np.array(ltt)[self._distance_rank]

        self._time_differences = [0.0]
        self._T0 = [0.0]
        T0 = 0.0
        for i in range(len(ltt) - 1):

            dt = ltt[i + 1] - ltt[i]
            assert (
                dt > 0
            ), "The time diferences should be positive if the ranking worked!"

            T0 += dt
            self._T0.append(T0)
            self._time_differences.append(dt)

        self._T0 = self._T0[unsort]
        self._time_differences[unsort]

    def _create_light_curves(self):
        """FIXME! briefly describe function

        :returns: 
        :rtype: 

        """

        for t0, (name, detector) in zip(self._T0, self._detectors.items()):

            detector.build_light_curve(self._grb, t0)
