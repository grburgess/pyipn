import numpy as np

class Detector(object):

    def __init__(self, location, pointing, effective_area):


        self._effective_area = effective_area
        self._pointing = pointing
        self._location = location
        

    @property
    def effective_area(self):

        return self._effective_area

    @property
    def pointing(self):

        return self._pointing

    @property
    def location(self):

        return self._pointing
