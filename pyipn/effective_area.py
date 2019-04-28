import numpy as np

class EffectiveArea(object):

    def __init__(self, total_area):
        """FIXME! briefly describe function

        :param total_area: 
        :returns: 
        :rtype: 

        """

        self._total_area = total_area
        self._seperation_angle = None
        
        
    @property
    def total_area(self):
        return self._total_area


    @property
    def effective_area(self):
        """FIXME! briefly describe function

        :returns: 
        :rtype: 

        """

        if self._seperation_angle is not None:

            return self.effective_area_at(self._seperation_angle)

        else:

            return self._total_area
        
        
    def set_seperation_angle(self, angle):
        """FIXME! briefly describe function

        :param angle: 
        :returns: 
        :rtype: 

        """

        self._seperation_angle = angle

        
    def effective_area_at(self, angle):
        """FIXME! briefly describe function

        :param angle: 
        :returns: 
        :rtype: 

        """

        if angle >= 0.5 * np.pi:

            return 0.

        else:

            return self._total_area * np.cos(angle)


    
