import numpy as np

__all__ = ('reference_angle', 'reference_angle_deg',
           'wrapped_angle', 'wrapped_angle_deg')


def reference_angle(a):
    """Convert an angle to a reference angle between -pi and pi."""
    a = np.mod(a, 2 * np.pi)
    return np.where(a <= np.pi, a, a - 2 * np.pi)


def reference_angle_deg(a):
    """Convert an angle to a reference angle between -180 and 180 degrees."""
    a = np.mod(a, 360)
    return np.where(a <= 180, a, a - 360)


def wrapped_angle(a):
    """Convert an angle to a reference angle between 0 and 2*pi."""
    return np.mod(a, 2 * np.pi)


def wrapped_angle_deg(a):
    """Convert an angle to a reference angle between 0 and 2*pi."""
    return np.mod(a, 360)
