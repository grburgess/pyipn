import numpy as np
import numba as nb


@nb.njit(fastmath=True)
def RFF(time, omega1, omega2, beta1, beta2, bw, scale):

    features1 = bw * np.outer(time, omega1)
    features2 = bw * np.outer(time, omega2)

    log_fhat = scale * (
        np.dot(np.cos(features1) + np.cos(features2), beta1)
        + np.dot(np.sin(features1) + np.sin(features2), beta2)
    )

    return np.exp(log_fhat)


@nb.njit(fastmath=True)
def RFF_multiscale(time, omega1, omega2, beta1, beta2, bw, scale1, scale2):

    features1 = bw * np.outer(time, omega1)
    features2 = bw * np.outer(time, omega2)

    log_fhat = np.dot(
        scale1 * np.cos(features1) + scale2 * np.cos(features2), beta1
    ) + np.dot(scale1 * np.sin(features1) + scale2 * np.sin(features2), beta2)

    return np.exp(log_fhat)
