import numpy as np
from numba import jit, njit
from pyipn.numba_array import VectorFloat64

@njit
def norris(x, K, t_start, t_rise, t_decay):
    if x > t_start:
        return (
            K
            * np.exp(2 * (t_rise / t_decay) ** (1 / 2))
            * np.exp(-t_rise / (x - t_start) - (x - t_start) / t_decay)
        )
    else:
        return 0.0


@njit
def source_poisson_generator(tstart, tstop, K, p_start, t_rise, t_decay):
    """
    Non-homogeneous poisson process generator
    for a given max rate and time range, this function
    generates time tags sampled from the energy integrated
    lightcurve.
    """
    if K == 0.:
        return np.empty(1)
    else:
        num_time_steps = 1000

        time_grid = np.linspace(tstart, tstop + 1.0, num_time_steps)

        tmp = np.zeros(num_time_steps)

        for i in range(num_time_steps):
            tmp[i] = norris(time_grid[i], K, p_start, t_rise, t_decay)

        fmax = tmp.max() #zeros if p_start > tstop!

        time = tstart

        arrival_times = VectorFloat64(0)
        
        if tstart >= p_start:
        
            arrival_times.append(tstart)
            
        while time < tstop:

            time = time - (1.0 / fmax) * np.log(np.random.rand())
            test = np.random.rand()

            p_test = norris(time, K, p_start, t_rise, t_decay) / fmax

            if test <= p_test:
                arrival_times.append(time)

        return arrival_times.arr


@njit
def background_poisson_generator(tstart, tstop, slope, intercept):
    """
    Non-homogeneous poisson process generator
    for a given max rate and time range, this function
    generates time tags sampled from the energy integrated
    lightcurve.
    """

    num_time_steps = 1000

    time_grid = np.linspace(tstart, tstop + 1.0, num_time_steps)

    tmp = intercept + slope * time_grid

    fmax = tmp.max()

    time = tstart
    arrival_times = VectorFloat64(0)
    arrival_times.append(tstart)

    while time < tstop:

        time = time - (1.0 / fmax) * np.log(np.random.rand())
        test = np.random.rand()

        p_test = (intercept + slope * time) / fmax

        if test <= p_test:
            arrival_times.append(time)

    return arrival_times.arr
