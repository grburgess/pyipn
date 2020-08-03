import numpy as np
from numba import jit, njit
from numba.extending import overload
from numba import types
import numba as nb
from pyipn.numba_array import VectorFloat64


@nb.njit(fastmath=True)
def norris(x, K, t_start, t_rise, t_decay):
    if x > t_start:
        return (
            K
            * np.exp(2.0 * np.sqrt(t_rise / t_decay))
            * np.exp((-t_rise / (x - t_start)) - (x - t_start) / t_decay)
        )
    else:
        return 0.0


@nb.njit(fastmath=True)
def mulit_pulse(x, Ks, t_starts, t_rises, t_decays):

    n_pulses = len(Ks)

    out = 0.0

    for n in range(n_pulses):

        out += norris(x, Ks[n], t_starts[n], t_rises[n], t_decays[n])

    return out


def _pulse(x, K, t_start, t_rise, t_decay):
    pass


@overload(_pulse)
def impl_pulse(x, K, t_start, t_rise, t_decay):

    if isinstance(K, types.Float) or isinstance(K, types.Integer):

        if isinstance(x, types.Float) or isinstance(x, types.Integer):

            def impl(x, K, t_start, t_rise, t_decay):
                return norris(x, K, t_start, t_rise, t_decay)

        elif isinstance(x, types.Array):

            def impl(x, K, t_start, t_rise, t_decay):
                N = len(x)
                out = np.empty(N)

                for i in range(N):

                    out[i] = norris(x[i], K, t_start, t_rise, t_decay)

                return out

    elif isinstance(K, types.Array):

        if isinstance(x, types.Float) or isinstance(x, types.Integer):

            def impl(x, K, t_start, t_rise, t_decay):
                return mulit_pulse(x, K, t_start, t_rise, t_decay)

        elif isinstance(x, types.Array):

            def impl(x, K, t_start, t_rise, t_decay):
                N = len(x)
                out = np.empty(N)

                for i in range(N):

                    out[i] = mulit_pulse(x[i], K, t_start, t_rise, t_decay)

                return out

    return impl


@nb.njit()
def pulse(x, K, t_start, t_rise, t_decay):
    return _pulse(x, K, t_start, t_rise, t_decay)


@nb.njit(fastmath=True)
def source_poisson_generator(tstart, tstop, K, p_start, t_rise, t_decay):
    """
    Non-homogeneous poisson process generator
    for a given max rate and time range, this function
    generates time tags sampled from the energy integrated
    lightcurve.
    """

    num_time_steps = 1000

    time_grid = np.linspace(tstart, tstop + 1.0, num_time_steps)

    # tmp = np.zeros(num_time_steps)

    # for i in range(num_time_steps):
    #     tmp[i] = norris(time_grid[i], K, p_start, t_rise, t_decay)

    tmp = pulse(time_grid, K, p_start, t_rise, t_decay)

    fmax = tmp.max()  # zeros if p_start > tstop!

    time = tstart

    arrival_times = VectorFloat64(0)

    # if tstart >= p_start:

    arrival_times.append(tstart)

    while time < tstop:

        time = time - (1.0 / fmax) * np.log(np.random.rand())
        test = np.random.rand()

        p_test = pulse(time, K, p_start, t_rise, t_decay) / fmax

        if test <= p_test:
            arrival_times.append(time)

    return arrival_times.arr


# @njit(fastmath=True)
# def source_poisson_generator_multi(tstart, tstop, K, p_start, t_rise, t_decay):
#     """
#     Non-homogeneous poisson process generator
#     for a given max rate and time range, this function
#     generates time tags sampled from the energy integrated
#     lightcurve.
#     """

#     num_time_steps = 1000

#     time_grid = np.linspace(tstart, tstop + 1.0, num_time_steps)

#     tmp = np.zeros(num_time_steps)

#     for i in range(num_time_steps):
#         tmp[i] = mulit_pulse(time_grid[i], K, p_start, t_rise, t_decay)

#     fmax = tmp.max()  # zeros if p_start > tstop!

#     time = tstart

#     arrival_times = VectorFloat64(0)

#     if tstart >= min(p_start):

#         arrival_times.append(tstart)

#     while time < tstop:

#         time = time - (1.0 / fmax) * np.log(np.random.rand())
#         test = np.random.rand()

#         p_test = mulit_pulse(time, K, p_start, t_rise, t_decay) / fmax

#         if test <= p_test:
#             arrival_times.append(time)

#     return arrival_times.arr


@nb.njit
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
