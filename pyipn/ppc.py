import numpy as np
import numba as nb





@nb.njit()
def ppc_generator(time,exposure ,omega1, omega2, beta1, beta2, bw, scale, amplitude, dt, N):

    out = np.empty((N, len(time)))
    
    
        for n in range(N):

        rate = amplitude[n] * RFF_multiscale(
            time - dt[n],
            omega1[:, n],
            omega2[:, n],
            beta1[:, n],
            beta2[:, n],
            bw[n],
            scale[0, n],
            scale[1, n],
        )


        for i in range(len(time)):
        
           out[n,i] =  np.random.poisson(rate[i] * exposure[i])


    return out


