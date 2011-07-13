#!/usr/bin/env python

__all__ = ['leastsquared', 'foaw', 'levant', 'median_filter']

from pylab import *
from scipy.signal import lfilter, butter
import os, subprocess, threading

cimport numpy as np
cimport cython

DTYPE = double
ctypedef np.double_t DTYPE_t
cdef inline DTYPE_t dmin(DTYPE_t x, DTYPE_t y): return x if x <= y else y
cdef inline int imin(int x, int y): return x if x <= y else y

# Least squared 15 (from Freedom6S API)
def leastsquared(n=15):
    dTemp = (n - 1) / 2.0;
    dTemp2 = 12.0/ (n * ((n*n) - 1));
    return (dTemp-arange(n))*dTemp2

# First-Order Adaptive Windowing (FOAW)
@cython.boundscheck(False)
@cython.wraparound(False)
def foaw(np.ndarray[DTYPE_t] pos not None, double sr, double noise_max,
                 int n=16, int best=False):
    cdef double T = 1/sr
    cdef np.ndarray[DTYPE_t] result = zeros(len(pos), dtype=DTYPE)
    cdef int k, i, j, outside
    cdef DTYPE_t b, ykj, velocity
    for k in range(len(pos)):
        velocity = 0
        for i in range(1,imin(n,k)):
            # Calculate slope over interval
            if (best):
                # least squared method (best-fit-FOAW)
                b = ( ( i*sum([pos[k-j] for j in range(i+1)])
                      - 2*sum([pos[k-j]*j for j in range(i+1)]) )
                    / (T*i*(i+1)*(i+2)/6) )
            else:
                # direct method (end-fit-FOAW)
                b = (pos[k]-pos[k-i]) / (i*T)

            # Check the linear estimate of each middle point
            outside = False
            for j in range(1,i):
                ykj = pos[k]-(b*j*T)

                # Compare to the measured value within the noise margin
                # If it's outside noise margin, return last estimate
                if ykj < (pos[k-j]-noise_max) or ykj > (pos[k-j]+noise_max):
                    outside = True
                    break
            if outside: break
            velocity = b

        result[k] = velocity

    return result

# No need to call out to compiled C version in Cython, this one is
# plenty fast.
fast_foaw = foaw

@cython.boundscheck(False)
@cython.wraparound(False)
def median_filter(np.ndarray[DTYPE_t] pos, int n=5):
    cdef np.ndarray[DTYPE_t] result = zeros(len(pos), dtype=DTYPE)
    cdef int k
    for k in range(1,len(pos)):
        result[k] = median(pos[k-imin(n,k):k])
    return result

# Levant's differentiator, from Levant A. (1998). "Robust exact
# differentiation via sliding mode technique." Automatica, 34(3),
# 379-384.  Suggested for use with force-feedback devices in Chawda et
# al., "Application of Levant's Differentiator for Velocity Estimation
# and Increased Z-Width in Haptic Interfaces", WHC 2011.

# Note that it's not very well-suited to the test data in this file
# because it is sensitive to an estimate of maximum acceleration,
# which in the case of this highly discontinuous velocity is very
# large.  On sinusoidal test data it fairs much better, and gets
# better as sampling rate increases (as opposed to the other
# techniques here).

# Moreover, the papers suggest that Lambda and alpha variables can be
# better-tuned.

# Lipschitz's constant 'C' = maximum absolute acceleration, must be
# provided.

@cython.boundscheck(False)
@cython.wraparound(False)
def levant(np.ndarray[DTYPE_t] pos not None, DTYPE_t sr,
           DTYPE_t C, _alpha=None, _Lambda=None, int n=1):
    cdef DTYPE_t T = 1/sr
    cdef np.ndarray[DTYPE_t] result = zeros(len(pos), dtype=DTYPE)
    cdef DTYPE_t alpha, Lambda
    # Coefficients derived from C
    alpha = 1.1 * C if _alpha==None else _alpha
    Lambda = sqrt(C) if _Lambda==None else _Lambda
    cdef DTYPE_t x = 0, u1 = 0, e, u
    cdef int k, i
    for k in xrange(len(pos)):
        for i in xrange(n):
            e = x - pos[k]
            u1 = u1 - alpha * sign(e) * T/n
            u = u1 - Lambda * sqrt(abs(e)) * sign(e)
            x = x + u * T/n
        result[k] = u
    return result

# Plotting, velocity curves and derivatives
def plotcurves(curves, titles, vel_yrange=None, dif_yrange=None):
    for n, v in enumerate(curves):
        acc = v-vel
        subplot(len(curves),2,n*2+1)
        plot(time, v)
        if (vel_yrange!=None):
            axis([time[0],time[-1],vel_yrange[0],vel_yrange[1]])
        title(titles[n]+': velocity')
        subplot(len(curves),2,n*2+2)
        plot(time, acc)
        if (dif_yrange!=None):
            axis([time[0],time[-1],dif_yrange[0],dif_yrange[1]])
        title(titles[n]+': ideal difference')
