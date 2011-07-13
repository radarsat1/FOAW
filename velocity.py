#!/usr/bin/env python

__all__ = ['leastsquared', 'foaw', 'fast_foaw', 'median_filter', 'levant']

from pylab import *
from scipy.signal import lfilter, butter
import os, subprocess, threading

# Least squared 15 (from Freedom6S API)
def leastsquared(n=15):
    dTemp = (n - 1) / 2.0;
    dTemp2 = 12.0/ (n * ((n*n) - 1));
    return (dTemp-arange(n))*dTemp2

# First-Order Adaptive Windowing (FOAW)
def foaw(pos, sr, noise_max, n=16, best=False):
    T = 1/sr
    result = zeros(len(pos))
    for k in range(len(pos)):
        velocity = 0
        for i in range(1,min(n,k)):
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

def fast_foaw(pos, sr, noise_max, n=16, best=False):
    """Run a faster version of FOAW by calling to C compiled code."""
    path = '.'.join(__file__.split('.')[:-1]+['py'])
    program = os.path.join(os.path.dirname(os.path.realpath(path)),'deriv')
    cmd = [program, str(sr), str(noise_max), str(n), str(best and 1 or 0)]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, close_fds=True)

    def writer():
        for i in pos:
            print >>p.stdin, i
        p.stdin.close()

    t = threading.Thread(target=writer)
    t.start()
    out = []
    for i in p.stdout:
        out.append(float(i))
    t.join()
    return array(out)

def median_filter(pos, n=5):
    result = zeros(len(pos))
    for k in range(1,len(pos)):
        result[k] = median(pos[k-min(n,k):k])
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

def f(alpha,Lambda,p,u1,x):
    e = x-p
    return array([ -alpha * sign(e),
                    u1-Lambda * sqrt(abs(e)) * sign(e) ])

def levant(pos, sr, C, alpha=None, Lambda=None, rk=1):
    T = 1/sr
    result = zeros(len(pos))
    # Coefficients derived from C
    if alpha == None:
        alpha = 1.1 * C
    if Lambda == None:
        Lambda = sqrt(C)
    x = 0
    u1 = 0
    if rk==4:
        for k in range(len(pos)):
            k1du1, k1dx = f(alpha,Lambda,pos[k], u1, x)
            k2du1, k2dx = f(alpha,Lambda,pos[k], u1+(T/2)*k1du1, x+(T/2)*k1dx)
            k3du1, k3dx = f(alpha,Lambda,pos[k], u1+(T/2)*k2du1, x+(T/2)*k2dx)
            k4du1, k4dx = f(alpha,Lambda,pos[k], u1+T*k3du1, x+T*k3dx)
            u1 = u1 + (T/6)*(k1du1 + 2*k2du1 + 2*k3du1 + k4du1)
            u = (1.0/6)*(k1dx + 2*k2dx + 2*k3dx + k4dx)
            x = x + u*T
            result[k] = u
    elif rk==2:
        for k in range(len(pos)):
            k1du1, k1dx = f(alpha,Lambda,pos[k],u1,x)
            tu1 = u1 + k1du1*(T/2)
            tx = x + k1dx*(T/2)
            k2du1, k2dx = f(alpha,Lambda,pos[k],tu1,tx)
            u1 = u1 + k2du1*T
            x = x + k2dx*T
            result[k] = k2dx
    elif rk==1:
        for k in range(len(pos)):
            k1du1, k1dx = f(alpha,Lambda,pos[k],u1,x)
            u1 = u1 + k1du1*T
            x = x + k1dx*T
            result[k] = k1dx
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
