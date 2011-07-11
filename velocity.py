#!/usr/bin/env python

__all__ = ['leastsquared', 'foaw', 'fast_foaw', 'median_filter']

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
# al., "Application of Levant’s Differentiator for Velocity Estimation
# and Increased Z-Width in Haptic Interfaces", WHC 2011.

# Note that it's not very well-suited to the test data in this file
# because it is sensitive to an estimate of maximum acceleration,
# which in the case of this highly discontinuous velocity is very
# large.  On sinusoidal test data it fairs much better, and gets
# better as sampling rate increases (as opposed to the other
# techniques here).

# Moreover, the papers suggest that Lambda and alpha variables can be
# better-tuned.

def levant(pos, sr, n=1, C=None, alpha=None, Lambda=None):
    T = 1/sr
    result = zeros(len(pos))
    # Lipschitz's constant 'C' = maximum absolute acceleration
    if C == None:
        C = max(abs(vel[1:]-vel[:-1]))/T
    # Coefficients derived from C
    if alpha == None:
        alpha = 1.1 * C
    if Lambda == None:
        Lambda = sqrt(C)
    x = 0
    u1 = 0
    for k in range(len(pos)):
        for i in range(n):
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

if __name__=="__main__":
    # Constants
    sr = 1000.0;
    T = 1/sr;
    r = int(sr/100);
    noise_max = 1e-05;  # This is ||e_k||inf

    # Define a velocity curve
    vel = array([0.]*(15*r) + [1.]*(4*r) + [2.]*(25*r) + [0.]*(5*r)
                + [-1.]*(3*r) + [-1.]*(20*r))
    time = arange(len(vel))/float(sr);

    # Another interesting test signal
    # vel = (((0.5+sin(time*50)*pow(2,-time*1))
    #         + (0.2+sin(time*500)*0.2*pow(2,-time*1)))
    #        *concatenate((ones(len(time)/2),
    #                      zeros(len(time)/2))))

    # Integrate it to get position
    pos = lfilter([1], [1,-1], vel)*T

    # Add some noise
    pos = pos + rand(len(pos))*noise_max

    # Finite difference
    fdvel = lfilter([1,-1],[1],pos)/T

    # Butterworth 300 Hz
    [B,A] = butter(2, 300/(sr/2))
    bwvel = lfilter(B,A,fdvel)

    # FD skip 3
    dist = 3
    fd3vel = lfilter(array([1]+[0]*(dist-1)+[-1])/float(dist),[1],pos)/T

    lsvel = lfilter(leastsquared(15), 1, pos)/T

    levantvel = levant(pos, sr)

    endfitfoawvel = foaw(pos, sr, noise_max, n=16, best=False)
    bestfitfoawvel = foaw(pos, sr, noise_max, n=16, best=True)
    mpos = median_filter(pos, n=3)
    endfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=False)
    bestfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=True)

    curves = [fdvel, fd3vel, bwvel, lsvel, levantvel]
    titles = ['Simple Finite Difference',
              'Finite difference 3',
              'Butterworth 300 Hz',
              'Least Squared',
              "Levant's Differentator"]

    figure(1)
    clf()
    plotcurves(curves, titles, vel_yrange = [-1.5, 2.5],
               dif_yrange = [-0.3, 0.3])

    curves = [endfitfoawvel,bestfitfoawvel,endfitfoawvelm,bestfitfoawvelm]
    titles = ['end-fit-FOAW','best-fit-FOAW','end-fit-FOAW w/ median',
              'best-fit-FOAW w/ median']

    figure(2)
    clf()
    plotcurves(curves, titles, vel_yrange = [-1.5, 2.5],
               dif_yrange = [-0.3, 0.3])
    show()

    figure(3)
    clf()
    plot(vel, label='ideal')
    plot(lsvel, label='ls')
    plot(bestfitfoawvel, label='bf-foaw')
    plot(levantvel, label='levant')
    legend()

    r = len(levantvel)/5
    print 'bf-foaw error (%d Hz) ='%sr, sqrt(sum((bestfitfoawvel[r:] - vel[r:])*(bestfitfoawvel[r:] - vel[r:])))
    print 'ef-foaw error (%d Hz) ='%sr, sqrt(sum((endfitfoawvel[r:] - vel[r:])*(endfitfoawvel[r:] - vel[r:])))
    print 'bw2-300 error (%d Hz) ='%sr, sqrt(sum((bwvel[r:] - vel[r:])*(bwvel[r:] - vel[r:])))
    print 'levant error (%d Hz) ='%sr, sqrt(sum((levantvel[r:] - vel[r:])*(levantvel[r:] - vel[r:])))
    print 'fd error (%d Hz) ='%sr, sqrt(sum((fdvel[r:] - vel[r:])*(fdvel[r:] - vel[r:])))
