#!/usr/bin/env python

from pylab import *
from scipy.signal import lfilter, butter

# Constants
sr = 1000.0;
T = 1/sr;
r = int(sr/100);
noise_max = 0.1*T;  # This is ||e_k||inf

# Define a velocity curve
vel = array([0.]*(15*r) + [1.]*(4*r) + [2.]*(25*r) + [0.]*(5*r)
            + [-1.]*(3*r) + [-1.]*(20*r))
time = arange(len(vel))/float(sr);

# Integrate it to get position
pos = lfilter([1], [1,-1], vel)*T;

# Add some noise
pos = pos + rand(len(pos))*noise_max

# Finite difference
fdvel = lfilter([1,-1],[1],pos)/T

# Butterworth 100 Hz
[B,A] = butter(1, 0.1)
bwvel = lfilter(B,A,fdvel)

# FD skip 3
dist = 3
fd3vel = lfilter(array([1]+[0]*(dist-1)+[-1])/float(dist),[1],pos)/T

# Least squared 15 (from Freedom6S API)
def leastsquared(n=15):
    dTemp = (n - 1) / 2.0;
    dTemp2 = 12.0/ (n * ((n*n) - 1));
    return (dTemp-arange(n))*dTemp2

lsvel = lfilter(leastsquared(15), 1, pos)/T

# First-Order Adaptive Windowing (FOAW)
def foaw(pos, n=16, best=False):
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

def median_filter(pos, n=5):
    result = zeros(len(pos))
    for k in range(1,len(pos)):
        result[k] = median(pos[k-min(n,k):k])
    return result

endfitfoawvel = foaw(pos, n=16, best=False)
bestfitfoawvel = foaw(pos, n=16, best=True)
mpos = median_filter(pos, n=3)
endfitfoawvelm = foaw(mpos, n=16, best=False)
bestfitfoawvelm = foaw(mpos, n=16, best=True)

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
    show()

curves = [fdvel, fd3vel, bwvel, lsvel]
titles = ['Simple Finite Difference',
          'Finite difference 3',
          'Butterworth %d Hz'%(sr*0.1),
          'Least Squared']

plotcurves(curves, titles, vel_yrange = [-1.5, 2.5], dif_yrange = [-0.3, 0.3])

curves = [endfitfoawvel,bestfitfoawvel,endfitfoawvelm,bestfitfoawvelm]
titles = ['end-fit-FOAW','best-fit-FOAW','end-fit-FOAW w/ median','best-fit-FOAW w/ median']

plotcurves(curves, titles, vel_yrange = [-1.5, 2.5], dif_yrange = [-0.3, 0.3])

