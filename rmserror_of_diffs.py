#!/usr/bin/env python

from pylab import *
from scipy.signal import *
from velocity import *
from collections import defaultdict
from pprint import pprint

def do_rms_for_methods(sr, vel, noise_max):
    # Constants
    T = 1/sr;
    r = int(sr/100);

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

    levantC = C=max(abs(vel[1:]-vel[:-1]))/T
    print 'Lipschitz\'s constant C = %f, sr = %f'%(levantC,sr)
    levantvel1 = levant(pos, sr, C=levantC, rk=1)
    levantvel2 = levant(pos, sr, C=levantC, rk=2)
    levantvel4 = levant(pos, sr, C=levantC, rk=4)

    endfitfoawvel = foaw(pos, sr, noise_max, n=16, best=False)
    bestfitfoawvel = foaw(pos, sr, noise_max, n=16, best=True)
    mpos = median_filter(pos, n=3)
    endfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=False)
    bestfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=True)

    def rms(x):
        return sqrt(sum((x[r:] - vel[r:])*(x[r:] - vel[r:])))

    r = len(levantvel1)/5
    return {'bf-foaw16': (sr, rms(bestfitfoawvel)),
            'ef-foaw16': (sr, rms(endfitfoawvel)),
            'bw2-300': (sr, rms(bwvel)),
            'levant1': (sr, rms(levantvel1)),
            'levant2': (sr, rms(levantvel2)),
            'levant4': (sr, rms(levantvel4)),
            'fd': (sr, rms(fdvel))}

def do_plot(results, name):
    marks = 'o^s+D'*10
    for i, algo in enumerate(sort(results.keys())):
        loglog(*zip(*results[algo]), marker=marks[i], label=algo)
        xlim(1000, 100000)
        xlabel('sampling rate (Hz)')
        ylabel('RMS error')
        legend()
        title(name)

    print
    print name
    print '-'*len(name)
    print
    pprint(dict(results))
    savefig(name.replace(' ','_')+'.png')

if __name__=="__main__":
    rc('legend',fontsize=8)

    rates = [1000, 2000, 4000, 6000, 8000,
             10000, 20000, 40000, 60000, 80000, 100000]

    # Two sinusoids without noise
    def genvel(sr,t):
        time = arange(sr*t)/float(sr);
        return time, (((0.5+sin(time*50)*pow(2,-time*10))
                       + (0.2+sin(time*500)*0.2*pow(2,-time*10))))
    figure(1)
    clf()
    results = defaultdict(lambda: [])
    for sr in rates:
        time, vel = genvel(sr,1)
        res = do_rms_for_methods(float(sr), vel, noise_max=0)
        for r in res:
            results[r].append(res[r])
    do_plot(results, 'two sinusoids without noise')

    # Two sinusoids with noise
    def genvel(sr,t):
        time = arange(sr*t)/float(sr);
        return time, ((0.5+sin(time*50)*pow(2,-time*10))
                      + (0.2+sin(time*500)*0.2*pow(2,-time*10)))
    figure(2)
    clf()
    results = defaultdict(lambda: [])
    for sr in rates:
        time, vel = genvel(sr,1)
        res = do_rms_for_methods(float(sr), vel, noise_max=1e-05)
        for r in res:
            results[r].append(res[r])
    do_plot(results, 'two sinusoids with noise = 1e-05')

    # Two sinusoids with low noise
    def genvel(sr,t):
        time = arange(sr*t)/float(sr);
        return time, ((0.5+sin(time*50)*pow(2,-time*10))
                      + (0.2+sin(time*500)*0.2*pow(2,-time*10)))
    figure(3)
    clf()
    results = defaultdict(lambda: [])
    for sr in rates:
        time, vel = genvel(sr,1)
        res = do_rms_for_methods(float(sr), vel, noise_max=1e-06)
        for r in res:
            results[r].append(res[r])
    do_plot(results, 'two sinusoids with noise = 1e-06')

    # Two sinusoids + noise + discontinuity
    def genvel(sr,t):
        time = arange(sr*t)/float(sr);
        return time, (((0.5+sin(time*50)*pow(2,-time*10))
                       + (0.2+sin(time*500)*0.2*pow(2,-time*10)))
                      *concatenate((ones(len(time)/2),
                                    zeros(len(time)/2))))
    figure(4)
    clf()
    results = defaultdict(lambda: [])
    for sr in rates:
        time, vel = genvel(sr,1)
        res = do_rms_for_methods(float(sr), vel, noise_max=1e-05)
        for r in res:
            results[r].append(res[r])
    do_plot(results, 'two sinusoids + noise=1e-05 + discontinuity')

    show()
