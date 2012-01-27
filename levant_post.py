#!/usr/bin/env python

from pylab import *
from scipy.signal import lfilter, butter
import velocity

# Test some Levant post-processing filters

sr = 10000.0
t = arange(int(sr))/sr
f0 = 10
p = cos(t*2*pi*f0)*exp(-t*0.5)
v = lfilter([1,-1],1,p)*sr[2:]
noise = uniform(-0.001,0.001,len(t))

C = max(lfilter([1,-1],1,v)*sr[2:])
n = 70

lv = velocity.levant(C=C,sr=sr,pos=p+noise)
v2 = velocity.avgfilt(lv, n=n)
v3 = velocity.maxmin(lv, n=n)
v4 = lfilter(*butter(2,0.02,'low'),x=lv)

clf()
subplot(2,1,1)
plot(t,v)
plot(t,lv)

subplot(2,1,2)
plot(t[4900:5100],v[4900:5100])
plot(t[4900:5100],lv[4900:5100])
plot(t[4900:5100],v2[4900:5100])
plot(t[4900:5100],v3[4900:5100])
plot(t[4900:5100],v4[4900:5100])
