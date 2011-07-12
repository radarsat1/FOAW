
from pylab import *
from scipy.signal import *
from velocity import *

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

    levantvel = levant(pos, sr, C=max(abs(vel[1:]-vel[:-1]))/T)

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

    show()
