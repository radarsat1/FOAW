
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

    levantvel1 = levant(pos, sr, C=max(abs(vel[1:]-vel[:-1]))/T, rk=1)
    levantvel2 = levant(pos, sr, C=max(abs(vel[1:]-vel[:-1]))/T, rk=2)
    levantvel4 = levant(pos, sr, C=max(abs(vel[1:]-vel[:-1]))/T, rk=4)

    endfitfoawvel = foaw(pos, sr, noise_max, n=16, best=False)
    bestfitfoawvel = foaw(pos, sr, noise_max, n=16, best=True)
    mpos = median_filter(pos, n=3)
    endfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=False)
    bestfitfoawvelm = foaw(mpos, sr, noise_max, n=16, best=True)

    curves = [fdvel, fd3vel, bwvel, lsvel]
    titles = ['Simple Finite Difference',
              'Finite difference 3',
              'Butterworth 300 Hz',
              'Least Squared']

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

    curves = [levantvel1, levantvel2, levantvel4]
    titles = ['Levant RK=1',
              'Levant RK=2',
              'Levant RK=4']

    figure(3)
    clf()
    plotcurves(curves, titles, vel_yrange = [-1.5, 2.5],
               dif_yrange = [-0.3, 0.3])

    figure(4)
    clf()
    plot(vel, label='ideal')
    plot(lsvel, label='ls')
    plot(bestfitfoawvel, label='bf-foaw')
    plot(levantvel1, label='levant1')
    plot(levantvel2, label='levant2')
    plot(levantvel4, label='levant4')
    legend()

    def rms(x):
        return sqrt(sum((x[r:] - vel[r:])*(x[r:] - vel[r:])))

    r = len(levantvel1)/5
    print 'bf-foaw error (%d Hz) ='%sr, rms(bestfitfoawvel)
    print 'ef-foaw error (%d Hz) ='%sr, rms(endfitfoawvel)
    print 'bw2-300 error (%d Hz) ='%sr, rms(bwvel)
    print 'levant1 error (%d Hz) ='%sr, rms(levantvel1)
    print 'levant2 error (%d Hz) ='%sr, rms(levantvel2)
    print 'levant4 error (%d Hz) ='%sr, rms(levantvel4)
    print 'fd error (%d Hz) ='%sr, rms(fdvel)

    show()
