
First-Order Adaptive Windowing velocity estimation in C
=======================================================

Perform the FOAW velocity estimation routine.

This algorithm is described here:

Janabi-Sharifi, F.; Hayward, V.; Chen, C.-S.J., "Discrete-time
adaptive windowing for velocity estimation," Control Systems
Technology, IEEE Transactions on , vol.8, no.6, pp.1003-1009, Nov
2000

http://www.cim.mcgill.ca/~haptic/pub/FS-VH-CSC-TCST-00.pdf 

This implementation (C)2008 Stephen Sinclair, IDMIL, McGill
University.  This work is covered by the GPL-compatible version of
the BSD license, please see the following URL for more information:

http://www.opensource.org/licenses/bsd-license.html

The exact license is listed in the file COPYING, which you should
have received with this source code.

Usage
=====

See the routine `do_foaw()` in the file `foaw.c` for an example of how
to call the velocity estimation routine on a per-step basis.

The algorithm is implemented in,

    float do_foaw_sample(float *posbuf, int size, int *k,
                         float current_pos, int best)

The `posbuf` array must be initialized to zero.  `k` should be
initialized to zero, and is a variable used by the algorithm to track
the current position within `posbuf`.  `size` is the maximum window
size used for estimation.  A value of 15 will give good results, but
you'll have to experiment to see what is needed for your application.
The `best` boolean switches between "best-fit" and "end-fit"
techniques described in the paper.  In general "best-fit" (`best=1`)
gives the best results but takes an order higher in time complexity.
The value of `current_pos` should be set to the latest position value.
The output of the function is the best estimate of the current
velocity.

This velocity estimation algorithm is appropriate for use whenever
velocity is derived from sampled position.  It provides the velocity
estimate with the widest window that stays within a given noise
margin.  In general it will chose a small window when velocity is high
and a large window when velocity is low, thus giving the best veloctiy
estimate while preserving frequency response.  It's particularly
useful for haptic simulation (as mentioned in the paper), since this
is an application area where bandwidth is of utmost importance.

The noise margin is specified by the `NOISE` macro defined at the top
of `foaw.c`, and should be based the amplitude of noise observed in
your position signal.  The `SR` macro should be set to your sample
rate in order to get an accurate velocity estimate.

The program `deriv` runs the FOAW routine on a stream of incoming
ASCII-encoded floating-point values with parameters given as
command-line arguments.
