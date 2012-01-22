
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef double TFLOAT;

/******** First-Order Adaptive Windowing *********/

/*
 * Perform the FOAW velocity estimation routine.
 * This algorithm is described here:
 *
 * Janabi-Sharifi, F.; Hayward, V.; Chen, C.-S.J., "Discrete-time
 * adaptive windowing for velocity estimation," Control Systems
 * Technology, IEEE Transactions on , vol.8, no.6, pp.1003-1009, Nov
 * 2000
 *
 * http://www.cim.mcgill.ca/~haptic/pub/FS-VH-CSC-TCST-00.pdf 
 *
 * This implementation (C)2008 Stephen Sinclair, IDMIL, McGill
 * University.  This work is covered by the GPL-compatible version of
 * the BSD license, please see the following URL for more information:
 *
 * http://www.opensource.org/licenses/bsd-license.html
 *
 * The exact license is listed in the file COPYING, which you should
 * have received with this source code.
 */

static inline TFLOAT do_foaw_sample(TFLOAT *posbuf, int size, int *k,
                                    TFLOAT current_pos, int best,
                                    TFLOAT noise_max, TFLOAT T)
{
    int i, j, l, bad;
    TFLOAT b, ykj;
    TFLOAT velocity = 0;

    /* circular buffer */
    *k = (*k+1)%size;
    posbuf[*k] = current_pos;

    for (i=1; i<size; i++)
    {
        if (best)
        {
            // best-fit-FOAW
            b = 0;
            for (l=0; l<(i+1); l++)
                b +=  i*posbuf[(*k-l+size)%size]
                    - 2*posbuf[(*k-l+size)%size]*l;
            b = b / (T*i*(i+1)*(i+2)/6);
        }
        else
            // end-fit-FOAW
            b = (posbuf[*k]-posbuf[(*k-i+size)%size]) / (i*T);
        bad = 0;
        for (j=1; j<i; j++)
        {
            ykj = posbuf[*k]-(b*j*T);
            if (   (ykj < (posbuf[(*k-j+size)%size]-noise_max))
                || (ykj > (posbuf[(*k-j+size)%size]+noise_max)))
            {
                bad = 1;
                break;
            }
        }
        if (bad) break;
        velocity = b;
    }

    return velocity;
}

void foaw_best_fit(TFLOAT SR, int N, TFLOAT noise,
                   TFLOAT *input, TFLOAT *output,
                   int size)
{
    TFLOAT T = 1.0/SR;
    TFLOAT posbuf[N];
    int k=0;
    memset(posbuf, input[0], N*sizeof(TFLOAT));
    while (size--) {
        *(output++) = do_foaw_sample(posbuf, N, &k, *(input++), 1,
                                     noise, T);
    }
}

void foaw_end_fit(TFLOAT SR, int N, TFLOAT noise,
                  TFLOAT *input, TFLOAT *output,
                  int size)
{
    TFLOAT T = 1.0/SR;
    TFLOAT posbuf[N];
    int k=0;
    int i=0;
    for (i=0; i<N; i++)
        posbuf[i] = input[0];
    while (size--)
        *(output++) = do_foaw_sample(posbuf, N, &k, *(input++), 0, noise, T);
}

/******** Levant's differentiator *********/

/*
 * Levant's differentiator, from Levant A. (1998). "Robust exact
 * differentiation via sliding mode technique." Automatica, 34(3),
 * 379-384.  Suggested for use with force-feedback devices in Chawda et
 * al., "Application of Levant's Differentiator for Velocity Estimation
 * and Increased Z-Width in Haptic Interfaces", WHC 2011.

 * Note that it's not very well-suited to the test data in this file
 * because it is sensitive to an estimate of maximum acceleration,
 * which in the case of this highly discontinuous velocity is very
 * large.  On sinusoidal test data it fairs much better, and gets
 * better as sampling rate increases (as opposed to the other
 * techniques here).

 * Moreover, the papers suggest that Lambda and alpha variables can be
 * better-tuned.

 * Lipschitz's constant 'C' = maximum absolute acceleration, must be
 * provided.
*/

/* Table lookup for a 2nd-order fit to sqrt between 0 and 1 */
static TFLOAT sqrttbl[100][3] = {
    { -5.97546696e+02,   1.40461896e+01,   1.62673996e-02},
    { -7.02157732e+01,   6.22495658e+00,   4.49043263e-02},
    { -3.19746276e+01,   4.77081197e+00,   5.88286283e-02},
    { -1.91974969e+01,   4.02063743e+00,   6.98777678e-02},
    { -1.31389771e+01,   3.54176558e+00,   7.93591255e-02},
    { -9.71290158e+00,   3.20177791e+00,   8.78045862e-02},
    { -7.55515321e+00,   2.94422075e+00,   9.54971804e-02},
    { -6.09320979e+00,   2.74034500e+00,   1.02609728e-01},
    { -5.04885132e+00,   2.57374507e+00,   1.09257226e-01},
    { -4.27221921e+00,   2.43427995e+00,   1.15520891e-01},
    { -3.67616738e+00,   2.31529653e+00,   1.21460651e-01},
    { -3.20691652e+00,   2.21222372e+00,   1.27122246e-01},
    { -2.82966427e+00,   2.12180283e+00,   1.32541544e-01},
    { -2.52099609e+00,   2.04163945e+00,   1.37747291e-01},
    { -2.26463978e+00,   1.96992936e+00,   1.42762952e-01},
    { -2.04897170e+00,   1.90528363e+00,   1.47607977e-01},
    { -1.86549063e+00,   1.84661332e+00,   1.52298706e-01},
    { -1.70785035e+00,   1.79305090e+00,   1.56849024e-01},
    { -1.57122761e+00,   1.74389558e+00,   1.61270847e-01},
    { -1.45189802e+00,   1.69857423e+00,   1.65574492e-01},
    { -1.34694484e+00,   1.65661291e+00,   1.69768967e-01},
    { -1.25405475e+00,   1.61761590e+00,   1.73862188e-01},
    { -1.17137198e+00,   1.58124978e+00,   1.77861157e-01},
    { -1.09739237e+00,   1.54723138e+00,   1.81772104e-01},
    { -1.03088520e+00,   1.51531848e+00,   1.85600603e-01},
    { -9.70834857e-01,   1.48530245e+00,   1.89351660e-01},
    { -9.16396543e-01,   1.45700249e+00,   1.93029796e-01},
    { -8.66862460e-01,   1.43026106e+00,   1.96639107e-01},
    { -8.21635589e-01,   1.40494015e+00,   2.00183320e-01},
    { -7.80209213e-01,   1.38091829e+00,   2.03665837e-01},
    { -7.42150767e-01,   1.35808804e+00,   2.07089773e-01},
    { -7.07088992e-01,   1.33635404e+00,   2.10457989e-01},
    { -6.74703641e-01,   1.31563127e+00,   2.13773122e-01},
    { -6.44717165e-01,   1.29584365e+00,   2.17037606e-01},
    { -6.16887956e-01,   1.27692290e+00,   2.20253695e-01},
    { -5.91004816e-01,   1.25880751e+00,   2.23423479e-01},
    { -5.66882409e-01,   1.24144193e+00,   2.26548904e-01},
    { -5.44357491e-01,   1.22477581e+00,   2.29631781e-01},
    { -5.23285776e-01,   1.20876341e+00,   2.32673804e-01},
    { -5.03539310e-01,   1.19336309e+00,   2.35676554e-01},
    { -4.85004269e-01,   1.17853682e+00,   2.38641516e-01},
    { -4.67579100e-01,   1.16424980e+00,   2.41570081e-01},
    { -4.51172935e-01,   1.15047011e+00,   2.44463557e-01},
    { -4.35704256e-01,   1.13716841e+00,   2.47323178e-01},
    { -4.21099742e-01,   1.12431770e+00,   2.50150104e-01},
    { -4.07293281e-01,   1.11189305e+00,   2.52945432e-01},
    { -3.94225127e-01,   1.09987143e+00,   2.55710198e-01},
    { -3.81841163e-01,   1.08823151e+00,   2.58445385e-01},
    { -3.70092270e-01,   1.07695350e+00,   2.61151921e-01},
    { -3.58933769e-01,   1.06601903e+00,   2.63830689e-01},
    { -3.48324944e-01,   1.05541102e+00,   2.66482526e-01},
    { -3.38228619e-01,   1.04511352e+00,   2.69108229e-01},
    { -3.28610790e-01,   1.03511168e+00,   2.71708556e-01},
    { -3.19440298e-01,   1.02539161e+00,   2.74284228e-01},
    { -3.10688548e-01,   1.01594034e+00,   2.76835934e-01},
    { -3.02329251e-01,   1.00674569e+00,   2.79364330e-01},
    { -2.94338202e-01,   9.97796259e-01,   2.81870045e-01},
    { -2.86693084e-01,   9.89081335e-01,   2.84353678e-01},
    { -2.79373289e-01,   9.80590852e-01,   2.86815802e-01},
    { -2.72359760e-01,   9.72315340e-01,   2.89256968e-01},
    { -2.65634852e-01,   9.64245876e-01,   2.91677701e-01},
    { -2.59182206e-01,   9.56374050e-01,   2.94078506e-01},
    { -2.52986636e-01,   9.48691924e-01,   2.96459867e-01},
    { -2.47034029e-01,   9.41191998e-01,   2.98822250e-01},
    { -2.41311251e-01,   9.33867182e-01,   3.01166101e-01},
    { -2.35806068e-01,   9.26710766e-01,   3.03491849e-01},
    { -2.30507071e-01,   9.19716396e-01,   3.05799908e-01},
    { -2.25403608e-01,   9.12878045e-01,   3.08090675e-01},
    { -2.20485725e-01,   9.06189999e-01,   3.10364533e-01},
    { -2.15744109e-01,   8.99646831e-01,   3.12621851e-01},
    { -2.11170041e-01,   8.93243384e-01,   3.14862985e-01},
    { -2.06755348e-01,   8.86974756e-01,   3.17088279e-01},
    { -2.02492362e-01,   8.80836281e-01,   3.19298062e-01},
    { -1.98373883e-01,   8.74823517e-01,   3.21492656e-01},
    { -1.94393146e-01,   8.68932230e-01,   3.23672370e-01},
    { -1.90543785e-01,   8.63158384e-01,   3.25837501e-01},
    { -1.86819811e-01,   8.57498130e-01,   3.27988339e-01},
    { -1.83215578e-01,   8.51947789e-01,   3.30125163e-01},
    { -1.79725765e-01,   8.46503851e-01,   3.32248244e-01},
    { -1.76345350e-01,   8.41162958e-01,   3.34357843e-01},
    { -1.73069592e-01,   8.35921900e-01,   3.36454215e-01},
    { -1.69894009e-01,   8.30777606e-01,   3.38537604e-01},
    { -1.66814366e-01,   8.25727133e-01,   3.40608249e-01},
    { -1.63826651e-01,   8.20767664e-01,   3.42666382e-01},
    { -1.60927070e-01,   8.15896498e-01,   3.44712226e-01},
    { -1.58112024e-01,   8.11111046e-01,   3.46745998e-01},
    { -1.55378103e-01,   8.06408824e-01,   3.48767911e-01},
    { -1.52722072e-01,   8.01787445e-01,   3.50778169e-01},
    { -1.50140858e-01,   7.97244620e-01,   3.52776971e-01},
    { -1.47631544e-01,   7.92778149e-01,   3.54764511e-01},
    { -1.45191357e-01,   7.88385916e-01,   3.56740977e-01},
    { -1.42817660e-01,   7.84065886e-01,   3.58706553e-01},
    { -1.40507944e-01,   7.79816104e-01,   3.60661417e-01},
    { -1.38259821e-01,   7.75634686e-01,   3.62605741e-01},
    { -1.36071014e-01,   7.71519818e-01,   3.64539694e-01},
    { -1.33939356e-01,   7.67469754e-01,   3.66463441e-01},
    { -1.31862780e-01,   7.63482810e-01,   3.68377141e-01},
    { -1.29839313e-01,   7.59557363e-01,   3.70280951e-01},
    { -1.27867073e-01,   7.55691849e-01,   3.72175022e-01},
    { -1.25944263e-01,   7.51884758e-01,   3.74059502e-01},
};

#define sign(x) ((x)==0.0?0.0:((x)<0.0?-1.0:1.0))

static inline TFLOAT tbl_sqrt(TFLOAT x)
{
    if (x>=0.0 && x<1.0) {
        int i = (int)(x*100);
        return x*x*sqrttbl[i][0]+x*sqrttbl[i][1]+sqrttbl[i][2];
    }
    else
        return sqrt(x);
}

static inline void f(TFLOAT alpha, TFLOAT Lambda,
                     TFLOAT p, TFLOAT u1, TFLOAT x,
                     TFLOAT *du1, TFLOAT *dx)
{
    TFLOAT e = x-p;
    *du1 = -alpha * sign(e);
    *dx = u1-Lambda * tbl_sqrt(fabs(e)) * sign(e);
}

void levant(TFLOAT sr, TFLOAT C, int rk,
            TFLOAT *pos, TFLOAT *result, int size)
{
    TFLOAT T = 1.0/sr;

    // Coefficients derived from C;
    TFLOAT alpha = 1.1 * C;
    TFLOAT Lambda = sqrt(C);

    TFLOAT x = 0, u1 = 0, u;
    TFLOAT k1du1, k1dx, k2du1, k2dx, k3du1, k3dx, k4du1, k4dx, tu1, tx;

    int i;

    if (rk==4) {
        for (i=0; i<size; i++) {
            f(alpha,Lambda,pos[i], u1, x, &k1du1, &k1dx);
            f(alpha,Lambda,pos[i], u1+(T/2)*k1du1, x+(T/2)*k1dx, &k2du1, &k2dx);
            f(alpha,Lambda,pos[i], u1+(T/2)*k2du1, x+(T/2)*k2dx, &k3du1, &k3dx);
            f(alpha,Lambda,pos[i], u1+T*k3du1, x+T*k3dx, &k4du1, &k4dx);
            u1 = u1 + (T/6)*(k1du1 + 2*k2du1 + 2*k3du1 + k4du1);
            u = (1.0/6)*(k1dx + 2*k2dx + 2*k3dx + k4dx);
            x = x + u*T;
            result[i] = u;
        }
    }
    else if (rk==2) {
        for (i=0; i<size; i++) {
            f(alpha,Lambda,pos[i],u1,x,&k1du1,&k1dx);
            tu1 = u1 + k1du1*(T/2);
            tx = x + k1dx*(T/2);
            f(alpha,Lambda,pos[i],tu1,tx, &k2du1, &k2dx);
            u1 = u1 + k2du1*T;
            x = x + k2dx*T;
            result[i] = k2dx;
        }
    }
    else if (rk==1) {
        for (i=0; i<size; i++) {
            f(alpha,Lambda,pos[i],u1,x,&k1du1,&k1dx);
            u1 = u1 + k1du1*T;
            x = x + k1dx*T;
            result[i] = k1dx;
        }
    }
    else {
        printf("[levant] Unknown rk==%d\n", rk);
    }
}

/******** Median filter *********/

static int tfloatcomp(const void *a, const void *b)
{
    if (*(const TFLOAT*)a < *(const TFLOAT*)b)
        return -1;
    else
        return *(const TFLOAT*)a > *(const TFLOAT*)b;
}

void median_filter(unsigned int n, TFLOAT *pos, TFLOAT *result, int size)
{
    int i;
    TFLOAT buf[n], med[n];
    memset(buf, 0, sizeof(buf));
    for (i=0; i<size; i++) {
        buf[i%n] = pos[i];
        memcpy(med, buf, sizeof(buf));
        qsort(med, n, sizeof(TFLOAT), tfloatcomp);
        result[i] = med[n/2];
    }
}
