
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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

float SR=1000;
#define T (1.0f/SR)
#define SIZE 1000
float NOISE=0;

float vel[SIZE];
float pos[SIZE];
float mpos[SIZE];
float fdvel[SIZE];
float foawvel[SIZE];
float dif[SIZE];

#ifndef min
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

void generate_velocity()
{
    int k;
    for (k=0; k<50; k++)
        vel[k] = 0;
    for (; k<150; k++)
        vel[k] = 1;
    for (; k<450; k++)
        vel[k] = 2;
    for (; k<600; k++)
        vel[k] = 0;
    for (; k<700; k++)
        vel[k] = -1;
    for (; k<1000; k++)
        vel[k] = -2;
}

void integrate_position()
{
    int k;
    pos[0] = 0;
    for (k=1; k<SIZE; k++)
        pos[k] = (vel[k]*T+pos[k-1]);
}

void add_position_noise()
{
    int k;
    for (k=0; k<SIZE; k++)
        pos[k] += rand()/(float)RAND_MAX * NOISE;
}

void finite_difference()
{
    int k;
    fdvel[0] = 0;
    for (k=1; k<SIZE; k++)
        fdvel[k] = (pos[k]-pos[k-1])/T;
}

float do_foaw_sample(float *posbuf, int size, int *k,
                     float current_pos, int best)
{
    int i, j, l, bad;
    float b, ykj;
    float velocity = 0;
    float noise_max = NOISE;

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

void subtract(float *result, float *from, float *what)
{
    int k;
    for (k=0; k<SIZE; k++)
        result[k] = from[k]-what[k];
}

void median_position(int n)
{
    int k, i, j, off=600;
    float *buf;
    buf = malloc(sizeof(float)*2*n);
    for (k=0; k<SIZE; k++)
    {
        int size = min(n,k);
        int a=n, b=n;
        buf[n] = pos[k];
        for (i=1; i<size; i++)
        {
            if (pos[k-i] > buf[b]) {
                buf[b+1] = pos[k-i];
                b++;
            }
            else if (pos[k-i] < buf[a]) {
                buf[a-1] = pos[k-i];
                a--;
            }
            else if (pos[k-i] > buf[b-1]) {
                buf[b+1] = buf[b];
                buf[b] = pos[k-i];
                b++;
            }
            else if (pos[k-i] < buf[a+1]) {
                buf[a-1] = buf[a];
                buf[a] = pos[k-i];
                a--;
            }
        }
        mpos[k] = buf[(b-a)/2+a];
    }
    free(buf);
}


int main(int argc, char *argv[])
{
    int n = 15;
    int best = 0;
    float *posbuf = 0;
    int k, i=0;
    float inp;

    if (argc > 1) {
        SR = (float)atoi(argv[1]);
        if (SR <= 0) {
            fprintf(stderr, "Bad sample rate.\n");
            return 1;
        }
    }
    else {
        printf("deriv <rate> <noise> <size> <bestfit=1,endfit=0>\n");
        exit(0);
    }

    if (argc > 2) {
        NOISE = atof(argv[2]);
        if (NOISE <= 0) {
            fprintf(stderr, "Bad noise margin.\n");
            return 1;
        }
    }
    else
        NOISE = (0.1*T);

    if (argc > 3) {
        n = atoi(argv[3]);
        if (n <= 0) {
            fprintf(stderr, "Bad N.\n");
            return 1;
        }
    }

    if (argc > 4) {
        best = atoi(argv[4]);
        if (best < 0 || best > 1) {
            fprintf(stderr, "Bad 'best'.\n");
            return 1;
        }
    }

    posbuf = malloc(sizeof(float)*n);
    memset(posbuf, 0, sizeof(float)*n);

    scanf("%g", &inp);
    do_foaw_sample(posbuf, n, &i, inp, best);
    while (scanf("%g", &inp)!=EOF) {
        printf("%g\n", do_foaw_sample(posbuf, n, &i, inp, best));
    }

    free(posbuf);
    return 0;
}
