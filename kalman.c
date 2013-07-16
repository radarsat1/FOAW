
// Translation of kalman.py

#include <stdio.h>
#include <string.h>
#include "matvec.h"
#include "kalman.h"

void init_kalman3_pa(kalman3_pa_t *k, KALMAN_FLOAT sampleRate,
                     KALMAN_FLOAT q, KALMAN_FLOAT pR, KALMAN_FLOAT aR)
{
    memset(k, 0, sizeof(kalman3_pa_t));

    KALMAN_FLOAT T = k->T = 1.0f/sampleRate;

    kalman3_pa_set_covariance(k, sampleRate, q, pR, aR);

    // Initial error covariance
    IDENTIFY_MATRIX_3X3(k->P);

    // Motion matrix
    IDENTIFY_MATRIX_3X3(k->A);
    k->A[0][1] = T;
    k->A[0][2] = 0.5*T*T;
    k->A[1][2] = T;

    // Initial state
    k->x_hat[0] = k->x_hat[1] = k->x_hat[2] = 0;
    k->x_hat_est[0] = k->x_hat_est[1] = k->x_hat_est[2] = 0;

    // Measurement matrix, measured = H*state + noise
    k->H[0][0] = 1.0f;
    k->H[0][1] = k->H[0][2] = 0;
    k->H[1][2] = 1.0f;
    k->H[1][0] = k->H[1][1] = 0;
}

void kalman3_pa_set_covariance(kalman3_pa_t *k, KALMAN_FLOAT sampleRate,
                               KALMAN_FLOAT q, KALMAN_FLOAT pR,
                               KALMAN_FLOAT aR)
{
    KALMAN_FLOAT T = k->T = 1.0f/sampleRate;

    // Process covariance
    k->Q[0][0]                           = T*T*T*T*T / 4 * q;
    k->Q[0][1] = k->Q[1][0]              = T*T*T*T / 2 * q;
    k->Q[0][2] = k->Q[1][1] = k->Q[2][0] = T*T*T * q;
    k->Q[0][2] /= 2;
    k->Q[1][2] = k->Q[2][1]              = T*T * q;
    k->Q[2][2]                           = T * q;

    // Measurement covariance
    k->R[0][0] = pR;
    k->R[0][1] = 0;
    k->R[1][0] = 0;
    k->R[1][1] = aR;
}

void step_kalman3_pa(kalman3_pa_t *k, KALMAN_FLOAT pos, KALMAN_FLOAT acc)
{
    KALMAN_FLOAT obs[2] = { pos, acc };
    KALMAN_FLOAT tmp2[3][3], tmp3[3][3];
    KALMAN_FLOAT tmpv[3], det;
    KALMAN_FLOAT I[3][3]={{1,0,0},{0,1,0},{0,0,1}};

    // Make prediction
    MAT_DOT_VEC_3X3(k->x_hat_est, k->A, k->x_hat);  // result 1x3

    MATRIX_PRODUCT_3X3_3X3T(tmp2, k->P, k->A);
    MATRIX_PRODUCT_3X3(tmp3, k->A, tmp2);
    MATRIX_ADD_3X3(k->P_est, tmp3, k->Q);

    // Update estimate
    MAT_DOT_VEC_3X2(tmpv, k->H, k->x_hat_est);  // result 

    VEC_DIFF_2(k->error_x, obs, tmpv);

    MATRIX_PRODUCT_3X3_3X2T(tmp2, k->P, k->H);   // result 2x3
    MATRIX_PRODUCT_3X2_2X3(tmp3, k->H, tmp2);   // result 2x2
    MATRIX_ADD_2X2(k->error_P, tmp3, k->R); // result 2x2

    MATRIX_PRODUCT_3X3_3X2T(tmp3, k->P_est, k->H); // result 2x3
    MAT_SOLVE_AX_EQ_B_2X3_2X2(k->K, tmp3, k->error_P); // result 2x3

    MAT_DOT_VEC_2X3_FULL(tmpv, k->K, k->error_x);  // result 1x3
    VEC_SUM(k->x_hat, k->x_hat_est, tmpv);         // result 1x3

    MATRIX_PRODUCT_2X3_3X2(tmp2, k->K, k->H);
    MATRIX_SUB_3X3(tmp3, I, tmp2);
    MATRIX_PRODUCT_3X3(k->P, tmp3, k->P_est);
}

void run_kalman3_pa(KALMAN_FLOAT sampleRate,
                    KALMAN_FLOAT q, KALMAN_FLOAT pR, KALMAN_FLOAT aR,
                    KALMAN_FLOAT *pos, KALMAN_FLOAT *acc,
                    KALMAN_FLOAT *result, int N)
{
    kalman3_pa_t k;
    init_kalman3_pa(&k, sampleRate, q, pR, aR);

    while (N-- > 0) {
        step_kalman3_pa(&k, *(pos++), *(acc++));
        *(result++) = k.x_hat[0];
        *(result++) = k.x_hat[1];
        *(result++) = k.x_hat[2];
    }
}


void init_kalman2_p(kalman2_p_t *k, KALMAN_FLOAT sampleRate,
                    KALMAN_FLOAT q, KALMAN_FLOAT pR)
{
    memset(k, 0, sizeof(kalman2_p_t));

    KALMAN_FLOAT T = k->T = 1.0f/sampleRate;

    // Process covariance
    k->Q[0][0] = T*T*T / 3 * q;
    k->Q[0][1] = k->Q[1][0] = T*T / 2 * q;
    k->Q[1][1] = T * q;

    // Measurement covariance
    k->R = pR;

    // Initial error covariance
    IDENTIFY_MATRIX_2X2(k->P);

    // Motion matrix
    IDENTIFY_MATRIX_2X2(k->A);
    k->A[0][0] = 1;
    k->A[0][1] = T;
    k->A[1][1] = 1;

    // Initial state
    k->x_hat[0] = k->x_hat[1] = 0;
    k->x_hat_est[0] = k->x_hat_est[1] = 0;

    // Measurement matrix, measured = H*state + noise
    k->H[0] = 1.0f;
    k->H[1] = 0;
}

void step_kalman2_p(kalman2_p_t *k, KALMAN_FLOAT pos)
{
    KALMAN_FLOAT obs = pos;
    KALMAN_FLOAT tmp2[2][2], tmp3[2][2];
    KALMAN_FLOAT tmp, tmpv[2], det;
    KALMAN_FLOAT I[2][2]={{1,0},{0,1}};

    // Make prediction
    MAT_DOT_VEC_2X2(k->x_hat_est, k->A, k->x_hat);  // result 2x1

    MATRIX_PRODUCT_2X2_2X2T(tmp2, k->P, k->A);
    MATRIX_PRODUCT_2X2(tmp3, k->A, tmp2);
    MATRIX_ADD_2X2(k->P_est, tmp3, k->Q);

    // Update estimate
    DOT_VEC_1X2_2X1(tmp, k->H, k->x_hat_est);  // result 1x1

    k->error_x = obs - tmp;

    MAT_DOT_VEC_2X2(tmpv, k->P, k->H);   // result 2x1
    DOT_VEC_1X2_2X1(tmp, k->H, tmpv);   // result 1
    k->error_P = tmp + k->R;

    MATRIX_PRODUCT_2X2_2X1T(tmpv, k->P_est, k->H); // result 2x1
    MAT_SOLVE_AX_EQ_B_1X2_1X1(k->K, tmpv, k->error_P); // result 2x1

    VEC_SCALE_2(tmpv, k->error_x, k->K);  // result 2x1
    VEC_SUM(k->x_hat, k->x_hat_est, tmpv);         // result 2x1

    MATRIX_PRODUCT_2X1_1X2(tmp2, k->K, k->H); // result 2x2
    MATRIX_SUB_2X2(tmp3, I, tmp2);
    MATRIX_PRODUCT_2X2(k->P, tmp3, k->P_est);
}

void run_kalman2_p(KALMAN_FLOAT sampleRate,
                    KALMAN_FLOAT q, KALMAN_FLOAT pR,
                    KALMAN_FLOAT *pos,
                    KALMAN_FLOAT *result, int N)
{
    kalman2_p_t k;
    init_kalman2_p(&k, sampleRate, q, pR);

    while (N-- > 0) {
        step_kalman2_p(&k, *(pos++));
        *(result++) = k.x_hat[0];
        *(result++) = k.x_hat[1];
    }
}
