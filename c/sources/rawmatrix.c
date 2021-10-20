/** @file Minimalistic library with basic matrix operations. */
#include <stdlib.h>
#include <math.h>
#include <rawmatrix.h>

#define mabs(x) (x >= 0.0 ? x : -x) 
#define TINY 1.0E-12
#define IsEqual(x, y, limit) ((x >= y ? (x - y) : (y - x)) <= limit ? 1 : 0)


int rawmatrix_check(const rawmatrix *const A)
{
    if (A == NULL) {
        return RAWMATRIX_NULL_HANDLE;
    }
}


int rawmatrix_lu_decomposition(rawmatrix_t *A,
                               unsigned int *idx,
							   double *s,
							   int *multiplier,
							   int *isnan)
{
    int rcode, i, j, k, q;
    double x, xmax, *a, *ai, *aj, *acol, *adiag, *sj;

    rcode = rawmatrix_check(A);
    if (rcode != RAWMATRIX_OK) {
        return rcode;
    }

    if (multiplier != NULL){ *multiplier = 1; }
    if (isnan != NULL){ *isnan = 0; }

    // Selection of the biggest element of each row
    if ((idx != NULL) && (s != NULL)){
        a = A->x;
        sj = s;
        for (i = A->n; i > 0; --i){
            *sj = 0.0;
            for (j = A->n; j > 0; --j){
                if ((x = mabs(*a)) > *sj) {
                    *sj = x;
                }
                ++a;
            }
            if ((*sj == 0.0) && (isnan != NULL) && (!*isnan)){
                *isnan = A->n - i + 1;
            }
            *sj = 1.0 / *sj;
            ++sj;
        }
    }

    acol = aj = A->x;
    aj -= A->n;

    // Col j
    for (j = 0; j < A->n; ++j){
        // Lines 1 <= i < j
        a = acol++;
        q = 0;
        for (i = j; i > 0; --i){
            for (k = q++; k > 0; --k){
                *a -= (*ai-- * *aj);
                aj -= A->n;
            }
            aj = a;
            a += A->n;
            ai = a - i ;
        }
        adiag = aj;
        if (s != NULL){ sj = s; }
        // Lines j <= i <= N
        xmax = 0.0;
        q = A->n - j;
        for (i = q; i > 0; --i){
            for (k = j; k > 0; --k){
                *a -= (*ai-- * *aj);
                aj -= A->n;
            }
            if (idx != NULL){
                if (s != NULL){
                    x = (*sj++ * mabs(*a));
                }
                else {
                    x = mabs(*a);
                }
                if (x > xmax){
                    xmax = x;
                    q = i;
                }
            }
            aj = adiag;
            a += A->n;
            ai = a - 1;
        }
        adiag += A->n;
        q = A->n - j - q;
        if (q){
            ai = adiag - j;
            aj = ai + (A->n * q);
            for (k = A->n; k > 0; --k) { // Exchanging line adiag (j) and amax
                x = *ai;
                *ai++ = *aj;
                *aj++ = x;
            }
            if (multiplier != NULL) {
                *multiplier = -*multiplier;
            }
            *(s + q) = *s;
        }
        *idx++ = j + q;
        if ((*adiag == 0.0) && (isnan != NULL) && (!*isnan)) {
            *isnan = j + 1;
        }
        x = 1.0 / *adiag;
        for (i = A->n - j - 1; i > 0; --i) {
            adiag += A->n; *adiag *= x; // Division by Bjj for j < i <= N
        }
        if (s != NULL){
            ++s;
        }
    }

    return(0);
}


int rawmatrix_lu_resolution(rawmatrix_t *LU, rawmatrix_t *B, unsigned int *idx)
{
    unsigned int i, j, k = 0;
    double *a, *ai, *bi, *bj, x;
    rawmatrix_t *A = LU;

    rcode = rawmatrix_check(LU);
    if (rcode != RAWMATRIX_OK) {
        return rcode;
    }
    rcode = rawmatrix_check(B);
    if (rcode != RAWMATRIX_OK) {
        return rcode;
    }
    if (A->n != A->m) {
        return RAWMATRIX_ERR_NOT_SQUARED;
    }
    if (idx == NULL) {
        return RAWMATRIX_ERR_NULL_INPUT;
    }
    if (B->m != 1) {
        return RAWMATRIX_ERR_NOT_A_VECTOR;
    }
    if (B->n != A->n) {
        return RAWMATRIX_ERR_SIZE_MISMATCH;
    }
    a = A->x;
    bi = B->x;

    if (idx != NULL){ // Permutation if idx is not NULL
        for (i = A->n; i > 0; --i){
            x = *(B->x + *idx);
            *(B->x + *idx++) = *bi;
            *bi++ = x;
        }
        bi = B->x;
    }

    for (i = A->n; i > 0; --i){
        if (k){ // k is used if the first elements of b are 0.0
            ai = a;
            bj = B->x + (A->n - k);
            for (j = k - i; j > 0; --j) {
                *bi -= (*ai++ * *bj++);
            }
        }
        else if (*bi != 0.0){ // The computation above will start when the first non null sum is found
            k = i;
            a += (A->n - k);
        }
        a += A->n;
        ++bi;
    }

    a -= (A->n - k);
    --a;

    for (i = A->n; i > 0; --i){
        ai = a + 1;
        bj = bi--;
        for (j = A->n - i; j > 0; --j){ *bi -= (*ai++ * *bj++); }
        *bi /= *a--;
        a -= A->n;
    }

    return(0);
}
