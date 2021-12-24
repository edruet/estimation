#include <stdlib.h>
#include <mem.h>
#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#include "Matrix.h"

#define mabs(x) (x >= 0.0 ? x : -x) 
#define TINY 1.0E-12
#define IsEqual(x,y,limit) ((x >= y ? (x - y) : (y - x)) <= limit ? 1 : 0)


int NumberIsEqual(double x,double y,double limit)
{
  return(IsEqual(x,y,limit));
}


// ----------------------------------------------------------------------------------------------------------------

int MatrixPrint(Matrix *A,char *format,char *msg)
{
  int i, j;
  double *a;
  char defformat[] = "%12.6lf";

  if (A == NULL){ return(1); }
  if ((a = A->x) == NULL){ return(2); }
  if (format == NULL){ format = defformat; }

  if (msg != NULL){ printf("\n%s\n",msg); }
  else { printf("\n"); }

  for (i = A->n ; i > 0; --i){
    for (j = A->m; j > 0; --j){
      printf(format,*a);
      ++a;
    }
    printf("\n");
  }

  return(0);
}


// ----------------------------------------------------------------------------------------------------------------

int MatrixSet(Matrix *A)
{
  if (A != NULL){ memset(A,0,sizeof(Matrix)); }
  return(0);
}


int MatrixAlloc(Matrix *A,unsigned int n,unsigned int m)
{
  if (A == NULL){ return(1); }
  if (n == 0){ return(2); }
  if (m == 0){ return(3); }

  A->n = n;
  A->m = m;

  if ((A->x = malloc(A->n * A->m * sizeof(double))) == NULL){ return(4); }

  return(0);
}


int MatrixRealloc(Matrix *A,unsigned int n,unsigned int m)
{
  if (A == NULL){ return(1); }
  if (n == 0){ return(2); }
  if (m == 0){ return(3); }

  A->n = n;
  A->m = m;

  if ((A->x = realloc(A->x,A->n * A->m * sizeof(double))) == NULL){ return(4); }

  return(0);
}


int MatrixFree(Matrix *A)
{
  if (A == NULL){ return(1); }
  if (A->x != NULL){ free(A->x); }
  memset(A,0,sizeof(Matrix));
  return(0);
}


int MatrixCopy(Matrix *A,Matrix *B)
{
  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if (A->x == NULL){ return(3); }
  if (B == A){ return(0); }
  if (B->n != A->n){ return(4); }
  if (B->m != A->m){ return(5); }

  memcpy(B->x,A->x,B->n * B->m * sizeof(double));

  return(0);
}

// ----------------------------------------------------------------------------------------------------------------

int MatrixEquality(Matrix *A,Matrix *B,double limit)
{
  unsigned int i;
  double *a, *b;

  if ((A == NULL) || (B == NULL) || (A->x == NULL) || (B->x == NULL) || (A->n != B->n) || (A->m != B->m)){ return(-1); }

  a = A->x;
  b = B->x;

  for (i = A->n * A->m; i > 0; --i){
    if (!IsEqual(*a,*b,limit)){ return((A->n * A->m) - i + 1); }
    ++a;
    ++b;
  }

  return(0);
}


int MatrixPermuteRows(Matrix *A,unsigned int n1,unsigned int n2)
{
  unsigned int i;
  double *a, *b, x;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (n1 == 0){ return(3);}
  if (n2 == 0){ return(4); }
  if (n1 > A->n){ return(5); }
  if (n2 > A->n){ return(6); }

  a = A->x + ((n1 - 1) * A->m);
  b = A->x + ((n2 - 1) * A->m);

  for (i = A->m; i > 0; --i){
    x = *a;
    *a++ = *b;
    *b++ = x;
  }

  return(0);
}


int MatrixPermuteCols(Matrix *A,unsigned int m1,unsigned int m2)
{
  unsigned int i;
  double *a, *b, x;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (m1 == 0){ return(3);}
  if (m2 == 0){ return(4); }
  if (m1 > A->m){ return(5); }
  if (m2 > A->m){ return(6); }

  a = A->x + (m1 - 1);
  b = A->x + (m2 - 1);

  for (i = A->m; i > 0; --i){
    x = *a;
    *a = *b;
    a += A->n;
    *b = x;
    b += A->n;
  }

  return(0);
}


int MatrixZero(Matrix *A)
{
  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (A->n == 0){ return(3); }
  if (A->m == 0){ return(4); }

  memset(A->x,0,A->n * A->m * sizeof(double));

  return(0);
}


int MatrixCst(Matrix *A,double x)
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if ((a = A->x) == NULL){ return(2); }

  for (i = A->n * A->m; i > 0; --i){ *a++ = x; }

  return(0);
}


int MatrixAdd(Matrix *A,double x)
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if ((a = A->x) == NULL){ return(2); }

  for (i = A->n * A->m; i > 0; --i){ *a++ += x; }

  return(0);
}


int MatrixMul(Matrix *A,double x)
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if ((a = A->x) == NULL){ return(2); }

  for (i = A->n * A->m; i > 0; --i){ *a++ *= x; }

  return(0);
}


int MatrixOp(Matrix *A,double (*Op)(double))
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if ((a = A->x) == NULL){ return(2); }
  if (Op == NULL){ return(3); }

  for (i = A->n * A->m; i > 0; --i){ *a++ = Op(*a); }

  return(0);
}


int MatrixSumElt(Matrix *A,Matrix *B)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (A->n != B->n){ return(5); }
  if (A->m != B->m){ return(6); }

  for (i = A->m * A->n; i > 0; --i){ *a++ += *b++; }

  return(0);
}


int MatrixSubElt(Matrix *A,Matrix *B)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (A->n != B->n){ return(5); }
  if (A->m != B->m){ return(6); }

  for (i = A->m * A->n; i > 0; --i){ *a++ -= *b++; }

  return(0);
}


int MatrixMulElt(Matrix *A,Matrix *B)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (A->n != B->n){ return(5); }
  if (A->m != B->m){ return(6); }

  for (i = A->m * A->n; i > 0; --i){ *a++ *= *b++; }

  return(0);
}


int MatrixDivElt(Matrix *A,Matrix *B)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (A->n != B->n){ return(5); }
  if (A->m != B->m){ return(6); }

  for (i = A->m * A->n; i > 0; --i){ *a++ /= *b++; }

  return(0);
}


int MatrixDivlElt(Matrix *A,Matrix *B,double limit,double x)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (A->n != B->n){ return(5); }
  if (A->m != B->m){ return(6); }

  for (i = A->m * A->n; i > 0; --i){
    if (((*b >= 0.0) && (*b < limit)) || ((*b < 0.0) && (-*b < limit))){ *a++ = x; }
    else { *a++ /= *b++; }
  }

  return(0);
}


int MatrixProd(Matrix *A,Matrix *B,Matrix *C)
{
  unsigned int i, j, k;
  double *a, *b, *c, *arow, *bcol;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if (C == NULL){ return(3); }
  if ((a = A->x) == NULL){ return(4); }
  if ((b = B->x) == NULL){ return(5); }
  if ((c = C->x) == NULL){ return(6); }
  if (A->m != B->n){ return(7); }
  if (C->n != A->n){ return(8); }
  if (C->m != B->m){ return(9); }

  arow = a;

  for (i = C->n; i > 0; --i){
    bcol = B->x;
    for (j = C->m; j > 0; --j){
      a = arow;
      b = bcol++;
      *c = 0.0;
      for (k = A->m; k > 0; --k){
	*c += (*a++ * *b);
	b += B->m;
      }
      ++c;
    }
    arow = a;
  }

  return(0);
}


// ----------------------------------------------------------------------------------------------------------------

double MatrixOpAssign(double a,double b)
{
  return(b);
}

double MatrixOpAdd(double a,double b)
{
  return(a + b);
}

double MatrixOpSub(double a,double b)
{
  return(a - b);
}

double MatrixOpMul(double a,double b)
{
  return(a * b);
}

double MatrixOpDiv(double a,double b)
{
  return(a / b);
}

double MatrixOpMod(double a,double b)
{
  return(fmod(a,b));
}

double MatrixOpPow(double a,double b)
{
  return(pow(a,b));
}

double MatrixOpMin(double a,double b)
{
  return(a <= b ? a : b);
}

double MatrixOpMax(double a,double b)
{
  return(a >= b ? a : b);
}


// ----------------------------------------------------------------------------------------------------------------

int MatrixGenOpUn(Matrix *A,double(*Operation)(double),unsigned int n,unsigned int m,int step,unsigned int nbr)
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (Operation == NULL){ return(3); }
  if (n == 0){ return(4); }
  if (m == 0){ return(5); }
  if (n > A->n){ return(6); }
  if (m > A->m){ return(7); }
  if (step == 0){ return(8); }
  if (nbr == 0){ return(9); }

  i = ((n - 1) * A->m) + (m - 1);

  if (step >= 0){
    if (((((unsigned int)step) * (nbr - 1)) + i) >= (A->n * A->m)){ return(10); }
  }
  else {
    if ((((unsigned int)(-step)) * (nbr - 1)) > i){ return(11); }
  }

  a = A->x + i;

  for (i = nbr; i > 0; --i){
    *a = Operation(*a);
    a += step;
  }

  return(0);
}


int MatrixGenOpCst(Matrix *A,double(*Operation)(double,double),unsigned int n,unsigned int m,int step,unsigned int nbr,double x)
{
  unsigned int i;
  double *a;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (Operation == NULL){ return(3); }
  if (n == 0){ return(4); }
  if (m == 0){ return(5); }
  if (n > A->n){ return(6); }
  if (m > A->m){ return(7); }
  if (step == 0){ return(8); }
  if (nbr == 0){ return(9); }

  i = ((n - 1) * A->m) + (m - 1);

  if (step >= 0){
    if (((((unsigned int)step) * (nbr - 1)) + i) >= (A->n * A->m)){ return(10); }
  }
  else {
    if ((((unsigned int)(-step)) * (nbr - 1)) > i){ return(11); }
  }

  a = A->x + i;

  for (i = nbr; i > 0; --i){
    *a = Operation(*a,x);
    a += step;
  }

  return(0);
}


int MatrixGenOpMulCst(Matrix *A,double(*Operation)(double,double),unsigned int n,unsigned int m,int step,unsigned int nbr,...)
{
  unsigned int i;
  double *a, x;
  va_list p;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }
  if (Operation == NULL){ return(3); }
  if (n == 0){ return(4); }
  if (m == 0){ return(5); }
  if (n > A->n){ return(6); }
  if (m > A->m){ return(7); }
  if (step == 0){ return(8); }
  if (nbr == 0){ return(9); }

  i = ((n - 1) * A->m) + (m - 1);

  if (step >= 0){
    if (((((unsigned int)step) * (nbr - 1)) + i) >= (A->n * A->m)){ return(10); }
  }
  else {
    if ((((unsigned int)(-step)) * (nbr - 1)) > i){ return(11); }
  }

  a = A->x + i;

  va_start(p,nbr);
  for (i = nbr; i > 0; --i){
    x = va_arg(p,double);
    *a = Operation(*a,x);
    a += step;
  }
  va_end(p);

  return(0);
}


int MatrixGenOpBin(Matrix *A,Matrix *B,double(*Operation)(double,double),
		   unsigned int na,unsigned int ma,int stepa,unsigned int nb,unsigned int mb,int stepb,unsigned int nbr)
{
  unsigned int i;
  double *a, *b;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if (A->x == NULL){ return(3); }
  if (B->x == NULL){ return(4); }
  if (Operation == NULL){ return(5); }
  if (na == 0){ return(6); }
  if (ma == 0){ return(7); }
  if (na > A->n){ return(8); }
  if (ma > A->m){ return(9); }
  if (nb == 0){ return(10); }
  if (mb == 0){ return(11); }
  if (nb > B->n){ return(12); }
  if (mb > B->m){ return(13); }
  if (stepa == 0){ return(14); }
  if (stepb == 0){ return(15); }
  if (nbr == 0){ return(16); }

  i = ((na - 1) * A->m) + (ma - 1);
  if (stepa >= 0){
    if ((((unsigned int)(stepa - 1)) * nbr) >= ((A->n * A->m) - i)){ return(17); }
  }
  else {
    if ((((unsigned int)(-stepa + 1)) * nbr) > i){ return(18); }
  }
  a = A->x + i;

  i = ((nb - 1) * B->m) + (mb - 1);
  if (stepb >= 0){
    if ((((unsigned int)(stepb - 1)) * nbr) >= ((B->n * B->m) - i)){ return(19); }
  }
  else {
    if ((((unsigned int)(-stepb + 1)) * nbr) > i){ return(20); }
  }
  b = B->x + i;

  for (i = nbr; i > 0; --i){
    *a = Operation(*a,*b);
    a += stepa;
    b += stepb;
  }

  return(0);
}


// ----------------------------------------------------------------------------------------------------------------

int MatrixTransp(Matrix *A)
{
  unsigned int i, j, k;
  double *a, *b, x;

  if (A == NULL){ return(1); }
  if (A->x == NULL){ return(2); }

  i = A->n;
  A->n = A->m;
  A->m = i;

  if ((A->n <= 1) || (A->m <= 1)){ return(0); }

  if (A->n == A->m){
    k = 0;
    a = A->x;
    for (i = A->n - 1; i > 0; --i){
      b = a++;
      for (j = i; j > 0; --j){
	b += A->n;
	x = *a;
	*a++ = *b;
	*b = x;
      }
      a += ++k;
    }
  }
  else {
    a = A->x + ((A->n * A->m) - 2);

    for (i = (A->n * A->m) - 2; i > 0; --i){
      j = i;
      do { j = A->n * (j % A->m) + j / A->m; } while (j > i);
      if (j != i){
        x = *a;
        *a = *(A->x + j);
        *(A->x + j) = x;
      }
      --a;
    }
  }

  return(0);
}


int MatrixCopyTransp(Matrix *A,Matrix *B)
{
  unsigned int i, j;
  double *a, *acol, *b;

  if (A == NULL){ return(1); }
  if ((acol = A->x) == NULL){ return(2); }
  if (B == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (B->m != A->n){ return(5); }
  if (B->n != A->m){ return(6); }
  if (B == A){ return(MatrixTransp(A)); }

  for (i = B->n; i > 0; --i){
    a = acol++;
    for (j = B->m; j > 0; --j){
      *b++ = *a;
      a += A->m;
    }
  }

  return(0);
}


int MatrixAATransp(Matrix *A,Matrix *B)
{
  unsigned int i, j, k;
  double *b, *a1, *a2;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((a1 = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (B->n != A->n){ return(5); }
  if (B->m != B->m){ return(6); }

  for (i = B->n; i > 0; --i){
    a2 = A->x;
    for (j = B->n; j > 0; --j){
      *b = 0.0;
      for (k = A->m; k > 0; --k){ *b += (*a1++ * *a2++); }
      ++b;
      a1 -= A->m;
    }
    a1 += A->m;
  }

  return(0);
}


int MatrixATranspA(Matrix *A,Matrix *B)
{
  unsigned int i, j, k;
  double *b, *a1, *a2, *col1, *col2;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if ((col1 = A->x) == NULL){ return(3); }
  if ((b = B->x) == NULL){ return(4); }
  if (B->n != A->m){ return(5); }
  if (B->m != B->m){ return(6); }

  for (i = B->n; i > 0; --i){
    col2 = A->x;
    for (j = B->n; j > 0; --j){
      *b = 0.0;
      a1 = col1;
      a2 = col2++;
      for (k = A->n; k > 0; --k){
	*b += (*a1 * *a2);
	a1 += A->m;
	a2 += A->m;
      }
      ++b;
    }
    ++col1;
  }

  return(0);
}


// ----------------------------------------------------------------------------------------------------------------

int MatrixLUDcmp(Matrix *A,unsigned int *idx,double *s,int *multiplier,int *isnan)
{
  int i, j, k, q;
  double x, xmax, *a, *ai, *aj, *acol, *adiag, *sj;

  if (A == NULL){ return(1); }
  if (A->n != A->m){ return(2); }
  if (A->x == NULL){ return(3); }

  if (multiplier != NULL){ *multiplier = 1; }
  if (isnan != NULL){ *isnan = 0; }

  // Selection of the biggest element of each row
  if ((idx != NULL) && (s != NULL)){
    a = A->x;
    sj = s;
    for (i = A->n; i > 0; --i){
      *sj = 0.0;
      for (j = A->n; j > 0; --j){
        if ((x = mabs(*a)) > *sj){ *sj = x; }
        ++a;
      }
      if ((*sj == 0.0) && (isnan != NULL) && (!*isnan)){ *isnan = A->n - i + 1; }
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
	if (s != NULL){ x = (*sj++ * mabs(*a)); }
        else { x = mabs(*a); }
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
      for (k = A->n; k > 0; --k){ // Exchanging line adiag (j) and amax
	x = *ai;
	*ai++ = *aj;
	*aj++ = x;
      }
      if (multiplier != NULL){ *multiplier = -*multiplier; }
      *(s + q) = *s;
    }
    *idx++ = j + q;
    if ((*adiag == 0.0) && (isnan != NULL) && (!*isnan)){ *isnan = j + 1; }
    x = 1.0 / *adiag;
    for (i = A->n - j - 1; i > 0; --i){ adiag += A->n; *adiag *= x; } // Division by Bjj for j < i <= N
    if (s != NULL){ ++s; }
  }

  return(0);
}


int MatrixLUSubst(Matrix *LU,Matrix *B,unsigned int *idx)
{
  unsigned int i, j, k = 0;
  double *a, *ai, *bi, *bj, x;
  Matrix *A = LU;

  if (A == NULL){ return(1); }
  if (B == NULL){ return(2); }
  if (A->n != A->m){ return(3); }
  if (B->m != 1){ return(4); }
  if (B->n != A->n){ return(5); }
  if ((a = A->x) == NULL){ return(6); }
  if ((bi = B->x) == NULL){ return(7); }

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
      for (j = k - i; j > 0; --j){ *bi -= (*ai++ * *bj++); }
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


int MatrixInv(Matrix *LU,Matrix *B,unsigned int *idx)
{
  int i;
  double *b;
  Matrix *A = LU;

  if (A == NULL){ return(1); }
  if (A->n != A->m){ return(2); }
  if (A->x == NULL){ return(3); }
  if (B == A){ return(4); }
  if (B == NULL){ return(5); }
  if (B->n != B->m){ return(6); }
  if ((b = B->x) == NULL){ return(7); }
  if (B->n != A->n){ return(8); }

  B->m = 1;

  memset(B->x,0,B->n * B->n * sizeof(double));

  for (i = B->n; i > 0; --i){
    *(B->x + (B->n - i)) = 1.0;
    MatrixLUSubst(A,B,idx);
    B->x += B->n;
  }

  B->x = b;
  B->m = B->n;

  MatrixTransp(B);

  return(0);
}


int MatrixDet(Matrix *LU,int multiplier,double *determinant)
{
  unsigned int i;
  double *a;
  Matrix *A = LU;

  if (A == NULL){ return(1); }
  if (A->n != A->m){ return(2); }
  if ((a = A->x) == NULL){ return(3); }
  if (determinant == NULL){ return(4); }

  *determinant = (double) multiplier;

  for (i = A->n; i > 0; --i){
    *determinant *= *a++;
    a += A->n;
  }

  return(0);
}


int MatrixLUSplit(Matrix *LU,Matrix *L,Matrix *U)
{
  unsigned int i, j;
  double *acol, *a, *bcol, *b;
  Matrix *A = LU;

  if (A == NULL){ return(1); }
  if (L == NULL){ return(2); }
  if (U == NULL){ return(2); }
  if ((acol = A->x) == NULL){ return(4); }
  if ((bcol = L->x) == NULL){ return(5); }
  if (A->n != A->m){ return(6); }
  if (L->n != L->m){ return(7); }
  if (U->n != U->m){ return(8); }
  if (L->n != A->n){ return(9); }
  if (U->n != A->n){ return(10); }

  memset(L->x,0,L->n * L->m * sizeof(double));
  memset(U->x,0,U->n * U->m * sizeof(double));

  for (i = A->n; i > 0; --i){
    acol += A->n;
    a = acol++;
    *bcol = 1.0;
    bcol += A->n;
    b = bcol++;
    for (j = i - 1; j > 0; --j){
      *b = *a;
      a += A->n;
      b += A->n;
    }
  }

  a = A->x;
  b = U->x;

  for (i = A->n; i > 0; --i){
    a += (A->n - i);
    b += (A->n - i);
    for (j = i; j > 0; --j){ *b++ = *a++;  }
  }

  return(0);
}


int MatrixReorderRows(Matrix *A,unsigned int *idx)
{
  unsigned int i, j;
  double *a, *b, x;

  if (A == NULL){ return(1); }
  if (idx == NULL){ return(2); }
  if ((a = A->x) == NULL){ return(3); }

  idx += A->n;

  for (i = A->n; i > 0; --i){
    --idx;
    if (*idx != (i - 1)){
      a = A->x + (A->n * (i - 1));
      b = A->x + (A->n * *idx);
      for (j = A->n; j > 0; --j){
        x = *a;
        *a++ = *b;
        *b++ = x;
      }
    }
  }

  return(0);
}

// ----------------------------------------------------------------------------------------------------------------
