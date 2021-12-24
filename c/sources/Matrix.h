#ifndef __MATHLIBRARY_MATRIX_H

#define __MATHLIBRARY_MATRIX_H

typedef struct struct_Matrix
{
  unsigned int n;
  unsigned int m;
  double *x;
}Matrix;


#ifdef __cplusplus
extern "C"
{
#endif

  int NumberIsEqual(double x,double y,double limit);

  // Display
  int MatrixPrint(Matrix *A,char *format,char *msg);

  // Allocation
  int MatrixSet(Matrix *A);
  int MatrixAlloc(Matrix *A,unsigned int n,unsigned int m);
  int MatrixRealloc(Matrix *A,unsigned int n,unsigned int m);
  int MatrixFree(Matrix *A);
  int MatrixCopy(Matrix *A,Matrix *B);

  // Basic operations
  int MatrixEquality(Matrix *A,Matrix *B,double limit);
  int MatrixPermuteRows(Matrix *A,unsigned int n1,unsigned int n2);
  int MatrixPermuteCols(Matrix *A,unsigned int m1,unsigned int m2);
  int MatrixZero(Matrix *A);
  int MatrixCst(Matrix *A,double x);
  int MatrixAdd(Matrix *A,double x);
  int MatrixMul(Matrix *A,double x);
  int MatrixOp(Matrix *A,double (*Op)(double));
  int MatrixSumElt(Matrix *A,Matrix *B);
  int MatrixSubElt(Matrix *A,Matrix *B);
  int MatrixMulElt(Matrix *A,Matrix *B);
  int MatrixDivElt(Matrix *A,Matrix *B);
  int MatrixDivlElt(Matrix *A,Matrix *B,double limit,double x);
  int MatrixProd(Matrix *A,Matrix *B,Matrix *C);

  // Generic operations
  double MatrixOpAssign(double a,double b);
  double MatrixOpAdd(double a,double b);
  double MatrixOpSub(double a,double b);
  double MatrixOpMul(double a,double b);
  double MatrixOpDiv(double a,double b);
  double MatrixOpMod(double a,double b);
  double MatrixOpPow(double a,double b);
  double MatrixOpMin(double a,double b);
  double MatrixOpMax(double a,double b);
  int MatrixGenOpUn(Matrix *A,double(*Operation)(double),unsigned int n,unsigned int m,int step,unsigned int nbr);
  int MatrixGenOpCst(Matrix *A,double(*Operation)(double,double),unsigned int n,unsigned int m,int step,unsigned int nbr,double x);
  int MatrixGenOpMulCst(Matrix *A,double(*Operation)(double,double),unsigned int n,unsigned int m,int step,unsigned int nbr,...);
  int MatrixGenOpBin(Matrix *A,Matrix *B,double(*Operation)(double,double),
		     unsigned int na,unsigned int ma,int stepa,unsigned int nb,unsigned int mb,int stepb,unsigned int nbr);

  // Tranpose
  int MatrixTransp(Matrix *A);
  int MatrixCopyTransp(Matrix *A,Matrix *B);
  int MatrixAATransp(Matrix *A,Matrix *B);
  int MatrixATranspA(Matrix *A,Matrix *B);

  // LU decomposition
  int MatrixLUDcmp(Matrix *A,unsigned int *idx,double *s,int *multiplier,int *isnan);
  int MatrixLUSubst(Matrix *LU,Matrix *B,unsigned int *idx);
  int MatirxInv(Matrix *LU,Matrix *B,unsigned int *idx);
  int MatrixDet(Matrix *LU,int multiplier,double *determinant);
  int MatrixLUSplit(Matrix *LU,Matrix *L,Matrix *U);
  int MatrixReorderRows(Matrix *A,unsigned int *idx);

#ifdef __cplusplus
}
#endif


#endif
