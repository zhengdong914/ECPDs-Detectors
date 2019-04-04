/*** mxutils.h ***/

#ifndef MXUTILS_H
#define MXUTILS_H

mxArray *mxGetColumn (const mxArray *mat, int col);
double   mxSumSquareArray (const mxArray *mat);
mxArray *scmult (const mxArray *mat, double s);
mxArray *dec (mxArray *A, const mxArray *B);
mxArray *minus (const mxArray *A, const mxArray *B);
mxArray *inc (mxArray *A, const mxArray *B);
mxArray *plus (const mxArray *A, const mxArray *B);
double   mxInnerProd (const mxArray *A, const mxArray *B);
mxArray *ones(int r, int c);
mxArray *eye(int r);
mxArray *mult (const mxArray *A, const mxArray *B);
mxArray *tr(const mxArray *A);
mxArray *subm(const mxArray *A, int k);
mxArray *subsubm(const mxArray *A, int k, int l);
void subcopy(mxArray *A, int k, const mxArray *B);
void subsubcopy(mxArray *A, int k, int l, const mxArray *B);
void subinc(mxArray *A, int k, const mxArray *B);
void subdec(mxArray *A, int k, const mxArray *B);
void copy (mxArray *A,const mxArray *B);
mxArray *diag(const mxArray *A);
mxArray *rdiv(const mxArray *M,const mxArray *v);
mxArray *cdiv(const mxArray *M,const mxArray *v);
int disp(const mxArray *A);
double trace (const mxArray *A);

typedef double	mrReal;

extern void  lubksb(mrReal **a, int n, int *indx, mrReal *b);
extern void  ludcmp(mrReal **a, int n, int *indx, mrReal *d);
 
extern void inv( const mxArray  *a,  /* in: matrix to be inverted */ mxArray *mxAInv, /* out: Inverse matrix */ mxArray *mxDet /* out: Determinant */ ) ; 
extern int getSquareDim( const mxArray *a ) ;
extern void copyMxToNR( const mxArray *mxA, mrReal **b ) ;
extern int getDimensions( const mxArray *a, int **dim ) ;
extern void mxResize( mxArray *a, int rows, int cols ) ;
extern double det(const mxArray *a);
mxArray *symm (mxArray *A);

#endif

