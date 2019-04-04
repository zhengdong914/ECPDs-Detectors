/*** mxutils.c ***/

#include <stdlib.h>
#include <math.h>

#include "mex.h"
#include "matrix.h"
#include "mxutils.h"
#include "nrutil.h"

/* mxGetColumn
 *
 * Returns the specified column from the matrix 
 */
mxArray *mxGetColumn (const mxArray *mat, int col)
{
  mxArray *ret;
  double *pr, *matPr;
  int idx,numRows;

  numRows = mxGetM(mat);
  ret = mxCreateDoubleMatrix(numRows,1,mxREAL);
  pr = mxGetPr(ret);
  matPr = mxGetPr(mat);
  for (idx=0; idx<numRows; idx++) {
    pr[idx] = matPr[numRows*col+idx];
  }
  return(ret);
}

/* mxSumSquareArray
 *
 * Computes the sum of the squares of the array elements.
 */
double mxSumSquareArray (const mxArray *mat)
{
  int idx, maxIdx;
  double *pr, ret=0.0;

  maxIdx = mxGetNumberOfElements(mat);
  pr = mxGetPr(mat);
  for (idx=0; idx<maxIdx;idx++)
    ret += pr[idx]*pr[idx];
  
  return ret;
}

/* scmult
 *
 * A multiplies a matlab array by a scalar value
 * and returns the result.
 */
mxArray *scmult (const mxArray *mat, double s)
{
  int idx, maxIdx;
  double *pr;
  double *prR;
  mxArray *R;

  R=mxCreateDoubleMatrix(mxGetM(mat),mxGetN(mat),mxREAL);
  maxIdx = mxGetNumberOfElements(mat);
  pr = mxGetPr(mat);
  prR = mxGetPr(R);
  for (idx=0; idx<maxIdx; idx++)
    prR[idx]=pr[idx] * s;

  return R;
}

/* dec - decrement
 *
 * Computes A minus B, stores the result in A and
 * returns a pointer to A. It assumes that A and 
 * B are the same size.
 */
mxArray *dec (mxArray *A, const mxArray *B)
{
  int idx, maxIdx;
  double *prA, *prB;

  prA = mxGetPr(A);
  prB = mxGetPr(B);
  maxIdx = mxGetNumberOfElements(A);
  for (idx=0; idx<maxIdx;idx++)
    prA[idx] -= prB[idx];

  return A;
}

/* minus
 *
 * Computes A minus B, returns new matrix.
 */

mxArray *minus (const mxArray *A, const mxArray *B)
{
  int idx, maxIdx;
  mxArray *R;
  double *prA, *prB, *prR;

  R=mxCreateDoubleMatrix(mxGetM(A),mxGetN(A),mxREAL);
 
  prA = mxGetPr(A);
  prB = mxGetPr(B);
  prR = mxGetPr(R);
  maxIdx = mxGetNumberOfElements(A);
  for (idx=0;idx<maxIdx;idx++)
    prR[idx] = prA[idx] - prB[idx];

  return R;
}

/* inc
 *
 * Computes A plus B, stores the result in A and
 * returns a pointer to A. It assumes that A and
 * B are the same size.
 */
mxArray *inc (mxArray *A, const mxArray *B)
{
  int idx, maxIdx, maxB;
  double *prA, *prB;

  prA = mxGetPr(A);
  prB = mxGetPr(B);
  maxIdx = mxGetNumberOfElements(A);
  maxB=mxGetNumberOfElements(B);
  if (maxIdx != maxB)
    mexErrMsgTxt("inc: sizes don't match");

  for (idx=0;idx<maxIdx;idx++)
    prA[idx] += prB[idx];

  return A;
}

/* plus
 *
 * Computes A plus B, stores the result in A and
 * returns a pointer to A. It assumes that A and
 * B are the same size.
 */
mxArray *plus (const mxArray *A, const mxArray *B)
{
  int idx, maxIdx;
  mxArray *R;
  double *prA, *prB, *prR;

  R=mxCreateDoubleMatrix(mxGetM(A),mxGetN(A),mxREAL);
 
  prA = mxGetPr(A);
  prB = mxGetPr(B);
  prR = mxGetPr(R);
  maxIdx = mxGetNumberOfElements(A);
  for (idx=0;idx<maxIdx;idx++)
    prR[idx] = prA[idx] + prB[idx];

  return R;
}

/* mxInnerProd
 * 
 * Computes A'*B where A and B are vectors.
 */
double mxInnerProd (const mxArray *A, const mxArray *B)
{
  int idx, maxIdx;
  double *prA, *prB;
  double ret  = 0.0;

  prA = mxGetPr(A);
  prB = mxGetPr(B);
  maxIdx = mxGetNumberOfElements(A);
  for(idx=0; idx<maxIdx;idx++)
    ret += prA[idx]*prB[idx];
  return ret;
}

/*******************************
 * mult
 *
 * Computes A*B where A and B are arbitrary matrices.
 * It creates a new matrix for storing the result.
 * !! You are responsible for freeing it. !!
 */
mxArray * mult (const mxArray *A, const mxArray *B)
{
  int i,j,k,maxI,maxJ,maxK,rowA,rowB,rowR;
  mxArray *R;
  double *prA,*prB,*prR;

  R=mxCreateDoubleMatrix(mxGetM(A),mxGetN(B),mxREAL);
  /*Note: R already contains zeros*/

  maxI=rowA=rowR=mxGetM(A);
  maxJ=mxGetN(A);
  rowB=mxGetM(B);
  maxK=mxGetN(B);

  prA = mxGetPr(A);
  prB = mxGetPr(B);
  prR = mxGetPr(R);

  for (i=0; i<maxI; i++)
    for(j=0; j<maxJ; j++)
      for(k=0; k<maxK; k++)
	prR[i+rowR*k] += prA[i+rowA*j]*prB[j+rowB*k];

  return R;	
}

/*******************************
 * ones 
 */

mxArray *ones(int r, int c)
{
  int i,j;
  mxArray *R;
  double *pR;

  R=mxCreateDoubleMatrix(r,c,mxREAL);
  
  pR=mxGetPr(R);
  
  for (i=0;i<r;i++)
    for (j=0;j<c;j++)
      pR[i+r*j]=1.0;

  return R;
}

/*******************************
 * identity matrix
 */

mxArray *eye(int r)
{
  int i,j;
  mxArray *R;
  double *pR;

  R=mxCreateDoubleMatrix(r,r,mxREAL);
  
  pR=mxGetPr(R);
  
  for (i=0;i<r;i++)
    pR[i+r*i]=1.0;
  
  return R;
}

/*******************************
 * transpose
 */

mxArray *tr(const mxArray *A)
{
  int r,c,i,j;
  mxArray *R;
  double *pR, *pA;

  r=mxGetM(A);
  c=mxGetN(A);
  R=mxCreateDoubleMatrix(c,r,mxREAL);
  pR=mxGetPr(R);
  pA=mxGetPr(A);

  for (i=0;i<r;i++)
    for (j=0;j<c;j++)
      pR[j+c*i]=pA[i+r*j];

  return R;
}
  
/*******************************
 * subsubmatrix operation from an ND array
 * creates and returns a new array 
 * A(:,..,k,l);
 */

mxArray *subsubm(const mxArray *A, int k, int l)
{
  int d,i,num;
  const int *dims;
  int *Rdims;
  mxArray *R;
  double *pR, *pA;

  d=mxGetNumberOfDimensions(A);
  dims=mxGetDimensions(A);
  
  if (l>dims[d-1] | k>dims[d-2])
    mexErrMsgTxt("subsubm: dimension out of bounds");
    
  if (d==2)
    {
      Rdims=(int *)mxCalloc(2,sizeof(int)); 
      Rdims[0]=1; Rdims[1]=1;
      R=mxCreateNumericArray(2,Rdims,mxDOUBLE_CLASS,mxREAL);
    }
  else if (d==3)
    {
      Rdims=(int *)mxCalloc(2,sizeof(int)); 
      Rdims[0]=dims[0]; Rdims[1]=1;
      R=mxCreateNumericArray(2,Rdims,mxDOUBLE_CLASS,mxREAL);
    }
  else /* d>3 */
    {
      Rdims=(int *)mxCalloc(d-2,sizeof(int)); 
      for (i=0;i<d-2;i++)
	Rdims[i]=dims[i];
      R=mxCreateNumericArray(d-2,Rdims,mxDOUBLE_CLASS,mxREAL);
    }

  pR=mxGetPr(R);
  pA=mxGetPr(A);

  num=mxGetNumberOfElements(R);
  
  for (i=0;i<num;i++)
    pR[i]=pA[((l-1)*dims[d-2]+k-1)*num+i];
  
  return R;
}
/*******************************
 * submatrix operation from an ND array
 * creates and returns a new array 
 */

mxArray *subm(const mxArray *A, int k)
{
  int d,i,num;
  const int *dims;
  int *Rdims;
  mxArray *R;
  double *pR, *pA;

  d=mxGetNumberOfDimensions(A);
  dims=mxGetDimensions(A);
  
  if (k>dims[d-1])
    mexErrMsgTxt("subm: dimension out of bounds");
    
  if (d==2)
    {
      Rdims=(int *)mxCalloc(2,sizeof(int)); 
      Rdims[0]=dims[0]; Rdims[1]=1;
      R=mxCreateNumericArray(2,Rdims,mxDOUBLE_CLASS,mxREAL);
    }
  else /* d>2 */
    {
      Rdims=(int *)mxCalloc(d-1,sizeof(int)); 
      for (i=0;i<d-1;i++)
	Rdims[i]=dims[i];
      R=mxCreateNumericArray(d-1,Rdims,mxDOUBLE_CLASS,mxREAL);
    }
      
  pR=mxGetPr(R);
  pA=mxGetPr(A);

  num=mxGetNumberOfElements(R);
  
  for (i=0;i<num;i++)
    pR[i]=pA[(k-1)*num+i];
  
  return R;
}
      
/*******************************
 * copy the N dim array B into kth major position of 
 * the N+1 dim array A 
 */

void subcopy(mxArray *A, int k, const mxArray *B)
{
  int dA, dB,sB,i,num;
  const int *dimsA;
  const int *dimsB;
  mxArray *R;
  double *pR, *pA, *pB;
  char buffer[100];
  
  dA=mxGetNumberOfDimensions(A);
  dimsA=mxGetDimensions(A);
  dB=mxGetNumberOfDimensions(B);
  dimsB=mxGetDimensions(B);

  if (dB != dA-1 && dA>2)
    mexErrMsgTxt("subcopy: B does not have one dimension less than A"); 
  
  if (dA==2)
    {
      if (dimsA[0] != dimsB[0] || dimsB[1] !=1)
	mexErrMsgTxt("subcopy: The sizes of B and A do not match");
    }
  else
    for (i=0;i<dB;i++)
      if (dimsA[i] != dimsB[i])
	mexErrMsgTxt("subcopy: The sizes of B and A do not match");
  
  sB=mxGetNumberOfElements(B);
  pB=mxGetPr(B);
  pA=mxGetPr(A);
  
  for (i=0;i<sB;i++)
    pA[(k-1)*sB+i]=pB[i];
}

/*******************************
 * copy the N dim array B into (k,l)th position of 
 * the N+2 dim array A 
 * A(:...,k,l)=B;
 */

void subsubcopy(mxArray *A, int k, int l, const mxArray *B)
{
  int dA, dB,sB,i,num,c;
  const int *dimsA;
  const int *dimsB;
  mxArray *R;
  double *pR, *pA, *pB;
  char buffer[100];
  
  dA=mxGetNumberOfDimensions(A);
  dimsA=mxGetDimensions(A);
  dB=mxGetNumberOfDimensions(B);
  dimsB=mxGetDimensions(B);
  
  if (dB != dA-2 && dA>3)
    mexErrMsgTxt("subsubcopy: A does not have two dimensions less than B");

  if ((dA==2 &&  (dimsB[0] != 1 || dimsB[1] != 1 || dB !=2)) ||
      (dA==3 &&  (dimsB[0] != dimsA[0] || dimsB[1] != 1 || dB !=2)))
      mexErrMsgTxt("subsubcopy: the sizes of B and A do not match");
  if (dA > 3)
    for (i=0;i<dA-2;i++)
      if (dimsB[i] != dimsA[i])
	mexErrMsgTxt("subsubcopy: the sizes of B and A do not match");

  if (l>dimsA[dA-1] || k>dimsA[dA-2])
    mexErrMsgTxt("subsubcopy: indexing outside the range of A");
    
  sB=mxGetNumberOfElements(B);
  pB=mxGetPr(B);
  pA=mxGetPr(A);

  for (i=0;i<sB;i++)
    pA[((l-1)*dimsA[dA-2]+k-1)*sB+i]=pB[i];
}

/*******************************
 * increment  the kth major position of 
 * the N+1 dim array A by  the N dim array B
 */

void subinc(mxArray *A, int k, const mxArray *B)
{
  int dA, dB,sB,i,num;
  const int *dimsA;
  const int *dimsB;
  mxArray *R;
  double *pR, *pA, *pB;
  char buffer[100];

  dA=mxGetNumberOfDimensions(A);
  dimsA=mxGetDimensions(A);
  dB=mxGetNumberOfDimensions(B);
  dimsB=mxGetDimensions(B);

  if (dB != dA-1 && dA>2)
    {
      sprintf(buffer, "%s: subinc: %s does not have one dimension less than %s",
	      __FILE__, mxGetName(B), mxGetName(A));
      mexErrMsgTxt(buffer);
    }
  if (dA>2)
    {
      for (i=0;i<dB;i++)
	if (dimsA[i] != dimsB[i])
	  mexErrMsgTxt("subinc: The sizes of B and A do not match");
    }
  else if (dA==2)
    {
      if (dimsA[0] != dimsB[0] && dimsB[1] != 1 )
	mexErrMsgTxt("subinc: The sizes of B and A do not match");
    }
    
  sB=mxGetNumberOfElements(B);
  pB=mxGetPr(B);
  pA=mxGetPr(A);

  for (i=0;i<sB;i++)
    pA[(k-1)*sB+i] += pB[i];
}
/*******************************
 * decrement the kth major position of 
 * the N+1 dim array A by  the N dim array B
 */

void subdec(mxArray *A, int k, const mxArray *B)
{
  int dA, dB,sB,i,num;
  const int *dimsA;
  const int *dimsB;
  mxArray *R;
  double *pR, *pA, *pB;
  char buffer[100];

  dA=mxGetNumberOfDimensions(A);
  dimsA=mxGetDimensions(A);
  dB=mxGetNumberOfDimensions(B);
  dimsB=mxGetDimensions(B);

  if (dB != dA-1 && dA>2)
    {
      sprintf(buffer, "%s: subdec: %s does not have one dimension less than %s",
	      __FILE__, mxGetName(B), mxGetName(A));
      mexErrMsgTxt(buffer);
    }
  for (i=0;i<dB;i++)
    if (dimsA[i] != dimsB[i])
     mexErrMsgTxt("subdec: The sizes of B and A do not match");

  sB=mxGetNumberOfElements(B);
  pB=mxGetPr(B);
  pA=mxGetPr(A);

  for (i=0;i<sB;i++)
    pA[(k-1)*sB+i] -= pB[i];
}

/*******************************
 * copy
 *
 * Copies B into A.
 */
void copy(mxArray *A,const mxArray *B)
{
  int i,maxI, dA, dB;
  const int *dimsA, *dimsB;
  double *prA,*prB;
  char buffer[100];
  
  dA=mxGetNumberOfDimensions(A);
  dimsA=mxGetDimensions(A);
  dB=mxGetNumberOfDimensions(B);
  dimsB=mxGetDimensions(B);

  if (dA !=dB)
    mexErrMsgTxt("copy: The dimensions of B and A do not match");

  for (i=0;i<dA;i++)
    if (dimsA[i] != dimsB[i])
      {
	sprintf(buffer, "copy: The sizes of B and A do not match on dimension %d",i+1);
	mexErrMsgTxt(buffer);
      }

  maxI=mxGetNumberOfElements(B);
  prA=mxGetPr(A);
  prB=mxGetPr(B);

  for(i=0; i<maxI; i++)
    prA[i]=prB[i];
}

/*******************************
 * diag -- create diagonal matrix from vector
 */

mxArray *diag(const mxArray *A)
{
  int r,i;
  mxArray *R;
  double *pR, *pA;

  r = mxGetNumberOfElements(A);
  R=mxCreateDoubleMatrix(r,r,mxREAL);
  pR=mxGetPr(R);
  pA=mxGetPr(A);

  for (i=0;i<r;i++)
    pR[i+r*i]=pA[i];

  return R;
}

/* Divide M by v row-wise*/
mxArray *rdiv(const mxArray *M,const mxArray *v)
{
  mxArray	*R;
  int m,n,i,j;
  double *pR, *pM, *pv;

  m=mxGetM(M);
  n=mxGetN(M);
  
  if (mxGetM(v) != m || mxGetN(v) !=1)
    mexErrMsgTxt("RDIV: row sizes don't match!\n");

  R = mxCreateDoubleMatrix(m,n,mxREAL);  

  pR = mxGetPr(R);
  pM = mxGetPr(M);
  pv = mxGetPr(v);

  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      pR[j*m+i] = pM[j*m + i]/pv[i];

  return R;

}

/* Divide M by v column-wise*/
mxArray *cdiv(const mxArray *M,const mxArray *v)
{
  mxArray	*R;
  int m,n,i,j;
  double *pR, *pM, *pv;

  m=mxGetM(M);
  n=mxGetN(M);
  
  if (mxGetN(v) != n || mxGetM(v) !=1)
    mexErrMsgTxt("CDIV: row sizes don't match!\n");

  R = mxCreateDoubleMatrix(m,n,mxREAL);  

  pR = mxGetPr(R);
  pM = mxGetPr(M);
  pv = mxGetPr(v);

  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      pR[j*m+i] = pM[j*m + i]/pv[j];

  return R;

}

/**************************************************************
 *
 * MODULE DESCRIPTION
 *
 * Routines for computing the inverse of a square matrix.
 * Can be called from a mex function
 *
 *  Calling sequence
 *
 *   inv(const mxArray *a, mxArray *inv, mxArray *det);
 *
 *  The caller should allocate the return arrays (inv and det),
 *  If they are allocated with the wrong dimensions, they are
 *  silently resized.
 *
 * $Author: Mike Revow$
 *
 *************************************************************/

void inv(
const mxArray 	*a,	  /* in: matrix to be inverted */
mxArray		*mxAInv,  /* out: Inverse matrix */
mxArray		*mxDet   /*  out: Determinant */
)
{ 
  int  		i, j, N, Ni;
  static int	*indx= NULL;
  static mrReal	*col = NULL;
  static mrReal **B = NULL;
  static int	size = 0;
  mrReal	d;
  int		dims[2];
  double 	*aInv, *pa, *det;

  N = getSquareDim(a);
  if (N <= 0)
    {
      mexErrMsgTxt("Error in inv: N<0\n");
      return;
    }
  
  if (N != getSquareDim((const mxArray *)mxAInv))
    mxResize(mxAInv, N, N);
  
  if (1 != getSquareDim((const mxArray *)mxDet))
    mxResize(mxDet, 1, 1);
  
  if (N > size)
  {
    if (B)  free_matrix(B, 1, size, 1, size);
    if (indx)  free_ivector(indx, 1, size);
    if (col) free_vector(col, 1, size);

    B = matrix(1, N, 1, N);
    col = vector(1, N);
    indx = ivector(1, N);
    size = N;
  }

  copyMxToNR(a,B);

  aInv = mxGetPr(mxAInv);
  pa = aInv;
  det = mxGetPr(mxDet);

  ludcmp(B, N, indx, &d);
  for (j=1; j <= N; j++)
  {
    for (i=1; i<= N; i++)
      col[i] = 0.0;
    
    col[j] = 1.0;

    lubksb(B, N, indx, col);

    for (i=1; i<= N; i++)
      *pa++ = col[i];
  }
  *det = 1.0; 
  for (j=1; j<= N; j++) 
    *det *= B[j][j];
  *det *= (double) d;
}    

double det(
const mxArray 	*a	  /* in: matrix  */
)
{ 
  int  		i, j, N;
  static int	*indx= NULL;
  static mrReal	*col = NULL;
  static mrReal **B = NULL;
  static int	size = 0;
  mrReal	d;
  double detv;

  N = getSquareDim(a);
  if (N <= 0)
    return -999.0;

  if (N > size)
    {
      if (B)      free_matrix(B, 1, size, 1, size);
      if (indx)      free_ivector(indx, 1, size);
      if (col)      free_vector(col, 1, size);
      B = matrix(1, N, 1, N);
      col = vector(1, N);
      indx = ivector(1, N);
      size = N;
    }

  copyMxToNR(a,B);

  ludcmp(B, N, indx, &d);

  detv = (double ) d; 
  for (j=1; j<= N; j++)     detv *= B[j][j];

  return detv;
}    


int getSquareDim(
const mxArray		*a
)
{
  int		N;
  int		*dims;
  char		buffer[100];

  if ( 2 != getDimensions(a, &dims))
  {
    sprintf(buffer, "%s: getSquareDim: Expected array %s to have 2 dimensions, got %d",
	    __FILE__, mxGetName(a), mxGetNumberOfDimensions(a));
    mexErrMsgTxt(buffer);
  }

  if (dims[0] != dims[1])
  {
    sprintf(buffer, "%s: Matrix %s is not square",
	    __FILE__, mxGetName(a));
    mexErrMsgTxt(buffer);
  }

  return (dims[0]);
}

void copyMxToNR(
const mxArray	*mxA,
mrReal		**b
)
{
  int		i,j, M,N;
  double	*a;
  int		*dims;
  char		buffer[100];

  if ( 2 != getDimensions(mxA, &dims))
  {
    sprintf(buffer, "%s: copyMxToNR: Expected array %s to have 2 dimensions, got %d",
	    __FILE__, mxGetName(mxA), mxGetNumberOfDimensions(mxA));
    mexErrMsgTxt(buffer);
  }

  M = dims[0];
  N = dims[1];
  a = mxGetPr(mxA);

  for (j = 1 ; j <= N ; ++j)
  {
    for (i = 1 ; i <= M ; ++i)
    {
      b[i][j] = *a++;
    }
  }

}

int getDimensions(
const mxArray		*a,
int			**dim
)
{
  *dim = (int *)mxGetDimensions(a);
  return (mxGetNumberOfDimensions(a));
}
  
void mxResize(
mxArray		*a,
int		rows,
int		cols
)
{
  int		dims[2];
  double	*pr;

  dims[0] = rows;
  dims[1] = cols;

  if (rows *  cols > 0)
  {
    mxSetDimensions(a, dims, 2); 
    pr = mxGetPr(a);

    if (a)
    {
      mxFree(pr);
    }

    pr = mxCalloc(rows * cols, sizeof(double));

    mxSetPr(a, pr); 
        
  }
}

/*******************************
 * disp
 */

int disp(const mxArray *A)
{
  int r,c,i,j;
  double *pA;

  r=mxGetM(A);
  c=mxGetN(A);
  pA=mxGetPr(A);

  for (i=0;i<r;i++)
    {
      for (j=0;j<c;j++)
	mexPrintf("\t%f",pA[i+r*j]);
      mexPrintf("\n");
    }
  mexPrintf("\n");
  return 1;
}

/* matrix trace 
 *
 */
double trace (const mxArray *A)
{
  int i, m, n;
  double *prA;
  double res=0.0;

  prA = mxGetPr(A);
  m=mxGetM(A);  n=mxGetN(A);
  if (n != m)
    {
      mexErrMsgTxt("Error in trace: Matrix must be square\n");
      return 0;
    }
  for (i=0; i<n;i++)
    res += prA[i*n+i];
  return res;
}

/* symm
 *
 * Symmetrize matrix A=(A+A')/2
 * modifies A and returns a pointer to it
 */
mxArray *symm (mxArray *A)
{
  int r,c,i,j;
  double *pA;

  r=mxGetM(A);
  c=mxGetN(A);
  pA=mxGetPr(A);
  
  if (r !=c)
    mexErrMsgTxt("symm: Matrix must be square");
  
  for (i=0;i<r;i++)
    for (j=0;j<c;j++)
      pA[j*r+i] = 0.5*(pA[j*r + i] + pA[i*r + j]);

  return A;
}
