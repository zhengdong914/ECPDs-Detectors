/**************************************************************
 * MODULE DESCRIPTION
 *
 *   Kalman Smoother
 *   (Z Ghahramani)
 *
 *************************************************************/
#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include "matrix.h"
#include "mxutils.h"

#define zeros(X,Y) mxCreateDoubleMatrix(X,Y,mxREAL);
#define sendup(X) mxSetName(X,#X); mexPutArray(X,"base")

double sumall (const mxArray *M, const mxArray *P);
/* Kalman smoother */
/* returns: [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3] */

/******************************************************/
void kalmansmooth( 
const mxArray *A, 
const mxArray *C, 
const mxArray *Q, 
const mxArray *R, 
const mxArray *x0, 
const mxArray *P0, 
const mxArray *Y, 
mxArray	*lik, 
mxArray *Xfin, 
mxArray *Pfin,
mxArray *Ptsum,
mxArray *YX,
mxArray *A1,
mxArray *A2,
mxArray *A3
)
{

  double tiny=exp(-700);
  int K,p,T,N;
  int Xdims[3], Pdims[3];
  const int *dims;
  mxArray *Xpre, *Xcur;
  mxArray *Ppre, *Pcur, *Pt, *Pcov, *Pc, *J;
  mxArray *Kcur, *invP, *invR;
  mxArray *temp1, *temp2, *temp3, *temp3b, *temp4, *CP, *KC;
  mxArray *Ydiff, *Xfint, *I, *detb;
  int t,n,k,i,j;
  double detiP;
  double knorm;
  double *plik, *pD;

  K=mxGetM(A);
  p=mxGetM(C); /* cols */
  dims=mxGetDimensions(Y); /* Y should be N x p x T */
  T=dims[2];
  N=dims[0];

  if (mxGetN(A) != K || dims[1] != p || mxGetN(C) != K)
    {
      mexPrintf("Input dimensions are inconsistent!\n");
      return;
    }
  
  knorm=pow(2.0*M_PI,-((double) p)/2.0);
  
  Xdims[0]=N;  Xdims[1]=K;  Xdims[2]=T;
  Pdims[0]=K;  Pdims[1]=K;  Pdims[2]=T;

  Xpre = zeros(N,K); 
  Xcur = mxCreateNumericArray(3,Xdims,mxDOUBLE_CLASS,mxREAL);
  Xfint=zeros(N,K);
  Ppre = mxCreateNumericArray(3,Pdims,mxDOUBLE_CLASS,mxREAL);
  Pcur = mxCreateNumericArray(3,Pdims,mxDOUBLE_CLASS,mxREAL);
  J = mxCreateNumericArray(3,Pdims,mxDOUBLE_CLASS,mxREAL);
  Pt = zeros(K,K); 
  Pcov = zeros(K,K); 
  Kcur = zeros(K,p); 
  invP = zeros(p,p);
  /* temp1 = zeros(p,K);   temp2 = zeros(p,K);   temp3 = zeros(K,K); */
  /* temp4 = zeros(K,p);  invP */
  temp3b = zeros(K,K);
  detb = zeros(1,1);
  CP = zeros(K,p); 
  KC = zeros(K,K); 
  Ydiff = zeros(N,p); 
  I = eye(K);

  /* FORWARD PASS */

  plik=mxGetPr(lik);
  *plik=0.0;
  t=1; 

  /* index is (t-1)*N*K+k*N+n for NxKxT matrix*/
  
  copy(Xpre,mult(ones(N,1),tr(x0)));
  
  subcopy(Ppre,1,P0);

  invR=rdiv(eye(p),R); /* diag(1./(R+(R==0)*tiny)); */

  /* BEGIN FORWARD */

  for (t=1;t<=T;t++)
    {
      if (K<p)
	{
	  temp1=rdiv(C,R);
	  temp2=mult(temp1,subm(Ppre,t));
	  temp3=mult(tr(C),temp2);
	  inv(plus(I,temp3),temp3b,detb); 
	  temp4=mult(temp3b,tr(temp1));
	  invP=minus(invR,mult(temp2,temp4));
	  CP= dec(tr(temp1),mult(temp3,temp4));
	  detiP=sqrt(det(invP));
	}
      else
	{
	  temp1=inc(diag(R),mult(C,mult(subm(Ppre,t),tr(C))));
	  invP=zeros(p,p);
	  inv(temp1,invP,detb);
	  CP=mult(tr(C),invP);
	  pD=mxGetPr(detb);
	  detiP=1.0/sqrt(pD[0]);
	}
      copy(Kcur, mult(subm(Ppre,t),CP));
      copy(KC,mult(Kcur,C));
      copy(Ydiff,dec(subm(Y,t),mult(Xpre,tr(C))));

      subcopy(Xcur, t, Xpre);
      subinc(Xcur, t, mult(Ydiff,tr(Kcur)));
      subcopy(Pcur, t, subm(Ppre,t));
      subdec(Pcur, t, mult(KC,subm(Ppre,t)));

      if (t<T)
	{
	  copy(Xpre, mult(subm(Xcur,t),tr(A)));
	  subcopy(Ppre, t+1, mult(A,mult(subm(Pcur,t),tr(A))));
	  subinc(Ppre, t+1, Q);
	}
      /*   CALCULATE LIKELIHOOD */

      /* detiP = sqrt(det(invP)); double */
      *plik += N*log(detiP) - 0.5*sumall(Ydiff,invP);
      
    }
  *plik += N*T*log(knorm);

  /* BACKWARD PASS */
  /* mexPrintf("backward pass!\n"); */
  
  t=T;
  subcopy(Xfin, t, subm(Xcur, t));
  subcopy(Pfin, t, subm(Pcur, t));
  copy(Pt, subm(Pfin, t));
  copy(Xfint,subm(Xfin,t));
  inc(Pt, scmult(mult(tr(Xfint),Xfint),(double)(1.0/ ((double)N))));
  copy(A2,Pt);
  scmult(A2,(double) -1.0);
  copy(Ptsum,Pt);
  
  copy(YX,mult(tr(subm(Y,t)),Xfint));
  temp3=zeros(K,K);

  for (t=T-1;t>0;t--)
    {
      inv(subm(Ppre,t+1),temp3,detb);
      subcopy(J, t, mult(mult(subm(Pcur,t),tr(A)),temp3));
      subcopy(Xfin, t, subm(Xcur,t));
      subinc(Xfin, t, mult(dec(subm(Xfin,t+1),mult(subm(Xcur,t),tr(A))),
			   tr(subm(J,t))));
      copy(Xfint,subm(Xfin,t));
      subcopy(Pfin, t, subm(Pcur,t));
      subinc(Pfin, t, mult( subm(J,t), 
			    mult(dec(subm(Pfin,t+1),subm(Ppre,t+1)),
				 tr(subm(J,t))))); 
      copy(Pt, inc(subm(Pfin,t), scmult( mult(tr(Xfint),Xfint), (double) (1.0/N) )));
      inc(Ptsum,Pt);
      inc(YX,mult(tr(subm(Y,t)),Xfint));
    }

  /* DEBUGGING */
  /*    sendup(Xfin); sendup(Pfin); sendup(YX); sendup(Pt); sendup(Ptsum); return; */
  /* END DEBUGGING */

  copy(A3, minus(Ptsum,Pt));
  inc(A2,Ptsum);
  
  t=T;
  copy(Pcov,mult(minus(I,KC),mult(A,subm(Pcur,t-1))));
  inc(A1,inc(scmult(mult(tr(subm(Xfin,t)),subm(Xfin,t-1)), (double) (1.0/(double)N)),Pcov));
  
  for (t=T-1;t>1;t--)
    {
      copy(Pcov,mult(inc(subm(Pcur,t),
			 mult(subm(J,t),minus(Pcov,mult(A,subm(Pcur,t))))),
		     tr(subm(J,t-1)))); 
      inc(A1,inc(scmult(mult(tr(subm(Xfin,t)),subm(Xfin,t-1)), (double) (1.0/(double)N)),Pcov));
    }
  
}


double sumall (const mxArray *M, const mxArray *P) /* sum of M*P*M' */
{
  int i, j, k, n, m;
  double *pM, *pP, ret=0.0;

  n=mxGetN(M);
  m=mxGetM(M);
  pM=mxGetPr(M);
  pP=mxGetPr(P);

  if (mxGetM(P) != n || mxGetN(P) != n)
    mexErrMsgTxt("SUMALL: sizes don't match!\n");

  for (i=0;i<n;i++)
    for (k=0;k<n;k++)
      for (j=0;j<m; j++)
	ret += pM[i*m+j]*pM[k*m+j]*pP[k*n+i];

  return ret;
}

/******************************************************/
/******************************************************/
/******************************************************/
void mexFunction(
int nlhs, mxArray *plhs[],
int nrhs, const mxArray *prhs[]
)
{
  const mxArray *A, *C, *Q, *R, *x0, *P0, *Y;
  mxArray	*Xfin, *Pfin, *Ptsum, *YX, *A1, *A2, *A3, *lik;

  int		i = 0;
  int           K, p, N, T;
  const int *dims;
  int Xdims[3], Pdims[3];
  double	*retX=0, *y;
  static int	expectRhs = 7;
  static int	expectLhs = 8;
  char		buffer[100];

  /* Check for proper number of arguments */
  
  if (nrhs != expectRhs) 
  {
    sprintf(buffer, "kalmansmooth requires %d input arguments. Got %d",
	    expectRhs, nrhs);
    mexErrMsgTxt(buffer);
  } else if (nlhs > expectLhs) 
  {
    sprintf(buffer, "kalmansmooth requires %d outputs arguments. Got %d",
	    expectLhs, nlhs);
    mexErrMsgTxt(buffer);
  }

  A = prhs[i++];
  C = prhs[i++];
  Q = prhs[i++];
  R = prhs[i++];
  x0 = prhs[i++];
  P0 = prhs[i++];
  Y = prhs[i++];
  
  K=mxGetM(A); /* rows */
  p=mxGetM(C); /* cols */
  dims=mxGetDimensions(Y); /* Y should be N x p x T */
  T=dims[2];
  N=dims[0];
  if (dims[1] != p)
    mexErrMsgTxt("Dimensions don't match in kalmansmooth\n");
  
  Xdims[0]=N;  Xdims[1]=K;  Xdims[2]=T;

  Xfin = mxCreateNumericArray(3,Xdims,mxDOUBLE_CLASS,mxREAL);
  
  Pdims[0]=K;  Pdims[1]=K;  Pdims[2]=T;

  Pfin = mxCreateNumericArray(3,Pdims,mxDOUBLE_CLASS,mxREAL);

  Ptsum = zeros(K,K); 
  YX = zeros(p,K);

  A1 = zeros(K,K); 
  A2 = zeros(K,K);
  A3 = zeros(K,K); 
  lik= zeros(1,1);

  kalmansmooth(A, C, Q, R, x0, P0, Y,lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3);

  i=0;
  
  plhs[i++]=lik;
  plhs[i++]=Xfin;
  plhs[i++]=Pfin;
  plhs[i++]=Ptsum;
  plhs[i++]=YX;
  plhs[i++]=A1;
  plhs[i++]=A2;
  plhs[i++]=A3;
}

  /*  bytes_to_copy = TOTAL_ELEMENTS * mxGetElementSize(array_ptr);
      memcpy(start_of_pr, real_data, bytes_to_copy); */

