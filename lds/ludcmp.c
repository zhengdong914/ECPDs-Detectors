#include <math.h>
#include <stdlib.h>

#define TINY 1.0e-20;
extern double *vector();
int ludcmp(a,n,indx,d)
int n,*indx;
double **a,*d;
{
	int i,imax=-1,j,k;
	double big,dum,sum,temp;
	static double *vv = (double *)NULL;
	static int size_vv = 0;
	void nrerror(),free_vector();

	if (n > size_vv)
	{
	    if (vv)
	    {
		free_vector(vv,1,n);
	    }
		
	    vv=vector(1,n);
	    size_vv = n;
	}
	
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) {
		  mexErrMsgTxt("Singular matrix in routine LUDCMP\n"); 
		}
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	/* free_vector(vv,1,n); */
	return (1);
}

#undef TINY
