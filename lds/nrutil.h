#ifndef H_NR_UTIL_H
#define H_NR_UTIL_H


#ifndef _COMPLEX_
#define _COMPLEX_
typedef struct RCOMPLEX {float r,i;} rcomplex;
typedef struct DCOMPLEX {double r,i;} dcomplex;
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern void nrerror(char error_text[]);
extern double *vector(int nl, int nh);
extern rcomplex *fcvector(int nl, int nh);
extern dcomplex *dcvector(int nl,  int nh);
extern int *ivector(int nl, int nh);
extern char *charvector(int nl, int nh);
extern double *dvector(int nl, int nh);
extern double **matrix(int nrl, int nrh, int ncl, int nch);
extern double **dmatrix(int nrl, int nrh, int ncl, int nch);
extern int **imatrix(int nrl, int nrh, int ncl, int nch);
extern char **cmatrix(int nrl, int nrh, int ncl, int nch);
extern double **submatrix(double **a, int oldrl,int oldrh,int oldcl,int oldch,int newrl,int newcl);
extern double* free_vector(double *v,int nl, int nh);
extern int* free_ivector(int *v,int nl, int nh);
extern double* free_dvector(double *v,int nl, int nh);
extern void free_fcvector(rcomplex *v,int nl, int nh);
extern void free_dcvector(dcomplex *v,int nl, int nh);
extern double** free_matrix(double **m,int nrl, int nrh, int ncl,int nch);
extern double** free_dmatrix(double **m,int nrl, int nrh, int ncl,int nch);
extern int** free_imatrix(int **m,int nrl, int nrh, int ncl,int  nch);
extern void free_submatrix(double **b,int nrl, int nrh, int ncl,int nch);
extern double **convert_matrix(double *a,int nrl, int nrh, int ncl, int nch);
extern void free_convert_matrix(double *b,int nrl, int nrh, int ncl,int nch);

#ifdef __cplusplus
  }
#endif

#endif
