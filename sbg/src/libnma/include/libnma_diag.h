#ifndef LIBNMA_DIAG_H_
#define LIBNMA_DIAG_H_

#include "nma.h"
//#include "libmrrr/mrrr.h"

//*  DSYEVR computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
//*  selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
int diag_dsyevr(int , int, double *, double **, double **);

// Mon (21/12/2010): This one does not allocate memory for eigenvectors (diag_dsygvx yes)
//                   ("eigvec" sould be a size*evec_size array of doubles)
//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx2(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu);

//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx(double *hess_matrix,double *mass_matrix,double *eigval, int size, int il, int iu);

//*  DSYGVD computes ALL the eigenvalues, and optionally, the eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
//*  B are assumed to be symmetric and B is also positive definite.
//*  If eigenvectors are desired, it uses a divide and conquer algorithm.
void diag_dsygvd(double *hess_matrix,double *mass_matrix,double *eigval, int size);



//*  DSPEVX computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//*  can be selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//   (SMALL MEMORY REQUIREMENTS) (default compute eigenvalues and eigenvectors)
//void diag_dspevx(double *hess_matrix,double **ppeigval, double **ppeigvec, int size, int il, int iu, char jobz='V');
void diag_dspevx(double *hess_matrix,double **ppeigval, double **ppeigvec, int size, int il, int iu, double vl=0.0, double vu=0.0, char jobz='V', char range='I',int *pemeF=NULL);

//*  SSPEVX computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//*  can be selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//   (SMALL MEMORY REQUIREMENTS) (default compute eigenvalues and eigenvectors)
void diag_sspevx(float *hess_matrix,float **ppeigval, float **ppeigvec, int size, int iu, char jobz='V');


//*  DSPGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
//   (SMALL MEMORY REQUIREMENTS) (by default compute eigenvalues and eigenvectors)
int diag_dspgvx(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu);

//// MRRR from P.Bientinesi (14/12/2010)
////######################################################################
////Symmetric-definite case (packed storage):                            #
////######################################################################
//int diag_dspgeig(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu);


//*  SSPGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
//   (SMALL MEMORY REQUIREMENTS)
void diag_sspgvx(float *hess_matrix,float *mass_matrix,float *eigval, float *eigvec, int size, int iu);

//   XSPGVX --> valid for single/double precision floating point data
//   WARNING: floating precision detected by "sizeof()" !!!
//*  DSPGVX and SSPGVX compute SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
void diag_xspgvx(floating *hess_matrix,floating *mass_matrix,floating *eigval, floating *eigvec, int size, int iu);

// ARPACK-based diagonalization routine (from: J.I.Aliaga & E.Quintana)
// Very fast if small number of eigenpairs requested. Not optimized for packed storage, thus it does not scale linear with number of threads.
// A --> Hessian matrix
// B --> Kinetic energy matrix
// d --> Eigenvalues array
// z --> Eigenvectors array
// nrows --> Problem dimension (size)
// nev --> Number of eigenpairs to be computed
// ncv --> Factor for the Number of columns of  matrix V (=ncv*nev). It should be optimized for a given problem (OPTIONAL)
void dsdrv1_AP_BP_W_mon(double *B, double *A, double *d, double *z, int nrows, int nev, float ncvf = 5);

// ARPACK-based diagonalization routine (from: J.I.Aliaga & E.Quintana)
// Very fast if small number of eigenpairs requested. It scales linear with number of threads.
// Although matrices should be provided in packed storage, optimized for packed storage, thus it does not
// A --> Hessian matrix
// B --> Kinetic energy matrix
// d --> Eigenvalues array
// z --> Eigenvectors array
// nrows --> Problem dimension (size)
// nev --> Number of eigenpairs to be computed
// ncv --> Factor for the Number of columns of matrix V (=ncv*nev). It should be optimized for a given problem (OPTIONAL)
void dsdrv1_A_B_W_mon(double *B, double *A, double *d, double *z, int nrows, int nev, float ncvf = 5);

// DSPGVD  computes  ALL the eigenvalues, and optionally, the eigenvectors
// of a real generalized  symmetric-definite  eigenproblem,  of  the  form
// A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
// are assumed to be symmetric, stored in PACKED FORMAT,  and  B  is  also
// positive definite.
// If eigenvectors are desired, it uses a divide and conquer algorithm.
// (BIG MEMORY REQUIREMENTS: If JOBZ = 'V' & N > 1, LWORK >= 1 + 6*N + 2*N**2.)
void diag_dspgvd(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size);

// DSPGV computes ALL the eigenvalues and, optionally, the eigenvectors of
// a  real  generalized  symmetric-definite  eigenproblem,  of  the   form
// A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
// are assumed to be symmetric, stored in packed format,  and  B  is  also
// positive definite.
// (MINIMAL MEMORY REQUIREMENTS)
void diag_dspgv(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size);

//*  SSPGV computes ALL the eigenvalues and, optionally, the eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
//*  Here A and B are assumed to be symmetric, stored in packed format,
//*  and B is also positive definite.
//   (MINIMAL MEMORY REQUIREMENTS)
void diag_sspgv(float *hess_matrix,float *mass_matrix,float *eigval,float *eigvec, int size);




// INVERSE MATRIX COMPUTATION sub-rutine (Double precision, double)
// Uses LAPACK's dsptri_. (SYMMETRIC, PACKED STORAGE, MINIMAL MEMORY USAGE)
// Uses LAPACK's dsptrf_. (Triangular factorization: LU, Cholesky)
//*  DSPTRI computes the inverse of a real symmetric indefinite matrix
//*  A in packed storage using the factorization A = U*D*U**T or
//*  A = L*D*L**T computed by DSPTRF.
void inv_dsptri(double *ap, int size);

// INVERSE MATRIX COMPUTATION sub-rutine (Double precision, double)
// Uses LAPACK's dtptri_. (TRIANGULAR, PACKED STORAGE, MINIMAL MEMORY USAGE)
// DSPTRI computes the inverse of a real symmetric indefinite matrix
void inv_dtptri(double *ap, int size);

//*  Reduces a real symmetric-definite generalized eigenproblem
//*  to standard form, using packed storage.
// Uses LAPACK's dspgst_. (SYMMETRIC, PACKED STORAGE, MINIMAL MEMORY USAGE)
// Uses LAPACK's dsptrf_. (Triangular factorization: LU, Cholesky)
int reduce_dspgst(double *ap, double *bp, int size);


// *******************************************************************************
// PABLO's OLD DIAGONALIZATION-FUNCTION (DEFPROT-back compatibility)
// *******************************************************************************
// lpack diagonalization
int nma (int , int, float *, float **, float **);

#endif /*LIBNMA_DIAG_H_*/
