/*************************************************************************
 *                   libnma's HEADER: nma.h                              *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   libnma's main header.                                               *
 *   (defines, common definitions, data-types, data-structures, externs) *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#ifndef NMA_H_
#define NMA_H_

#include <math.h>
#include <iostream>
//#include <stdlib.h>
//#include "libpdb/ResIni.h"
//#include <libpdb/world.h>
#include <libpdb/include/Macromolecule.h>
#include <libpdb/include/pdbIter.h>
#include <libtools/include/timer.h>
#include <libtools/include/Mersenne.h> // Mersenne Twister Random Number generator

#define FILE_NAME 300

#ifndef ZERO
#define ZERO 0.0
#endif

#ifndef LINE_LENGTH
#define LINE_LENGTH 256
#endif

typedef struct
{
  int nat;
  /* # of pseudo-atoms */
  int nan;
  /* # of dihedral angles */
  int k1;
  /* index of first atom */
} tri;

// New "twid"
typedef struct
{
	int k; // atom indices
	int l;
	int i; // residue indices
	int j;
	double d;
// Mon (4/2/2011) "Not needed anymore!"
//	double Eel;
//	double Ivw;
//	double Ihb;
	double C;
}
twid;

// "derivatives" structure (K-matrix)
typedef struct
{
  double x;
  double y;
  double z;
}
  trd;

typedef struct // needed by nmaview (single precision version of "trd")
{
    float x;
    float y;
    float z;
} trs;

typedef struct // needed by read_TSfunc()
{
	char i;
	char j;
	int t; // topology
	float a;
	float b;
	float c;
} TSfunc;

typedef struct // needed by move_vwMFAx() (parallel)
{
	int model; // atomic model: CA/CA3/C5/HA
	pdbIter *iter;
	double step;
	double *evec;
	int index;
	int size;
	double ***V;
	double ***W;
	int **body1;
	int first;
	int last;
} moveVW_data;

typedef struct // needed by hessianMCAxHD_par() (parallel)
{
	// 0. Get Input: "a0","b0","boxsize","asipas", "isready", "Ut", "Ur", "hess_matrix", "undh", "size".
	long int k; // Top-right row box.
	long int l; // Top-right column box.
	int icbox; // Boxsize in internal coordinates.
	twid ***asipas; // Array with pointers to an array of "sorted ipas" belonging only to a single Hessian box.
	int *nasipas; // number of sipas whithin each box
	int nipa; // number of ipas
	bool *isready; // This bool array tells if box(k,l) is ready to be computed.
	bool *isdone; // This bool array tells if box(k,l) computation is over.
	double ****Ut,****Ur; // Ut = U-top, Ur = U-right. Both matrices will store corresponding U interaction elements fore every Hessian-Box.
	double *hess_matrix; // The Hessian matrix.
	int *undh; // This returns the first unit-index on the left side of the dihedral.
	double **erx; // "erx" array
	long int size; // Total number of variables.
	int *jobs_done; // The number of jobs done.
	float *coord; // Current atomic model array of cartesian coordinates.
	pthread_mutex_t *mtx_isready; // Mutex to update "isready" matrix.
	pthread_cond_t *cond_isready; // Condition to check "isready".
//	pthread_mutex_t *mtx_control; // Mutex for control.
//	pthread_cond_t *cond_control; // Condition for control.
} hessianCApar_data;

typedef struct // needed by di2cartVW()
{
	double *CCevec;
	double *ICevec;
	float *coord;
	int nmodes;
	int natoms;
	int model;
	int size;
	double ***V;
	double ***W;
	int **body1;
	int first;
	int last;
} di2cartVW_data;

// USED in "nmamon" and "libnma_def.cpp" (18/3/2010)
// This allow to perform diagonalization using either single or double precision
// After changing floating point number definition you have to change the
// ssyevr_(LAPACK)[single] or dsyevr_(LAPACK)[double] functions inside nma.cpp.
typedef float myfloat;

// This allow to perform diagonalization using either single or double precision
// After changing floating point number definition you have to change the
// ssyevr_(LAPACK)[single] or dsyevr_(LAPACK)[double] functions inside nma.cpp.
//typedef float floating;
typedef double floating;


// EXTERN FUNCTIONS DECLARATION:
extern "C"
{
	void dsyevr_(char *jobz, char * range, char *uplo, int *n,
			double *a, int *lda, double *vl, double *vu, int *il,
			int *iu, double *abstol, int* m, double *w,
			double *z__, int *ldz, int *isuppz, double *work,
			int *lwork, int *iwork, int *liwork, int *info);
	// Use the following function instead to force double precision diagonalization
	// int dsyevr_(char *jobz, char * range, char *uplo, int *n,
	// Use the following function instead to force single precision diagonalization
	//int ssyevr_(char *jobz, char * range, char *uplo, int *n,
	//                myfloat *a, int *lda, myfloat *vl, myfloat *vu, int *il,
	//                int *iu, myfloat *abstol, int* m, myfloat *w,
	//                myfloat *z__, int *ldz, int *isuppz, myfloat *work,
	//                int *lwork, int *iwork, int *liwork, int *info);

	// (D,S)LAMCH - determine (double,single) precision machine parameters
	double dlamch_(char *cmach);
	double slamch_(char *cmach);

	void    dsygvd_(int *itype,char *jobz,char *uplo,int *n,double *a,int *lda,double *b,
			int *ldb,double *w,double *work,int *lwork,int *iwork,int *liwork,int *info);

	//       SUBROUTINE DSPGVD( ITYPE,  JOBZ,  UPLO,  N,  AP,  BP,  W, Z, LDZ, WORK,
	//                          LWORK, IWORK, LIWORK, INFO )
	//           CHARACTER      JOBZ, UPLO
	//           INTEGER        INFO, ITYPE, LDZ, LIWORK, LWORK, N
	//           INTEGER        IWORK( * )
	//           DOUBLE         PRECISION AP( * ), BP( * ), W( * ), WORK(  *  ),  Z(
	//                          LDZ, * )
	void    dspgvd_(int *itype,char *jobz,char *uplo,int *n,double *a,double *b,double *w,
			double *z,int *ldz,double *work,int *lwork,int *iwork,int *liwork,int *info);

	//SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
	//$                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
	//$                   LWORK, IWORK, IFAIL, INFO )

	void    dsygvx_(int *itype,char *jobz,char *range,char *uplo,int *n,double *a,int *lda,double *b,int *ldb,
			int *vl,int *vl2,int *il,int *iu,double *abstol,int *m,double *w,double *z,int *ldz,double *work,
			int *lwork,int *iwork,int *liwork,int *info);

	// SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, INFO )
	//           CHARACTER     JOBZ, UPLO
	//           INTEGER       INFO, ITYPE, LDZ, N
	//           DOUBLE        PRECISION AP( * ), BP( * ), W( * ),  WORK(  *  ),  Z(
	//                         LDZ, * )
	void    dspgv_(int *itype,char *jobz,char *uplo,int *n,double *a,double *b,double *w,
			double *z,int *ldz,double *work,int *info);

	//SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
	//$                  INFO )
	//CHARACTER          JOBZ, UPLO
	//INTEGER            INFO, ITYPE, LDZ, N
	//REAL               AP( * ), BP( * ), W( * ), WORK( * ),
	//$                   Z( LDZ, * )
	void sspgv_(int *itype,char *jobz,char *uplo,int *n,float *a,float *b,float *w,
			float *z,int *ldz,float *work,int *info);

	//SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
	//$                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,
	//$                   INFO )
	//CHARACTER          JOBZ, RANGE, UPLO
	//INTEGER            IL, INFO, IU, LDZ, M, N
	//DOUBLE PRECISION   ABSTOL, VL, VU
	//INTEGER            IFAIL( * ), IWORK( * )
	//DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
	void dspevx_( char *jobz, char *range, char *uplo, int *size, double *AP, double *nulo, double *nulo2,
			int *il, int *iu, double *toler, int *eme, double *eigval, double *eigvect, int *size2,
			double *work, int *iwork, int *ifail, int *info );

	//SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
	//$                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,
	//$                   INFO )
	//CHARACTER          JOBZ, RANGE, UPLO
	//INTEGER            IL, INFO, IU, LDZ, M, N
	//REAL               ABSTOL, VL, VU
	//INTEGER            IFAIL( * ), IWORK( * )
	//REAL               AP( * ), W( * ), WORK( * ), Z( LDZ, * )
	void sspevx_( char *jobz, char *range, char *uplo, int *size, float *AP, float *nulo, float *nulo2,
			int *il, int *iu, float *toler, int *eme, float *eigval, float *eigvect, int *size2,
			float *work, int *iwork, int *ifail, int *info );

	//SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
	// CHARACTER          UPLO
	// INTEGER            INFO, N
	// INTEGER            IPIV( * )
	// DOUBLE PRECISION   AP( * ), WORK( * )
	void dsptrf_(char *uplo, int *size, double *ap, int *ipiv, int *info); // computes IPIV

	//SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
	// CHARACTER          UPLO
	// INTEGER            INFO, N
	// INTEGER            IPIV( * )
	// DOUBLE PRECISION   AP( * ), WORK( * )
	void dsptri_(char *uplo, int *size, double *ap, int *ipiv, double *work, int *info); // computes Inverse-Matrix

	//SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO )
	//LAPACK routine (version 3.2) --
	//CHARACTER          DIAG, UPLO
	//INTEGER            INFO, N
	//DOUBLE PRECISION   AP( * )
	void dtptri_(char *diag, char *uplo, int *size, double *ap, int *info); // computes Inverse-Matrix

	//SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
	//*  -- LAPACK routine (version 3.0) --
	//CHARACTER          UPLO
	//INTEGER            INFO, ITYPE, N
	//DOUBLE PRECISION   AP( * ), BP( * )
	void dspgst_(int *itype, char *uplo, int *size, double *ap, double *bp, int *info);

	//SUBROUTINE DSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,
	//$                   IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,
	//$                   IFAIL, INFO )
	//*  -- LAPACK driver routine (version 3.1) --
	//CHARACTER          JOBZ, RANGE, UPLO
	//INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N
	//DOUBLE PRECISION   ABSTOL, VL, VU
	//INTEGER            IFAIL( * ), IWORK( * )
	//DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ),
	//$                   Z( LDZ, * )
	void dspgvx_( int *itype, char *jobz, char *range, char *uplo, int *size,
			double *AP, double *BP, double *nulo, double *nulo2, int *il, int *iu,
			double *tol, int *eme, double *eigval, double *eigvect, int *ldz,
			double *work, int *iwork, int *ifail, int *info );
	void sspgvx_( int *itype, char *jobz, char *range, char *uplo, int *size,
			float *AP, float *BP, float *nulo, float *nulo2, int *il, int *iu,
			float *tol, int *eme, float *eigval, float *eigvect, int *ldz,
			float *work, int *iwork, int *ifail, int *info );

	// ARPACK needed functions, from J.I.Aliaga (17/5/2013)
	int dpptrf_(char *uplo, int *n, double *ap, int *info);
	int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
	void dsaupd_(int *ido, char *bmat, int *n, char *which,
				int *nev, double *tol, double *resid, int *ncv,
				double *v, int *ldv, int *iparam, int *ipntr,
				double *workd, double *workl, int *lworkl,
				int *info);
	int dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);
	int dtptrs_(char *uplo, char *trans, char *diag, int *n,
		int *nrhs, double *ap, double *b, int *ldb, int *info);
	int dtrtrs_(char *uplo, char *trans, char *diag, int *n,
			int *nrhs, double *a, int *lda, double *b,
			int *ldb, int *info );
	int dspmv_(char *uplo, int *n, double *alpha,
		double *ap, double *x, int *incx, double *beta,
		double *y, int *incy);
	int dgemv_(char *trans, int *m, int *n, double *
		alpha, double *a, int *lda, double *x, int *incx,
		double *beta, double *y, int *incy);
	void dseupd_(int *rvec, char *All, int *select, double *d,
				double *z, int *ldz, double *sigma,
				char *bmat, int *n, char *which, int *nev,
				double *tol, double *resid, int *ncv, double *v,
				int *ldv, int *iparam, int *ipntr, double *workd,
				double *workl, int *lworkl, int *ierr);


	// MRRR from P.Bientinesi (14/12/2010)
	//######################################################################
	//Symmetric-definite case (packed storage):                            #
	//######################################################################
	 int dspgeig(int *itype, char *jobz, char *range, char *uplo,
	             int *n, double *AP, double *BP, double *vl,
	             double *vu, int *il, int *iu, int *m, double *W,
	             double *Z, int *ldz);
}

/*========================================================================================*/


#endif /*NMA_H_*/
