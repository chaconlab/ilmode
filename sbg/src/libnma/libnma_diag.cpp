/************************************************************************
 *                     LIBRARY: libnma_deriv                             *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Library with "in-house" functions for calling LAPACK's              *
 *   diagonalization rutines, and related ones.                          *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

//#include <nma.h>
#include <libnma_diag.h>
#include <libnma_io.h>

// LAPACK's Symmetric Eigenproblems (SEP)
// http://www.netlib.org/lapack/lug/node30.html


// CHECK THIS!!!
// Uses ARPACK's dsaupd_.???? It's faster???
// DSAUPD has been designed to compute approximations to a
// few eigenpairs of a linear operator OP that is real and symmetric
// with respect to a real positive semi-definite symmetric matrix B:
//*	Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite


//*  DSYEVR computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
//*  selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
int diag_dsyevr(int size, int neigval, double *hessian, double ** eigvalp, double ** eigvectp)
{
	bool debug=true;
	int i;
	double * eigval, * eigvect;
	double toler_zero_mode;

	toler_zero_mode=0.0001;

	/* Lapack input/output */
	char jobz, range, uplo;
	int lda, il, iu, ldz, * isuppz, lwork, liwork, info, * iwork;
	double vl, vu, toler, * work;

	/* allocate memory */
	eigval = ( double * ) malloc( (size+1) * sizeof( double ) );
	for ( i = 0; i < (size+1); i++ ) eigval[i] = 0.0;
	eigvect = ( double * ) malloc( (size+1) * neigval * sizeof( double ) );
	for ( i = 0; i < (size+1) * neigval; i++ ) eigvect[i] = 0.0;

	/**  */
	/* diagonalize matrix here */
	/* USING SSYEVR/DSYEVR FORTRAN SUBRUTINE FROM LAPACK */
	/**  */

	jobz = 'V';	// 'V':  Compute eigenvalues and eigenvectors
	// 'N':  Compute eigenvalues only;
	if ( neigval == size ) range = 'V'; // 'A': all eigenvalues will be found.
	// 'V': all eigenvalues in the half-open interval (VL,VU] will be found.
	else range = 'I'; // 'I': the IL-th through IU-th eigenvalues will be found

	// If  RANGE='I', the indices (in float ascending order) of the smallest and largest eigenvalues to be returned.
	//  1 <= IL <= IU <= N, if N > 0;
	//	IL = 1 and IU = 0 if N = 0.
	//	Not referenced if RANGE = 'A' or 'V'
	if ( range == 'I' )
	{
		il = 1;
		iu = neigval;
	}
	uplo = 'U';	// 'U':  Upper triangle of A is stored;
				// 'L':  Lower triangle of A is stored.

	if(range=='V')
	{
		// vl = -0.0001;
		vl = -9E29;
		vu = 9E29;	// If RANGE='V', the lower and upper bounds of the  interval  to  be  searched  for eigenvalues.
		// VL < VU.  Not referenced if RANGE = 'A' or 'I'.
	}

	toler = 0.0;

	lda = size;
	ldz = size;

	if ( range == 'I' )
	{
		isuppz = ( int * ) malloc( 2 * neigval * sizeof( int ) );
		for ( i = 0; i < 2 * neigval; i++ ) isuppz[i] = 0;
	}
	else
	{
		isuppz = ( int * ) malloc( 2 * size * sizeof( int ) );
		for ( i = 0; i < 2 * size; i++ ) isuppz[i] = 0;
	}

	// Constant workspace!!!!
	lwork = size * 26;
	work = ( double * ) malloc( lwork * sizeof( double ) );
	for ( i = 0; i < lwork; i++ ) work[i] = 0.0;

	liwork = size * 10;
	iwork = ( int * ) malloc( liwork * sizeof( int ) );
	for ( i = 0; i < liwork; i++ ) iwork[i] = 0;

	if(debug) printf("Msg(diag_dsyevr): Diagonalization:  size=%d, neigval=%d\n",size,neigval);

	dsyevr_( & jobz, & range, & uplo, & size, hessian, & lda,
			& vl, & vu, & il, & iu, & toler,
			& neigval, eigval, eigvect,
			& ldz, isuppz, work, & lwork, iwork, & liwork, & info );

	double dump;
	for ( i = 6; i < neigval - 1; i++ )
	{
		if(eigval[i]>1e-11)
			break;
	}

	if ( info )
	{
		printf( "nmac> An error occured in the matrix diagonalization: %d\n", info );
		exit( 1 );
	}
	else
		printf( "nmac> %d eigvalues found (first %d null)\n", neigval, i );

	free(isuppz);
	free(work);
	free(iwork);

	*eigvalp=eigval;
	*eigvectp=eigvect;

	return i; // Number of Null modes.
}
//SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
//$                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
//$                   IWORK, LIWORK, INFO )
//*
//*  -- LAPACK driver routine (version 3.2) --
//*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
//*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
//DOUBLE PRECISION   ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            ISUPPZ( * ), IWORK( * )
//DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
//*  selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//*
//*  DSYEVR first reduces the matrix A to tridiagonal form T with a call
//*  to DSYTRD.  Then, whenever possible, DSYEVR calls DSTEMR to compute
//*  the eigenspectrum using Relatively Robust Representations.  DSTEMR
//*  computes eigenvalues by the dqds algorithm, while orthogonal
//*  eigenvectors are computed from various "good" L D L^T representations
//*  (also known as Relatively Robust Representations). Gram-Schmidt
//*  orthogonalization is avoided as far as possible. More specifically,
//*  the various steps of the algorithm are as follows.
//*
//*  For each unreduced block (submatrix) of T,
//*     (a) Compute T - sigma I  = L D L^T, so that L and D
//*         define all the wanted eigenvalues to high relative accuracy.
//*         This means that small relative changes in the entries of D and L
//*         cause only small relative changes in the eigenvalues and
//*         eigenvectors. The standard (unfactored) representation of the
//*         tridiagonal matrix T does not have this property in general.
//*     (b) Compute the eigenvalues to suitable accuracy.
//*         If the eigenvectors are desired, the algorithm attains full
//*         accuracy of the computed eigenvalues only right before
//*         the corresponding vectors have to be computed, see steps c) and d).
//*     (c) For each cluster of close eigenvalues, select a new
//*         shift close to the cluster, find a new factorization, and refine
//*         the shifted eigenvalues to suitable accuracy.
//*     (d) For each eigenvalue with a large enough relative separation compute
//*         the corresponding eigenvector by forming a rank revealing twisted
//*         factorization. Go back to (c) for any clusters that remain.
//*
//*  The desired accuracy of the output can be specified by the input
//*  parameter ABSTOL.
//*
//*  For more details, see DSTEMR's documentation and:
//*  - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
//*    to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
//*    Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
//*  - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
//*    Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
//*    2004.  Also LAPACK Working Note 154.
//*  - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
//*    tridiagonal eigenvalue/eigenvector problem",
//*    Computer Science Division Technical Report No. UCB/CSD-97-971,
//*    UC Berkeley, May 1997.
//*
//*
//*  Note 1 : DSYEVR calls DSTEMR when the full spectrum is requested
//*  on machines which conform to the ieee-754 floating point standard.
//*  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
//*  when partial spectrum requests are made.
//*
//*  Normal execution of DSTEMR may create NaNs and infinities and
//*  hence may abort due to a floating point exception in environments
//*  which do not handle NaNs and infinities in the ieee standard default
//*  manner.
//*
//*  Arguments
//*  =========
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found.
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found.
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
//********** DSTEIN are called
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A is stored;
//*          = 'L':  Lower triangle of A is stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix A.  N >= 0.
//*
//*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
//*          On entry, the symmetric matrix A.  If UPLO = 'U', the
//*          leading N-by-N upper triangular part of A contains the
//*          upper triangular part of the matrix A.  If UPLO = 'L',
//*          the leading N-by-N lower triangular part of A contains
//*          the lower triangular part of the matrix A.
//*          On exit, the lower triangle (if UPLO='L') or the upper
//*          triangle (if UPLO='U') of A, including the diagonal, is
//*          destroyed.
//*
//*  LDA     (input) INTEGER
//*          The leading dimension of the array A.  LDA >= max(1,N).
//*
//*  VL      (input) DOUBLE PRECISION
//*  VU      (input) DOUBLE PRECISION
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) DOUBLE PRECISION
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing A to tridiagonal form.
//*
//*          See "Computing Small Singular Values of Bidiagonal Matrices
//*          with Guaranteed High Relative Accuracy," by Demmel and
//*          Kahan, LAPACK Working Note #3.
//*
//*          If high relative accuracy is important, set ABSTOL to
//*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
//*          eigenvalues are computed to high relative accuracy when
//*          possible in future releases.  The current code does not
//*          make any guarantees about high relative accuracy, but
//*          future releases will. See J. Barlow and J. Demmel,
//*          "Computing Accurate Eigensystems of Scaled Diagonally
//*          Dominant Matrices", LAPACK Working Note #7, for a discussion
//*          of which matrices define their eigenvalues to high relative
//*          accuracy.
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) DOUBLE PRECISION array, dimension (N)
//*          The first M elements contain the selected eigenvalues in
//*          ascending order.
//*
//*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          If JOBZ = 'N', then Z is not referenced.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*          Supplying N columns is always safe.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
//*          The support of the eigenvectors in Z, i.e., the indices
//*          indicating the nonzero elements in Z. The i-th eigenvector
//*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
//*          ISUPPZ( 2*i ).
//********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
//*
//*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
//*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
//*
//*  LWORK   (input) INTEGER
//*          The dimension of the array WORK.  LWORK >= max(1,26*N).
//*          For optimal efficiency, LWORK >= (NB+6)*N,
//*          where NB is the max of the blocksize for DSYTRD and DORMTR
//*          returned by ILAENV.
//*
//*          If LWORK = -1, then a workspace query is assumed; the routine
//*          only calculates the optimal size of the WORK array, returns
//*          this value as the first entry of the WORK array, and no error
//*          message related to LWORK is issued by XERBLA.
//*
//*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
//*          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
//*
//*  LIWORK  (input) INTEGER
//*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
//*
//*          If LIWORK = -1, then a workspace query is assumed; the
//*          routine only calculates the optimal size of the IWORK array,
//*          returns this value as the first entry of the IWORK array, and
//*          no error message related to LIWORK is issued by XERBLA.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  Internal error
//*
//*  Further Details
//*  ===============
//*
//*  Based on contributions by
//*     Inderjit Dhillon, IBM Almaden, USA
//*     Osni Marques, LBNL/NERSC, USA
//*     Ken Stanley, Computer Science Division, University of
//*       California at Berkeley, USA
//*     Jason Riedy, Computer Science Division, University of
//*       California at Berkeley, USA
//*
//* =====================================================================

// Mon (21/12/2010): This one does not allocate memory for eigenvectors (diag_dsygvx yes)
//                   ("eigvec" sould be a size*evec_size array of doubles)
//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx2(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu)
{
  /* Lapack input/output */
  char jobz, uplo, range;
  int  itype, lda, ldb, ldz, lwork, info, *iwork, *ifail, nulo, m, il;
  double *work,abstol;
//  double *p_evec;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
  uplo='U';  /* Upper triangle of A is stored */
  lda=size;
  ldb=size;
  ldz=size;
  il=1;
  abstol=0;
  lwork = 2*size*size+6*size+1;

  // Allocating workspace & initialization
  work = (double *) malloc(lwork * sizeof(double));
  iwork = (int *) malloc(5*size * sizeof(int));
  ifail = (int *) malloc(size * sizeof(int));
  if(!work || !iwork || !ifail)
  {
	  printf("Msg(diag_dsygvx2): Unable to allocate memory! Forcing exit!\n");
	  exit(1);
  }
  for(int i=0; i<size; i++)
	  ifail[i]=0;
  for(int i=0; i<lwork; i++)
	  work[i]=0.0;
  for(int i=0; i<5*size; i++)
	  iwork[i]=0;

//  evec_size = (iu-il+1);
//  p_evec = (double *) malloc( size * evec_size * sizeof(double) );
//  for(int i=0; i<size*evec_size; i++)
//	  p_evec[i]=0;

  // Diagonalization
//  dsygvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, &lda,
//		  mass_matrix, &ldb, &nulo, &nulo, &il, &iu, &abstol, &m,
//		  eigval, p_evec, &ldz, work, &lwork, iwork,
//          ifail, &info);
  dsygvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, &lda,
		  mass_matrix, &ldb, &nulo, &nulo, &il, &iu, &abstol, &m,
		  eigval, eigvec, &ldz, work, &lwork, iwork,
          ifail, &info);
//  printf("Msg(diag_dsygvx): %d eigenvectors found!",m);

//  // Outputting as usual... (eigenvectors inside hess_matrix)
//  for(int i=0; i<size*evec_size; i++)
//	  hess_matrix[i] = p_evec[i]; // copying eigenvectors into hess_matrix

//  free(p_evec);
  free(ifail);
  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork); // <-- "iwork" could be allocated only once outside!!!
}



//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx(double *hess_matrix,double *mass_matrix,double *eigval, int size, int il, int iu)
{
  /* Lapack input/output */
  char jobz, uplo, range;
  int  itype, lda, ldb, ldz, lwork, info, *iwork, *ifail, nulo, m, evec_size;
  double *work,abstol,*p_evec;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
  uplo='U';  /* Upper triangle of A is stored */
  lda=size;
  ldb=size;
  ldz=size;
  abstol=0;
  lwork = 2*size*size+6*size+1;
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++) work[i]=0.0;

  iwork = (int *) malloc(5*size * sizeof(int));
  for(int i=0; i<5*size; i++) iwork[i]=0;

  ifail = (int *) malloc(size * sizeof(int));
  for(int i=0; i<size; i++) ifail[i]=0;

  evec_size = (iu-il+1);
  p_evec = (double *) malloc( size * evec_size * sizeof(double) );
  for(int i=0; i<size*evec_size; i++) p_evec[i]=0;

  // Diagonalization
  dsygvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, &lda,
		  mass_matrix, &ldb, &nulo, &nulo, &il, &iu, &abstol, &m,
		  eigval, p_evec, &ldz, work, &lwork, iwork,
          ifail, &info);
//  printf("Msg(diag_dsygvx): %d eigenvectors found!",m);

  // Outputting as usual... (eigenvectors inside hess_matrix)
  for(int i=0; i<size*evec_size; i++)
	  hess_matrix[i] = p_evec[i]; // copying eigenvectors into hess_matrix

  free(p_evec);
  free(ifail);
  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork); // <-- "iwork" could be allocated only once outside!!!
}

//SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
//$                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
//$                   LWORK, IWORK, IFAIL, INFO )
//*
//*  -- LAPACK driver routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
//DOUBLE PRECISION   ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            IFAIL( * ), IWORK( * )
//DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),
//$                   Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSYGVX computes selected eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          Specifies the problem type to be solved:
//*          = 1:  A*x = (lambda)*B*x
//*          = 2:  A*B*x = (lambda)*x
//*          = 3:  B*A*x = (lambda)*x
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found.
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found.
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A and B are stored;
//*          = 'L':  Lower triangle of A and B are stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix pencil (A,B).  N >= 0.
//*
//*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
//*          On entry, the symmetric matrix A.  If UPLO = 'U', the
//*          leading N-by-N upper triangular part of A contains the
//*          upper triangular part of the matrix A.  If UPLO = 'L',
//*          the leading N-by-N lower triangular part of A contains
//*          the lower triangular part of the matrix A.
//*
//*          On exit, the lower triangle (if UPLO='L') or the upper
//*          triangle (if UPLO='U') of A, including the diagonal, is
//*          destroyed.
//*
//*  LDA     (input) INTEGER
//*          The leading dimension of the array A.  LDA >= max(1,N).
//*
//*  B       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
//*          On entry, the symmetric matrix B.  If UPLO = 'U', the
//*          leading N-by-N upper triangular part of B contains the
//*          upper triangular part of the matrix B.  If UPLO = 'L',
//*          the leading N-by-N lower triangular part of B contains
//*          the lower triangular part of the matrix B.
//*
//*          On exit, if INFO <= N, the part of B containing the matrix is
//*          overwritten by the triangular factor U or L from the Cholesky
//*          factorization B = U**T*U or B = L*L**T.
//*
//*  LDB     (input) INTEGER
//*          The leading dimension of the array B.  LDB >= max(1,N).
//*
//*  VL      (input) DOUBLE PRECISION
//*  VU      (input) DOUBLE PRECISION
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) DOUBLE PRECISION
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing A to tridiagonal form.
//*
//*          Eigenvalues will be computed most accurately when ABSTOL is
//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
//*          If this routine returns with INFO>0, indicating that some
//*          eigenvectors did not converge, try setting ABSTOL to
//*          2*DLAMCH('S').
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) DOUBLE PRECISION array, dimension (N)
//*          On normal exit, the first M elements contain the selected
//*          eigenvalues in ascending order.
//*
//*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'N', then Z is not referenced.
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          The eigenvectors are normalized as follows:
//*          if ITYPE = 1 or 2, Z**T*B*Z = I;
//*          if ITYPE = 3, Z**T*inv(B)*Z = I.
//*
//*          If an eigenvector fails to converge, then that column of Z
//*          contains the latest approximation to the eigenvector, and the
//*          index of the eigenvector is returned in IFAIL.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
//*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
//*
//*  LWORK   (input) INTEGER
//*          The length of the array WORK.  LWORK >= max(1,8*N).
//*          For optimal efficiency, LWORK >= (NB+3)*N,
//*          where NB is the blocksize for DSYTRD returned by ILAENV.
//*
//*          If LWORK = -1, then a workspace query is assumed; the routine
//*          only calculates the optimal size of the WORK array, returns
//*          this value as the first entry of the WORK array, and no error
//*          message related to LWORK is issued by XERBLA.
//*
//*  IWORK   (workspace) INTEGER array, dimension (5*N)
//*
//*  IFAIL   (output) INTEGER array, dimension (N)
//*          If JOBZ = 'V', then if INFO = 0, the first M elements of
//*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
//*          indices of the eigenvectors that failed to converge.
//*          If JOBZ = 'N', then IFAIL is not referenced.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  DPOTRF or DSYEVX returned an error code:
//*             <= N:  if INFO = i, DSYEVX failed to converge;
//*                    i eigenvectors failed to converge.  Their indices
//*                    are stored in array IFAIL.
//*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
//*                    minor of order i of B is not positive definite.
//*                    The factorization of B could not be completed and
//*                    no eigenvalues or eigenvectors were computed.
//*
//*  Further Details
//*  ===============
//*
//*  Based on contributions by
//*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
//*
//* =====================================================================



//*  DSYGVD computes ALL the eigenvalues, and optionally, the eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
//*  B are assumed to be symmetric and B is also positive definite.
//*  If eigenvectors are desired, it uses a divide and conquer algorithm.
void diag_dsygvd(double *hess_matrix,double *mass_matrix,double *eigval, int size)
{
  /* Lapack input/output */
  char jobz, uplo;
  int  itype, lda, ldb, lwork, liwork, info, *iwork;
  double *work;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  uplo='U';  /* Upper triangle of A is stored */
  lda=size;
  ldb=size;
  lwork = 2*size*size+6*size+1;
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++) work[i]=0.0;

  liwork = 5*size+3;
  iwork = (int *) malloc(liwork * sizeof(int));
  for(int i=0; i<liwork; i++) iwork[i]=0;

  // Diagonalization
  dsygvd_(&itype, &jobz, &uplo, &size, hess_matrix, &lda,
          mass_matrix, &ldb, eigval, work, &lwork, iwork,
          &liwork, &info);

  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork); // <-- "iwork" could be allocated only once outside!!!
}

/*
      SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     $                   LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
*  of a real generalized symmetric-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
*  B are assumed to be symmetric and B is also positive definite.
*  If eigenvectors are desired, it uses a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          matrix Z of eigenvectors.  The eigenvectors are normalized
*          as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
*          or the lower triangle (if UPLO='L') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the symmetric matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK >= 1.
*          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
*          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK >= 1.
*          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  DPOTRF or DSYEVD returned an error code:
*             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
*                    failed to converge; i off-diagonal elements of an
*                    intermediate tridiagonal form did not converge to
*                    zero;
*                    if INFO = i and JOBZ = 'V', then the algorithm
*                    failed to compute an eigenvalue while working on
*                    the submatrix lying in rows and columns INFO/(N+1)
*                    through mod(INFO,N+1);
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*/


//*  DSPEVX computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//*  can be selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//   (SMALL MEMORY REQUIREMENTS) (default compute eigenvalues and eigenvectors)
// - When using "vl" and "vu" to define the range, set [iu,il] range to indicate
// how many elements will be allocated: "emeI" (iu-il+1)
// - If using range='V', *pemeF="number of eigenpairs obtained"
void diag_dspevx(double *hess_matrix,double **ppeigval, double **ppeigvec, int size, int il, int iu, double vl, double vu, char jobz, char range, int *pemeF)
{
	// Lapack input/output
//	char jobz;
//	char uplo, range;
	char uplo;
	int  lwork, liwork, lifail, info;
	double *work;
	int *iwork,*ifail;
	double nulo=0.0;
//	int il = 1;
	int emeI,emeF;
	double *eigval;
	double *eigvec;

	eigval = *ppeigval;
	eigvec = *ppeigvec;

	emeI = iu-il+1; // M = IU-IL+1.
//	jobz='V';  // Compute eigenvalues and eigenvectors
//	range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
	uplo='U';  // Upper triangle of A is stored

	//	double tol = 0.0; // Abs tol
	//*          Eigenvalues will be computed most accurately when ABSTOL is
	//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
	//*          If this routine returns with INFO>0, indicating that some
	//*          eigenvectors did not converge, try setting ABSTOL to
	//*          2*DLAMCH('S').
	char cmach='S';
	double tol = 2 * dlamch_(&cmach);
//	printf("Msg(diag_dspevx): ABSTOL = %f\n",tol);

	//*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
	lwork = 8*size;
	if(!(work = (double *) malloc(lwork * sizeof(double))))
	{
		printf("Msg(diag_dspevx): I'm sorry, WORK memory allocation failed!\n"
				"Forcing exit!\n");
		exit(1);
	}

	for(int i=0; i<lwork; i++)
		work[i]=0.0;
	//*  IWORK   (workspace) INTEGER array, dimension (5*N)
	liwork = 5*size;
	iwork = (int *) malloc( liwork * sizeof(int) );
	for(int i=0; i<liwork; i++)
		iwork[i]=0.0;
	//*  IFAIL   (output) INTEGER array, dimension (N)
	lifail = size;
	ifail = (int *) malloc( lifail * sizeof(int) );
	for(int i=0; i<lifail; i++)
		ifail[i]=0.0;

	// Allocating requested eigenvalues (if requested by NULL)
	if(*ppeigval==NULL)
	{
		if( !(eigval = (double *) malloc( sizeof(double) * emeI) ) )
		{
			printf("Msg(diag_dspevx): I'm sorry, Eigenvalues memory allocation failed!\n"
					"Forcing exit!\n");
			exit(1);
		}
		*ppeigval = eigval; // output
	}

	if(jobz=='V')
	{
		// Allocating eigenvectors (if requested by NULL)
		if(*ppeigvec==NULL)
		{
			if( !(eigvec = (double *) malloc( sizeof(double) * size * emeI ) ) )
			{
				printf("Msg(diag_dspevx): I'm sorry, Eigenvectors memory allocation failed!\n"
						"Forcing exit!\n");
				exit(1);
			}
			*ppeigvec = eigvec; // output
		}
	}

	// Diagonalization
	dspevx_(&jobz, &range, &uplo, &size, hess_matrix, &vl, &vu, &il, &iu,
			&tol, &emeF, eigval, eigvec, &size, work, iwork, ifail, &info);

//printf("eigenvalues computed emeF= %d\n",emeF);
//for(int i=0; i<emeF; i++)
//	printf("%i %f",i,eigval[i]);
//printf("eigenvalues:\n");

	//SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
	//$                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,
	//$                   INFO )
	if(range=='V')
		*pemeF=emeF; // Output the number of eigenpairs obtained

	// Final checking
	if(range=='I')
		if(emeI!=emeF)
		{
			printf("Msg(diag_dspevx): I'm sorry, requested eigenvectors (%d) different from computed ones (%d)\n",emeI,emeF);
			exit(3);
		}
	if(info != 0)
	{
		printf("Msg(diag_dspevx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}

	free(work); // <-- "work" could be allocated only once outside!!!
	free(iwork);
	free(ifail);
	*ppeigval = eigval;
	*ppeigvec = eigvec;
}
//SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
//$                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,
//$                   INFO )
//*
//*  -- LAPACK driver routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, IU, LDZ, M, N
//DOUBLE PRECISION   ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            IFAIL( * ), IWORK( * )
//DOUBLE PRECISION   AP( * ), W( * ), WORK( * ), Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSPEVX computes selected eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//*  can be selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//*
//*  Arguments
//*  =========
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found;
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found;
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A is stored;
//*          = 'L':  Lower triangle of A is stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix A.  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, AP is overwritten by values generated during the
//*          reduction to tridiagonal form.  If UPLO = 'U', the diagonal
//*          and first superdiagonal of the tridiagonal matrix T overwrite
//*          the corresponding elements of A, and if UPLO = 'L', the
//*          diagonal and first subdiagonal of T overwrite the
//*          corresponding elements of A.
//*
//*  VL      (input) DOUBLE PRECISION
//*  VU      (input) DOUBLE PRECISION
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) DOUBLE PRECISION
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing AP to tridiagonal form.
//*
//*          Eigenvalues will be computed most accurately when ABSTOL is
//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
//*          If this routine returns with INFO>0, indicating that some
//*          eigenvectors did not converge, try setting ABSTOL to
//*          2*DLAMCH('S').
//*
//*          See "Computing Small Singular Values of Bidiagonal Matrices
//*          with Guaranteed High Relative Accuracy," by Demmel and
//*          Kahan, LAPACK Working Note #3.
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) DOUBLE PRECISION array, dimension (N)
//*          If INFO = 0, the selected eigenvalues in ascending order.
//*
//*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          If an eigenvector fails to converge, then that column of Z
//*          contains the latest approximation to the eigenvector, and the
//*          index of the eigenvector is returned in IFAIL.
//*          If JOBZ = 'N', then Z is not referenced.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
//*
//*  IWORK   (workspace) INTEGER array, dimension (5*N)
//*
//*  IFAIL   (output) INTEGER array, dimension (N)
//*          If JOBZ = 'V', then if INFO = 0, the first M elements of
//*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
//*          indices of the eigenvectors that failed to converge.
//*          If JOBZ = 'N', then IFAIL is not referenced.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  if INFO = i, then i eigenvectors failed to converge.
//*                Their indices are stored in array IFAIL.
//*
//*  =====================================================================

//*  SSPEVX computes SELECTED eigenvalues and, optionally, eigenvectors
//*  of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//*  can be selected by specifying either a range of values or a range of
//*  indices for the desired eigenvalues.
//   (SMALL MEMORY REQUIREMENTS) (default compute eigenvalues and eigenvectors)
void diag_sspevx(float *hess_matrix,float **ppeigval, float **ppeigvec, int size, int iu, char jobz)
{
	// Lapack input/output
//	char jobz;
	char uplo, range;
	int  lwork, liwork, lifail, info;
	float *work;
	int *iwork,*ifail;
	float nulo=0.0;
	int il = 1;
	int eme;
	float *eigval;
	float *eigvec;

	eigval = *ppeigval;
	eigvec = *ppeigvec;

	eme = iu-il+1; // M = IU-IL+1.
//	jobz='V';  // Compute eigenvalues and eigenvectors
	range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
	uplo='U';  // Upper triangle of A is stored

	//	double tol = 0.0; // Abs tol
	//*          Eigenvalues will be computed most accurately when ABSTOL is
	//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
	//*          If this routine returns with INFO>0, indicating that some
	//*          eigenvectors did not converge, try setting ABSTOL to
	//*          2*DLAMCH('S').
	char cmach='S';
	float tol = 2 * dlamch_(&cmach);
	printf("Msg(diag_sspevx): ABSTOL = %f\n",tol);

	//*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
	lwork = 8*size;
	if(!(work = (float *) malloc(lwork * sizeof(float))))
	{
		printf("Msg(diag_sspevx): I'm sorry, WORK memory allocation failed!\n"
				"Forcing exit!\n");
		exit(1);
	}

	for(int i=0; i<lwork; i++)
		work[i]=0.0;
	//*  IWORK   (workspace) INTEGER array, dimension (5*N)
	liwork = 5*size;
	iwork = (int *) malloc( liwork * sizeof(int) );
	for(int i=0; i<liwork; i++)
		iwork[i]=0.0;
	//*  IFAIL   (output) INTEGER array, dimension (N)
	lifail = size;
	ifail = (int *) malloc( lifail * sizeof(int) );
	for(int i=0; i<lifail; i++)
		ifail[i]=0.0;

	// Allocating eigenvalues (if requested by NULL)
	if(*ppeigval==NULL)
	{
		if( !(eigval = (float *) malloc( sizeof(float) * iu) ) )
		{
			printf("Msg(diag_dspevx): I'm sorry, Eigenvalues memory allocation failed!\n"
					"Forcing exit!\n");
			exit(1);
		}
		*ppeigval = eigval; // output
	}

	// Allocating eigenvectors (if requested by NULL)
	if(*ppeigvec==NULL)
	{
		if( !(eigvec = (float *) malloc( sizeof(float) * size * iu ) ) )
		{
			printf("Msg(diag_sspevx): I'm sorry, Eigenvectors memory allocation failed!\n"
					"Forcing exit!\n");
			exit(1);
		}
		*ppeigvec = eigvec; // output
	}

	// Diagonalization
	sspevx_(&jobz, &range, &uplo, &size, hess_matrix, &nulo, &nulo, &il, &iu,
			&tol, &eme, eigval, eigvec, &size, work, iwork, ifail, &info);
	//SUBROUTINE SSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,
	//$                   ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,
	//$                   INFO )
	if(info != 0)
	{
		printf("Msg(diag_sspevx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
	free(work); // <-- "work" could be allocated only once outside!!!
	free(iwork);
	free(ifail);
	*ppeigval = eigval;
	*ppeigvec = eigvec;
}
//  -- LAPACK driver routine (version 3.1) --
//       Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//       November 2006
//
//
//    Purpose
//    =======
//
//    SSPEVX computes selected eigenvalues and, optionally, eigenvectors
//    of a real symmetric matrix A in packed storage.  Eigenvalues/vectors
//    can be selected by specifying either a range of values or a range of
//    indices for the desired eigenvalues.
//
//    Arguments
//    =========
//
//    JOBZ    (input) CHARACTER*1
//            = 'N':  Compute eigenvalues only;
//            = 'V':  Compute eigenvalues and eigenvectors.
//
//    RANGE   (input) CHARACTER*1
//            = 'A': all eigenvalues will be found;
//            = 'V': all eigenvalues in the half-open interval (VL,VU]
//                   will be found;
//            = 'I': the IL-th through IU-th eigenvalues will be found.
//
//    UPLO    (input) CHARACTER*1
//            = 'U':  Upper triangle of A is stored;
//            = 'L':  Lower triangle of A is stored.
//
//    N       (input) INTEGER
//            The order of the matrix A.  N >= 0.
//
//    AP      (input/output) REAL array, dimension (N*(N+1)/2)
//            On entry, the upper or lower triangle of the symmetric matrix
//            A, packed columnwise in a linear array.  The j-th column of A
//            is stored in the array AP as follows:
//            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//
//            On exit, AP is overwritten by values generated during the
//            reduction to tridiagonal form.  If UPLO = 'U', the diagonal
//            and first superdiagonal of the tridiagonal matrix T overwrite
//            the corresponding elements of A, and if UPLO = 'L', the
//            diagonal and first subdiagonal of T overwrite the
//            corresponding elements of A.
//
//    VL      (input) REAL
//    VU      (input) REAL
//            If RANGE='V', the lower and upper bounds of the interval to
//            be searched for eigenvalues. VL < VU.
//            Not referenced if RANGE = 'A' or 'I'.
//
//    IL      (input) INTEGER
//    IU      (input) INTEGER
//            If RANGE='I', the indices (in ascending order) of the
//            smallest and largest eigenvalues to be returned.
//            1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//            Not referenced if RANGE = 'A' or 'V'.
//
//    ABSTOL  (input) REAL
//            The absolute error tolerance for the eigenvalues.
//            An approximate eigenvalue is accepted as converged
//            when it is determined to lie in an interval [a,b]
//            of width less than or equal to
//
//                    ABSTOL + EPS *   max( |a|,|b| ) ,
//
//            where EPS is the machine precision.  If ABSTOL is less than
//            or equal to zero, then  EPS*|T|  will be used in its place,
//            where |T| is the 1-norm of the tridiagonal matrix obtained
//            by reducing AP to tridiagonal form.
//
//            Eigenvalues will be computed most accurately when ABSTOL is
//            set to twice the underflow threshold 2*SLAMCH('S'), not zero.
//            If this routine returns with INFO>0, indicating that some
//            eigenvectors did not converge, try setting ABSTOL to
//            2*SLAMCH('S').
//
//            See "Computing Small Singular Values of Bidiagonal Matrices
//            with Guaranteed High Relative Accuracy," by Demmel and
//            Kahan, LAPACK Working Note #3.
//
//    M       (output) INTEGER
//            The total number of eigenvalues found.  0 <= M <= N.
//            If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//
//    W       (output) REAL array, dimension (N)
//            If INFO = 0, the selected eigenvalues in ascending order.
//
//    Z       (output) REAL array, dimension (LDZ, max(1,M))
//            If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//            contain the orthonormal eigenvectors of the matrix A
//            corresponding to the selected eigenvalues, with the i-th
//            column of Z holding the eigenvector associated with W(i).
//            If an eigenvector fails to converge, then that column of Z
//            contains the latest approximation to the eigenvector, and the
//            index of the eigenvector is returned in IFAIL.
//            If JOBZ = 'N', then Z is not referenced.
//            Note: the user must ensure that at least max(1,M) columns are
//            supplied in the array Z; if RANGE = 'V', the exact value of M
//            is not known in advance and an upper bound must be used.
//
//    LDZ     (input) INTEGER
//            The leading dimension of the array Z.  LDZ >= 1, and if
//            JOBZ = 'V', LDZ >= max(1,N).
//
//    WORK    (workspace) REAL array, dimension (8*N)
//
//    IWORK   (workspace) INTEGER array, dimension (5*N)
//
//    IFAIL   (output) INTEGER array, dimension (N)
//            If JOBZ = 'V', then if INFO = 0, the first M elements of
//            IFAIL are zero.  If INFO > 0, then IFAIL contains the
//            indices of the eigenvectors that failed to converge.
//            If JOBZ = 'N', then IFAIL is not referenced.
//
//    INFO    (output) INTEGER
//            = 0:  successful exit
//            < 0:  if INFO = -i, the i-th argument had an illegal value
//            > 0:  if INFO = i, then i eigenvectors failed to converge.
//                  Their indices are stored in array IFAIL.
//
//    =====================================================================



//*  DSPGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
//   (SMALL MEMORY REQUIREMENTS) (by default compute eigenvalues and eigenvectors)
int diag_dspgvx(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu)
{
  bool debug = false;
  /* Lapack input/output */
  char jobz, uplo, range;
  int  itype, lwork, liwork, lifail, info;
  double *work;
  int *iwork,*ifail;
  double nulo=0.0;
  int il = 1;
  int eme = iu;

  itype=1;   // specifies A*x = lambda*B*x
  jobz='V';  // Compute eigenvalues and eigenvectors
  uplo='U';  // Upper triangle of A is stored
  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
//  ldz=size;

//  double tol = 0.0; // Abs tol
	//*          Eigenvalues will be computed most accurately when ABSTOL is
	//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
	//*          If this routine returns with INFO>0, indicating that some
	//*          eigenvectors did not converge, try setting ABSTOL to
	//*          2*DLAMCH('S').
 	char cmach='S';
  	double tol = 2 * dlamch_(&cmach);
  	if(debug)
  		printf("Msg(diag_dspgvx): ABSTOL = %E\n",tol);
//	tol = 1e-17;

  //*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
  lwork = 8*size;
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++)
	  work[i]=0.0;
  //*  IWORK   (workspace) INTEGER array, dimension (5*N)
  liwork = 5*size;
  iwork = (int *) malloc( liwork * sizeof(int) );
  for(int i=0; i<liwork; i++)
	  iwork[i]=0.0;
  //*  IFAIL   (output) INTEGER array, dimension (N)
  lifail = size;
  ifail = (int *) malloc( lifail * sizeof(int) );
  for(int i=0; i<lifail; i++)
	  ifail[i]=0.0;

  // Diagonalization
  dspgvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, mass_matrix, &nulo, &nulo, &il, &iu,
		  &tol, &eme, eigval, eigvec, &size, work, iwork, ifail, &info);
  if(info != 0)
  {
//	  printf("Msg(diag_dspgvx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
//	  exit(3);
	  return(info); // Errors (info != 0) must be handled outside!
  }

  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork);
  free(ifail);
  return 0; // successful exit!
}
//SUBROUTINE DSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,
//$                   IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,
//$                   IFAIL, INFO )
//*
//*  -- LAPACK driver routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N
//DOUBLE PRECISION   ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            IFAIL( * ), IWORK( * )
//DOUBLE PRECISION   AP( * ), BP( * ), W( * ), WORK( * ),
//$                   Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSPGVX computes selected eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.  Eigenvalues and eigenvectors can be
//*  selected by specifying either a range of values or a range of indices
//*  for the desired eigenvalues.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          Specifies the problem type to be solved:
//*          = 1:  A*x = (lambda)*B*x
//*          = 2:  A*B*x = (lambda)*x
//*          = 3:  B*A*x = (lambda)*x
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found.
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found.
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A and B are stored;
//*          = 'L':  Lower triangle of A and B are stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix pencil (A,B).  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, the contents of AP are destroyed.
//*
//*  BP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          B, packed columnwise in a linear array.  The j-th column of B
//*          is stored in the array BP as follows:
//*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
//*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
//*
//*          On exit, the triangular factor U or L from the Cholesky
//*          factorization B = U**T*U or B = L*L**T, in the same storage
//*          format as B.
//*
//*  VL      (input) DOUBLE PRECISION
//*  VU      (input) DOUBLE PRECISION
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) DOUBLE PRECISION
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing A to tridiagonal form.
//*
//*          Eigenvalues will be computed most accurately when ABSTOL is
//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
//*          If this routine returns with INFO>0, indicating that some
//*          eigenvectors did not converge, try setting ABSTOL to
//*          2*DLAMCH('S').
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) DOUBLE PRECISION array, dimension (N)
//*          On normal exit, the first M elements contain the selected
//*          eigenvalues in ascending order.
//*
//*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'N', then Z is not referenced.
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          The eigenvectors are normalized as follows:
//*          if ITYPE = 1 or 2, Z**T*B*Z = I;
//*          if ITYPE = 3, Z**T*inv(B)*Z = I.
//*
//*          If an eigenvector fails to converge, then that column of Z
//*          contains the latest approximation to the eigenvector, and the
//*          index of the eigenvector is returned in IFAIL.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
//*
//*  IWORK   (workspace) INTEGER array, dimension (5*N)
//*
//*  IFAIL   (output) INTEGER array, dimension (N)
//*          If JOBZ = 'V', then if INFO = 0, the first M elements of
//*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
//*          indices of the eigenvectors that failed to converge.
//*          If JOBZ = 'N', then IFAIL is not referenced.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  DPPTRF or DSPEVX returned an error code:
//*             <= N:  if INFO = i, DSPEVX failed to converge;
//*                    i eigenvectors failed to converge.  Their indices
//*                    are stored in array IFAIL.
//*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
//*                    minor of order i of B is not positive definite.
//*                    The factorization of B could not be completed and
//*                    no eigenvalues or eigenvectors were computed.
//*
//*  Further Details
//*  ===============
//*
//*  Based on contributions by
//*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
//*
//* =====================================================================


//// MRRR from P.Bientinesi (14/12/2010)
////######################################################################
////Symmetric-definite case (packed storage):                            #
////######################################################################
////# int dspgeig(int *itype, char *jobz, char *range, char *uplo,       #
////#             int *n, double *AP, double *BP, double *vl,            #
////#             double *vu, int *il, int *iu, int *m, double *W,       #
////#             double *Z, int *ldz);                                  #
////######################################################################
//int diag_dspgeig(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size, int iu)
//{
//  bool debug = false;
//  /* Lapack input/output */
//  char jobz, uplo, range;
//  int  itype, lwork, liwork, lifail, info;
//  double *work;
//  int *iwork,*ifail;
//  double nulo=0.0;
//  int il = 1;
//  int eme = iu;
//
//  itype=1;   // specifies A*x = lambda*B*x
//  jobz='V';  // Compute eigenvalues and eigenvectors
//  uplo='U';  // Upper triangle of A is stored
//  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
////  ldz=size;
//
////  double tol = 0.0; // Abs tol
//	//*          Eigenvalues will be computed most accurately when ABSTOL is
//	//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
//	//*          If this routine returns with INFO>0, indicating that some
//	//*          eigenvectors did not converge, try setting ABSTOL to
//	//*          2*DLAMCH('S').
// 	char cmach='S';
//  	double tol = 2 * dlamch_(&cmach);
//  	if(debug)
//  		printf("Msg(diag_dspgvx): ABSTOL = %E\n",tol);
////	tol = 1e-17;
//
//  //*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
//  lwork = 8*size;
//  work = (double *) malloc(lwork * sizeof(double));
//  for(int i=0; i<lwork; i++)
//	  work[i]=0.0;
//  //*  IWORK   (workspace) INTEGER array, dimension (5*N)
//  liwork = 5*size;
//  iwork = (int *) malloc( liwork * sizeof(int) );
//  for(int i=0; i<liwork; i++)
//	  iwork[i]=0.0;
//  //*  IFAIL   (output) INTEGER array, dimension (N)
//  lifail = size;
//  ifail = (int *) malloc( lifail * sizeof(int) );
//  for(int i=0; i<lifail; i++)
//	  ifail[i]=0.0;
//
//  // Diagonalization
////  dspgvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, mass_matrix, &nulo, &nulo, &il, &iu,
////		  &tol, &eme, eigval, eigvec, &size, work, iwork, ifail, &info);
//
//  //# int dspgeig(int *itype, char *jobz, char *range, char *uplo,       #
//  //#             int *n, double *AP, double *BP, double *vl,            #
//  //#             double *vu, int *il, int *iu, int *m, double *W,       #
//  //#             double *Z, int *ldz);                                  #
//  info = dspgeig(&itype, &jobz, &range, &uplo, &size, hess_matrix, mass_matrix, &nulo, &nulo, &il, &iu,
//		  &eme, eigval, eigvec, &size);
//
//  if(info != 0)
//  {
////	  printf("Msg(diag_dspgvx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
////	  exit(3);
//	  return(info); // Errors (info != 0) must be handled outside!
//  }
//
//  free(work); // <-- "work" could be allocated only once outside!!!
//  free(iwork);
//  free(ifail);
//  return 0; // successful exit!
//}

//######################################################################
//# C function prototype:                                              #
//######################################################################
//#                                                                    #
//# int dspgeig(int *itype, char *jobz, char *range, char *uplo,       #
//#             int *n, double *AP, double *BP, double *vl,            #
//#             double *vu, int *il, int *iu, int *m, double *W,       #
//#             double *Z, int *ldz);                                  #
//#                                                                    #
//# Arguments:                                                         #
//# ----------                                                         #
//#                                                                    #
//# INPUTS:                                                            #
//# -------                                                            #
//# itype              1  - A*x = lambda*B*x                           #
//#                    2  - A*B*x = lambda*x                           #
//#                    3  - B*A*x = lambda*x                           #
//# jobz              "N" - compute only eigenvalues                   #
//#                   "V" - compute also eigenvectors                  #
//# range             "A" - all                                        #
//#                   "V" - by interval: (VL,VU]                       #
//#                   "I" - by index:     IL-IU                        #
//# uplo              "L" - Upper triangle of A and B stored           #
//#                   "U" - Lower triangle of A and B stored           #
//# n                 Order of the matrix A and B                      #
//# ldz               Leading dimension of matrix Z;                   #
//#                   often equal to matrix size n                     #
//#                                                                    #
//# INPUT + OUTPUT:                                                    #
//# ---------------                                                    #
//# AP                On entry symmetric input matrix, stored in       #
//# (double[s])       packed format by columns. Depending on the       #
//# s = (n*(n+1))/2   value of 'uplo' only the upper or lower          #
//#                   triangular part is stored                        #
//#                   On output the array will contain the             #
//#                   'm' computed eigenvectors                        #
//# BP                On entry symmetric definite input matrix, stored #
//# (double[s])       in packed format by columns.                     #
//# s = (n*(n+1))/2   Depending on the value of 'uplo' only the upper  #
//#                   or lower triangular part is stored               #
//#                   On output overwritten                            #
//# vl                If range="V", lower bound of interval            #
//#                   (vl,vu], on output refined                       #
//#                   If range="A" or "I" not referenced as input      #
//#                   On output the interval (vl,vu] contains ALL      #
//#                   the computed eigenvalues                         #
//# vu                If range="V", upper bound of interval            #
//#                   (vl,vu], on output refined                       #
//#                   If range="A" or "I" not referenced as input.     #
//#                   On output the interval (vl,vu] contains ALL      #
//#                   the computed eigenvalues                         #
//# il                If range="I", lower index (1-based indexing) of  #
//#                   the subset 'il' to 'iu'                          #
//#                   If range="A" or "V" not referenced as input.     #
//#                   On output the eigenvalues with index il to iu    #
//#                   are computed                                     #
//# iu                If range="I", upper index (1-based indexing) of  #
//#                   the subset 'il' to 'iu'                          #
//#                   If range="A" or "V" not referenced as input      #
//#                   On output the eigenvalues with index il to iu    #
//#                   are computed                                     #
//#                                                                    #
//# OUTPUT:                                                            #
//# -------                                                            #
//# m                 Number of eigenvalues and eigenvectors computed  #
//# W (double[n])     Eigenvalues                                      #
//#                   The first 'm' entries contain the eigenvalues    #
//# Z (double[m*n])   Eigenvectors                                     #
//#                   Enough space must be provided to store the       #
//#                   vectors. 'm' should be bigger or equal           #
//#                   to 'n' for range="A" or "V" and 'iu-il+1' for    #
//#                   range="I".                                       #
//#                                                                    #
//# NOTICE:           The routine will allocate work space of size     #
//#                   double[n*n] for range="A" or "V" and double[m*n] #
//#                   for range="I"                                    #
//#                                                                    #
//# RETURN VALUE:                                                      #
//# -------------                                                      #
//#                 0 - Success                                        #
//#                 1 - Wrong input parameter                          #
//#                 2 - Misc errors                                    #
//#                                                                    #
//######################################################################



//*  SSPGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
//   (SMALL MEMORY REQUIREMENTS)
void diag_sspgvx(float *hess_matrix,float *mass_matrix,float *eigval, float *eigvec, int size, int iu)
{
  /* Lapack input/output */
  char jobz, uplo, range;
  int  itype, lwork, liwork, lifail, info;
  float *work;
  int *iwork,*ifail;
  float nulo=0.0;
  int il = 1;
  int eme = iu;

  itype=1;   // specifies A*x = lambda*B*x
  jobz='V';  // Compute eigenvalues and eigenvectors
  uplo='U';  // Upper triangle of A is stored
  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
  float tol = 0.0; // Abs tol

  //*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
  lwork = 8*size;
  work = (float *) malloc(lwork * sizeof(float));
  for(int i=0; i<lwork; i++)
	  work[i]=0.0;
  //*  IWORK   (workspace) INTEGER array, dimension (5*N)
  liwork = 5*size;
  iwork = (int *) malloc( liwork * sizeof(int) );
  for(int i=0; i<liwork; i++)
	  iwork[i]=0.0;
  //*  IFAIL   (output) INTEGER array, dimension (N)
  lifail = size;
  ifail = (int *) malloc( lifail * sizeof(int) );
  for(int i=0; i<lifail; i++)
	  ifail[i]=0.0;

  // Diagonalization
  sspgvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, mass_matrix, &nulo, &nulo, &il, &iu,
		  &tol, &eme, eigval, eigvec, &size, work, iwork, ifail, &info);
  if(info != 0)
  {
	  printf("Msg(diag_sspgvx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }

  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork);
  free(ifail);
}
//SUBROUTINE SSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,
//$                   IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,
//$                   IFAIL, INFO )
//*
//*  -- LAPACK driver routine (version 3.2) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, ITYPE, IU, LDZ, M, N
//REAL               ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            IFAIL( * ), IWORK( * )
//REAL               AP( * ), BP( * ), W( * ), WORK( * ),
//$                   Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  SSPGVX computes selected eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.  Eigenvalues and eigenvectors can be
//*  selected by specifying either a range of values or a range of indices
//*  for the desired eigenvalues.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          Specifies the problem type to be solved:
//*          = 1:  A*x = (lambda)*B*x
//*          = 2:  A*B*x = (lambda)*x
//*          = 3:  B*A*x = (lambda)*x
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found.
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found.
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A and B are stored;
//*          = 'L':  Lower triangle of A and B are stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix pencil (A,B).  N >= 0.
//*
//*  AP      (input/output) REAL array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, the contents of AP are destroyed.
//*
//*  BP      (input/output) REAL array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          B, packed columnwise in a linear array.  The j-th column of B
//*          is stored in the array BP as follows:
//*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
//*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
//*
//*          On exit, the triangular factor U or L from the Cholesky
//*          factorization B = U**T*U or B = L*L**T, in the same storage
//*          format as B.
//*
//*  VL      (input) REAL
//*  VU      (input) REAL
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) REAL
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing A to tridiagonal form.
//*
//*          Eigenvalues will be computed most accurately when ABSTOL is
//*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
//*          If this routine returns with INFO>0, indicating that some
//*          eigenvectors did not converge, try setting ABSTOL to
//*          2*SLAMCH('S').
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) REAL array, dimension (N)
//*          On normal exit, the first M elements contain the selected
//*          eigenvalues in ascending order.
//*
//*  Z       (output) REAL array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'N', then Z is not referenced.
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          The eigenvectors are normalized as follows:
//*          if ITYPE = 1 or 2, Z**T*B*Z = I;
//*          if ITYPE = 3, Z**T*inv(B)*Z = I.
//*
//*          If an eigenvector fails to converge, then that column of Z
//*          contains the latest approximation to the eigenvector, and the
//*          index of the eigenvector is returned in IFAIL.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace) REAL array, dimension (8*N)
//*
//*  IWORK   (workspace) INTEGER array, dimension (5*N)
//*
//*  IFAIL   (output) INTEGER array, dimension (N)
//*          If JOBZ = 'V', then if INFO = 0, the first M elements of
//*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
//*          indices of the eigenvectors that failed to converge.
//*          If JOBZ = 'N', then IFAIL is not referenced.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  SPPTRF or SSPEVX returned an error code:
//*             <= N:  if INFO = i, SSPEVX failed to converge;
//*                    i eigenvectors failed to converge.  Their indices
//*                    are stored in array IFAIL.
//*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
//*                    minor of order i of B is not positive definite.
//*                    The factorization of B could not be completed and
//*                    no eigenvalues or eigenvectors were computed.
//*
//*  Further Details
//*  ===============
//*
//*  Based on contributions by
//*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
//*
//* =====================================================================


//   XSPGVX --> valid for single/double precision floating point data
//   WARNING: floating precision detected by "sizeof()" !!!
//*  DSPGVX and SSPGVX compute SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric, stored in packed storage, and B
//*  is also positive definite.
void diag_xspgvx(floating *hess_matrix,floating *mass_matrix,floating *eigval, floating *eigvec, int size, int iu)
{
  /* Lapack input/output */
  char jobz, uplo, range;
  int  itype, lwork, liwork, lifail, info;
  floating *work;
  float *swork;
  double *dwork;
  int *iwork,*ifail;
  floating nulo=0.0;
  float snulo;
  double dnulo;
  int il = 1;
  int eme = iu;
  float *shess,*smass,*seval,*sevec;
  double *dhess,*dmass,*deval,*devec;

  itype=1;   // specifies A*x = lambda*B*x
  jobz='V';  // Compute eigenvalues and eigenvectors
  uplo='U';  // Upper triangle of A is stored
  range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
//  ldz=size;

//  double tol = 0.0; // Abs tol
	//*          Eigenvalues will be computed most accurately when ABSTOL is
	//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
	//*          If this routine returns with INFO>0, indicating that some
	//*          eigenvectors did not converge, try setting ABSTOL to
	//*          2*DLAMCH('S').
 	char cmach='S';
// 	floating tol = 2 * dlamch_(&cmach);
 	float stol;
 	double dtol;
//	printf("Msg(diag_xspevx): ABSTOL = %E\n",tol);
//	tol = 1e-17;

  //*  WORK    (workspace) DOUBLE PRECISION array, dimension (8*N)
  lwork = 8*size;
  work = (floating *) malloc(lwork * sizeof(floating));
  for(int i=0; i<lwork; i++)
	  work[i]=0.0;
  //*  IWORK   (workspace) INTEGER array, dimension (5*N)
  liwork = 5*size;
  iwork = (int *) malloc( liwork * sizeof(int) );
  for(int i=0; i<liwork; i++)
	  iwork[i]=0.0;
  //*  IFAIL   (output) INTEGER array, dimension (N)
  lifail = size;
  ifail = (int *) malloc( lifail * sizeof(int) );
  for(int i=0; i<lifail; i++)
	  ifail[i]=0.0;

  // Diagonalization
  if( sizeof(floating) == 4 ) // SINGLE precision
  {
	  stol = (float) 2 * slamch_(&cmach);
	  printf("Msg(diag_xspevx): SINGLE precision! ABSTOL  stol= %E\n",stol);
	  swork = (float *) work;
	  snulo = (float) nulo;
	  shess = (float *) hess_matrix;
	  smass = (float *) mass_matrix;
	  seval = (float *) eigval;
	  sevec = (float *) eigvec;
	  sspgvx_(&itype, &jobz, &range, &uplo, &size, shess, smass, &snulo, &snulo, &il, &iu,
			  &stol, &eme, seval, sevec, &size, swork, iwork, ifail, &info);
  }
  else // DOUBLE precision
  {
	  dtol = (double) 2 * dlamch_(&cmach);
//	  printf("Msg(diag_xspevx): DOUBLE precision! ABSTOL dtol= %E\n",dtol);
	  dwork = (double *) work;
	  dnulo = (double) nulo;
	  dhess = (double *) hess_matrix;
	  dmass = (double *) mass_matrix;
	  deval = (double *) eigval;
	  devec = (double *) eigvec;
	  dspgvx_(&itype, &jobz, &range, &uplo, &size, dhess, dmass, &dnulo, &dnulo, &il, &iu,
		  &dtol, &eme, deval, devec, &size, dwork, iwork, ifail, &info);
  }

  if(info != 0)
  {
	  printf("Msg(diag_xspgvx): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }

  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork);
  free(ifail);
}

// ARPACK-based diagonalization routine (from: J.I.Aliaga & E.Quintana)
// Very fast if small number of eigenpairs requested. Not optimized for packed storage, thus it does not scale linear with number of threads.
// A --> Hessian matrix
// B --> Kinetic energy matrix
// d --> Eigenvalues array
// z --> Eigenvectors array
// nrows --> Problem dimension (size)
// nev --> Number of eigenpairs to be computed
// ncv --> Factor for the Number of columns of matrix V (=ncv*nev). It should be optimized for a given problem (OPTIONAL)
void dsdrv1_AP_BP_W_mon(double *B, double *A, double *d, double *z, int nrows, int nev, float ncvf)
{
//	int nrows;

	char *bmat = "I";

	//  WHICH   Character*2.  (INPUT)
	//		  Specify which of the Ritz values of OP to compute.
	//		  'LM' - compute the NEV largest (in magnitude) eigenvalues.
	char *which = "LM";

	//  IPARAM(3) = MXITER
	//  On INPUT:  maximum number of Arnoldi update iterations allowed.
	//  On OUTPUT: actual number of Arnoldi update iterations taken.
	int maxit = 50; //

	int i,j, inc = 1, ido = 0, ldv = nrows, info = 0, nsteps = 0;

	// Extracted from DSAUPD remarks:
	//	  4. At present there is no a-priori analysis to guide the selection
	//	     of NCV relative to NEV.  The only formal requirement is that NCV > NEV.
	//	     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
	//	     the same type are to be solved, one should experiment with increasing
	//	     NCV while keeping NEV fixed for a given test problem.  This will
	//	     usually decrease the required number of OP*x operations but it
	//	     also increases the work and storage required to maintain the orthogonal
	//	     basis vectors.   The optimal "cross-over" with respect to CPU time
	//	     is problem dependent and must be determined empirically.
	int ncv = (int) (ncvf * (float) nev);

	int iparam[11], ipntr[11], lworkl = (ncv * ncv + 8 * ncv);
	double tol = 0.0, alpha = 1.0, beta = 0.0;
	double *resid = NULL, *v = NULL, *workd = NULL, *workl = NULL, *x = NULL,
			*y = NULL, *xy = NULL;
//	double elap1, ucpu1, elap2, ucpu2, time[10][3];

	char *howny = "A";
	int rvec = 1, *select = NULL, ldz = nrows;
	double sigma = 0.0, normA = 0.0, normB = 0.0, norm =	0.0;


	// Compute Cholesky factorization of B
	//	print_band (B, nrows, "B", 0, 10, 1, 10);
	info = 0;
	//  DPPTRF computes the Cholesky factorization of a real symmetric
	//  positive definite matrix A stored in packed format.
	//  The factorization has the form
	//	 A = U**T * U,  if UPLO = 'U', or
	//	 A = L  * L**T,  if UPLO = 'L',
	//  where U is an upper triangular matrix and L is lower triangular.

	//	show_matrix(B, nrows, "B-matrix from Mon");
//	fprintf(stderr, "BEFORE nrows %d, info %d\n B[0]= ", nrows, info);
//	for(i=0; i<10; i++)
//		fprintf(stderr, "%f ", B[i]);
//	fprintf(stderr,"\n");

	dpptrf_("U", &nrows, B, &info);

//	fprintf(stderr, "AFTER nrows %d, info %d\n B[0]= ", nrows, info);
//	for(i=0; i<10; i++)
//		fprintf(stderr, "%f ", B[i]);
//	fprintf(stderr,"\n");

	//	print_band (B, nrows, "B_G", 0, 10, 1, 10);

	// Create ARPACK structures for ARPACK loop
//	v = create_rect_matrix(ncv, ldv);
	if ( !(v = (double *) malloc (ncv * ldv * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

//	resid = create_vector(nrows);
	if ( !(resid = (double *) malloc (nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

//	workd = create_vector(3 * nrows);
	if ( !(workd = (double *) malloc (3 * nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

//	workl = create_vector(lworkl);
	if ( !(workl = (double *) malloc (lworkl * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

//	xy = create_vector(nrows);
	if ( !(xy = (double *) malloc (nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

	for (i = 0; i < 11; i++)
	{
		iparam[i] = 0;
		ipntr[i] = 0;
	}
	iparam[0] = 1;
	iparam[2] = maxit;
	iparam[6] = 1;

	// ARPACK loop
	nsteps = 0;

	//  Reverse communication interface for the Implicitly Restarted Arnoldi
	//  Iteration.  For symmetric problems this reduces to a variant of the Lanczos
	//  method.  This method has been designed to compute approximations to a
	//  few eigenpairs of a linear operator OP that is real and symmetric
	//  with respect to a real positive semi-definite symmetric matrix B,
	//  i.e.:  B*OP = (OP')*B
	dsaupd_(&ido, bmat, &nrows, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
//	fprintf(stderr,"ido %d, bmat %s, nrows %d, which %s, nev %d, tol %f, resid %016lx, ncv %d, v %016lx, ldv %d, iparam %016lx, ipntr %016lx, workd %016lx, workl %016lx, lworkl %d, info %d\n",
//			        ido, bmat, nrows, which, nev, tol, (unsigned long)resid, ncv, (unsigned long)v, ldv, (unsigned long)iparam, (unsigned long)ipntr,
//			        (unsigned long)workd, (unsigned long)workl, lworkl, info);

	while ((ido == 1) || (ido == -1))
	{
		x = workd + ipntr[0] - 1;
		y = workd + ipntr[1] - 1;
		nsteps++;

		// dcopy_ copies a vector, x, to a vector, y.
		// uses unrolled loops for increments equal to one.
		dcopy_(&nrows, x, &inc, xy, &inc); // memcpy (xy, x, nrows*sizeof(double));

		//	DTPTRS solves a triangular system of the form
		//	   A * X = B  or  A**T * X = B,
		//	where A is a triangular matrix of order N stored in packed format,
		//	and B is an N-by-NRHS matrix.  A check is made to verify that A is
		//	nonsingular.
		dtptrs_("U", "N", "N", &nrows, &inc, B, xy, &nrows, &info);

		//  DSPMV  performs the matrix-vector operation */
		//     y := alpha*A*x + beta*y, */
		//  where alpha and beta are scalars, x and y are n element vectors and */
		//  A is an n by n symmetric matrix, supplied in packed form. */
		dspmv_("U", &nrows, &alpha, A, xy, &inc, &beta, y, &inc);

		//	DTPTRS solves a triangular system of the form
		//	   A * X = B  or  A**T * X = B,
		//	where A is a triangular matrix of order N stored in packed format,
		//	and B is an N-by-NRHS matrix.  A check is made to verify that A is
		//	nonsingular.
		dtptrs_("U", "T", "N", &nrows, &inc, B, y, &nrows, &info);

		dsaupd_(&ido, bmat, &nrows, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

	}
//	printf("IDO = %d , INFO = %d \n", ido, info);

	// Create ARPACK structures to extract eigenvalues and eigenvectors
//	select = create_indices(ncv);
	if( !(select = (int *) malloc (ncv * sizeof(int))) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}
	for (i = 0; i < ncv; i++)
		select[i] = 1;
//	d = create_vector(nev); // EIGENVALUES
//	if ( !(d = (double *) malloc (nev * sizeof(double)) ) ) {
//		printf ("Memory allocation error\n"); exit(-1);
//	}
	ldz = nrows;
//	z = create_rect_matrix(nev, nrows); // EIGENVECTORS
//	if ( !(z = (double *) malloc (nev * nrows * sizeof(double)) ) ) {
//		printf ("Memory allocation error\n"); exit(-1);
//	}


	//    DSEUPD returns the converged approximations to eigenvalues
	//    of A*z = lambda*B*z and (optionally):
	//
	//        (1) the corresponding approximate eigenvectors,
	//        (2) an orthonormal (Lanczos) basis for the associated approximate
	//            invariant subspace,
	//        (3) Both.
	//
	//    There is negligible additional cost to obtain eigenvectors.  An orthonormal
	//    (Lanczos) basis is always computed.  There is an additional storage cost
	//    of n*nev if both are requested (in this case a separate array Z must be
	//    supplied).
	//
	//    These quantities are obtained from the Lanczos factorization computed
	//    by DSAUPD for the linear operator OP prescribed by the MODE selection
	//    (see IPARAM(7) in DSAUPD documentation.)  DSAUPD must be called before
	//    this routine is called. These approximate eigenvalues and vectors are
	//    commonly called Ritz values and Ritz vectors respectively.  They are
	//    referred to as such in the comments that follow.   The computed orthonormal
	//    basis for the invariant subspace corresponding to these Ritz values is
	//    referred to as a Lanczos basis.
	//
	//    See documentation in the header of the subroutine DSAUPD for a definition
	//    of OP as well as other terms and the relation of computed Ritz values
	//    and vectors of OP with respect to the given problem  A*z = lambda*B*z.
	//
	//    The approximate eigenvalues of the original problem are returned in
	//    ascending algebraic order.  The user may elect to call this routine
	//    once for each desired Ritz vector and store it peripherally if desired.
	//    There is also the option of computing a selected set of these vectors
	//    with a single call.
//	fprintf(stderr,"\nBEFORE: rvec %d, howny %s, select %016lx, d %016lx, z %016lx, ldz %d, sigma %f, bmat %s, nrows %d, which %s, nev %d, tol %f, resid %016lx, ncv %d, v %016lx, ldv %d, iparam %016lx, ipntr %016lx, workd %016lx, workl %016lx, lworkl %d, info %d\n",
//			        rvec, howny, (unsigned long)select, (unsigned long)d, (unsigned long)z, ldz, sigma, bmat, nrows, which, nev, tol, (unsigned long)resid, ncv, (unsigned long)v, ldv, (unsigned long)iparam, (unsigned long)ipntr, (unsigned long)workd, (unsigned long)workl, lworkl, info);
    dseupd_(&rvec, howny, select, d, z, &ldz, &sigma, bmat, &nrows, which,
			&nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, &info);
//	fprintf(stderr,"\nAFTER: rvec %d, howny %s, select %016lx, d %016lx, z %016lx, ldz %d, sigma %f, bmat %s, nrows %d, which %s, nev %d, tol %f, resid %016lx, ncv %d, v %016lx, ldv %d, iparam %016lx, ipntr %016lx, workd %016lx, workl %016lx, lworkl %d, info %d\n",
//			        rvec, howny, (unsigned long)select, (unsigned long)d, (unsigned long)z, ldz, sigma, bmat, nrows, which, nev, tol, (unsigned long)resid, ncv, (unsigned long)v, ldv, (unsigned long)iparam, (unsigned long)ipntr, (unsigned long)workd, (unsigned long)workl, lworkl, info);

    // Extract eigenvalues and eigenvectors
//	printf("Extract eigenvalues and eigenvectors\n");
//	for(int i = 0; i<nev; i++)
//		printf ("ev= %d  = %f\n",i,d[i]);

	// MON: Revert eigenvector order... by swapping.
	// From the LARGEST EIGENVALUES/VECTORS OF (B*x=sigma*A*x)
	// to SMALLEST EIGENVALUES/VECTORS OF (A*x=lamba*B*x)
//		printf("Using DTPTRS_\n");
		info = 0;
		dtptrs_("U", "N", "N", &nrows, &nev, B, z, &nrows, &info);

//		printf("Reverting into standard order by swapping\n");
		double swap;
		for (i = 0; i < nev/2; i++)
		{
			swap = d[i];
			d[i] = 1/d[nev-1-i];
			d[nev-1-i] = 1/swap;
			for(j = 0; j<nrows; j++)
			{
//				swap = z[i][j];
				swap = z[i*nrows + j];
//				z[i][j] = z[nev-1-i][j] * sqrt(d[i]);
				z[i*nrows + j] = z[(nev-1-i)*nrows + j] * sqrt(d[i]);
//				z[nev-1-i][j] = swap * sqrt(d[nev-1-i]);
				z[(nev-1-i)*nrows + j] = swap * sqrt(d[nev-1-i]);
			}
		}
		if(nev % 2 == 1)
		{
			i=nev/2;
			d[i] = 1/d[i];
			for(j = 0; j<nrows; j++)
//				z[i][j] = z[i][j] * sqrt(d[i]);
				z[i*nrows + j] *= sqrt(d[i]);
		}

		//	// MON: Saving eigenvectors & eigenvalues (in ptraj)
		// Saves & Normalizes eigenpairs according to "ptraj" format
//		save_ptraj_modes (nomResult, nrows, 0, nev, d, z,0);
//		printf("Write into %s:  %d rows  %d neigs\n",nomResult,nrows,nev);
		// exit(0);

	// printf("REVERSE_STEPS = %d   NCV= %d\n", nsteps, ncv);

	free(v);
	free(select);
	free(workl);
	free(workd);
	free(resid);
	free(xy);
}

// ARPACK-based diagonalization routine (from: J.I.Aliaga & E.Quintana)
// Very fast if small number of eigenpairs requested. It scales linear with number of threads.
// Although matrices should be provided in packed storage, optimized for packed storage, thus it does not
// A --> Hessian matrix (squared, no-packed)
// B --> Kinetic energy matrix (squared, no-packed)
// d --> Eigenvalues array
// z --> Eigenvectors array
// nrows --> Problem dimension (size)
// nev --> Number of eigenpairs to be computed
// ncv --> Factor for the Number of columns of matrix V (=ncv*nev). It should be optimized for a given problem (OPTIONAL)
void dsdrv1_A_B_W_mon(double *B, double *A, double *d, double *z, int nrows, int nev, float ncvf)
{
	int i,j, inc = 1, ido = 0, ldv = nrows, info = 0, nsteps = 0;
	char *bmat = "I";
	//  WHICH   Character*2.  (INPUT)
	//		  Specify which of the Ritz values of OP to compute.
	//		  'LM' - compute the NEV largest (in magnitude) eigenvalues.
	char *which = "LM";
	//  IPARAM(3) = MXITER
	//  On INPUT:  maximum number of Arnoldi update iterations allowed.
	//  On OUTPUT: actual number of Arnoldi update iterations taken.
	int maxit = 50; //

	// Extracted from DSAUPD remarks:
	//	  4. At present there is no a-priori analysis to guide the selection
	//	     of NCV relative to NEV.  The only formal requirement is that NCV > NEV.
	//	     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
	//	     the same type are to be solved, one should experiment with increasing
	//	     NCV while keeping NEV fixed for a given test problem.  This will
	//	     usually decrease the required number of OP*x operations but it
	//	     also increases the work and storage required to maintain the orthogonal
	//	     basis vectors.   The optimal "cross-over" with respect to CPU time
	//	     is problem dependent and must be determined empirically.
	int ncv = (int) (ncvf * (float) nev);

	int iparam[11], ipntr[11], lworkl = (ncv * ncv + 8 * ncv);
	double tol = 0.0, alpha = 1.0, beta = 0.0;
	double *resid = NULL, *v = NULL, *workd = NULL, *workl = NULL, *x = NULL, *y = NULL, *xy = NULL;
	char *howny = "A";
	int rvec = 1, *select = NULL, ldz = nrows;
	double sigma = 0.0;

	// Compute Cholesky factorization of B
	info = 0;
	dpotrf_("U", &nrows, B, &nrows, &info);

	// Create ARPACK structures for ARPACK loop
	ldv = nrows;

//	v = create_rect_matrix(ncv, ldv);
//	resid = create_vector(nrows);
//	workd = create_vector(3 * nrows);
//	workl = create_vector(lworkl);
//	xy = create_vector(nrows);

	if ( !(v = (double *) malloc (ncv * ldv * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

	if ( !(resid = (double *) malloc (nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

	if ( !(workd = (double *) malloc (3 * nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

	if ( !(workl = (double *) malloc (lworkl * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}

	if ( !(xy = (double *) malloc (nrows * sizeof(double)) ) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}


	for (i = 0; i < 11; i++)
	{
		iparam[i] = 0;
		ipntr[i] = 0;
	}
	iparam[0] = 1;
	iparam[2] = maxit;
	iparam[6] = 1;

	// ARPACK loop
	nsteps = 0;
	dsaupd_(&ido, bmat, &nrows, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	while ((ido == 1) || (ido == -1))
	{
		x = workd + ipntr[0] - 1;
		y = workd + ipntr[1] - 1;
		nsteps++;
		dcopy_(&nrows, x, &inc, xy, &inc); // memcpy (xy, x, nrows*sizeof(double));
		dtrtrs_("U", "N", "N", &nrows, &inc, B, &nrows, xy, &nrows, &info);
		dgemv_("N", &nrows, &nrows, &alpha, A, &nrows, xy, &inc, &beta, y, &inc); // (24'73, 27'97)
		dtrtrs_("U", "T", "N", &nrows, &inc, B, &nrows, y, &nrows, &info);
		dsaupd_(&ido, bmat, &nrows, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
	}
	// printf("IDO = %d , INFO = %d \n", ido, info);

	// Create ARPACK structures to extract eigenvalues and eigenvectors
//	select = create_indices(ncv);
	if( !(select = (int *) malloc (ncv * sizeof(int))) ) {
		printf ("Memory allocation error\n"); exit(-1);
	}
	for (i = 0; i < ncv; i++)
		select[i] = 1;
//	d = create_vector(nev);
	ldz = nrows;
//	z = create_rect_matrix(nev, nrows);

	// Extract eigenvalues and eigenvectors
	dseupd_(&rvec, howny, select, d, z, &ldz, &sigma, bmat, &nrows, which,
			&nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl,
			&lworkl, &info);

	// MON: Revert eigenvector order... by swapping.
	// From the LARGEST EIGENVALUES/VECTORS OF (B*x=sigma*A*x)
	// to SMALLEST EIGENVALUES/VECTORS OF (A*x=lamba*B*x)

		// DTRTRS solves a triangular system of the form
		//*     A * X = B  or  A**T * X = B,
		//*  where A is a triangular matrix of order N, and B is an N-by-NRHS
		//*  matrix.  A check is made to verify that A is nonsingular.
		// printf("Using DTPTRS_\n");
		info = 0;
		dtrtrs_("U", "N", "N", &nrows, &nev, B, &nrows, z, &nrows, &info);

		// printf("Reverting into standard order by swapping\n");
		double swap;
		for (i = 0; i < nev/2; i++)
		{
			swap = d[i];
			d[i] = 1/d[nev-1-i];
			d[nev-1-i] = 1/swap;
			for(j = 0; j<nrows; j++)
			{
				swap = z[i*nrows + j];
				z[i*nrows + j] = z[(nev-1-i)*nrows + j] * sqrt(d[i]);
				z[(nev-1-i)*nrows + j] = swap * sqrt(d[nev-1-i]);
			}
		}
		if(nev % 2 == 1)
		{
			i=nev/2;
			d[i] = 1/d[i];
			for(j = 0; j<nrows; j++)
				z[i*nrows + j] *= sqrt(d[i]);
		}

//		//	// MON: Saving eigenvectors & eigenvalues (in ptraj)
//		// Saves & Normalizes eigenpairs according to "ptraj" format
//		save_ptraj_modes (nomResult, nrows, 0, nev, d, z,0);
//		printf("Write into %s:  %d rows  %d neigs\n",nomResult,nrows,nev);
//		// exit(0);

//	// Print results
//	x = (double *) malloc(nrows * sizeof(double));
//	for (i = 0; i < nev; i++) {
//		alpha = 1.0, beta = 0.0;
//		dcopy_(&nrows, z[i], &inc, x, &inc); // memcpy (x, z[i], nrows*sizeof(double));
//		dtrtrs_("U", "N", "N", &nrows, &inc, B[0], &nrows, x, &nrows, &info);
//		dgemv_("N", &nrows, &nrows, &alpha, A[0], &nrows, x, &inc, &beta, xy,
//				&inc);
//		dtrtrs_("U", "T", "N", &nrows, &inc, B[0], &nrows, xy, &nrows, &info);
//		alpha = -d[i];
//		daxpy_(&nrows, &alpha, z[i], &inc, xy, &inc);
//		beta = ddot_(&nrows, xy, &inc, xy, &inc);
//		beta = dnrm2_(&nrows, xy, &inc);
//		tol = dnrm2_(&nrows, z[i], &inc);
//		printf(
//				"EIG(%d) = %10.5e , NORM(%d) = %10.5e , ERROR(%d) = (%10.5e/%10.5e) = %10.5g\n",
//				i, d[i], i, tol, i, beta, norm, beta / norm);
//	}
//	free(x);
//	x = NULL;

	// printf("REVERSE_STEPS = %d\n", nsteps);

//	// Write eigenvalues and eigenvectors
//	if (nomResult != NULL) {
//		info = 0;
//		dtrtrs_("U", "N", "N", &nrows, &nev, B[0], &nrows, z[0], &nrows, &info);
//		write_eigs(nomResult, z, nrows, d, nev);
//	}

	free(select);
	free(xy);
	free(workl);
	free(workd);
	free(resid);
	free(v);
}

// DSPGVD  computes  ALL the eigenvalues, and optionally, the eigenvectors
// of a real generalized  symmetric-definite  eigenproblem,  of  the  form
// A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
// are assumed to be symmetric, stored in PACKED FORMAT,  and  B  is  also
// positive definite.
// If eigenvectors are desired, it uses a divide and conquer algorithm.
// (BIG MEMORY REQUIREMENTS: If JOBZ = 'V' & N > 1, LWORK >= 1 + 6*N + 2*N**2.)
void diag_dspgvd(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size)
{
  bool debug = false;
  /* Lapack input/output */
  char jobz, uplo;
  int  itype, ldz, lwork, liwork, info, *iwork;
  double *work;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  uplo='U';  /* Upper triangle of A is stored */
  ldz=size;
  lwork = 1; // query
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++) work[i]=0.0;
  lwork = -1; // query

  liwork = 1; // query
  iwork = (int *) malloc(liwork * sizeof(int));
  for(int i=0; i<liwork; i++) iwork[i]=0;
  liwork = -1;

  // Memory Query
  dspgvd_(&itype, &jobz, &uplo, &size, hess_matrix, mass_matrix,
		  eigval, eigvec, &ldz, work, &lwork, iwork,
          &liwork, &info);
  if(info != 0)
  {
	  printf("Msg(diag_dspgvd): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }
  lwork = (int) work[0];
  liwork = iwork[0];
  if(debug)
	  printf("\nMsg(diag_dspgvd):  Optimal LWORK= %d (%d Mb) (std_LWORK=%d) and LIWORK= %d\n",
		  lwork,lwork*sizeof(double)/1000000,2*size*size+6*size+1,liwork);

  // Allocating LAPACK working memory
  free(work);
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++)
	  work[i]=0.0;
  free(iwork);
  iwork = (int *) malloc(liwork * sizeof(int));
  for(int i=0; i<liwork; i++)
	  iwork[i]=0;

  // Diagonalization
  dspgvd_(&itype, &jobz, &uplo, &size, hess_matrix, mass_matrix,
		  eigval, eigvec, &ldz, work, &lwork, iwork,
          &liwork, &info);
  if(info != 0)
  {
	  printf("Msg(diag_dspgvd): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }

  free(work); // <-- "work" could be allocated only once outside!!!
  free(iwork); // <-- "iwork" could be allocated only once outside!!!
}
//DSPGVD(3)             LAPACK driver routine (version 3.1)            DSPGVD(3)
//
//
//
//NAME
//       DSPGVD  -  all  the  eigenvalues, and optionally, the eigenvectors of a
//       real  generalized  symmetric-definite   eigenproblem,   of   the   form
//       A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x
//
//SYNOPSIS
//       SUBROUTINE DSPGVD( ITYPE,  JOBZ,  UPLO,  N,  AP,  BP,  W, Z, LDZ, WORK,
//                          LWORK, IWORK, LIWORK, INFO )
//
//           CHARACTER      JOBZ, UPLO
//
//           INTEGER        INFO, ITYPE, LDZ, LIWORK, LWORK, N
//
//           INTEGER        IWORK( * )
//
//           DOUBLE         PRECISION AP( * ), BP( * ), W( * ), WORK(  *  ),  Z(
//                          LDZ, * )
//
//PURPOSE
//       DSPGVD  computes  all the eigenvalues, and optionally, the eigenvectors
//       of a real generalized  symmetric-definite  eigenproblem,  of  the  form
//       A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
//       are assumed to be symmetric, stored in packed format,  and  B  is  also
//       positive definite.
//       If eigenvectors are desired, it uses a divide and conquer algorithm.
//
//       The  divide  and  conquer  algorithm  makes very mild assumptions about
//       floating point arithmetic. It will work on machines with a guard  digit
//       in add/subtract, or on those binary machines without guard digits which
//       subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2. It  could
//       conceivably  fail on hexadecimal or decimal machines without guard dig-
//       its, but we know of none.
//
//
//ARGUMENTS
//       ITYPE   (input) INTEGER
//               Specifies the problem type to be solved:
//               = 1:  A*x = (lambda)*B*x
//               = 2:  A*B*x = (lambda)*x
//               = 3:  B*A*x = (lambda)*x
//
//       JOBZ    (input) CHARACTER*1
//               = 'N':  Compute eigenvalues only;
//               = 'V':  Compute eigenvalues and eigenvectors.
//
//       UPLO    (input) CHARACTER*1
//               = 'U':  Upper triangles of A and B are stored;
//               = 'L':  Lower triangles of A and B are stored.
//
//       N       (input) INTEGER
//               The order of the matrices A and B.  N >= 0.
//
//       AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//               On entry, the upper or lower triangle of the  symmetric  matrix
//               A,  packed  columnwise in a linear array.  The j-th column of A
//               is stored in the array AP as follows: if UPLO  =  'U',  AP(i  +
//               (j-1)*j/2)  =  A(i,j)  for  1<=i<=j;  if  UPLO  =  'L',  AP(i +
//               (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//               On exit, the contents of AP are destroyed.
//
//       BP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//               On entry, the upper or lower triangle of the  symmetric  matrix
//               B,  packed  columnwise in a linear array.  The j-th column of B
//               is stored in the array BP as follows: if UPLO  =  'U',  BP(i  +
//               (j-1)*j/2)  =  B(i,j)  for  1<=i<=j;  if  UPLO  =  'L',  BP(i +
//               (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
//
//               On exit, the triangular factor U or L from the Cholesky factor-
//               ization B = U**T*U or B = L*L**T, in the same storage format as
//               B.
//
//       W       (output) DOUBLE PRECISION array, dimension (N)
//               If INFO = 0, the eigenvalues in ascending order.
//
//       Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
//               If JOBZ = 'V', then if INFO = 0, Z contains  the  matrix  Z  of
//               eigenvectors.   The  eigenvectors are normalized as follows: if
//               ITYPE = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z =  I.
//               If JOBZ = 'N', then Z is not referenced.
//
//       LDZ     (input) INTEGER
//               The  leading dimension of the array Z.  LDZ >= 1, and if JOBZ =
//               'V', LDZ >= max(1,N).
//
//       WORK      (workspace/output)   DOUBLE   PRECISION   array,    dimension
//       (MAX(1,LWORK))
//               On exit, if INFO = 0, WORK(1) returns the required LWORK.
//
//       LWORK   (input) INTEGER
//               The   dimension   of   the   array   WORK.    If   N   <=    1,
//               LWORK  >= 1.  If JOBZ = 'N' and N > 1, LWORK >= 2*N.  If JOBZ =
//               'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
//
//               If LWORK = -1, then a workspace query is assumed;  the  routine
//               only  calculates  the  required  sizes  of  the  WORK and IWORK
//               arrays, returns these values as the first entries of  the  WORK
//               and  IWORK  arrays,  and  no  error message related to LWORK or
//               LIWORK is issued by XERBLA.
//
//       IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
//               On exit, if INFO = 0, IWORK(1) returns the required LIWORK.
//
//       LIWORK  (input) INTEGER
//               The dimension of the array IWORK.  If JOBZ  = 'N' or  N  <=  1,
//               LIWORK >= 1.  If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
//
//               If  LIWORK = -1, then a workspace query is assumed; the routine
//               only calculates the  required  sizes  of  the  WORK  and  IWORK
//               arrays,  returns  these values as the first entries of the WORK
//               and IWORK arrays, and no error  message  related  to  LWORK  or
//               LIWORK is issued by XERBLA.
//
//       INFO    (output) INTEGER
//               = 0:  successful exit
//               < 0:  if INFO = -i, the i-th argument had an illegal value
//               > 0:  DPPTRF or DSPEVD returned an error code:
//               <=  N:   if INFO = i, DSPEVD failed to converge; i off-diagonal
//               elements of an intermediate tridiagonal form did  not  converge
//               to  zero;  >  N:    if  INFO = N + i, for 1 <= i <= N, then the
//               leading minor of order i of B is not  positive  definite.   The
//               factorization of B could not be completed and no eigenvalues or
//               eigenvectors were computed.
//
//FURTHER DETAILS
//       Based on contributions by
//          Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA



// DSPGV computes ALL the eigenvalues and, optionally, the eigenvectors of
// a  real  generalized  symmetric-definite  eigenproblem,  of  the   form
// A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
// are assumed to be symmetric, stored in packed format,  and  B  is  also
// positive definite.
// (MINIMAL MEMORY REQUIREMENTS)
void diag_dspgv(double *hess_matrix,double *mass_matrix,double *eigval, double *eigvec, int size)
{
  /* Lapack input/output */
  char jobz, uplo;
  int  itype, ldz, lwork, liwork, info;
  double *work;
  bool debug=true;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  uplo='U';  /* Upper triangle of A is stored */
  ldz=size;
  lwork = 3*size;
  work = (double *) malloc(lwork * sizeof(double));
  for(int i=0; i<lwork; i++) work[i]=0.0;
  if(debug) printf("\nMsg(diag_dspgv):  Constant LWORK= %d (%d Mb) \n",lwork,lwork*sizeof(double)/1000000);

  // Diagonalization
  dspgv_(&itype, &jobz, &uplo, &size, hess_matrix, mass_matrix, eigval, eigvec, &ldz, work, &info);
  if(info != 0)
  {
	  printf("Msg(diag_dspgv): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }

  free(work); // <-- "work" could be allocated only once outside!!!
}
//DSPGV(3)              LAPACK driver routine (version 3.1)             DSPGV(3)
//
//
//
//NAME
//       DSPGV - all the eigenvalues and, optionally, the eigenvectors of a real
//       generalized    symmetric-definite    eigenproblem,    of    the    form
//       A*x=(lambda)*B*x, A*Bx=(lambda)*x, or B*A*x=(lambda)*x
//
//SYNOPSIS
//       SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, INFO )
//
//           CHARACTER     JOBZ, UPLO
//
//           INTEGER       INFO, ITYPE, LDZ, N
//
//           DOUBLE        PRECISION AP( * ), BP( * ), W( * ),  WORK(  *  ),  Z(
//                         LDZ, * )
//
//PURPOSE
//       DSPGV computes all the eigenvalues and, optionally, the eigenvectors of
//       a  real  generalized  symmetric-definite  eigenproblem,  of  the   form
//       A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and B
//       are assumed to be symmetric, stored in packed format,  and  B  is  also
//       positive definite.
//
//
//ARGUMENTS
//       ITYPE   (input) INTEGER
//               Specifies the problem type to be solved:
//               = 1:  A*x = (lambda)*B*x
//               = 2:  A*B*x = (lambda)*x
//               = 3:  B*A*x = (lambda)*x
//
//       JOBZ    (input) CHARACTER*1
//               = 'N':  Compute eigenvalues only;
//               = 'V':  Compute eigenvalues and eigenvectors.
//
//       UPLO    (input) CHARACTER*1
//               = 'U':  Upper triangles of A and B are stored;
//               = 'L':  Lower triangles of A and B are stored.
//
//       N       (input) INTEGER
//               The order of the matrices A and B.  N >= 0.
//
//       AP      (input/output) DOUBLE PRECISION array, dimension
//               (N*(N+1)/2)  On  entry, the upper or lower triangle of the sym-
//               metric matrix A, packed columnwise in a linear array.  The j-th
//               column  of  A  is  stored in the array AP as follows: if UPLO =
//               'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;  if  UPLO  =  'L',
//               AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//
//               On exit, the contents of AP are destroyed.
//
//       BP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//               On  entry,  the upper or lower triangle of the symmetric matrix
//               B, packed columnwise in a linear array.  The j-th column  of  B
//               is  stored  in  the  array BP as follows: if UPLO = 'U', BP(i +
//               (j-1)*j/2) =  B(i,j)  for  1<=i<=j;  if  UPLO  =  'L',  BP(i  +
//               (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
//
//               On exit, the triangular factor U or L from the Cholesky factor-
//               ization B = U**T*U or B = L*L**T, in the same storage format as
//               B.
//
//       W       (output) DOUBLE PRECISION array, dimension (N)
//               If INFO = 0, the eigenvalues in ascending order.
//
//       Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
//               If  JOBZ  =  'V',  then if INFO = 0, Z contains the matrix Z of
//               eigenvectors.  The eigenvectors are normalized as  follows:  if
//               ITYPE  = 1 or 2, Z**T*B*Z = I; if ITYPE = 3, Z**T*inv(B)*Z = I.
//               If JOBZ = 'N', then Z is not referenced.
//
//       LDZ     (input) INTEGER
//               The leading dimension of the array Z.  LDZ >= 1, and if JOBZ  =
//               'V', LDZ >= max(1,N).
//
//       WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
//
//       INFO    (output) INTEGER
//               = 0:  successful exit
//               < 0:  if INFO = -i, the i-th argument had an illegal value
//               > 0:  DPPTRF or DSPEV returned an error code:
//               <=  N:   if  INFO = i, DSPEV failed to converge; i off-diagonal
//               elements of an intermediate tridiagonal form did  not  converge
//               to  zero.   >  N:    if INFO = n + i, for 1 <= i <= n, then the
//               leading minor of order i of B is not  positive  definite.   The
//               factorization of B could not be completed and no eigenvalues or
//               eigenvectors were computed.



//*  SSPGV computes ALL the eigenvalues and, optionally, the eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
//*  Here A and B are assumed to be symmetric, stored in packed format,
//*  and B is also positive definite.
//   (MINIMAL MEMORY REQUIREMENTS)
void diag_sspgv(float *hess_matrix,float *mass_matrix,float *eigval,float *eigvec, int size)
{
  /* Lapack input/output */
  char jobz, uplo;
  int  itype, ldz, lwork, liwork, info;
  float *work;
  bool debug=true;

  itype=1;   /* specifies A*x = lambda*B*x */
  jobz='V';  /* Compute eigenvalues and eigenvectors */
  uplo='U';  /* Upper triangle of A is stored */
  ldz=size;
  lwork = 3*size;
  work = (float *) malloc(lwork * sizeof(float));
  for(int i=0; i<lwork; i++) work[i]=0.0;
  if(debug) printf("\nMsg(diag_sspgv):  Constant LWORK= %d (%d Mb) \n",lwork,lwork*sizeof(float)/1000000);

  // Diagonalization
  sspgv_(&itype, &jobz, &uplo, &size, hess_matrix, mass_matrix, eigval, eigvec, &ldz, work, &info);
  if(info != 0)
  {
	  printf("Msg(diag_sspgv): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
	  exit(3);
  }

  free(work); // <-- "work" could be allocated only once outside!!!
}
//SUBROUTINE SSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,
//$                  INFO )
//*
//*  -- LAPACK driver routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, UPLO
//INTEGER            INFO, ITYPE, LDZ, N
//*     ..
//*     .. Array Arguments ..
//REAL               AP( * ), BP( * ), W( * ), WORK( * ),
//$                   Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  SSPGV computes all the eigenvalues and, optionally, the eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
//*  Here A and B are assumed to be symmetric, stored in packed format,
//*  and B is also positive definite.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          Specifies the problem type to be solved:
//*          = 1:  A*x = (lambda)*B*x
//*          = 2:  A*B*x = (lambda)*x
//*          = 3:  B*A*x = (lambda)*x
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangles of A and B are stored;
//*          = 'L':  Lower triangles of A and B are stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrices A and B.  N >= 0.
//*
//*  AP      (input/output) REAL array, dimension
//*                            (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, the contents of AP are destroyed.
//*
//*  BP      (input/output) REAL array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          B, packed columnwise in a linear array.  The j-th column of B
//*          is stored in the array BP as follows:
//*          if UPLO = 'U', BP(i + (j-1)*j/2) = B(i,j) for 1<=i<=j;
//*          if UPLO = 'L', BP(i + (j-1)*(2*n-j)/2) = B(i,j) for j<=i<=n.
//*
//*          On exit, the triangular factor U or L from the Cholesky
//*          factorization B = U**T*U or B = L*L**T, in the same storage
//*          format as B.
//*
//*  W       (output) REAL array, dimension (N)
//*          If INFO = 0, the eigenvalues in ascending order.
//*
//*  Z       (output) REAL array, dimension (LDZ, N)
//*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
//*          eigenvectors.  The eigenvectors are normalized as follows:
//*          if ITYPE = 1 or 2, Z**T*B*Z = I;
//*          if ITYPE = 3, Z**T*inv(B)*Z = I.
//*          If JOBZ = 'N', then Z is not referenced.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace) REAL array, dimension (3*N)
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  SPPTRF or SSPEV returned an error code:
//*             <= N:  if INFO = i, SSPEV failed to converge;
//*                    i off-diagonal elements of an intermediate
//*                    tridiagonal form did not converge to zero.
//*             > N:   if INFO = n + i, for 1 <= i <= n, then the leading
//*                    minor of order i of B is not positive definite.
//*                    The factorization of B could not be completed and
//*                    no eigenvalues or eigenvectors were computed.
//*
//*  =====================================================================




// INVERSE MATRIX COMPUTATION sub-rutine (Double precision, double)
// Uses LAPACK's dsptri_. (SYMMETRIC, PACKED STORAGE, MINIMAL MEMORY USAGE)
// Uses LAPACK's dsptrf_. (Triangular factorization: LU, Cholesky)
//*  DSPTRI computes the inverse of a real symmetric indefinite matrix
//*  A in packed storage using the factorization A = U*D*U**T or
//*  A = L*D*L**T computed by DSPTRF.
void inv_dsptri(double *ap, int size)
{
	char uplo='U';  // Upper triangle of A is stored
	int *ipiv;
	int info;
	double *matrix; // temporal matrix storage
	double *work; // working array

//	matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 ); // (DOUBLE) array, dimension (N*(N+1)/2)
//	for(int i=0; i < size*(size+1)/2; i++)
//		matrix[i] = ap[i]; // needed to avoid AP-matrix lost!

	ipiv = (int *) malloc( sizeof(int) * size );
	//SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
//	dsptrf_(&uplo, &size, matrix, ipiv, &info); // computes IPIV
	dsptrf_(&uplo, &size, ap, ipiv, &info); // computes IPIV
	if(info != 0)
	{
		printf("Msg(lapack's DSPTRF): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
//	free( matrix );

	work = (double *) malloc( sizeof(double) * size );
	//SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
	dsptri_(&uplo, &size, ap, ipiv, work, &info);
	if(info != 0)
	{
		printf("Msg(lapack's DSPTRI): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
	free( ipiv );
	free( work );
}
//SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
//*
//*  -- LAPACK routine (version 3.0) --
//*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//*     Courant Institute, Argonne National Lab, and Rice University
//*     March 31, 1993
//*
//*     .. Scalar Arguments ..
// CHARACTER          UPLO
// INTEGER            INFO, N
//*     ..
//*     .. Array Arguments ..
// INTEGER            IPIV( * )
// DOUBLE PRECISION   AP( * ), WORK( * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSPTRI computes the inverse of a real symmetric indefinite matrix
//*  A in packed storage using the factorization A = U*D*U**T or
//*  A = L*D*L**T computed by DSPTRF.
//*
//*  Arguments
//*  =========
//*
//*  UPLO    (input) CHARACTER*1
//*          Specifies whether the details of the factorization are stored
//*          as an upper or lower triangular matrix.
//*          = 'U':  Upper triangular, form is A = U*D*U**T;
//*          = 'L':  Lower triangular, form is A = L*D*L**T.
//*
//*  N       (input) INTEGER
//*          The order of the matrix A.  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the block diagonal matrix D and the multipliers
//*          used to obtain the factor U or L as computed by DSPTRF,
//*          stored as a packed triangular matrix.
//*
//*          On exit, if INFO = 0, the (symmetric) inverse of the original
//*          matrix, stored as a packed triangular matrix. The j-th column
//*          of inv(A) is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
//*
//*  IPIV    (input) INTEGER array, dimension (N)
//*          Details of the interchanges and the block structure of D
//*          as determined by DSPTRF.
//*
//*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
//*
//*  INFO    (output) INTEGER
//*          = 0: successful exit
//*          < 0: if INFO = -i, the i-th argument had an illegal value
//*          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
//*               inverse could not be computed.
//*
//*  =====================================================================


//SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
//*
//*  -- LAPACK routine (version 3.0) --
//*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//*     Courant Institute, Argonne National Lab, and Rice University
//*     June 30, 1999
//*
//*     .. Scalar Arguments ..
//CHARACTER          UPLO
//INTEGER            INFO, N
//*     ..
//*     .. Array Arguments ..
//INTEGER            IPIV( * )
//DOUBLE PRECISION   AP( * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSPTRF computes the factorization of a real symmetric matrix A stored
//*  in packed format using the Bunch-Kaufman diagonal pivoting method:
//*
//*     A = U*D*U**T  or  A = L*D*L**T
//*
//*  where U (or L) is a product of permutation and unit upper (lower)
//*  triangular matrices, and D is symmetric and block diagonal with
//*  1-by-1 and 2-by-2 diagonal blocks.
//*
//*  Arguments
//*  =========
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A is stored;
//*          = 'L':  Lower triangle of A is stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix A.  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, the block diagonal matrix D and the multipliers used
//*          to obtain the factor U or L, stored as a packed triangular
//*          matrix overwriting A (see below for further details).
//*
//*  IPIV    (output) INTEGER array, dimension (N)
//*          Details of the interchanges and the block structure of D.
//*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
//*          interchanged and D(k,k) is a 1-by-1 diagonal block.
//*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
//*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
//*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
//*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
//*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
//*
//*  INFO    (output) INTEGER
//*          = 0: successful exit
//*          < 0: if INFO = -i, the i-th argument had an illegal value
//*          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
//*               has been completed, but the block diagonal matrix D is
//*               exactly singular, and division by zero will occur if it
//*               is used to solve a system of equations.
//*
//*  Further Details
//*  ===============
//*
//*  5-96 - Based on modifications by J. Lewis, Boeing Computer Services
//*         Company
//*
//*  If UPLO = 'U', then A = U*D*U', where
//*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
//*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
//*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
//*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
//*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
//*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
//*
//*             (   I    v    0   )   k-s
//*     U(k) =  (   0    I    0   )   s
//*             (   0    0    I   )   n-k
//*                k-s   s   n-k
//*
//*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
//*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
//*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
//*
//*  If UPLO = 'L', then A = L*D*L', where
//*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
//*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
//*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
//*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
//*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
//*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
//*
//*             (   I    0     0   )  k-1
//*     L(k) =  (   0    I     0   )  s
//*             (   0    v     I   )  n-k-s+1
//*                k-1   s  n-k-s+1
//*
//*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
//*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
//*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
//*
//*  =====================================================================



// INVERSE MATRIX COMPUTATION sub-rutine (Double precision, double)
// Uses LAPACK's dtptri_. (TRIANGULAR, PACKED STORAGE, MINIMAL MEMORY USAGE)
// DSPTRI computes the inverse of a real symmetric indefinite matrix
void inv_dtptri(double *ap, int size)
{
	char uplo='U';  // Upper triangle of A is stored
	char diag='N'; // Non-unit triangular (the diagonal elements are != 1)
	int info;

	// SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO )
	dtptri_(&uplo, &diag, &size, ap, &info);

	if(info != 0)
	{
		printf("Msg(lapack's DTPTRI): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
}
//SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO )
//*
//*  -- LAPACK routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
// CHARACTER          DIAG, UPLO
// INTEGER            INFO, N
//*     ..
//*     .. Array Arguments ..
// DOUBLE PRECISION   AP( * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DTPTRI computes the inverse of a real upper or lower triangular
//*  matrix A stored in packed format.
//*
//*  Arguments
//*  =========
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  A is upper triangular;
//*          = 'L':  A is lower triangular.
//*
//*  DIAG    (input) CHARACTER*1
//*          = 'N':  A is non-unit triangular;
//*          = 'U':  A is unit triangular.
//*
//*  N       (input) INTEGER
//*          The order of the matrix A.  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangular matrix A, stored
//*          columnwise in a linear array.  The j-th column of A is stored
//*          in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n.
//*          See below for further details.
//*          On exit, the (triangular) inverse of the original matrix, in
//*          the same packed storage format.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular
//*                matrix is singular and its inverse can not be computed.
//*
//*  Further Details
//*  ===============
//*
//*  A triangular matrix A can be transferred to packed storage using one
//*  of the following program segments:
//*
//*  UPLO = 'U':                      UPLO = 'L':
//*
//*        JC = 1                           JC = 1
//*        DO 2 J = 1, N                    DO 2 J = 1, N
//*           DO 1 I = 1, J                    DO 1 I = J, N
//*              AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)
//*      1    CONTINUE                    1    CONTINUE
//*           JC = JC + J                      JC = JC + N - J + 1
//*      2 CONTINUE                       2 CONTINUE
//*
//*  =====================================================================



//*  Reduces a real symmetric-definite generalized eigenproblem
//*  to standard form, using packed storage.
// Uses LAPACK's dspgst_. (SYMMETRIC, PACKED STORAGE, MINIMAL MEMORY USAGE)
// Uses LAPACK's dsptrf_. (Triangular factorization: LU, Cholesky)
int reduce_dspgst(double *ap, double *bp, int size)
{
	char uplo='U';  // Upper triangle of A is stored
	int *ipiv;
	int info;
	int itype=1; // If ITYPE = 1, the problem is A*x = lambda*B*x,

	ipiv = (int *) malloc( sizeof(int) * size );
	//SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
	dsptrf_(&uplo, &size, bp, ipiv, &info); // computes IPIV and BP Cholesky factorization
	if(info != 0)
	{
		printf("Msg(lapack's DSPTRF): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
	free( ipiv );

	//SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
	dspgst_(&itype, &uplo, &size, ap, bp, &info);
	if(info != 0)
	{
		printf("Msg(lapack's DSPGST): I'm sorry, INFO=%d, please check LAPACK documentation!\n",info);
		exit(3);
	}
	return 0;
}
//SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
//*
//*  -- LAPACK routine (version 3.0) --
//*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
//*     Courant Institute, Argonne National Lab, and Rice University
//*     March 31, 1993
//*
//*     .. Scalar Arguments ..
//CHARACTER          UPLO
//INTEGER            INFO, ITYPE, N
//*     ..
//*     .. Array Arguments ..
//DOUBLE PRECISION   AP( * ), BP( * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSPGST reduces a real symmetric-definite generalized eigenproblem
//*  to standard form, using packed storage.
//*
//*  If ITYPE = 1, the problem is A*x = lambda*B*x,
//*  and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)
//*
//*  If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
//*  B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.
//*
//*  B must have been previously factorized as U**T*U or L*L**T by DPPTRF.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);
//*          = 2 or 3: compute U*A*U**T or L**T*A*L.
//*
//*  UPLO    (input) CHARACTER
//*          = 'U':  Upper triangle of A is stored and B is factored as
//*                  U**T*U;
//*          = 'L':  Lower triangle of A is stored and B is factored as
//*                  L*L**T.
//*
//*  N       (input) INTEGER
//*          The order of the matrices A and B.  N >= 0.
//*
//*  AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          On entry, the upper or lower triangle of the symmetric matrix
//*          A, packed columnwise in a linear array.  The j-th column of A
//*          is stored in the array AP as follows:
//*          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
//*          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
//*
//*          On exit, if INFO = 0, the transformed matrix, stored in the
//*          same format as A.
//*
//*  BP      (input) DOUBLE PRECISION array, dimension (N*(N+1)/2)
//*          The triangular factor from the Cholesky factorization of B,
//*          stored in the same format as A, as returned by DPPTRF.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*
//*  =====================================================================
//



// Test to state that a matrix is positive definite (like the Kinetic Energy Matrix must be)
int is_posdefinite(double *matrix, int size)
{
	return 0;
}
// http://www.ece.uwaterloo.ca/~ece204/TheBook/04LinearAlgebra/posdef/complete.html
// A symmetric matrix is positive definite if:
//   1. all the diagonal entries are positive, and
//   2. each diagonal entry is greater than the sum of the absolute values of all other entries in the same row.
// An arbitrary matrix is positive definite if and only if each of its principal submatrices has a positive determinant.
// Note that only the last case does the implication go both ways.




////DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
//double prec_dlamch( char *cmach )
//{
//
//}

//DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
//*
//*  -- LAPACK auxiliary routine (version 3.2) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          CMACH
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DLAMCH determines double precision machine parameters.
//*
//*  Arguments
//*  =========
//*
//*  CMACH   (input) CHARACTER*1
//*          Specifies the value to be returned by DLAMCH:
//*          = 'E' or 'e',   DLAMCH := eps
//*          = 'S' or 's ,   DLAMCH := sfmin
//*          = 'B' or 'b',   DLAMCH := base
//*          = 'P' or 'p',   DLAMCH := eps*base
//*          = 'N' or 'n',   DLAMCH := t
//*          = 'R' or 'r',   DLAMCH := rnd
//*          = 'M' or 'm',   DLAMCH := emin
//*          = 'U' or 'u',   DLAMCH := rmin
//*          = 'L' or 'l',   DLAMCH := emax
//*          = 'O' or 'o',   DLAMCH := rmax
//*
//*          where
//*
//*          eps   = relative machine precision
//*          sfmin = safe minimum, such that 1/sfmin does not overflow
//*          base  = base of the machine
//*          prec  = eps*base
//*          t     = number of (base) digits in the mantissa
//*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
//*          emin  = minimum exponent before (gradual) underflow
//*          rmin  = underflow threshold - base**(emin-1)
//*          emax  = largest exponent before overflow
//*          rmax  = overflow threshold  - (base**emax)*(1-eps)
//*
//* =====================================================================


// *******************************************************************************
// PABLO's OLD DIAGONALIZATION-FUNCTION (DEFPROT-back compatibility)
// *******************************************************************************
/*
         SUBROUTINE SSYEVR( JOBZ, RANGE, UPLO,  N,  A,  LDA,  VL,  VU,  IL,  IU,
                          ABSTOL,  M,  W,  Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,
                          LIWORK, INFO )
           CHARACTER      JOBZ, RANGE, UPLO
           INTEGER        IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
           REAL           ABSTOL, VL, VU
           INTEGER        ISUPPZ( * ), IWORK( * )
           REAL           A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
        */
extern "C"
{
int ssyevr_(char *jobz, char *range, char *uplo, int *n,
                float *a, int *lda, float *vl, float *vu, int *il,
                int *iu, float *abstol, int* m, float *w,
                float *z__, int *ldz, int *isuppz, float *work,
                int *lwork, int *iwork, int *liwork, int *info);

/* void ssyevr_(char jobz, char range, char uplo, int  n,
						    float *a,  int  lda, float vl, float vu, int il,
						    int iu, float abstol, int m, float *w,
						    float *z, int ldz,int *isuppz, int *info);
*/
}

int nma (int size, int neigval, float *hessian, float ** eigvalp, float ** eigvectp)
{
int i;
float * eigval, * eigvect;
float toler_zero_mode;

toler_zero_mode=0.0001;

/* Lapack input/output */
 char jobz, range, uplo;
 int lda, il, iu, ldz, * isuppz, lwork, liwork, info, * iwork;
 float vl, vu, toler, * work;

 /* allocate memory */
 eigval = ( float * ) malloc( size * sizeof( float ) );
 for ( i = 0; i < size; i++ ) eigval[i] = 0.0;

 eigvect = ( float * ) malloc( size * neigval * sizeof( float ) );
 for ( i = 0; i < size * neigval; i++ ) eigvect[i] = 0.0;

 /**  */
 /* diagonalize matrix here */
 /* USING SSYEVR/DSYEVR FORTRAN SUBRUTINE FROM LAPACK */
 /**  */

 jobz = 'V';
 /* 'N':  Compute eigenvalues only; 'V':  Compute eigenvalues and eigenvectors */

 if ( neigval == size )
   range = 'A';
 else
   range = 'I';
 /* 'A': all eigenvalues will be found. 'V': all eigenvalues in the half-open interval (VL,VU] will be found.
 'I': the IL-th through IU-th eigenvalues will be found */

 if ( range == 'I' )
 {
   /* If  RANGE='I', the indices (in ascending order) of the smallest and largest eigenvalues to be returned.
   1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  Not referenced if RANGE = 'A' or 'V' */

   il = 1;
   iu = neigval;
 }
 else
 {
   il = 0;
   iu = 0;
 }

 uplo = 'U';
 /* 'U':  Upper triangle of A is stored; 'L':  Lower triangle of A is stored.
 On  entry, the symmetric matrix If UPLO = 'U', the leading N-by-N upper triangular part of the matrix
 input contains the upper triangular part of the matrix A. If UPLO = 'L', the leading N-by-N  lower  triangular
 part  of  A  contains  the  lower  triangular  part of the matrix A.
 On exit, the lower triangle (if UPLO='L') or the upper triangle (if UPLO='U') of A, including the diagonal, is destroyed.
 */

 vl = 0.0;
 vu = 9E29;
 /* If RANGE='V', the lower and upper bounds of the  interval  to  be  searched  for eigenvalues.
 VL < VU.  Not referenced if RANGE = 'A' or 'I'. */

 toler = 0;

 /**  */
 /* OUTPUT */
 /**  */

 /*   neigval      The total number of eigenvalues found. 0 <= neigval <= N.  If RANGE = 'A', neigval = N,
 and if RANGE = 'I',  neigval = IU-IL+1. */


 /* The first neigval elements of eigval will contain the selected eigenvalues in ascending order. */

 /* If  JOBZ = 'V', then if INFO = 0, the first neigval columns of
 eigvect contain the orthonormal eigenvectors of the matrix A corresponding to the selected eigenvalues, with the i-th column
 of eigvect holding the eigenvector associated with eigval[i].

 If JOBZ = 'N', then eigvect is not referenced. Note: the user must
 ensure that at least max(1,neigval) columns are supplied in the array Z;

 if RANGE = 'V', the exact value of neigval is not known in advance and an upper bound must be used. */

 lda = size;
 ldz = size;

 if ( range == 'I' )
 {
   isuppz = ( int * ) malloc( 2 * neigval * sizeof( int ) );
   for ( i = 0; i < 2 * neigval; i++ ) isuppz[i] = 0;
 }
 else
 {
   isuppz = ( int * ) malloc( 2 * size * sizeof( int ) );
   for ( i = 0; i < 2 * size; i++ ) isuppz[i] = 0;
 }


 lwork = size * 26;
 work = ( float * ) malloc( lwork * sizeof( float ) );
 for ( i = 0; i < lwork; i++ ) work[i] = 0.0;

 liwork = size * 10;
 iwork = ( int * ) malloc( liwork * sizeof( int ) );
 for ( i = 0; i < liwork; i++ ) iwork[i] = 0;


 /*ssyevr( &jobz,  &range,  &uplo,  &size, hessian,  &lda,
          &vl,  &vu,  &il,  &iu,  &toler,
          &neigval, eigval, eigvect,
         & ldz, isuppz, work, & lwork, iwork, & liwork, & info );
 */
/* ssyevr_( jobz,  range,  uplo,  size, hessian,  lda,
          vl,  vu,  il,  iu,  toler,
          neigval, eigval, eigvect,
         ldz, isuppz, &info );
*/

ssyevr_( & jobz, & range, & uplo, & size, hessian, & lda,
        & vl, & vu, & il, & iu, & toler,
        & neigval, eigval, eigvect,
         & ldz, isuppz, work, & lwork, iwork, & liwork, & info );

float dump;
for ( i = 6; i < neigval - 1; i++ )
  {

  dump=fabs(eigval[i]) - 0.0001*eigval[i + 1];
  if (dump>0.0001)  break;

  }


 if ( info )
 {
   printf( "nmac> An error occured in the matrix diagonalization %d\n", info );
   exit( 1 );
 }
 else
   printf( "nmac> %d eigvalues found (first %d null)\n", neigval, i );




free(isuppz);
free(work);
free(iwork);


*eigvalp=eigval;
*eigvectp=eigvect;

return i;

}
