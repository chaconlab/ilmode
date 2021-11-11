#ifndef LIBNMA_MATRIX_H_
#define LIBNMA_MATRIX_H_

#include "nma.h"

// Mon made (1/10/2008)
// Computes the rotation matrix corresponding to a
// rotation angle "a" (in radians) around the vector v[3].
void rotmat(double a, float v[3], double **rot);

// Infinitesimal rotation matrix (commutative properties holds!)
// Mon made (4/1/2010)
void d_rotmat(double adeg, float v[3], double **rot);

// Mon made (21/5/2008)
// Multiply matrices by accumulation
// a x b = b' (Warning: "b" will be overwritten)
void matacc(double a[3][3], double b[3][3]);

// Computes the inverse of the 3x3 matrix A, returning it in X (obtained from internet)
int inverse( double A[3][3], double X[3][3] );

// Mon made (1/10/2008)
// 3x3 matrix multiplication: M1 x M2 = M
void mult3(double M1[3][3], double M2[3][3], double M[3][3]);
void mult3(double **M1, double **M2, double **M);

// Mon made (1/10/2008)
// (3x3) x (3x1) matrix-vector multiplication:
// R x t = tf (Applies a rotation to a position vector)
void multvec3(double **R, Tcoor t, Tcoor tf);

// Mon made (1/10/2008)
// Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
void rmot(double **R, Tcoor t, Tcoor pos, Tcoor posf);

// Rounds floating point variables according to the desired precission
double rounder( double in, int prec );

#endif /*LIBNMA_MATRIX_H_*/
