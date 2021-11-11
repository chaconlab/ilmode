/*************************************************************************
 *                     LIBRARY: libnma_matrix                            *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Library with some Matrix operation rutines.                         *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

//#include <nma.h>
#include <libnma_matrix.h>

// Mon made (1/10/2008)
// Computes the rotation matrix corresponding to a
// rotation angle "a" (in radians) around the vector v[3].
void rotmat(double a, float v[3], double **rot)
{
	double sa, ca, uca, t, v1, v2, v3;
	double v11, v12, v13, v22, v23, v33;
	int m;
	// "a" is in radians (Mon: 20/1/2009)

//	t=sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
	t=sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
//	for(m=0;m<3;m++)     /* normalize v */
//		v[m] /= t;
	v[0] /= t; v[1] /= t; v[2] /= t; // normalize v
	v1=v[0]; v2=v[1]; v3=v[2];
	v11=v1*v1; v12=v1*v2; v13=v1*v3;
	v22=v2*v2; v23=v2*v3; v33=v3*v3;
	ca=cos(a); sa=sin(a);
//	ca = 1 - a*a/2; sa = a - a*a*a/6; // Taylor expansion... (2nd not null terms)
	uca=1.0-ca;

	rot[0][0] = ca+v11*uca;
	rot[0][1] = v12*uca-v3*sa;
	rot[0][2] = v13*uca+v2*sa;
	rot[1][0] = v12*uca+v3*sa;
	rot[1][1] = ca+v22*uca;
	rot[1][2] = v23*uca-v1*sa;
	rot[2][0] = v13*uca-v2*sa;
	rot[2][1] = v23*uca+v1*sa;
	rot[2][2] = ca+v33*uca;
}

// Mon made (21/5/2008)
// Multiply matrices by accumulation
// a x b = b' (Warning: "b" will be overwritten)
void matacc(double a[3][3], double b[3][3])
{
	// matrix accumulation: M1 * M2 = M2 (3x3)
	//    c{i,j} = SUM( a{i,r}b{r,j} )
	double c[3][3];
	int i,j,k;


	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			c[i][j] = 0.0; // initialization
			for(k=0;k<3;k++)
				c[i][j] += a[i][k] * b[k][j];
		}

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			b[i][j] = c[i][j]; // output
}

// Mon made (1/10/2008)
// 3x3 matrix multiplication: M1 x M2 = M
void mult3(double **M1, double **M2, double **M)
{
	int i,j,k;
	double temp;
	for(j=0; j<3; j++) // rows
	{
		for(i=0; i<3; i++) // cols
		{
			temp = 0.0;
			for(k=0; k<3; k++)
				temp += M1[j][k] * M2[k][i]; // row_1 x col_2, ...
			M[j][i] = temp;
		}
	}
}

// Mon made (1/10/2008)
// (3x3) x (3x1) matrix-vector multiplication:
// R x t = tf (Applies a rotation to a position vector)
void multvec3(double **R, Tcoor t, Tcoor tf)
{
	int i,j,k;
	double temp;
	for(j=0; j<3; j++) // rows
	{
		temp = 0.0;
		for(i=0; i<3; i++) // cols
			temp += R[j][i] * t[i]; // row_1 x col_2, ...
		tf[j] = temp;
	}
}

// Mon made (1/10/2008)
// Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
void rmot(double **R, Tcoor t, Tcoor pos, Tcoor posf)
{
	posf[0] = R[0][0]*pos[0] + R[0][1]*pos[1] + R[0][2]*pos[2] + t[0];
	posf[1] = R[1][0]*pos[0] + R[1][1]*pos[1] + R[1][2]*pos[2] + t[1];
	posf[2] = R[2][0]*pos[0] + R[2][1]*pos[1] + R[2][2]*pos[2] + t[2];
}

// Computes 3x3 inverse matrix (obtained from internet)
// PennyBoki @ </dream.in.code>
// #include <stdio.h>
int inverse( double A[3][3], double X[3][3])
{
	double B[3][3];//the transpose of a matrix A
	double C[3][3];//the adjunct matrix of transpose of a matrix A not adjunct of A

	int i,j;
	double x,n=0;//n is the determinant of A

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			B[i][j]=0;
			C[i][j]=0;
		}
	}

	for(i=0,j=0;j<3;j++)
	{
		if(j==2)
			n+=A[i][j]*A[i+1][0]*A[i+2][1];
		else if(j==1)
			n+=A[i][j]*A[i+1][j+1]*A[i+2][0];
		else
			n+=A[i][j]*A[i+1][j+1]*A[i+2][j+2];
	}
	for(i=2,j=0;j<3;j++)
	{
		if(j==2)
			n-=A[i][j]*A[i-1][0]*A[i-2][1];
		else if(j==1)
			n-=A[i][j]*A[i-1][j+1]*A[i-2][0];
		else
			n-=A[i][j]*A[i-1][j+1]*A[i-2][j+2];
	}

	if(n!=0) x=1.0/n;
	else
	{
		printf("Division by 0, not good!\n");
		printf("=====================================================================\n\n");
		return 0;
	}

	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			B[i][j]=A[j][i];

	C[0][0]=B[1][1]*B[2][2]-(B[2][1]*B[1][2]);
	C[0][1]=(-1)*(B[1][0]*B[2][2]-(B[2][0]*B[1][2]));
	C[0][2]=B[1][0]*B[2][1]-(B[2][0]*B[1][1]);

	C[1][0]=(-1)*(B[0][1]*B[2][2]-B[2][1]*B[0][2]);
	C[1][1]=B[0][0]*B[2][2]-B[2][0]*B[0][2];
	C[1][2]=(-1)*(B[0][0]*B[2][1]-B[2][0]*B[0][1]);

	C[2][0]=B[0][1]*B[1][2]-B[1][1]*B[0][2];
	C[2][1]=(-1)*(B[0][0]*B[1][2]-B[1][0]*B[0][2]);
	C[2][2]=B[0][0]*B[1][1]-B[1][0]*B[0][1];

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			X[i][j]=C[i][j]*x;

		}
	}

	return 0;
}

// Rounds floating point variables according to the desired precission
double rounder( double in, int prec )
{
	double temp = pow( (double)10, prec );
	in *= temp;
	in = round( in );
	in /= temp;
	return( in );
}
