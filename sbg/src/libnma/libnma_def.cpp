// Jose Ramon Lopez-Blanco.
// Pablo Chacon's Structural Bioinformatics Group
// (url: sbg.cib.csic.es)
// (CIB-CSIC)

// ######################################################################################
//
// NMAC related library
// (Deformability)
//
// ######################################################################################

#include <nma.h>
#include <libnma_def.h>

#define TINY 1.0e-20
#define TOL 1.0e-04
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
     a[k][l]=h+s*(g-h*tau);


#define ptr_check(p) { \
    if ((p) == NULL) { \
        fprintf(stderr, "Memory allocation failure. Exiting.\n"); \
        exit(1); \
    } \
}


/*
extern "C" {
void ludcmp(double **, int , int *, double *, double * );
void hardy(int, double *, double **, double **, int *,double *);
void grad_hmq(float *, int, double **, double **, int *, double *,
              double *, double *, double *, double * );
double lnorm(double, double, double, double, double, double);
void djacobi3( double [4][4], int , double [4], double [4][4], int * );
void lubksb(double **, int, int *, double []);
void compute_def (int , float *, double *, float *, float *, int, int, double **, double **);
void save_ascii_defmob (int, double *, double *);
void norm_defmob (int , double *);
}
*/

  /* computes eigenvalues and eigenvectors of a real symmetric matrix a */
  /* taken from Numerical Recipes */



/* ================================================================== */
double lnorm( double a11, double a12, double a13, double a22, double a23, double a33 )
{
  /* computes the lambda-norm (i.e., the maximum of the absolute values of
  the eigenvalues) of a 3x3 symmetric matrix a, whose upper triangular part is given as input */

  double mev, temp;
  int nrot;


  //double eigval[4], eigvect[4][4], a[4][4];
  double a[4][4], v[4][4], b[4], z[4], d[4];

  a[1] [1] = a11; a[1] [2] = a12; a[1] [3] = a13;
  a[2] [1] = a[1] [2]; a[2] [2] = a22; a[2] [3] = a23;
  a[3] [1] = a[1] [3]; a[3] [2] = a[2] [3]; a[3] [3] = a33;


  /* computes eigenvalues and eigenvectors of a real symmetric matrix a */
  /* taken from Numerical Recipes */

int j, iq, ip, i,n;
double tresh, theta, tau, t, sm, s, h, g, c;


n=3;

for ( ip = 1; ip <= n; ip++ )
{
for ( iq = 1; iq <= n; iq++ ) v[ip] [iq] = 0.0;
v[ip] [ip] = 1.0;
}
for ( ip = 1; ip <= n; ip++ )
{
b[ip] = d[ip] = a[ip] [ip];
z[ip] = 0.0;
}
nrot = 0;
for ( i = 1; i <= 50; i++ )
{
sm = 0.0;
for ( ip = 1; ip <= n - 1; ip++ )
{
 for ( iq = ip + 1; iq <= n; iq++ )
   sm += fabs( a[ip] [iq] );
}
if ( sm == 0.0 )
{
 i=500;
 break;
}
if ( i < 4 )
 tresh = 0.2 * sm / ( n * n );
else
 tresh = 0.0;
for ( ip = 1; ip <= n - 1; ip++ )
{
 for ( iq = ip + 1; iq <= n; iq++ )
 {
   g = 100.0 * fabs( a[ip] [iq] );
   if ( i > 4 && ( float )( fabs( d[ip] ) + g ) == ( float )fabs( d[ip] )
        && ( float )( fabs( d[iq] ) + g ) == ( float )fabs( d[iq] ) )
          a[ip] [iq] = 0.0;
   else if ( fabs( a[ip] [iq] ) > tresh )
   {
     h = d[iq] - d[ip];
     if ( ( float )( fabs( h ) + g ) == ( float )fabs( h ) )
       t = ( a[ip] [iq] ) / h;
     else
     {
       theta = 0.5 * h / ( a[ip] [iq] );
       t = 1.0 / ( fabs( theta ) + sqrt( 1.0 + theta * theta ) );
       if ( theta < 0.0 ) t = -t;
     }
     c = 1.0 / sqrt( 1 + t * t );
     s = t * c;
     tau = s / ( 1.0 + c );
     h = t * a[ip] [iq];
     z[ip] -= h;
     z[iq] += h;
     d[ip] -= h;
     d[iq] += h;
     a[ip] [iq] = 0.0;
     for ( j = 1; j <= ip - 1; j++ )
     {
       ROTATE( a, j, ip, j, iq )
     }
     for ( j = ip + 1; j <= iq - 1; j++ )
     {
       ROTATE( a, ip, j, j, iq )
     }
     for ( j = iq + 1; j <= n; j++ )
     {
       ROTATE( a, ip, j, iq, j )
     }
     for ( j = 1; j <= n; j++ )
     {
       ROTATE( v, j, ip, j, iq )
     }
     ++nrot;
   }
 }
}
for ( ip = 1; ip <= n; ip++ )
{
 b[ip] += z[ip];
 d[ip] = b[ip];
 z[ip] = 0.0;
}
}



  mev = fabs( d[1] );
  if ( ( temp = fabs( d[2] ) ) > mev ) mev = temp;
  if ( ( temp = fabs( d[3] ) ) > mev ) mev = temp;


  return mev;
}




/* ============================================================= */
void lubksb( double * * a, int n, int * indx, double b[] )
{
  /* performs the backsubstitution using the output from ludcmp */
  /* taken from Numerical recipes */

  int i, ii = 0, ip, j;
  double sum;

  for ( i = 1; i <= n; i++ )
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if ( ii )
      for ( j = ii; j <= i - 1; j++ ) sum -= a[i] [j] * b[j];
    else if ( sum ) ii = i;
    b[i] = sum;
  }
  for ( i = n; i >= 1; i-- )
  {
    sum = b[i];
    for ( j = i + 1; j <= n; j++ ) sum -= a[i] [j] * b[j];
    b[i] = sum / a[i] [i];
  }
}
/*=================================================================*/
void ludcmp(double **a, int n, int *indx, double *d, double *vv)
{
/* performs the LU decomposition of matrix a */
/* taken from Numerical Recipes */
/* warning ALWAYS DOUBLE */

        int i,imax,j,k;
        double big,dum,temp;
        long double sum, dummy;

        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs( a[i][j])) > big) big=temp;
                if (big == 0.0) fprintf(stderr,"Singular matrix in routine ludcmp\n");
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum = sum -  a[i][k]*a[k][j];
                        a[i][j]=(double) sum;
                }
                big=0.0;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum = sum - a[i][k]*a[k][j];
                        a[i][j]=double(sum);
                        if ( (dum=vv[i]*fabsl(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum= a[imax][k];
                                a[imax][k]= a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                indx[j]=imax;
                if (a[j][j] == 0.0 ) a[j][j]=TINY;
                if (j != n) {
                        dummy=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] = (double) dummy*a[i][j];
                }
        }

}



/* ======================================================================= */
void hardy( int num_points, double * dist_matrix, double ** a, double ** b, int * indx, double *vv, double dc )
{
	//  computes the matrix for the hardy hyperbolic multiquadric interpolant
	//  corresponding to the set of points in pdb_model, and performs its LU
	//  decomposition. Returns the original matrix in b, the factorized matrix in a and the permutations in indx.
	//  These have to be used as input to grad_hmq.

	// dc is a (squared) 'characteristic distance' that specifies the shape of the hyperboloids. dc = [10-20] works fine!
	// Previously recommended value: 0.1*diam(pdb_model)^2, was system size dependent!!!

//	double d, diam, dist, temp;
	double d, dist, temp;
	int i, j;
	long l;

//	/* compute diameter of pdb_model */
//	diam = 0.0;
//	for ( l = 0; l < num_points * ( num_points - 1 ) / 2; l++ )
//		if ( ( temp = dist_matrix[l] ) > diam ) diam = temp;

	for ( i = 0; i < num_points+1; i++ )
		indx[i]=0;


	for(i=0;i<num_points+1;i++)
		for(j=0;j<num_points+1;j++)
			b[i][j]=a[i][j]=0;

	/* compute matrix, copying it to b because a will be destroyed by ludcmp */
	for(i=0;i<num_points;i++)
	{
		for(j=0;j<num_points;j++)
		{
			if(j==i)
				dist=0.0;
			else if(i<j)
				dist=dist_matrix[num_points*i+j-1-i*(i+3)/2];
			else
				dist=dist_matrix[num_points*j+i-1-j*(j+3)/2];
			b[i+1][j+1]=a[i+1][j+1]=sqrt(dist*dist+dc);
		}
	}

	/* do the LU decomposition of a */
	ludcmp( a, num_points, indx, & d , vv);
}


/* ======================================================================= */
void grad_hmq( float *pdb_model, int num_points, double * * a, double * * b, int * indx, double * f, double * f1,
     double * f2, double * f3 , double * ft)
     {
       /* computes the gradient of a function given by its values
       on an arbitrary (scattered) set of points in R^3, by computing the Hardy
       hyperbolic multiquadric that interpolates the given function values.
       Unlike grad_hess_poly, this function evaluates the derivatives at all points of the set. */
       /* a, b and indx are the output of hardy, which has to be called once before
       making calls to the present function with various f */
       /* f is the input array of function values; f1,f2,f3 are the first partial derivatives (output). */

       double  c1;
       double fx, fy, fz;
       int i, j;


       /* move f to ft (shifted by 1) (lubksb overwrites ft) */
       for ( j = 0; j < num_points; j++ )
       {
         ft[j + 1] = f[j];
       }

       /* do the backsubstitution using ft */
       lubksb( a, num_points, indx, ft );

       /* now ft contains the c's of the Hardy multiquadric */

       /* compute the partials */
       for ( j = 0; j < num_points; j++ )
       {
         fx = fy = fz = 0.0;
         for ( i = 0; i < num_points; i++ )
         {
           c1 = ft[i + 1] / b[i + 1] [j + 1];
           fx += c1 * ( pdb_model[j*3    ] - pdb_model[i*3    ] );
           fy += c1 * ( pdb_model[j*3 + 1] - pdb_model[i*3 + 1] );
           fz += c1 * ( pdb_model[j*3 + 2] - pdb_model[i*3 + 2] );
         }
         f1[j] = fx; f2[j] = fy; f3[j] = fz;
       }

}


/* ================================================================== */

void compute_def (int num_atoms, float *pdb_model, double *dist_matrix,
                  myfloat *eigval, myfloat *eigvect, int neigval, int nf1,
                  double **defp, double  **mobp) {
int i,j,size;
// FILE *out;
double scale_factor;
char file_name[8];

size=num_atoms*3;

/* compute and factorize the interpolant matrix */
/* amat is the factorized matrix and bmat is the original matrix */

double **amat, **bmat, *work;
int *indx;
double *u1, *u2, *u3;
double *u11, *u12, *u13, *u21, *u22, *u23, *u31, *u32, *u33;
double tr3;
double current_mob, maximum_mob, current_def, maximum_def;


work = ( double *) malloc( sizeof( double ) * (num_atoms+1) );  ptr_check(work);
indx = ( int *) malloc( sizeof( int ) * (num_atoms+1) );  ptr_check(indx);

amat = ( double  ** ) malloc( sizeof( double* ) * (num_atoms+1) );  ptr_check(amat);
for ( i = 0; i < num_atoms+1; i++ ) {
       amat[i] = ( double * ) malloc( sizeof( double ) * (num_atoms+1) );
       ptr_check(amat[i]);
}
bmat = ( double  ** ) malloc( sizeof( double* ) * (num_atoms+1) ); ptr_check(amat);
for ( i = 0; i < num_atoms+1; i++ ) {
          bmat[i] = ( double * ) malloc( sizeof( double ) * (num_atoms+1) );
          ptr_check(bmat[i]);
}
for(i=0;i<num_atoms;i++)
  for(j=0;j<num_atoms;j++)
     amat[i][j]=bmat[i][j]=0;


u1 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u1);
u2 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u2);
u3 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u3);

u11 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u11);
u12 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u12);
u13 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u13);
u21 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u21);
u22 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u22);
u23 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u23);
u31 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u31);
u32 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u32);
u33 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u32);




/* initialize deformability and mobility arrays */
double *deform, *mob;
deform = (double *) malloc(num_atoms * sizeof(double));
mob = (double *) malloc(num_atoms * sizeof(double));

  for(i=0;i<num_atoms;i++){
    deform[i]=0.0;
    mob[i]=0.0;
  }


hardy(num_atoms, dist_matrix, amat, bmat, indx, work);

scale_factor = sqrt(eigval[nf1]);


//printf("nmat>scale_factor %d %d %d %f\n", nf1, neigval, size, scale_factor);
printf("nmac>\nnmac>                           Computing deformability    \nnmac>\n" );
printf("nmac>\nnmac>          MODE    EIGENVALUE      MAX_MOB          MAX_DEF\nnmac>\n");

/* compute deform and mob */
for(i=nf1; i<neigval; i++){
  for(j=0;j<num_atoms;j++){
    u1[j]=eigvect[i*size+3*j+0];
    u2[j]=eigvect[i*size+3*j+1];
    u3[j]=eigvect[i*size+3*j+2];
  }
  sprintf(file_name,"mod%d",i+1);

//  out = fopen(file_name, "w");

//  fprintf(out,"#>T Eigenvector_%d\n",i+1);
//  fprintf(out,"#>  Res         VX            VY            VZ            Mob            Def\n");

  grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u1,
                u11, u12, u13, work);
  grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u2,
                u21, u22, u23, work);
  grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u3,
                u31, u32, u33, work);

  maximum_mob=0.0;
  maximum_def=0.0;
  for(j=0;j<num_atoms;j++){
    tr3 = (u11[j]+u22[j]+u33[j])/3.0;     /* div(u)/3 */
    current_def = lnorm(u11[j]-tr3,0.5*(u12[j]+u21[j]),
                        0.5*(u13[j]+u31[j]),u22[j]-tr3,
                        0.5*(u23[j]+u32[j]),u33[j]-tr3)/eigval[i];
    current_def *= scale_factor;  /* to make deform indep. of scaling of mass and contacts */
    deform[j] += pow(current_def,2);
    current_mob = (u1[j]*u1[j]+u2[j]*u2[j]+u3[j]*u3[j])/eigval[i];
    mob[j] += current_mob;
    current_mob = sqrt(current_mob);
//           fprintf(out,"%6d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n",
//                j+1, u1[j], u2[j], u3[j], current_mob, current_def);
          if (current_mob > maximum_mob)
      maximum_mob = current_mob;
          if (current_def > maximum_def)
      maximum_def = current_def;
  }

//  fclose(out);
   if (i<27) printf("nmat>   Mode %4d %15.3E %15.3E %15.3E\n",
               i+1, eigval[i], maximum_mob, maximum_def);
//  fprintf(log,"  %4d %15.5E %15.5E %15.5E     %s\n", i+1, eigval[i], maximum_mob, maximum_def, file_name);

}



for(j=0;j<num_atoms;j++){
  deform[j]=sqrt(deform[j]);
  mob[j] = sqrt(mob[j]);
}

// fclose(log);

*mobp   = mob;
*defp = deform;

free(u1);free(u2);free(u3);
free(u11);free(u21);free(u31);
free(u12);free(u22);free(u32);
free(u13);free(u23);free(u33);

}

// Mon modified (25/3/2010) (just changed "myfloat" by "floating")
void compute_def(int num_atoms, float *pdb_model, double *dist_matrix,
                  floating *eigval, floating *eigvect, int neigval, int nf1,
                  double **defp, double  **mobp, double dc)
{
	int i,j,size;
	// FILE *out;
	double scale_factor;
	char file_name[8];

	size=num_atoms*3;

	/* compute and factorize the interpolant matrix */
	/* amat is the factorized matrix and bmat is the original matrix */

	double **amat, **bmat, *work;
	int *indx;
	double *u1, *u2, *u3;
	double *u11, *u12, *u13, *u21, *u22, *u23, *u31, *u32, *u33;
	double tr3;
	double current_mob, maximum_mob, current_def, maximum_def;


	work = ( double *) malloc( sizeof( double ) * (num_atoms+1) );  ptr_check(work);
	indx = ( int *) malloc( sizeof( int ) * (num_atoms+1) );  ptr_check(indx);

	amat = ( double  ** ) malloc( sizeof( double* ) * (num_atoms+1) );  ptr_check(amat);
	for ( i = 0; i < num_atoms+1; i++ ) {
		amat[i] = ( double * ) malloc( sizeof( double ) * (num_atoms+1) );
		ptr_check(amat[i]);
	}
	bmat = ( double  ** ) malloc( sizeof( double* ) * (num_atoms+1) ); ptr_check(amat);
	for ( i = 0; i < num_atoms+1; i++ ) {
		bmat[i] = ( double * ) malloc( sizeof( double ) * (num_atoms+1) );
		ptr_check(bmat[i]);
	}
	for(i=0;i<num_atoms;i++)
		for(j=0;j<num_atoms;j++)
			amat[i][j]=bmat[i][j]=0;

	u1 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u1);
	u2 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u2);
	u3 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u3);

	u11 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u11);
	u12 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u12);
	u13 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u13);
	u21 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u21);
	u22 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u22);
	u23 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u23);
	u31 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u31);
	u32 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u32);
	u33 = (double *) malloc(num_atoms * sizeof(double)); ptr_check(u32);

	/* initialize deformability and mobility arrays */
	double *deform, *mob;
	deform = (double *) malloc(num_atoms * sizeof(double));
	mob = (double *) malloc(num_atoms * sizeof(double));

	for(i=0;i<num_atoms;i++)
	{
		deform[i]=0.0;
		mob[i]=0.0;
	}

	hardy(num_atoms, dist_matrix, amat, bmat, indx, work, dc);

	scale_factor = sqrt(eigval[nf1]);

	//printf("nmat>scale_factor %d %d %d %f\n", nf1, neigval, size, scale_factor);
	printf("nmac>\nnmac>                           Computing deformability    \nnmac>\n" );
	printf("nmac>\nnmac>          MODE    EIGENVALUE      MAX_MOB          MAX_DEF\nnmac>\n");

	/* compute deform and mob */
	for(i=nf1; i<neigval; i++)
	{
		for(j=0;j<num_atoms;j++)
		{
			u1[j]=eigvect[i*size+3*j+0];
			u2[j]=eigvect[i*size+3*j+1];
			u3[j]=eigvect[i*size+3*j+2];
		}
		// sprintf(file_name,"mod%d",i+1);

		//  out = fopen(file_name, "w");

		//  fprintf(out,"#>T Eigenvector_%d\n",i+1);
		//  fprintf(out,"#>  Res         VX            VY            VZ            Mob            Def\n");

		grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u1,
				u11, u12, u13, work);
		grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u2,
				u21, u22, u23, work);
		grad_hmq(pdb_model, num_atoms, amat, bmat, indx, u3,
				u31, u32, u33, work);

		// fprintf(stderr,"\nDeformability of mode %d\n",i);
		maximum_mob=0.0;
		maximum_def=0.0;
		for(j=0;j<num_atoms;j++)
		{
			tr3 = (u11[j]+u22[j]+u33[j])/3.0;     /* div(u)/3 */
//			current_def = lnorm(u11[j]-tr3, 0.5*(u12[j]+u21[j]), 0.5*(u13[j]+u31[j]),
//								u22[j]-tr3,	0.5*(u23[j]+u32[j]), u33[j]-tr3) / eigval[i];
			current_def = lnorm(u11[j]-tr3, 0.5*(u12[j]+u21[j]), 0.5*(u13[j]+u31[j]),
								u22[j]-tr3,	0.5*(u23[j]+u32[j]), u33[j]-tr3);
//			current_def *= scale_factor;  /* to make deform indep. of scaling of mass and contacts */
			current_def *= sqrt(scale_factor);  /* to make deform indep. of scaling of mass and contacts */
			deform[j] += pow(current_def,2) / eigval[i];
			// fprintf(stderr,"%5d %10f %10f\n",j,deform[j],pow(current_def,2) / eigval[i]);
			current_mob = (u1[j]*u1[j]+u2[j]*u2[j]+u3[j]*u3[j]) / eigval[i];
			mob[j] += current_mob;
			current_mob = sqrt(current_mob);
			//           fprintf(out,"%6d  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n",
			//                j+1, u1[j], u2[j], u3[j], current_mob, current_def);
			if (current_mob > maximum_mob)
				maximum_mob = current_mob;
			if (current_def > maximum_def)
				maximum_def = current_def;
		}

		//  fclose(out);
		if (i<20) printf("nmat>   Mode %4d %15.3E %15.3E %15.3E\n",	i+1, eigval[i], maximum_mob, maximum_def);
		//  fprintf(log,"  %4d %15.5E %15.5E %15.5E     %s\n", i+1, eigval[i], maximum_mob, maximum_def, file_name);
	}

	for(j=0;j<num_atoms;j++)
	{
		deform[j]=sqrt(deform[j]);
		mob[j] = sqrt(mob[j]);
	}

	// fclose(log);

	*mobp   = mob;
	*defp = deform;

	free(u1);free(u2);free(u3);
	free(u11);free(u21);free(u31);
	free(u12);free(u22);free(u32);
	free(u13);free(u23);free(u33);

}

void save_ascii_defmob (int num_atoms, double *deform, double *mob,  double *bf, char *list, char *file)
{

FILE *mf;
double maxd, mind, maxb, minb,maxm, minm;
int j;
char file_name[12];
char amino[4];
amino[3]='\0';

// compute max min
maxd=-1000; mind=1000;
maxm=-1000; minm=1000;
maxb=-1000; minb=1000;
for(j=0;j<num_atoms;j++) {
  if (maxd<deform[j])  maxd=deform[j];
  if (mind>deform[j])  mind=deform[j];
  if (maxm<mob[j])     maxm=mob[j];
  if (minm>mob[j])     minm=mob[j];
  if (maxb<bf[j])     maxb=bf[j];
  if (minb>bf[j])     minb=bf[j];
}



/* output */

//sprintf(file_name,"defmob.tab");
//mf = fopen(file_name, "w");

// Modified by Mon (26/3/2010)
if(file==NULL)
{
	sprintf(file_name,"defmob.tab");
	mf = fopen(file_name, "w");
}
else
	mf = fopen(file, "w");

fprintf(mf,"#>T t   \n");
fprintf(mf,"#>---r---s----d-----------dn----------m-----------mn----------bf----------bfn-------\n");

for(j=0;j<num_atoms;j++)
{
	if(list!=NULL)
	{
		amino[0]=list[j*3];
		amino[1]=list[(j*3)+1];
		amino[2]=list[(j*3)+2];
	}
	else
		amino[0] = amino[1] = amino[2] = ' ';
	if (maxb<0.0000001)
	{


		fprintf(mf,"%6d %s %10.5E %10.5E %10.5E %10.5E %10.5E %10.5E\n",
				j+1,amino, deform[j], deform[j]/maxd, mob[j], mob[j]/maxm, bf[j], bf[j] );
		//              (mob[j]-minm)/(maxm-minm),


	} else
		fprintf(mf,"%6d %s %10.5E %10.5E %10.5E %10.5E %10.5E %10.5E\n",
				j+1, amino, deform[j], deform[j]/maxd, mob[j], mob[j]/maxm, bf[j], bf[j]/maxb );
	//              (mob[j]-minm)/(maxm-minm),
	//              (deform[j]-mind)/(maxd-mind));
}
fclose(mf);



//for(j=0;j<num_atoms;j++)
//  pdb_model[j].tempFactor=
//   (pdb_model[j].tempFactor-minx)/(maxx-minx);

//writepdb("def.pdb",num_atoms,pdb_model);

}

// Just for back-compatibility with "nmac"
void save_ascii_defmob (int num_atoms, double *deform, double *mob,  double *bf) {

FILE *mf;
double maxd, mind, maxb, minb,maxm, minm;
int j;
char file_name[12];

// compute max min
maxd=-1000; mind=1000;
maxm=-1000; minm=1000;
maxb=-1000; minb=1000;
for(j=0;j<num_atoms;j++) {
  if (maxd<deform[j])  maxd=deform[j];
  if (mind>deform[j])  mind=deform[j];
  if (maxm<mob[j])     maxm=mob[j];
  if (minm>mob[j])     minm=mob[j];
  if (maxb<bf[j])     maxb=bf[j];
  if (minb>bf[j])     minb=bf[j];
}



/* output */
sprintf(file_name,"defmob.tab");
mf = fopen(file_name, "w");
fprintf(mf,"#>T t   \n");
fprintf(mf,"#>---r---d-----------dn----------m-----------mn----------bf----------bfn-------\n");

for(j=0;j<num_atoms;j++) {

  if (maxb<0.0000001) {

    fprintf(mf,"%6d  %10.5E %10.5E %10.5E %10.5E %10.5E %10.5E\n",
              j+1, deform[j], deform[j]/maxd, mob[j], mob[j]/maxm, bf[j], bf[j] );
//              (mob[j]-minm)/(maxm-minm),


  } else
  fprintf(mf,"%6d  %10.5E %10.5E %10.5E %10.5E %10.5E %10.5E\n",
              j+1, deform[j], deform[j]/maxd, mob[j], mob[j]/maxm, bf[j], bf[j]/maxb );
//              (mob[j]-minm)/(maxm-minm),
//              (deform[j]-mind)/(maxd-mind));
}
fclose(mf);



//for(j=0;j<num_atoms;j++)
//  pdb_model[j].tempFactor=
//   (pdb_model[j].tempFactor-minx)/(maxx-minx);

//writepdb("def.pdb",num_atoms,pdb_model);

}


void norm_defmob (int num_atoms, double *deform)
{

double maxd, mind;
int j,posmax;

// compute max min
maxd=-1000; mind=1000; posmax=0;
for(j=0;j<num_atoms;j++) {
  if (maxd<deform[j])  {maxd=deform[j]; posmax=j;};
  if (mind>deform[j])  mind=deform[j];
}

if (posmax==0) maxd=(deform[0]+deform[1])/2.0;
else if (posmax==num_atoms-1) maxd=(deform[posmax]+deform[posmax-1])/2.0;
	else if (deform[posmax-1]>deform[posmax+1]) maxd=(deform[posmax]+deform[posmax-1])/2.0;
		else maxd=(deform[posmax]+deform[posmax+1])/2.0;

for(j=0;j<num_atoms;j++) {
if (deform[j]>=maxd) deform[j]=100;
else deform[j]= 100*(deform[j]-mind)/(maxd-mind);
// deform[j]= 1/ (1+exp(-(deform[j]-mind)/(maxd-mind)));
}
}
