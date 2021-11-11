/*************************************************************************
 *                     LIBRARY: libnma_kinetic                           *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Library to compute Kinetic Energy Matrices using several methods:   *
 *   	-In Internal Coordinate Space (ICS):                             *
 *      -naive K-matrix (O^3),                                           *
 *      -naive V/W-arrays (memmory efficient) (O^3),                     *
 *      -fast and memmory efficient Go's method (O^2).                   *
 *      -In CCS --> Diagonal matrix with atomic masses (not implemented) *
 *   See Go's classic papers about the so called "Naive" method.         *
 *   (It takes into account Muliple-Chains and different CG-models)      *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

//#include <nma.h>
#include <libnma_kinetic.h>
#include <libnma_misc.h>  // this should be removed when seg_props() is passed as an argument
#include <libnma_matrix.h>
#define DEVEL

// Calculates: elMY, emMY, el, em vectors. Needed for Fast Kinetic Energy computation.
// Table I. and (ec. 44) from Braun et al. (1984)
void calcKine2(double mb1, double *Y1, double mb3, double *Y3, double **I1, double **I3, double *Phi, double *Psi, double *elMY, double *emMY, double *el, double *em)
{
	double N[3],T[3],L[3]; // auxiliary vectors
	int m,n;

	// Table I. from Braun et al. (1984)
	// Phi = e_lambda  || e_mu
	// Psi = e_lambda  x  y_lambda  ||  e_mu x y_mu

	// Some terms of (ec. 44) from Braun et al. (1984)
	// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	// elMY[size] = M1 * Y1  x  (Psi) - I1 * e_lambda)

	for(m=0; m<3; m++)
		N[m] = mb1 * Y1[m]; // M1 * Y1
	//printf("PHI P= %f,%f,%f\n",Psi[0],Psi[1],Psi[2]);
	T[0] = N[1]*Psi[2] - N[2]*Psi[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
	T[1] = N[2]*Psi[0] - N[0]*Psi[2];
	T[2] = N[0]*Psi[1] - N[1]*Psi[0];
	//printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
	for(m=0; m<3; m++)
	{
		L[m] = 0.0;
		for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
			L[m] += I1[m][n] * Phi[n]; // I1 * e_lambda
		elMY[m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
		//  ==> (3x1)
	}
	//printf("PHI elMY= %f,%f,%f\n",elMY[j][0],elMY[j][1],elMY[j][2]);

	// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
	for(m=0; m<3; m++)
		N[m] = mb3 * Y3[m]; // M3 * Y3
	T[0] = N[1]*Psi[2] - N[2]*Psi[1]; // M3 * Y3  x  (e_mu  x  y_mu)
	T[1] = N[2]*Psi[0] - N[0]*Psi[2];
	T[2] = N[0]*Psi[1] - N[1]*Psi[0];
	//printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
	for(m=0; m<3; m++)
	{
		L[m] = 0.0;
		for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
			L[m] += I3[m][n] * Phi[n]; // I3 * e_mu
		emMY[m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
	}
	//printf("PHI emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

	// Remaining terms of (ec. 44) from Braun et al. (1984)
	// el[size]; // Psi - Phi
	// em[size]; // Psi - Phi
//	for( m = 0; m < 3; m++)
//	{
//		L[m] = Phi[m] - Y1[m]; // (y_lambda - Y1)
//		T[m] = Phi[m] - Y3[m]; // (y_mu - Y3)
//	}
	// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
//	el[0] = Phi[1]*L[2] - Phi[2]*L[1]; // e_lambda  x  (y_lambda - Y1)
//	el[1] = Phi[2]*L[0] - Phi[0]*L[2];
//	el[2] = Phi[0]*L[1] - Phi[1]*L[0];
	el[0] = Psi[0] - (Phi[1]*Y1[2] - Phi[2]*Y1[1]); // Psi - (Phi x Y1)
	el[1] = Psi[1] - (Phi[2]*Y1[0] - Phi[0]*Y1[2]);
	el[2] = Psi[2] - (Phi[0]*Y1[1] - Phi[1]*Y1[0]);
	//printf("PHI el= %f,%f,%f\n",el[j][0],el[j][1],el[j][2]);
//	em[0] = Phi[1]*T[2] - Phi[2]*T[1]; // e_mu  x  (y_mu - Y3)
//	em[1] = Phi[2]*T[0] - Phi[0]*T[2];
//	em[2] = Phi[0]*T[1] - Phi[1]*T[0];
	em[0] = Psi[0] - (Phi[1]*Y3[2] - Phi[2]*Y3[1]); // Psi - (Phi x Y3)
	em[1] = Psi[1] - (Phi[2]*Y3[0] - Phi[0]*Y3[2]);
	em[2] = Psi[2] - (Phi[0]*Y3[1] - Phi[1]*Y3[0]);

	//printf("PHI em= %f,%f,%f\n",em[j][0],em[j][1],em[j][2]);

}


// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain)
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naive( Macromolecule *mol, floating **p_mass_matrix, trd *der, int size)
{
    int ks,isi,jsi,ind2,ind3;
    double mak,mdx,mdy,mdz;

    floating *mass_matrix;

 	if(*p_mass_matrix == NULL)
 	{
 		if( !(mass_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
 		{
 			printf("Msg(kinetic_naive): I'm sorry, Kinetic-Energy Matrix memory allocation failed!\n"
 					"Forcing exit!\n");
 			exit(1);
 		}
 		*p_mass_matrix = mass_matrix; // outputs Hessian Matrix
 	}
 	else
 		mass_matrix = *p_mass_matrix;

	// Kinetic energy matrix initialization
	for(int i=0; i<size*(size+1)/2; i++)
		mass_matrix[i]=0.0;

  //Iterador para recorrer atomos
  pdbIter *iter1 = new pdbIter( mol );

  // Compute Mass matrix
  // Ec. (6) from Noguti & Go, 1983
  for( iter1->pos_atom = 0,ks=0; !iter1->gend_atom(); iter1->next_atom(),ks+=size )
  {
  	// screens ALL atoms
    mak= ( iter1->get_atom() )->getPdbocc(); // atom mass
    for(int i=0; i<size; i++)
    { // screens lamda-dihedrals' "dy/dTheta(lamda)" (K-matrix)
      ind2=ks+i;
      mdx=mak*der[ind2].x;
      mdy=mak*der[ind2].y;
      mdz=mak*der[ind2].z;
      for(int j=i;j<size;j++) // (including diagonal)
      { // screens mu-dihedrals' "dy/dTheta(mu)" (K-matrix)
        ind3=ks+j;
        // upper triangular part (including diagonal)
        mass_matrix[i + j*(j+1)/2] += mdx*der[ind3].x + mdy*der[ind3].y + mdz*der[ind3].z;// H-matrix
      }
    }
  }
  delete iter1;
}

// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain)
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naive_double( Macromolecule *mol, double **p_mass_matrix, trd *der, int size)
{
    int ks,isi,jsi,ind2,ind3;
    double mak,mdx,mdy,mdz;

    double *mass_matrix;

 	if(*p_mass_matrix == NULL)
 	{
 		if( !(mass_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
 		{
 			printf("Msg(kinetic_naive): I'm sorry, Kinetic-Energy Matrix memory allocation failed!\n"
 					"Forcing exit!\n");
 			exit(1);
 		}
 		*p_mass_matrix = mass_matrix; // outputs Hessian Matrix
 	}
 	else
 		mass_matrix = *p_mass_matrix;

	// Kinetic energy matrix initialization
	for(int i=0; i<size*(size+1)/2; i++)
		mass_matrix[i]=0.0;

  //Iterador para recorrer atomos
  pdbIter *iter1 = new pdbIter( mol );

  // Compute Mass matrix
  // Ec. (6) from Noguti & Go, 1983
  for( iter1->pos_atom = 0,ks=0; !iter1->gend_atom(); iter1->next_atom(),ks+=size )
  {
  	// screens ALL atoms
    mak= ( iter1->get_atom() )->getPdbocc(); // atom mass
    for(int i=0; i<size; i++)
    { // screens lamda-dihedrals' "dy/dTheta(lamda)" (K-matrix)
      ind2=ks+i;
      mdx=mak*der[ind2].x;
      mdy=mak*der[ind2].y;
      mdz=mak*der[ind2].z;
      for(int j=i;j<size;j++) // (including diagonal)
      { // screens mu-dihedrals' "dy/dTheta(mu)" (K-matrix)
        ind3=ks+j;
        // upper triangular part (including diagonal)
        mass_matrix[i + j*(j+1)/2] += mdx*der[ind3].x + mdy*der[ind3].y + mdz*der[ind3].z;// H-matrix
      }
    }
  }
  delete iter1;
}

// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain) (for all models and types)
// The same as kinetic_naive(), but computing derivatives "on the fly" via V/W arrays
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naiveVW( Macromolecule *mol, float *coord, floating **p_mass_matrix,double ***V,double ***W,bool **body1, int size)
{
    int i,j;
    double mak;
    float rk[3];
    floating *mass_matrix;

 	if(*p_mass_matrix == NULL)
 	{
 		if( !(mass_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
 		{
 			printf("Msg(kinetic_naiveVW): I'm sorry, Kinetic-Energy Matrix memory allocation failed!\n"
 					"Forcing exit!\n");
 			exit(1);
 		}
 		*p_mass_matrix = mass_matrix; // outputs Hessian Matrix
 	}
 	else
 		mass_matrix = *p_mass_matrix;

	// Kinetic energy matrix initialization
	for(i=0; i<size*(size+1)/2; i++)
		mass_matrix[i]=0.0;

	// Dummy array of vectors (stores derivatives for current atom)
	double *deri[3];
	for(i=0; i<3; i++)
		deri[i] = (double *) malloc( sizeof(double) * size);

  //Iterador para recorrer atomos
  pdbIter *iter1 = new pdbIter( mol );

  // Compute Mass matrix
  // Ec. (6) from Noguti & Go, 1983
  for( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() )
  { // screens ALL atoms
    mak = sqrt( ( iter1->get_atom() )->getPdbocc() ); // atom mass
    rk[0] = coord[iter1->pos_atom*3]; // k-atom position
    rk[1] = coord[iter1->pos_atom*3+1];
    rk[2] = coord[iter1->pos_atom*3+2];

    // Loading dummy array of derivatives for current atom (like in the Hessian case)
    for(i=0; i<size; i++)
    	if( body1[i][iter1->pos_atom] ) // if body 1 atom
    	{
    		deri[0][i] = mak * (V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1]);
    		deri[1][i] = mak * (V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2]);
    		deri[2][i] = mak * (V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0]);
    	}
    	else
    	{
    		deri[0][i] = mak * (V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1]);
    		deri[1][i] = mak * (V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2]);
    		deri[2][i] = mak * (V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0]);
    	}

    for(i=0; i<size; i++) // screens lamda-dihedrals' "dy/dTheta(lamda)" (K-matrix)
      for(j=i;j<size;j++) // screens mu-dihedrals' "dy/dTheta(mu)" (K-matrix)
        // upper triangular part (including diagonal)
        mass_matrix[i + j*(j+1)/2] += deri[0][i]*deri[0][j]
                                    + deri[1][i]*deri[1][j]
                                    + deri[2][i]*deri[2][j];// H-matrix
  }
  delete iter1;

  for(i=0; i<3; i++)
	free( deri[i] );
}

// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain) (for all models and types)
// The same as kinetic_naive(), but computing derivatives "on the fly" via V/W arrays
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naiveVW_double( Macromolecule *mol, float *coord, double **p_mass_matrix,double ***V,double ***W,bool **body1, int size)
{
    int i,j;
    double mak;
    float rk[3];
    double *mass_matrix;

 	if(*p_mass_matrix == NULL)
 	{
 		if( !(mass_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
 		{
 			printf("Msg(kinetic_naive): I'm sorry, Kinetic-Energy Matrix memory allocation failed!\n"
 					"Forcing exit!\n");
 			exit(1);
 		}
 		*p_mass_matrix = mass_matrix; // outputs Hessian Matrix
 	}
 	else
 		mass_matrix = *p_mass_matrix;

	// Kinetic energy matrix initialization
	for(i=0; i<size*(size+1)/2; i++)
		mass_matrix[i]=0.0;

	// Dummy array of vectors (stores derivatives for current atom)
	double *deri[3];
	for(i=0; i<3; i++)
		deri[i] = (double *) malloc( sizeof(double) * size);

  //Iterador para recorrer atomos
  pdbIter *iter1 = new pdbIter( mol );

  // Compute Mass matrix
  // Ec. (6) from Noguti & Go, 1983
  for( iter1->pos_atom = 0; !iter1->gend_atom(); iter1->next_atom() )
  { // screens ALL atoms
    mak = sqrt( ( iter1->get_atom() )->getPdbocc() ); // atom mass
    rk[0] = coord[iter1->pos_atom*3]; // k-atom position
    rk[1] = coord[iter1->pos_atom*3+1];
    rk[2] = coord[iter1->pos_atom*3+2];

    // Loading dummy array of derivatives for current atom (like in the Hessian case)
    for(i=0; i<size; i++)
    	if( body1[i][iter1->pos_atom] ) // if body 1 atom
    	{
    		deri[0][i] = mak * (V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1]);
    		deri[1][i] = mak * (V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2]);
    		deri[2][i] = mak * (V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0]);
    	}
    	else
    	{
    		deri[0][i] = mak * (V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1]);
    		deri[1][i] = mak * (V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2]);
    		deri[2][i] = mak * (V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0]);
    	}

    for(i=0; i<size; i++) // screens lamda-dihedrals' "dy/dTheta(lamda)" (K-matrix)
      for(j=i;j<size;j++) // screens mu-dihedrals' "dy/dTheta(mu)" (K-matrix)
        // upper triangular part (including diagonal)
        mass_matrix[i + j*(j+1)/2] += deri[0][i]*deri[0][j]
                                    + deri[1][i]*deri[1][j]
                                    + deri[2][i]*deri[2][j];// H-matrix
  }
  delete iter1;

  for(i=0; i<3; i++)
	free( deri[i] );
}

// Computes the Kinetic Energy Matrix (H) (Multi-Chain & FULL-ATOM & Protein/RNA/DNA/SMOL)
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMFAx( Macromolecule *mol, float *coord2, tri *props, int size, floating **p_Hmatrix, int type, int model, bool *fix, bool *addrot)
{
	bool debug=false;
	bool debug2=false;
	bool debug3=false;
	double mtot,mta,mpr,mb1,mb3,nr2,mfa,fm;
	double r[3],my[3],y[3],pos[3],e[3],Y1[3],Y3[3],myt[3],T[3],P[3],L[3],N[3],EL[3];
	double I[3][3],I1[3][3],I3[3][3],J[3][3];
	double temp,temp1,temp2;
	int m, n, num_atoms, i, k, j, k1, k2, isi;
	int resn,num_res,num_seg,index_res; //,chis;
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	index_res = 0;
	floating *Hmatrix;
	int CBindex;
	TMOL fragtype;
	int indexbase;
	int seg_atoms;

	// (DOUBLE) array, dimension (N*(N+1)/2)
	if( *p_Hmatrix == NULL ) // if not NULL
	{
		if( !(Hmatrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 ) ) )
		{
			printf("Msg(kineticMFAx): Memory allocation failed!\nForcing exit\n");
			exit(1);
		}
		*p_Hmatrix = Hmatrix;
	}
	else
		Hmatrix = *p_Hmatrix;

	if(type == 1)
	{
		printf("Msg(kineticMFAx): I'm sorry, invalid Dihedral angle ordering type: Phi,Psi,Chi.\n");
		printf("Msg(kineticMFAx): (Only type=0,2 are valid!) Forcing exit!\n");
		exit(1);
	}

	j = 0; // current dihedral index
	int j2 = 0; // active dihedral index (screens all mobile dihedrals)
	mb1 = 0.0; // Mass of Body 1
	mpr = 0.0; // mpr = mass of previous residues
	my[0] = my[1] = my[2] = 0.0; // Sum(mass*coord) of previous residues
	k1 = 0;
	k2 = props[0].nat - 1; // 1st and last indices of atoms of residue i // Residue 0 ?

	// Centering with double precision (necessary with huge systems!!! <-- the fucking bug!!!)
	double *coord; // Now, "coord" is "double" instead of "single" precision. (and "coord2" is single)
	if( !(coord = (double *) malloc(sizeof(double)*3*num_atoms)) )
	{
		printf("Msg(kineticMCA_disc): Unable to allocate memory! Forcing exit!\n");
		exit(1);
	}

	// Computing the PDB's Center of Mass (CoM)
	/* total mass */
	r[0] = r[1] = r[2] = 0.0;
	mtot = 0.0;
	i = 0;
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )  // screens all-atoms
	{
		mta = ( iter->get_atom() )->getPdbocc(); // Load mass...
		mtot += mta;
		// Sum(mass*coord) before putting the CoM at 0
		r[0] += mta * coord2[i * 3];
		r[1] += mta * coord2[i * 3 + 1];
		r[2] += mta * coord2[i * 3 + 2];
		i++;
	}
	for(m=0;m<3;m++)
		r[m] /= mtot; // the Center of Mass
	if(debug)
		printf( "Msg(kineticMFAx): Mass %8.8f Center %8.8f %8.8f %8.8f \n",mtot,r[0],r[1],r[2]);

	// Setting "double" centered coordinates
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() )  // screens all-atoms
	{
		coord[iter->pos_atom*3]     = coord2[iter->pos_atom*3]     - r[0];
		coord[iter->pos_atom*3 + 1] = coord2[iter->pos_atom*3 + 1] - r[1];
		coord[iter->pos_atom*3 + 2] = coord2[iter->pos_atom*3 + 2] - r[2];
	}

	// Compute Inertia Matrix (Inertia Momentum) (coords CoM must be centered 0,0,0)
	// Ec. 22 (Noguti & Go, 1983)
	for ( int i = 0; i < 3; i++ )
		for ( int k = 0; k < 3; k++ )
			I[i][k] = 0.0; // Initializing Inertia-Matrix
	int ii = 0;
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screens all-atoms
	{
		mta = ( iter->get_atom() )->getPdbocc(); // mass
		r[0] = coord[ii * 3];
		r[1] = coord[ii * 3 + 1];
		r[2] = coord[ii * 3 + 2];
		ii++;
		temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) )
		for (int  m = 0; m < 3; m++ )
			for (int n = 0; n < 3; n++ )
			{
				temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
				if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
				I[m] [n] += mta * temp2;
			}
	}

	// compute the inverse of I
	inverse( I, J ); // J = Inverse of I (3x3 matrix)
	if(debug)
	{
		printf("TOTAL I-Inertia-Matrix:\n");
		for (int k = 0; k < 3; k++ )
		{
			for (int l = 0; l < 3; l++ )
				printf("%10.4f ",I[k][l]);
			printf("\n");
		}
	}

	double *M1; // body-1 mass
	double *M3; // body-3 mass
	double **el; // e_lambda  x  (y_lambda - Y1)
	double **em; // e_mu  x  (y_mu - Y3)
	double **elMY; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	double **emMY; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
	bool *ischi; // CHIs are different...
	bool has_oxt=false; // true = current segment has OXT
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

	M1 = (double *) malloc( sizeof(double) * size ); // body-1 mass
	M3 = (double *) malloc( sizeof(double) * size ); // body-3 mass
	el = (double **) malloc( sizeof(double *) * size ); // e_lambda  x  (y_lambda - Y1)
	em = (double **) malloc( sizeof(double *) * size ); // e_mu  x  (y_mu - Y3)
	elMY = (double **) malloc( sizeof(double *) * size ); // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	emMY = (double **) malloc( sizeof(double *) * size ); // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
	for(int i=0; i<size; i++)
	{
		el[i] = (double *) malloc( sizeof(double) * 3);
		em[i] = (double *) malloc( sizeof(double) * 3);
		elMY[i] = (double *) malloc( sizeof(double) * 3);
		emMY[i] = (double *) malloc( sizeof(double) * 3);
	}
	ischi = (bool *) malloc( sizeof(bool) * size ); // CHIs are different...

	// Pre-Computing M1, M3, el, em, elMY, emMY
	// Screen ALL residues
	Segment * seg;
	Atom * atom;
	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();

	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(debug)
			printf ("\nProcessing segment %d  num_res= %d  k1 %d k2 %d ***************\n"
					,iter_seg->pos_segment,num_res,k1,k2);

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Computing MASSES (updating)
			// M1 and M2
			mb1 = mpr; // Left_half = previous
			mb3 = mtot - mb1; // Right-half = Tot - Left_half
			if(debug3)
				printf("ROT-TRANS pos_segment= %d  j= %d  mpr= %f  mb1= %f  mb3= %f\n",iter_seg->pos_segment,j,mpr,mb1,mb3);

			// y_lambda (NH pos 0) --> NOT EXIST FOR TRASLATION
			// e_lambda --> NOT EXIST FOR TRASLATION

			// Computing CoM
			for ( m = 0; m < 3; m++ )
			{
				// [ (Previous CoM) * (Previous mass) + (NH mass)*(NH pos) ] / mb1 = body 1 CoM
				Y1[m] = my[m] / mb1;
				// [ (Remaining CoM) * (remaining mass) - (NH mass)*(NH pos) ] / mb3 = body 3 CoM
				Y3[m] = -my[m] / mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
			}
			//				printf("TRANS y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
			//						,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

			// Adding 3 TRASLATIONS
			for(int axis=0; axis<3; axis++) // 3 translations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{

					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3; // Right-half = Tot - Left_half

					if(debug)
						printf("TRANS pos_segment= %d  j= %d  mpr= %f  mb1= %f  mb3= %f  k1= %d  k2= %d\n"
								,iter_seg->pos_segment,j,mpr,mb1,mb3,k1,k2);

					// Table I. --> Braun et al. (1984)
					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda
					// elMY[size] = M1 * Y1  x  gamma_v
					for(m=0; m<3; m++)
					{
						P[m] = 0.0; // gamma_v initialization (Psi)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					}
					P[axis] = -1.0; // x,y,z unit vector (Psi)

					//				printf("j= %d   TRANS P= %f,%f,%f\n",j,P[0],P[1],P[2]);
					elMY[j2][0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  gamma_v,lambda
					elMY[j2][1] = N[2]*P[0] - N[0]*P[2];
					elMY[j2][2] = N[0]*P[1] - N[1]*P[0];
					//				printf("TRANS elMY= %f,%f,%f\n",elMY[j][0],elMY[j][1],elMY[j][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					// emMY[size] = M3 * Y3  x  gamma_v,mu
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					emMY[j2][0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  gamma_v,mu
					emMY[j2][1] = N[2]*P[0] - N[0]*P[2];
					emMY[j2][2] = N[0]*P[1] - N[1]*P[0];
					//	 			printf("TRANS emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

					// el[size]; // gamma_v,lambda
					// em[size]; // gamma_v,mu
					for(m=0; m<3; m++)
					{
						el[j2][m] = P[m];
						em[j2][m] = P[m];
					}
					//				printf("TRANS el= %f,%f,%f\n",el[j][0],el[j][1],el[j][2]);
					//				printf("TRANS em= %f,%f,%f\n",em[j][0],em[j][1],em[j][2]);
					ischi[j2] = false;

					j2++;
				}
				j++; // adding traslational degree of freedom
			} // 3 trans added

			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
			{
				// Calculate matrix I1 (Inertia Matrix 1)
				for(m=0;m<3;m++)
					for(n=0;n<3;n++)
						I1[m][n] = 0.0;

				if(debug)
					printf("ROT pos_segment= %d  j= %d  mpr= %f  mb1= %f  mb3= %f  k1= %d  k2= %d\n"
							,iter_seg->pos_segment,j,mpr,mb1,mb3,k1,k2);

				for ( iter->pos_atom = 0; iter->pos_atom < k1; iter->next_atom() )
				{
					mta = ( iter->get_atom() )->getPdbocc(); // mass
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
					{
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
					}
				}
				// now matrix I3
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I3[m][n] = I[m][n] - I1[m][n];
				if(debug)
				{
					printf("ROT  I1-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%10.4f ",I1[k][l]);
						printf("\n");
					}
					printf("ROT  I3-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%9.6f ",I3[k][l]);
						printf("\n");
					}
				}
				// End Inertia matrix computation

				// Adding 3 ROTATIONS
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{

						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3; // Right-half = Tot - Left_half

						if(debug)
							printf("ROT pos_segment= %d  j= %d  mpr= %f  mb1= %f  mb3= %f  k1= %d  k2= %d\n"
									,iter_seg->pos_segment,j,mpr,mb1,mb3,k1,k2);

						// Centers of Mass
						for ( m = 0; m < 3; m++ )
						{
							e[m] = 0.0; // delta_v initialization
							y[m] = 0.0;
						}
						e[axis] = 1.0; // delta_v

						//			printf("PHI y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
						//					,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// elMY[size] = M1 * Y1  x  (delta_v  x  u ) - I1 * delta_v
						P[0] = e[1]*y[2] - e[2]*y[1]; // delta_v  x  u  (Table I.Psi)
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						// printf("j= %d   ROT  P= %f,%f,%f\n",j,P[0],P[1],P[2]);
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (delta_v  x  u)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						//			printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * delta_v
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (delta_v  x  u) - I1 * delta_v)
							//  ==> (3x1)
						}
						if(debug2)
							printf("ROT  elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						// emMY[size]; // M3 * Y3  x  (delta_v  x  u ) - I3 * delta_v
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (delta_v  x  u)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						//			printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I3[m][n] * e[n]; // I3 * delta_v
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (delta_v  x  u ) - I3 * delta_v
							// ==> (3x1)
						}
						if(debug2)
							printf("ROT  emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

						// el[size]; // e_lambda  x  (y_lambda - Y1)
						// el[size]; // (delta_v  x  u)  -  (delta_v  x  Y1)
						T[0] = e[1]*Y1[2] - e[2]*Y1[1]; // (delta_v  x  Y1)
						T[1] = e[2]*Y1[0] - e[0]*Y1[2];
						T[2] = e[0]*Y1[1] - e[1]*Y1[0];
						for(m=0;m<3;m++)
							el[j2][m] = P[m] - T[m]; // (delta_v  x  u)  -  (delta_v  x  Y1)
						if(debug2)
							printf("ROT el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);

						// em[size]; // e_mu  x  (y_mu - Y3)
						// el[size]; // (delta_v  x  u)  -  (delta_v  x  Y3)
						T[0] = e[1]*Y3[2] - e[2]*Y3[1]; // (delta_v  x  Y3)
						T[1] = e[2]*Y3[0] - e[0]*Y3[2];
						T[2] = e[0]*Y3[1] - e[1]*Y3[0];
						for(m=0;m<3;m++)
							em[j2][m] = P[m] - T[m]; // (delta_v  x  u)  -  (delta_v  x  Y3)
						if(debug2)
							printf("ROT em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);
						ischi[j2] = false;
						j2++;
					}
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // adding rotational degree of freedom
				} // 3 rots added
			}
		}

		// Screening residues (fragments)
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666; // SMOL's can not be tabulated...

			if(debug)
				printf("FRAG index_res= %d  pos_fragment= %d (tmol= %d) %s  k1= %d  k2= %d\n"
						,index_res,iter_frag->pos_fragment,fragtype,res->getName(),k1,k2);

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{

				// NOT FIRST, NOT PRO -> PHI
				// ("PHI-bond")
				if ( iter_frag->pos_fragment != 0 && resn != PRO )
				{
					// printf ("%d Not first not PRO -> PHI \n",j);
					//  q_j = phi_i
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						//					printf ("PHI  j= %d  j2= %d\n",j,j2);
						// Calculate matrix I1 (Inertia Matrix 1)
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0;
						for ( iter->pos_atom = 0; iter->pos_atom <= k1; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc(); // mass
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for ( m = 0; m < 3; m++ )
							{
								for ( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
							}
						}

						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m][n] = I[m][n] - I1[m][n];
						//			printf("PHI I1-Inertia-Matrix:\n");
						//			for (int k = 0; k < 3; k++ )
						//			{
						//				for (int l = 0; l < 3; l++ )
						//					printf("%10.4f ",I1[k][l]);
						//				printf("\n");
						//			}
						//			printf("PHI I3-Inertia-Matrix:\n");
						//			for (int k = 0; k < 3; k++ )
						//			{
						//				for (int l = 0; l < 3; l++ )
						//					printf("%9.6f ",I3[k][l]);
						//				printf("\n");
						//			}

						iter->pos_atom = k1; // first residue atomic index (NH pseudo-atom)
						// ("k2" is updated below (for the 1st not PRO -> PHI)

						// Computing MASSES (updating)
						mta= ( iter->get_atom() )->getPdbocc(); // "mta" --> actual mass (just NH !)
						/* M1 and M2 */
						mb1 = mpr + mta; // Left_half = previous + actual (NH)
						mb3 = mtot - mb1; // Right-half = Tot - Left_half
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3; // Right-half = Tot - Left_half
						//					printf("PHI  mb1= %f  mb3= %f\n",mb1,mb3);

						// y_lambda (NH pos 0)
						y[0] = coord[k1 * 3];
						y[1] = coord[k1 * 3 + 1];
						y[2] = coord[k1 * 3 + 2];
						// CA pos 1
						e[0] = coord[(k1+1) * 3];
						e[1] = coord[(k1+1) * 3 + 1];
						e[2] = coord[(k1+1) * 3 + 2];
						// e_lambda
						e[0] -= y[0]; // NH --> CA
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							// [ (Previous CoM) * (Previous mass) + (NH mass)*(NH pos) ] / mb1 = body 1 CoM
							Y1[m] = ( my[m] +  mta * y[m] ) / mb1;
							// [ (Remaining CoM) * (remaining mass) - (NH mass)*(NH pos) ] / mb3 = body 3 CoM
							Y3[m] = -(mb1*Y1[m])/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						}
						//			printf("PHI y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
						//					,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						if(debug3)
							printf("PHI  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb3,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);
//						if(debug3)
//							for(int x=0;x<3;x++)
//								printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n",
//										I[x][0],I[x][1],I[x][2],I1[x][0],I1[x][1],I1[x][2],I3[x][0],I3[x][1],I3[x][2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						//printf("PHI P= %f,%f,%f\n",P[0],P[1],P[2]);
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						//printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							//  ==> (3x1)
						}
						//printf("PHI elMY= %f,%f,%f\n",elMY[j][0],elMY[j][1],elMY[j][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						//printf("PHI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						//printf("PHI emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

						// el[size]; // e_lambda  x  (y_lambda - Y1)
						// em[size]; // e_mu  x  (y_mu - Y3)
						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1]; // e_lambda  x  (y_lambda - Y1)
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						//printf("PHI el= %f,%f,%f\n",el[j][0],el[j][1],el[j][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1]; // e_mu  x  (y_mu - Y3)
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						//printf("PHI em= %f,%f,%f\n",em[j][0],em[j][1],em[j][2]);

						ischi[j2] = false; // it's not chi

						j2++;
					}

					j++; // der[] index
				}  // NOT FIRST NOT PRO

				// Checking if last residue has OXT (Ot)
				if(iter_frag->pos_fragment == num_res-1 && model == 2)
				{
					iter->pos_atom = k2;
					if( strcmp(iter->get_atom()->getName()," OXT ")==0 )
					{
						has_oxt = true;
						k2--;
					}
				}

				if(model==1)
					CBindex = 3;
				else
					CBindex = 4;

				// ********************************************************************************************
				//  LATERAL CHAIN-->CHI
				//
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if(type==2)
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frag->pos_fragment==0 ||
											(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						//	printf ("%d Lateral Chain --> CHI\n", j);
						//  q_j = chi_i
						if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
						{
							//						printf ("CHI  j= %d  j2= %d\n",j,j2);

							// Calculate matrix I3
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I3[m][n] = 0;
							// this should be expressed different!
							// screening [k1+2,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
								//				if ( iter->pos_atom > k1+2 ) // CB & R
							{
								mta = ( iter->get_atom() )->getPdbocc();
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];
								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for ( m = 0; m < 3; m++ )
								{
									for ( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I3[m] [n] += mta * mfa;
									}
								}
							}
							// now matrix I1
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I1[m][n] = I[m][n]-I3[m][n];

							// MASSES
							// screening CB & -R atoms (CB: "k1+3")
							mb3=0;
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
								mb3 += ( iter->get_atom() )->getPdbocc(); // atom mass
							// Now CB & -R is the whole 2nd body, the remaining atoms will be the 1st one!
							mb1 = mtot-mb3;
							// M1 & M3 arrays
							M1[j2] = mb1;
							M3[j2] = mb3; // Right-half = Tot - Left_half
							//						printf("CHI  mb1= %f  mb3= %f\n",mb1,mb3);

							// CA pos 1
							y[0] = coord[(k1+1) * 3 + 0];
							y[1] = coord[(k1+1) * 3 + 1];
							y[2] = coord[(k1+1) * 3 + 2];
							// CB pos 3
							e[0] = coord[(k1+CBindex) * 3 + 0];
							e[1] = coord[(k1+CBindex) * 3 + 1];
							e[2] = coord[(k1+CBindex) * 3 + 2];
							// e_lambda ==> CA-->CB unit vector
							e[0] -= y[0];
							e[1] -= y[1];
							e[2] -= y[2];

							temp = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
							for(m=0;m<3;m++)
							{
								e[m] /= temp; // unit vector normalization
								Y1[m] = 0.0;    /* to initialize */
								// Mon:
								// Not initialization!
								// The CoM of the whole protein is = 0.0 !!!
								// We'll substract to the CoM (total) the contribution
								// of the 2nd body (CB-R). (See below)
							}

							// this should be expressed different!
							// screening [k1+2,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc(); // atom mass
								// Substracting to the CoM (total) the body 2 (CB-R) contribution
								for(m=0; m<3; m++)
									Y1[m] -= mta * coord[iter->pos_atom*3 + m];
							}

							for(m=0;m<3;m++)
							{
								Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
								Y1[m] /= mb1;
							}

							if(debug3)
								printf("CHI  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb3,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);
//							if(debug3)
//								for(int x=0;x<3;x++)
//									printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n",
//											I[x][0],I[x][1],I[x][2],I1[x][0],I1[x][1],I1[x][2],I3[x][0],I3[x][1],I3[x][2]);

							// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							for(m=0; m<3; m++)
								N[m] = mb1 * Y1[m]; // M1 * Y1
							P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
							P[1] = e[2]*y[0] - e[0]*y[2];
							P[2] = e[0]*y[1] - e[1]*y[0];
							T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
								elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
								// ==> (3x1)
							}

							// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
							for(m=0; m<3; m++)
								N[m] = mb3 * Y3[m]; // M3 * Y3
							T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
								emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
							}

							for( m = 0; m < 3; m++)
							{
								Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
								Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
							}
							// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
							el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
							el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
							el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
							em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
							em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
							em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
							ischi[j2] = true; // it is chi
							j2++;
						}
						j++; // der[] index
					}

				// NOT LAST RESIDUE--> PSI (but Full-Atom allways has PSI)
				// ("PSI-bond")
				if ( iter_frag->pos_fragment != num_res - 1 || model == 2) // Full-Atom allways has PSI!
				{
					// printf("%d Not last residue -> PSI\n",j);
					// q_j = psi_i
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// calculate matrix I1
						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
							if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								// Do this with iterators: **********************
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];

								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for( m = 0; m < 3; m++ )
								{
									for( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I1[m][n] +=  mta * mfa;
									}
								}
							}

						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];
						//			printf("PSI I1-Inertia-Matrix:\n");
						//			for (int k = 0; k < 3; k++ )
						//			{
						//				for (int l = 0; l < 3; l++ )
						//					printf("%10.4f ",I1[k][l]);
						//				printf("\n");
						//			}
						//			printf("PSI I3-Inertia-Matrix:\n");
						//			for (int k = 0; k < 3; k++ )
						//			{
						//				for (int l = 0; l < 3; l++ )
						//					printf("%9.6f ",I3[k][l]);
						//				printf("\n");
						//			}

						// Updating M1
						// Actual body 1 mass (M1)
						mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
						// Screening atoms [k1,k2] of the actual residue (except C-atom)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
								//					if ( iter->pos_atom != k1+2 && iter->pos_atom != k1+3 ) // NOT C and NOT O
								mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

						mb3 = mtot-mb1;
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3;
						//					printf("PSI  mb1= %f  mb3= %f\n",mb1,mb3);
						// get CA pos 1
						y[0] = coord[(k1+1) * 3 + 0];
						y[1] = coord[(k1+1) * 3 + 1];
						y[2] = coord[(k1+1) * 3 + 2];
						// get C pos 2 (C=O) or C in full-atom
						e[0] = coord[(k1+2) * 3 + 0];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= y[0];
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m]; /* to initialize */ // with previous "masses*positions"
						}

						// Updating Body 1 Center of Mass: "Y1[m]"
						// Screening atoms [k1,k2] of the actual residue (except C-atom)
						for(iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
								//					if ( iter->pos_atom != k1+2 && iter->pos_atom != k1+3 ) // NOT C and NOT O
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								Y1[0] += mta * coord[iter->pos_atom*3];
								Y1[1] += mta * coord[iter->pos_atom*3+1];
								Y1[2] += mta * coord[iter->pos_atom*3+2];
							}

						for(m=0;m<3;m++) // 3D coords
						{
							//						Y1[m] /= mb1;
							//						Y3[m] = -(mb1*Y1[m])/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y1[m] /= mb1;
						}
						//			printf("PSI y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
						//					,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						if(debug3)
							printf("PSI  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb3,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);
//						if(debug3)
//							for(int x=0;x<3;x++)
//								printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n",
//										I[x][0],I[x][1],I[x][2],I1[x][0],I1[x][1],I1[x][2],I3[x][0],I3[x][1],I3[x][2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						// printf("PSI P= %f,%f,%f\n",P[0],P[1],P[2]);
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							// ==> (3x1)
						}
						// printf("PSI elMY= %f,%f,%f\n",elMY[j][0],elMY[j][1],elMY[j][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						// printf("PSI emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						// printf("PSI el= %f,%f,%f\n",el[j][0],el[j][1],el[j][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						// printf("PSI em= %f,%f,%f\n",em[j][0],em[j][1],em[j][2]);
						ischi[j2] = false; // it's not chi
						j2++;
					}
					j++; // der[] dihedral index
				}  // NOT LAST

				if(has_oxt)
				{
					k2++; // set k2 into its correct value (last atom index)
					iter->pos_atom = k2;
				}
				has_oxt = false; // initialice for the next segment

			} // END PROTEIN FRAGMENT
			// RNA FRAGMENT
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				// ALPHA (bond between P and O5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (P, O1P, O2P)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (P, O1P, O2P)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;
					//					printf("PSI  mb1= %f  mb3= %f\n",mb1,mb3);

					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = coord[k1*3];
					y[1] = coord[k1*3 + 1];
					y[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (P, O1P, O2P)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("ALPHA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

					if(debug2)
						printf("ALPHA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					// printf("PSI P= %f,%f,%f\n",P[0],P[1],P[2]);
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("ALPHA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					// printf("PSI emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("ALPHA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("ALPHA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END ALPHA

				// BETA (bond between P and O5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (O5*)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (O5*)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+3; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;
					//					printf("BETA  mb1= %f  mb3= %f\n",mb1,mb3);

					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = coord[(k1+3)*3];
					y[1] = coord[(k1+3)*3 + 1];
					y[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (O5*)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+3; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("BETA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

					if(debug2)
						printf("BETA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					// printf("BETA P= %f,%f,%f\n",P[0],P[1],P[2]);
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("BETA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("BETA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("BETA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					if(debug2)
						printf("BETA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("BETA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("BETA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END BETA

				// GAMMA (bond between C5* and C4*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (P, O1P, O2P, O5*, C5*)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (P, O1P, O2P, O5*, C5*)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+4; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;

					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = coord[(k1+4)*3];
					y[1] = coord[(k1+4)*3 + 1];
					y[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (P, O1P, O2P, O5*, C5*)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+4; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("GAMMA M1= %f  M2= %f  after Y1= %f %f %f  mpr= %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2],mpr);
					if(debug2)
						printf("GAMMA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("GAMMA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("GAMMA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					if(debug2)
						printf("GAMMA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("GAMMA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("GAMMA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END GAMMA

				//  LATERAL CHAIN-->CHI
				if(type==2)
					{
						if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
						{
							//						printf ("CHI  j= %d  j2= %d\n",j,j2);

							// Calculate matrix I3
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I3[m][n] = 0;
							// Screens ALL body 3 (base)
							for ( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc();
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];
								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for ( m = 0; m < 3; m++ )
								{
									for ( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I3[m] [n] += mta * mfa;
									}
								}
							}
							// now matrix I1
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I1[m][n] = I[m][n]-I3[m][n];

							// MASSES
							// Screens ALL body 3 (base)
							mb3 = 0.0;
							for ( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
								mb3 += ( iter->get_atom() )->getPdbocc(); // atom mass
							// Now CB & -R is the whole 2nd body, the remaining atoms will be the 1st one!
							mb1 = mtot-mb3;
							// M1 & M3 arrays
							M1[j2] = mb1;
							M3[j2] = mb3; // Right-half = Tot - Left_half
							if(debug)
								printf("CHI  mb1= %f  mb3= %f\n",mb1,mb3);

							// get N1/N9 positions (N1,pyrimidines= +12) (N9,purines= +22 or +21)
							if(resn==ADE) // purine
								indexbase = 21;
							else if(resn==GUA) // purine
								indexbase = 22;
							else if(resn==CYT || resn==URA || resn == DCYT || resn == DTHY ) // pyrimidine
								indexbase = 12;
							else
							{
								printf("Unknown indexbase (N1/N9 index)\n");
								exit(3);
							}
							// get N1/N9 position
							e[0] = coord[(k1+indexbase)*3];
							e[1] = coord[(k1+indexbase)*3 + 1];
							e[2] = coord[(k1+indexbase)*3 + 2];
							// get C1* position (+11)
							y[0] = coord[(k1+11)*3];
							y[1] = coord[(k1+11)*3 + 1];
							y[2] = coord[(k1+11)*3 + 2];
							e[0] -= y[0]; // C1* --> N1/N9 unit vector
							e[1] -= y[1];
							e[2] -= y[2];

							temp = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
							for(m=0;m<3;m++)
							{
								e[m] /= temp; // unit vector normalization
								Y1[m] = 0.0;    /* to initialize */
								// Mon:
								// Not initialization!
								// The CoM of the whole protein is = 0.0 !!!
								// We'll do substract to the CoM (total) the contribution
								// of the 2nd body (CB-R). (See below)
							}

							// Screens ALL body 3 (base)
							for ( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc(); // atom mass
								// Substracting to the CoM (total) the body 2 (CB-R) contribution
								for(m=0; m<3; m++)
									Y1[m] -= mta * coord[iter->pos_atom*3 + m];
							}

							for(m=0;m<3;m++)
							{
								Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
								Y1[m] /= mb1;
							}

							if(debug2)
								printf("CHI y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
									,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

							// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							for(m=0; m<3; m++)
								N[m] = mb1 * Y1[m]; // M1 * Y1
							P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
							P[1] = e[2]*y[0] - e[0]*y[2];
							P[2] = e[0]*y[1] - e[1]*y[0];
							T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
								elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
								// ==> (3x1)
							}

							// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
							for(m=0; m<3; m++)
								N[m] = mb3 * Y3[m]; // M3 * Y3
							T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
								emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
							}

							for( m = 0; m < 3; m++)
							{
								Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
								Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
							}
							// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
							el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
							el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
							el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
							em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
							em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
							em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];

							ischi[j2] = true; // it is chi

							j2++;
						}

						j++; // der[] index
					}
				// END CHI

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{
					// EPSILON (bond between C5* and C4*)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1 (all except O3*)
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								// Do this with iterators: **********************
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];
								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for( m = 0; m < 3; m++ )
								{
									for( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I1[m][n] +=  mta * mfa;
									}
								}
							}
						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						// Updating M1 - body 1 mass (M1) (all except O3*)
						mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
								mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

						mb3 = mtot-mb1;
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3;

						// get O3* position (+8) (y_lambda)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = coord[(k1+7)*3];
						y[1] = coord[(k1+7)*3 + 1];
						y[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m]; // to initialize with previous "masses*positions"
						}

						// Updating Body 1 Center of Mass: "Y1[m]" (all except O3*)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								Y1[0] += mta * coord[iter->pos_atom*3];
								Y1[1] += mta * coord[iter->pos_atom*3+1];
								Y1[2] += mta * coord[iter->pos_atom*3+2];
							}

						for(m=0;m<3;m++) // 3D coords
						{
							Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y1[m] /= mb1;
						}
						if(debug)
							printf("EPSILON M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

						if(debug2)
							printf("EPSILON y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
								,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							// ==> (3x1)
						}
						if(debug2)
							printf("EPSILON elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("GAMMA T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						if(debug2)
							printf("EPSILON emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						if(debug2)
							printf("EPSILON el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						if(debug2)
							printf("EPSILON em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

						ischi[j2] = false; // it's not chi

						j2++;
					}
					j++; // der[] dihedral index
					// END EPSILON

					// ZETA (bond between C5* and C4*)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1 (all)
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc(); // mass
							// Do this with iterators: **********************
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for( m = 0; m < 3; m++ )
							{
								for( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
							}
						}
						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						// Updating M1 - body 1 mass (M1) (all)
						mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

						mb3 = mtot-mb1;
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3;

						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = coord[(k1+8)*3];
						y[1] = coord[(k1+8)*3 + 1];
						y[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m]; // to initialize with previous "masses*positions"
						}

						// Updating Body 1 Center of Mass: "Y1[m]" (all)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc(); // mass
							Y1[0] += mta * coord[iter->pos_atom*3];
							Y1[1] += mta * coord[iter->pos_atom*3+1];
							Y1[2] += mta * coord[iter->pos_atom*3+2];
						}
						if(debug2)
							printf("ZETA M1= %f  M3= %f  before Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

						for(m=0;m<3;m++) // 3D coords
						{
							Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y1[m] /= mb1;
						}
						if(debug)
							printf("ZETA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);
						if(debug2)
							printf("ZETA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
								,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							// ==> (3x1)
						}
						if(debug2)
							printf("ZETA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("ZETA T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						if(debug2)
							printf("ZETA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						if(debug2)
							printf("ZETA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						if(debug2)
							printf("ZETA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

						ischi[j2] = false; // it's not chi

						j2++;
					}
					j++; // der[] dihedral index
					// END ZETA
				}
			} // END RNA
			// DNA FRAGMENT
			else if( resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// ALPHA (bond between P and O5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (P, O1P, O2P)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (P, O1P, O2P)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;
					//					printf("PSI  mb1= %f  mb3= %f\n",mb1,mb3);

					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = coord[k1*3];
					y[1] = coord[k1*3 + 1];
					y[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (P, O1P, O2P)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("ALPHA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

					if(debug2)
						printf("ALPHA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					// printf("PSI P= %f,%f,%f\n",P[0],P[1],P[2]);
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("ALPHA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("PSI T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					// printf("PSI emMY= %f,%f,%f\n",emMY[j][0],emMY[j][1],emMY[j][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("ALPHA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("ALPHA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END ALPHA

				// BETA (bond between P and O5*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (O5*)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (O5*)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+3; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;
					//					printf("BETA  mb1= %f  mb3= %f\n",mb1,mb3);

					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = coord[(k1+3)*3];
					y[1] = coord[(k1+3)*3 + 1];
					y[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (O5*)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+3; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("BETA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

					if(debug2)
						printf("BETA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					// printf("BETA P= %f,%f,%f\n",P[0],P[1],P[2]);
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("BETA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("BETA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("BETA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					if(debug2)
						printf("BETA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("BETA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("BETA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END BETA

				// GAMMA (bond between C5* and C4*)
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0; // Initialization
					// Screens ALL body 1 (P, O1P, O2P, O5*, C5*)
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						// Do this with iterators: **********************
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
						{
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
						}
					}
					// now matrix I3
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I3[m] [n] = I[m] [n] - I1[m] [n];

					// Updating M1 - body 1 mass (M1) (P, O1P, O2P, O5*, C5*)
					mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+4; iter->next_atom() )
						mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

					mb3 = mtot-mb1;
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3;

					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = coord[(k1+4)*3];
					y[1] = coord[(k1+4)*3 + 1];
					y[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]; // to initialize with previous "masses*positions"
					}

					// Updating Body 1 Center of Mass: "Y1[m]" (P, O1P, O2P, O5*, C5*)
					for ( iter->pos_atom = k1; iter->pos_atom <= k1+4; iter->next_atom() )
					{
						mta = ( iter->get_atom() )->getPdbocc(); // mass
						Y1[0] += mta * coord[iter->pos_atom*3];
						Y1[1] += mta * coord[iter->pos_atom*3+1];
						Y1[2] += mta * coord[iter->pos_atom*3+2];
					}

					for(m=0;m<3;m++) // 3D coords
					{
						Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
						Y1[m] /= mb1;
					}
					if(debug)
						printf("GAMMA M1= %f  M2= %f  after Y1= %f %f %f  mpr= %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2],mpr);
					if(debug2)
						printf("GAMMA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
							,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
					for(m=0; m<3; m++)
						N[m] = mb1 * Y1[m]; // M1 * Y1
					P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
					P[1] = e[2]*y[0] - e[0]*y[2];
					P[2] = e[0]*y[1] - e[1]*y[0];
					T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
							L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
						elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						// ==> (3x1)
					}
					if(debug2)
						printf("GAMMA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

					// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
					for(m=0; m<3; m++)
						N[m] = mb3 * Y3[m]; // M3 * Y3
					T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
					T[1] = N[2]*P[0] - N[0]*P[2];
					T[2] = N[0]*P[1] - N[1]*P[0];
					// printf("GAMMA T= %f,%f,%f\n",T[0],T[1],T[2]);
					for(m=0; m<3; m++)
					{
						L[m] = 0;
						for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
							L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
						emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
					}
					if(debug2)
						printf("GAMMA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

					for( m = 0; m < 3; m++)
					{
						Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
						Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
					}
					// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
					el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
					el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
					el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
					if(debug2)
						printf("GAMMA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
					em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
					em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
					em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
					if(debug2)
						printf("GAMMA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

					ischi[j2] = false; // it's not chi

					j2++;
				}
				j++; // der[] dihedral index
				// END GAMMA

				//  LATERAL CHAIN-->CHI
				if(type==2)
					{
						if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
						{
							//						printf ("CHI  j= %d  j2= %d\n",j,j2);

							// Calculate matrix I3
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I3[m][n] = 0;
							// Screens ALL body 3 (base)
							for ( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() ) // DNA
							{
								mta = ( iter->get_atom() )->getPdbocc();
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];
								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for ( m = 0; m < 3; m++ )
								{
									for ( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I3[m] [n] += mta * mfa;
									}
								}
							}
							// now matrix I1
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I1[m][n] = I[m][n]-I3[m][n];

							// MASSES
							// Screens ALL body 3 (base)
							mb3 = 0.0;
							for ( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() ) // DNA
								mb3 += ( iter->get_atom() )->getPdbocc(); // atom mass
							// Now CB & -R is the whole 2nd body, the remaining atoms will be the 1st one!
							mb1 = mtot-mb3;
							// M1 & M3 arrays
							M1[j2] = mb1;
							M3[j2] = mb3; // Right-half = Tot - Left_half
							if(debug)
								printf("CHI  mb1= %f  mb3= %f\n",mb1,mb3);

							// get N1/N9 positions (N1,pyrimidines= +11) (N9,purines= +21 or +20)
							if(resn==DADE) // purine
								indexbase = 20; // 21 in RNA
							else if(resn==DGUA) // purine
								indexbase = 21; // 22 in RNA
							else if(resn == DCYT || resn == DTHY ) // pyrimidine
								indexbase = 11; // 12 in RNA
							else
							{
								printf("Unknown indexbase (N1/N9 index)\n");
								exit(3);
							}
							// get N1/N9 position
							e[0] = coord[(k1+indexbase)*3];
							e[1] = coord[(k1+indexbase)*3 + 1];
							e[2] = coord[(k1+indexbase)*3 + 2];
							// get C1* position (+10) // DNA
							y[0] = coord[(k1+10)*3];
							y[1] = coord[(k1+10)*3 + 1];
							y[2] = coord[(k1+10)*3 + 2];
							e[0] -= y[0]; // C1* --> N1/N9 unit vector
							e[1] -= y[1];
							e[2] -= y[2];

							temp = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
							for(m=0;m<3;m++)
							{
								e[m] /= temp; // unit vector normalization
								Y1[m] = 0.0;    /* to initialize */
								// Mon:
								// Not initialization!
								// The CoM of the whole protein is = 0.0 !!!
								// We'll do substract to the CoM (total) the contribution
								// of the 2nd body (CB-R). (See below)
							}

							// Screens ALL body 3 (base)
							for ( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc(); // atom mass
								// Substracting to the CoM (total) the body 2 (CB-R) contribution
								for(m=0; m<3; m++)
									Y1[m] -= mta * coord[iter->pos_atom*3 + m];
							}

							for(m=0;m<3;m++)
							{
								Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
								Y1[m] /= mb1;
							}

							if(debug2)
								printf("CHI y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
									,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

							// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							for(m=0; m<3; m++)
								N[m] = mb1 * Y1[m]; // M1 * Y1
							P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
							P[1] = e[2]*y[0] - e[0]*y[2];
							P[2] = e[0]*y[1] - e[1]*y[0];
							T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
								elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
								// ==> (3x1)
							}

							// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
							for(m=0; m<3; m++)
								N[m] = mb3 * Y3[m]; // M3 * Y3
							T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
							T[1] = N[2]*P[0] - N[0]*P[2];
							T[2] = N[0]*P[1] - N[1]*P[0];
							for(m=0; m<3; m++)
							{
								L[m] = 0;
								for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
									L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
								emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
							}

							for( m = 0; m < 3; m++)
							{
								Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
								Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
							}
							// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
							el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
							el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
							el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
							em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
							em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
							em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];

							ischi[j2] = true; // it is chi

							j2++;
						}

						j++; // der[] index
					}
				// END CHI

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{
					// EPSILON (bond between C5* and C4*)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1 (all except O3*)
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								// Do this with iterators: **********************
								r[0] = coord[iter->pos_atom*3];
								r[1] = coord[iter->pos_atom*3+1];
								r[2] = coord[iter->pos_atom*3+2];
								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
								for( m = 0; m < 3; m++ )
								{
									for( n = 0; n < 3; n++ )
									{
										mfa = -r[m] * r[n];
										if ( m == n ) mfa += nr2;
										I1[m][n] +=  mta * mfa;
									}
								}
							}
						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						// Updating M1 - body 1 mass (M1) (all except O3*)
						mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
								mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

						mb3 = mtot-mb1;
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3;

						// get O3* position (+8) (y_lambda)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = coord[(k1+7)*3];
						y[1] = coord[(k1+7)*3 + 1];
						y[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m]; // to initialize with previous "masses*positions"
						}

						// Updating Body 1 Center of Mass: "Y1[m]" (all except O3*)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							if(iter->pos_atom != k1+8) // not O3*
							{
								mta = ( iter->get_atom() )->getPdbocc(); // mass
								Y1[0] += mta * coord[iter->pos_atom*3];
								Y1[1] += mta * coord[iter->pos_atom*3+1];
								Y1[2] += mta * coord[iter->pos_atom*3+2];
							}

						for(m=0;m<3;m++) // 3D coords
						{
							Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y1[m] /= mb1;
						}
						if(debug)
							printf("EPSILON M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

						if(debug2)
							printf("EPSILON y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
								,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							// ==> (3x1)
						}
						if(debug2)
							printf("EPSILON elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("GAMMA T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						if(debug2)
							printf("EPSILON emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						if(debug2)
							printf("EPSILON el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						if(debug2)
							printf("EPSILON em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

						ischi[j2] = false; // it's not chi

						j2++;
					}
					j++; // der[] dihedral index
					// END EPSILON

					// ZETA (bond between C5* and C4*)
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1 (all)
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc(); // mass
							// Do this with iterators: **********************
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for( m = 0; m < 3; m++ )
							{
								for( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
							}
						}
						// now matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						// Updating M1 - body 1 mass (M1) (all)
						mb1 = mpr; // "mpr" --> mass of previous residue (=0 the 1st time)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
							mb1 += ( iter->get_atom() )->getPdbocc(); // adding non-C atoms mass

						mb3 = mtot-mb1;
						// M1 & M3 arrays
						M1[j2] = mb1;
						M3[j2] = mb3;

						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = coord[(k1+8)*3];
						y[1] = coord[(k1+8)*3 + 1];
						y[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m]; // to initialize with previous "masses*positions"
						}

						// Updating Body 1 Center of Mass: "Y1[m]" (all)
						for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc(); // mass
							Y1[0] += mta * coord[iter->pos_atom*3];
							Y1[1] += mta * coord[iter->pos_atom*3+1];
							Y1[2] += mta * coord[iter->pos_atom*3+2];
						}
						if(debug2)
							printf("ZETA M1= %f  M3= %f  before Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);

						for(m=0;m<3;m++) // 3D coords
						{
							Y3[m] = -Y1[m]/mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
							Y1[m] /= mb1;
						}
						if(debug)
							printf("ZETA M1= %f  M3= %f  after Y1= %f %f %f\n",mb1,mb3,Y1[0],Y1[1],Y1[2]);
						if(debug2)
							printf("ZETA y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
								,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

						// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
						for(m=0; m<3; m++)
							N[m] = mb1 * Y1[m]; // M1 * Y1
						P[0] = e[1]*y[2] - e[2]*y[1]; // e_lambda  x  y_lambda
						P[1] = e[2]*y[0] - e[0]*y[2];
						P[2] = e[0]*y[1] - e[1]*y[0];
						T[0] = N[1]*P[2] - N[2]*P[1]; // M1 * Y1  x  (e_lambda  x  y_lambda)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (3x1) = (3x1)
								L[m] += I1[m][n] * e[n]; // I1 * e_lambda ==> (3x1)
							elMY[j2][m] = T[m] - L[m]; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
							// ==> (3x1)
						}
						if(debug2)
							printf("ZETA elMY= %f,%f,%f\n",elMY[j2][0],elMY[j2][1],elMY[j2][2]);

						// emMY[size]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
						for(m=0; m<3; m++)
							N[m] = mb3 * Y3[m]; // M3 * Y3
						T[0] = N[1]*P[2] - N[2]*P[1]; // M3 * Y3  x  (e_mu  x  y_mu)
						T[1] = N[2]*P[0] - N[0]*P[2];
						T[2] = N[0]*P[1] - N[1]*P[0];
						// printf("ZETA T= %f,%f,%f\n",T[0],T[1],T[2]);
						for(m=0; m<3; m++)
						{
							L[m] = 0;
							for(n=0; n<3; n++) // (3x3) * (1x3) = (1x3)
								L[m] += I3[m][n] * e[n]; // I3 * e_mu ==> (3x1)
							emMY[j2][m] = T[m] - L[m]; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu) ==> (3x1)
						}
						if(debug2)
							printf("ZETA emMY= %f,%f,%f\n",emMY[j2][0],emMY[j2][1],emMY[j2][2]);

						for( m = 0; m < 3; m++)
						{
							Y1[m] = y[m] - Y1[m]; // e_lambda  x  (y_lambda - Y1)
							Y3[m] = y[m] - Y3[m]; // e_mu  x  (y_mu - Y3)
						}
						// a x b = (a2b3-a3b2)i + (a3b1-a1b3)j + (a1b2-a2b1)k
						el[j2][0] = e[1]*Y1[2] - e[2]*Y1[1];
						el[j2][1] = e[2]*Y1[0] - e[0]*Y1[2];
						el[j2][2] = e[0]*Y1[1] - e[1]*Y1[0];
						if(debug2)
							printf("ZETA el= %f,%f,%f\n",el[j2][0],el[j2][1],el[j2][2]);
						em[j2][0] = e[1]*Y3[2] - e[2]*Y3[1];
						em[j2][1] = e[2]*Y3[0] - e[0]*Y3[2];
						em[j2][2] = e[0]*Y3[1] - e[1]*Y3[0];
						if(debug2)
							printf("ZETA em= %f,%f,%f\n",em[j2][0],em[j2][1],em[j2][2]);

						ischi[j2] = false; // it's not chi

						j2++;
					}
					j++; // der[] dihedral index
					// END ZETA
				}
			} // END DNA
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype != tmol_smol )
			{
			    // NOT FOUND MOL-TYPE
				printf("Msg(kineticMFAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			// Initialize for the next residue
			for ( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
			{
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mpr += mta;
				// "my[]" will calculate 1st body's CoM
				my[0] += mta*coord[iter->pos_atom*3];
				my[1] += mta*coord[iter->pos_atom*3+1];
				my[2] += mta*coord[iter->pos_atom*3+2];
			}
			if(debug)
			{
				printf("UPDATE  mpr= %f  my= %f %f %f   mb1= %f  mb3= %f\n",mpr,my[0],my[1],my[2],mb1,mb3);
				printf("\tY1= %f,%f,%f  Y3= %f,%f,%f\n"
						,my[0]/mpr,my[1]/mpr,my[2]/mpr,-my[0]/(mtot-mpr),-my[1]/(mtot-mpr),-my[2]/(mtot-mpr));
			}

			if( !( (iter_frag->pos_fragment == num_res-1) && (iter_seg->pos_segment == num_seg - 1) ) )
			{
				k1 += props[index_res].nat;
				k2 += props[index_res + 1].nat;
			}
			index_res++;
		} // end fragment iterator
		delete(iter_frag);
	} // end segment iterator
	iter_seg->clean_virtual();
	delete(iter_seg);
	delete(iter);
	free(coord); // free double precision centered coordinates

//	printf("\nEverything was pre-computed..........\n");
	for(i=0; i<size; i++)
	{
//		printf("%2d Mon's elMY= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",elMY[i][k]);
//		printf("\n");
//		printf("%2d Mon's emMY= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",emMY[i][k]);
//		printf("\n");
//		printf("%2d Mon's el= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",el[i][k]);
//		printf("\n");
//		printf("%2d Mon's em= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",em[i][k]);
//		printf("\n");

		// T = elMY * I^-1    ==> (1x3) * (3x3) = (1x3)
		// T = I^-1 * emMY   ==> (3x3) * (3x1) = (3x1)
		if(ischi[i]) // if it is CHI
		{
			for(m=0; m<3; m++)
			{
				T[m] = 0.0;
				for(n=0; n<3; n++)
					T[m] -= emMY[i][n] * J[m][n]; // "-=" is due to the sign reversal of "e" when Lambda=CHI
				EL[m] = em[i][m]; // stores the used "el" (now "em")
			}
			fm = -M3[i]; // stores the used "M1" (First Mass) (now "M3")
			// "-" is due to the sign reversal of "e" when Lambda=CHI
			//		printf("{I^-1 * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			//		printf("It's CHI, then -->  elMY = emMY\n");
			//		printf("{elMY * I^-1}= %f,%f,%f\n",T[0],T[1],T[2]);
		}
		else  // If it's not CHI
		{
			for(m=0; m<3; m++)
			{
				T[m] = 0.0;
				for(n=0; n<3; n++)
					T[m] += elMY[i][n] * J[m][n];
				EL[m] = el[i][m]; // stores used "el"
			}
			fm = M1[i]; // stores the used "M1" (First Mass)
			//		printf("{I^-1 * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			//	printf("{elMY * I^-1}= %f,%f,%f\n",T[0],T[1],T[2]);
		}

		isi = i * size;
		for(j=i+1; j<size; j++) // upper triangle only!
		{
			// P = T * emMY   ==> (1x3) * (3x1) = (1x1)
			temp = 0.0;
			temp1 = 0.0;
			for(m=0; m<3; m++)
			{
				//				temp += elMY[i][m] * T[m];
				temp += T[m] * emMY[j][m];
				temp1 += EL[m] * em[j][m];
				//			printf("{elMY * J * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			}
			Hmatrix[ i + j*(j+1)/2 ] = (fm*M3[j]/mtot)*temp1 + temp;
			//			mass_tr[ i + j*(j+1)/2 ] = mass_matrix[i*size + j];
			//			printf("i= %d  j= %d  M1= %f  M3= %f  mtot= %f  temp1= %f  temp= %f\n",i,j,M1[i],M3[j],mtot,temp1,temp);
		}

		j=i; // Diagonal case... (it's as usual...)
		temp = 0.0;
		temp1 = 0.0;
		for(m=0; m<3; m++)
		{
			T[m] = 0.0;
			for(n=0; n<3; n++)
				//				T[m] += J[m][n] * emMY[j][n];
				//				T[m] += elMY[i][n] * J[n][m];
				T[m] += elMY[i][n] * J[m][n];
			temp += T[m] * emMY[j][m];
			temp1 += el[i][m] * em[j][m];

		}

		Hmatrix[ i + j*(j+1)/2 ] = (M1[i]*M3[j]/mtot)*temp1 + temp;
		//			printf("i= %d  j= %d  M1= %f  M3= %f  mtot= %f  temp1= %f  temp= %f\n",i,j,M1[i],M3[j],mtot,temp1,temp);
	}

	free( M1 );
	free( M3 );
	for(int i=0; i<size; i++)
	{
		free( el[i] );
		free( em[i] );
		free( elMY[i] );
		free( emMY[i] );
	}
	free( el );
	free( em );
	free( elMY );
	free( emMY );
	free( ischi );
}

// Computes the Kinetic Energy Matrix (H) (Multi-Chain & CA-model) (Huge-systems)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMCAx( Macromolecule *mol, float *coord, tri *props, int size, floating **p_Hmatrix, int model, bool *fix)
{
	bool debug=false;
	double mtot,mta,mb1,mb3,nr2,mfa,fm;
	double r[3],my[3],y[3],e[3],Y1[3],Y3[3],T[3],EL[3];
	double I[3][3],J[3][3];
	double temp,temp1,temp2;
	int m, n, num_atoms, i, j2, k1, k0;
	long int j; // Needed for huge systems!
	long int sizel = size; // Needed for huge systems!
	int resn,num_res,num_seg;
	floating *Hmatrix;
	double **I1,**I3;
	TMOL fragtype;
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments

	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	I1 = (double **) malloc( sizeof(double *) * 3 );
	I3 = (double **) malloc( sizeof(double *) * 3 );
	for(i=0; i<3; i++)
	{
		I1[i] = (double *) malloc( sizeof(double) * 3);
		I3[i] = (double *) malloc( sizeof(double) * 3);
	}

	// (DOUBLE) array, size: (N*(N+1)/2)
	if( *p_Hmatrix == NULL ) // if not NULL
	{
		if( !(Hmatrix = (floating *) malloc( sizeof(floating) * sizel*(sizel+1)/2 ) ) )
		{
			printf("Msg(kineticM): Memory allocation failed!\nForcing exit\n");
			exit(1);
		}
		*p_Hmatrix = Hmatrix;
	}
	else
		Hmatrix = *p_Hmatrix;

	j = 0; // current dihedral index
	j2 = 0; // active dihedral index (screens all mobile dihedrals)
	mb1 = 0; // Mass of Body 1
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */

	// Computing the PDB's Center of Mass (CoM) (in Inertia Matrix computation, see below)
	float com[3];
	com[0] = com[1] = com[2] = 0.0; // center of mass

	// Initializing Inertia-Matrix
	for ( int i = 0; i < 3; i++ )
		for ( int k = 0; k < 3; k++ )
			I[i][k] = 0.0;

	mtot = 0.0;
	k1 = 0; // residue index
	k0 = 0; // iters CA-model atoms (first-NH and last-CO included)
	iter_seg = new pdbIter(mol);

	// Compute Inertia Matrix (Inertia Momentum) (coords CoM must be centered 0,0,0)
	// Ec. 22 (Noguti & Go, 1983)
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		iter_frag = new pdbIter( iter_seg->get_segment() );
		num_res = iter_frag->num_fragment();
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			// first segment fragment has NH
			if(iter_frag->pos_fragment == 0)
			{
				iter->pos_atom = k0; // NH index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k1*3) * 3]; // NH position in coord
				r[1] = coord[(k1*3) * 3 + 1];
				r[2] = coord[(k1*3) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}

				// computing CoM
				if(debug)
				{
					com[0] += mta * r[0];
					com[1] += mta * r[1];
					com[2] += mta * r[2];
				}

				k0++;
			}

			// CA atom
			iter->pos_atom = k0; // CA index
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			mtot += mta; // total mass
			r[0] = coord[(k1*3+1) * 3]; // NH position in coord
			r[1] = coord[(k1*3+1) * 3 + 1];
			r[2] = coord[(k1*3+1) * 3 + 2];
			temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
			for ( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
				{
					temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
					if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
					I[m][n] += mta * temp2;
				}

			// computing CoM
			if(debug)
			{
				com[0] += mta * r[0];
				com[1] += mta * r[1];
				com[2] += mta * r[2];
			}

			k0++;

			// Last segment fragment has CO
			if(iter_frag->pos_fragment == num_res-1)
			{
				iter->pos_atom = k0; // CO index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k1*3+2) * 3]; // CO position in coord
				r[1] = coord[(k1*3+2) * 3 + 1];
				r[2] = coord[(k1*3+2) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}

				// computing CoM
				if(debug)
				{
					com[0] += mta * r[0];
					com[1] += mta * r[1];
					com[2] += mta * r[2];
				}

				k0++;
			}

			k1++;
		}
		delete iter_frag;
	}
	delete iter_seg;

	// Compute CoM
	if(debug)
	{
		com[0] /= mtot;
		com[1] /= mtot;
		com[2] /= mtot;
		printf( "Msg(kineticMCAx): Mass %8.8f Center %8.8f %8.8f %8.8f \n",mtot,com[0],com[1],com[2]);
	}

	// Compute the inverse of I
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	if(debug)
	{
		printf("TOTAL I-Inertia-Matrix:   mtot= %f\n",mtot);
		for (int k = 0; k < 3; k++ )
		{
			for (int l = 0; l < 3; l++ )
				printf("%10.4f ",I[k][l]);
			printf("\n");
		}
	}

	double *M1; // body-1 mass
	double *M3; // body-3 mass
	double **el; // e_lambda  x  (y_lambda - Y1)
	double **em; // e_mu  x  (y_mu - Y3)
	double **elMY; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	double **emMY; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)

	M1 = (double *) malloc( sizeof(double) * size ); // body-1 mass
	M3 = (double *) malloc( sizeof(double) * size ); // body-3 mass
	el = (double **) malloc( sizeof(double *) * size ); // e_lambda  x  (y_lambda - Y1)
	em = (double **) malloc( sizeof(double *) * size ); // e_mu  x  (y_mu - Y3)
	elMY = (double **) malloc( sizeof(double *) * size ); // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	emMY = (double **) malloc( sizeof(double *) * size ); // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
	for(int i=0; i<size; i++)
	{
		el[i] = (double *) malloc( sizeof(double) * 3);
		em[i] = (double *) malloc( sizeof(double) * 3);
		elMY[i] = (double *) malloc( sizeof(double) * 3);
		emMY[i] = (double *) malloc( sizeof(double) * 3);
	}

	// Calculate matrix I1 (Inertia Matrix 1)
	for(m=0;m<3;m++)
		for(n=0;n<3;n++)
			I1[m][n] = 0.0;

	// Pre-Computing M1, M3, el, em, elMY, emMY
	// Screen ALL residues
	k1 = 0; // residue index
	k0 = 0; // iters CA-model atoms (first-NH and last-CO included)
	Segment * seg;
	Atom * atom;
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(debug)
			printf ("\nProcessing segment %d  num_res= %d  k1 %d k0 %d ***************\n"
					,iter_seg->pos_segment,num_res,k1,k0);

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			//			printf("Entering segment %d\n",iter_seg->pos_segment);

			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				// Computing MASSES (updating)
				// M1 and M2
				mb3 = mtot - mb1; // Right-half = Tot - Left_half
				if(debug)
					printf("ROT-TRANS pos_segment= %d  j= %d  mb1= %f  mb3= %f\n",iter_seg->pos_segment,j,mb1,mb3);

				// y_lambda (NH pos 0) --> NOT EXIST FOR TRASLATION
				// e_lambda --> NOT EXIST FOR TRASLATION

				// Computing CoM
				for ( m = 0; m < 3; m++ )
				{
					// [ (Previous CoM) * (Previous mass) + (NH mass)*(NH pos) ] / mb1 = body 1 CoM
					Y1[m] = my[m] / mb1;
					// [ (Remaining CoM) * (remaining mass) - (NH mass)*(NH pos) ] / mb3 = body 3 CoM
					Y3[m] = -my[m] / mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
				}
				//				printf("TRANS y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
				//						,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

				// Calculate Inertia Matrix I3
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I3[m][n] = I[m][n] - I1[m][n];

				if(debug)
				{
					printf("ROT  I1-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%10.4f ",I1[k][l]);
						printf("\n");
					}
					printf("ROT  I3-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%9.6f ",I3[k][l]);
						printf("\n");
					}
				}
				// End Inertia matrix computation
			}

			// Adding 3 TRASLATIONS
			for(int axis=0; axis<3; axis++) // 3 translations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3; // Right-half = Tot - Left_half

					if(debug)
						printf("TRANS pos_segment= %d  j= %d  j2= %d  mb1= %f  mb3= %f  k1= %d  k0= %d\n"
								,iter_seg->pos_segment,j,j2,mb1,mb3,k1,k0);

					// Table I. --> Braun et al. (1984)
					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda
					// elMY[size] = M1 * Y1  x  gamma_v
					for(m=0; m<3; m++)
					{
						y[m] = 0.0; // y_lambda (Phi)
						e[m] = 0.0; // e_lambda (Psi)
					}
					e[axis] = -1.0; // e_lambda (unit vector) (Psi)

					calcKine2(mb1,Y1,mb3,Y3,I1,I3,y,e,elMY[j2],emMY[j2],el[j2],em[j2]);
					j2++;
				}
				j++; // adding traslational degree of freedom
			} // 3 trans added

			// Adding 3 ROTATIONS
			for(int axis=0; axis<3; axis++) // 3 rotations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3; // Right-half = Tot - Left_half

					if(debug)
						printf("ROT pos_segment= %d  j= %d  j2= %d  mb1= %f  mb3= %f  k1= %d  k0= %d\n"
								,iter_seg->pos_segment,j,j2,mb1,mb3,k1,k0);

					// Centers of Mass
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // delta_v initialization
//						y[0] = 0.0;
						y[m] = 0.0;
					}
					e[axis] = 1.0; // delta_v

					calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,y,elMY[j2],emMY[j2],el[j2],em[j2]);
					j2++;
				}
				j++; // adding rotational degree of freedom
			} // 3 rots added
		}

		// Screening residues (fragments)
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

//			if(debug)
//				printf("FRAG pos_fragment= %d (tmol= %d) %s  k1= %d  k0= %d\n"
//						,iter_frag->pos_fragment,fragtype,res->getName(),k1,k0);

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{
				// First NH
				if( iter_frag->pos_fragment == 0  )
				{
					// Iter first segment NH
					// NH position
					iter->pos_atom = k0;
					mta = ( iter->get_atom() )->getPdbocc(); // get mass
					mb1 += mta; // accumulating body 1 mass
					// Calculate Inertia Matrix (I1)
					r[0] = coord[k1*3*3];
					r[1] = coord[k1*3*3+1];
					r[2] = coord[k1*3*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
					{
						my[m] += r[m] * mta; // Accumulating body 1 momentum
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
					}
					if(debug)
					{
						printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mtot-mb1,my[0],my[1],my[2]);
						printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
					}
					k0++;
				}

				// ("PHI-bond")
				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
				if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
				{
					// printf ("%d Not first not PRO -> PHI \n",j);
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						//					printf ("PHI  j= %d  j2= %d\n",j,j2);

						// Calculate Inertia Matrix (I3) from I1 (previously accumulated)
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m][n] = I[m][n] - I1[m][n];

						if(debug)
							printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I3[0][0],I3[0][1],I3[0][2],I3[1][0],I3[1][1],I3[1][2],I3[2][0],I3[2][1],I3[2][2]);

						// M1 & M3 arrays
						mb3 = mtot - mb1;
						M1[j2] = mb1;
						M3[j2] = mb3;
						// printf("PHI  mb1= %f  mb3= %f\n",mb1,mb3);

						// y_lambda (NH pos 0)
						y[0] = coord[k1*3*3];
						y[1] = coord[k1*3*3 + 1];
						y[2] = coord[k1*3*3 + 2];
						// CA pos 1
						e[0] = coord[(k1*3+1) * 3]; // CA position
						e[1] = coord[(k1*3+1) * 3 + 1];
						e[2] = coord[(k1*3+1) * 3 + 2];
						// e_lambda
						e[0] -= y[0]; // NH --> CA
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1;
							Y3[m] = -my[m] / mb3; // thanks to the pre-centering of the coordinates!
						}
//						if(debug)
//							printf("PHI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
//								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);
						if(debug)
							printf("PHI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);

						// r = Psi = e_lambda  x  y_lambda
						r[0] = e[1]*y[2] - e[2]*y[1];
						r[1] = e[2]*y[0] - e[0]*y[2];
						r[2] = e[0]*y[1] - e[1]*y[0];

						calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,r,elMY[j2],emMY[j2],el[j2],em[j2]);
						j2++;
					}
					j++; // der[] index
				}  // NOT FIRST NOT PRO

				// Current CA
				iter->pos_atom = k0;
				mta = ( iter->get_atom() )->getPdbocc(); // get mass
				mb1 += mta; // accumulating body 1 mass
				// Calculate Inertia Matrix (I1)
				r[0] = coord[(k1*3+1)*3];
				r[1] = coord[(k1*3+1)*3+1];
				r[2] = coord[(k1*3+1)*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
				{
					my[m] += r[m] * mta; // Accumulating body 1 momentum
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}
				}
				if(debug)
				{
					printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mtot-mb1,my[0],my[1],my[2]);
					printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
				}
				k0++; // CA

				// "PSI-bond" (In CA-only model, non-first and non-last)
				if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
				{
					//				printf("%d Not last residue -> PSI\n",j);

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						//					printf ("PSI  j= %d  j2= %d\n",j,j2);
						// Calculate matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						if(debug)
							printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I3[0][0],I3[0][1],I3[0][2],I3[1][0],I3[1][1],I3[1][2],I3[2][0],I3[2][1],I3[2][2]);

						// M1 & M3 arrays
						mb3 = mtot - mb1;
						M1[j2] = mb1;
						M3[j2] = mb3;

						//					printf("PSI  mb1= %f  mb3= %f\n",mb1,mb3);
						// get C pos 2
						e[0] = coord[(k1*3+2) * 3];
						e[1] = coord[(k1*3+2) * 3 + 1];
						e[2] = coord[(k1*3+2) * 3 + 2];
						// get CA pos 1
						y[0] = coord[(k1*3+1) * 3]; // CA position
						y[1] = coord[(k1*3+1) * 3 + 1];
						y[2] = coord[(k1*3+1) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= y[0];
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1;
							Y3[m] = -my[m] / mb3;
						}
//						if(debug)
//							printf("PSI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
//								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);
						if(debug)
							printf("PSI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);

						// r = Psi = e_lambda  x  y_lambda
						r[0] = e[1]*y[2] - e[2]*y[1];
						r[1] = e[2]*y[0] - e[0]*y[2];
						r[2] = e[0]*y[1] - e[1]*y[0];

						calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,r,elMY[j2],emMY[j2],el[j2],em[j2]);

						j2++;
					}
					j++; // dihedral index
				}

				// Last CO
				if( iter_frag->pos_fragment == num_res-1 )
				{
					// CO position
					iter->pos_atom = k0;
					mta = ( iter->get_atom() )->getPdbocc(); // get mass
					mb1 += mta; // accumulating body 1 mass
					// Calculate Inertia Matrix (I1)
					r[0] = coord[(k1*3+2)*3];
					r[1] = coord[(k1*3+2)*3+1];
					r[2] = coord[(k1*3+2)*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
					{
						my[m] += r[m] * mta; // Accumulating body 1 momentum
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
					}
					if(debug)
					{
						printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mtot-mb1,my[0],my[1],my[2]);
						printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
					}
					k0++;
				}

			} // END PROTEIN FRAGMENT
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(kineticMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			// Initialize for the next residue
			k1++; // Residue index (CA index)
		} // end fragment iterator
		delete(iter_frag);
	} // end segment iterator
	delete(iter_seg);
	delete(iter);

	if(debug)
		printf("\nEverything was pre-computed..........\n");
	for(i=0; i<size; i++)
	{
//		printf("%2d Mon's elMY= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",elMY[i][k]);
//		printf("\n");
//		printf("%2d Mon's emMY= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",emMY[i][k]);
//		printf("\n");
//		printf("%2d Mon's el= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",el[i][k]);
//		printf("\n");
//		printf("%2d Mon's em= ",i);
//		for (int k = 0; k < 3; k++ )
//			printf("%10.4f ",em[i][k]);
//		printf("\n");

		// T = elMY * I^-1    ==> (1x3) * (3x3) = (1x3)
		// T = I^-1 * emMY   ==> (3x3) * (3x1) = (3x1)
		// If it's not CHI
		for(m=0; m<3; m++)
		{
			T[m] = 0.0;
			for(n=0; n<3; n++)
				T[m] += elMY[i][n] * J[m][n];
			EL[m] = el[i][m]; // stores used "el"
		}
		fm = M1[i]; // stores the used "M1" (First Mass)
		if(debug)
		{
//			printf("{I^-1 * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			printf("{elMY * I^-1}= %f,%f,%f\n",T[0],T[1],T[2]);
		}

		for(j=i+1; j<size; j++) // upper triangle only!
		{
			// P = T * emMY   ==> (1x3) * (3x1) = (1x1)
			temp = 0.0;
			temp1 = 0.0;
			for(m=0; m<3; m++)
			{
				//				temp += elMY[i][m] * T[m];
				temp += T[m] * emMY[j][m];
				temp1 += EL[m] * em[j][m];
				//			printf("{elMY * J * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			}
			Hmatrix[ i + j*(j+1)/2 ] = (fm*M3[j]/mtot)*temp1 + temp;
			if(debug)
				printf("i= %d  j= %d  M1= %f  M3= %f  mtot= %f  temp1= %f  temp= %f\n",i,j,M1[i],M3[j],mtot,temp1,temp);
		}

		j=i; // Diagonal case... (it's as usual...)
		temp = 0.0;
		temp1 = 0.0;
		for(m=0; m<3; m++)
		{
			T[m] = 0.0;
			for(n=0; n<3; n++)
				T[m] += elMY[i][n] * J[m][n];
			temp += T[m] * emMY[j][m];
			temp1 += el[i][m] * em[j][m];

		}
		Hmatrix[ i + j*(j+1)/2 ] = (M1[i]*M3[j]/mtot)*temp1 + temp;
		if(debug)
			printf("i= %d  j= %d  M1= %f  M3= %f  mtot= %f  temp1= %f  temp= %f\n",i,j,M1[i],M3[j],mtot,temp1,temp);
	}

	for(int i=0; i<3; i++)
	{
		free( I1[i] );
		free( I3[i] );
	}
	free( I1 );
	free( I3 );
	free( M1 );
	free( M3 );
	for(int i=0; i<size; i++)
	{
		free( el[i] );
		free( em[i] );
		free( elMY[i] );
		free( emMY[i] );
	}
	free( el );
	free( em );
	free( elMY );
	free( emMY );
}


// Computes the Kinetic Energy Matrix (H) (Multi-Chain & CA-model) (Huge-systems)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMCAxHD( Macromolecule *mol, float *coord2, tri *props, int size, floating **p_Hmatrix, int model, bool *fix, char *file, bool justfile)
{
	bool debug = false;
	double mtot,mta,mb1,mb3,nr2,mfa,fm;
	double r[3],my[3],y[3],e[3],Y1[3],Y3[3],T[3],EL[3];
	double I[3][3],J[3][3];
	double temp,temp1,temp2,kine;
	int m, n, num_atoms, i, j2, k1, k0;
	long int j; // Needed for huge systems!
	long int sizel = size; // Needed for huge systems!
	int resn,num_res,num_seg;
	floating *Hmatrix; // the Kinetic Energy matrix
	double **I1,**I3;
	TMOL fragtype;
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments

	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	I1 = (double **) malloc( sizeof(double *) * 3 );
	I3 = (double **) malloc( sizeof(double *) * 3 );
	for(i=0; i<3; i++)
	{
		I1[i] = (double *) malloc( sizeof(double) * 3);
		I3[i] = (double *) malloc( sizeof(double) * 3);
	}

	FILE *f_file = NULL;
	if(file != NULL) // the Hessian will be dumped into a HD disk file
	{
		if( !(f_file = fopen(file, "w")) ) // if file creation error...
		{
			printf("Msg(kineticMCAxHD): Sorry, Kinetic-Energy Matrix file creation failed!\nForcing exit!\n");
			exit(1);
		}
		double dsize = (double) size;
		fwrite(&dsize,sizeof(double),1,f_file); // first number is the matrix rank
	}

	if(!justfile) // the matrix will be dumped into RAM.
	{
		// (DOUBLE) array, size: (N*(N+1)/2)
		if( *p_Hmatrix == NULL ) // if not NULL
		{
			if( !(Hmatrix = (floating *) malloc( sizeof(floating) * sizel*(sizel+1)/2 ) ) )
			{
				printf("Msg(kineticM): Memory allocation failed!\nForcing exit\n");
				exit(1);
			}
			*p_Hmatrix = Hmatrix;
		}
		else
			Hmatrix = *p_Hmatrix;
	}

	j = 0; // current dihedral index
	j2 = 0; // active dihedral index (screens all mobile dihedrals)
	mb1 = 0; // Mass of Body 1
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */

	// Computing the PDB's Center of Mass (CoM) (in Intertia Matrix computation, see below)
	double com[3];
	com[0] = com[1] = com[2] = 0.0; // center of mass

	// Initializing Inertia-Matrix
	for ( int i = 0; i < 3; i++ )
		for ( int k = 0; k < 3; k++ )
			I[i][k] = 0.0;

	mtot = 0.0;
	k1 = 0; // residue index
	k0 = 0; // iters CA-model atoms (first-NH and last-CO included)
	iter_seg = new pdbIter(mol);

	// Centering with double precision (necessary with huge systems!!! <-- the fucking bug!!!)
	double com2[3];
	double *coord;
	com2[0] = com2[1] = com2[2] = 0.0; // center of mass
	if( !(coord = (double *) malloc(sizeof(double)*3*num_res*3)) )
	{
		printf("Msg(kineticMCA_disc): Unable to allocate memory! Forcing exit!\n");
		exit(1);
	}

	// Computing "double precision" center of mass (coords CoM must be nicely centered at 0,0,0)
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		iter_frag = new pdbIter( iter_seg->get_segment() );
		num_res = iter_frag->num_fragment();
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			// first segment fragment has NH
			if(iter_frag->pos_fragment == 0)
			{
				iter->pos_atom = k0; // NH index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				com2[0] += mta * coord2[(k1*3) * 3]; // NH position in coord
				com2[1] += mta * coord2[(k1*3) * 3 + 1];
				com2[2] += mta * coord2[(k1*3) * 3 + 2];

				k0++;
			}

			// CA atom
			iter->pos_atom = k0; // CA index
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			mtot += mta; // total mass
			com2[0] += mta * coord2[(k1*3+1) * 3]; // CA position in coord
			com2[1] += mta * coord2[(k1*3+1) * 3 + 1];
			com2[2] += mta * coord2[(k1*3+1) * 3 + 2];
			k0++;

			// Last segment fragment has CO
			if(iter_frag->pos_fragment == num_res-1)
			{
				iter->pos_atom = k0; // CO index
// fprintf(stderr,"Last segment fragment has CO, seg= %d  frag= %d  atom= %d\n",iter_seg->pos_segment,iter_frag->pos_fragment,iter->pos_atom);
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				com2[0] += mta * coord2[(k1*3+2) * 3]; // CO position in coord
				com2[1] += mta * coord2[(k1*3+2) * 3 + 1];
				com2[2] += mta * coord2[(k1*3+2) * 3 + 2];
				k0++;
			}

			k1++;
		}
		delete iter_frag;
	}
	com2[0] /= mtot;
	com2[1] /= mtot;
	com2[2] /= mtot;
	if(debug)
		printf( "Msg(kineticMCAx): Mass1 %18.15e  Center1 %18.15e %18.15e %18.15e \n",mtot,com2[0],com2[1],com2[2]);

	// Setting "double" centered coordinates
	num_res = mol->get_num_fragments();
	for(int x = 0; x < num_res; x++)
	{
		for(int y = 0; y < 3; y++)
		{
			coord[(x*3+y)*3] =   coord2[(x*3+y)*3] - com2[0];
			coord[(x*3+y)*3+1] =   coord2[(x*3+y)*3+1] - com2[1];
			coord[(x*3+y)*3+2] =   coord2[(x*3+y)*3+2] - com2[2];
		}
	}
	mtot = 0.0;
	k1 = 0; // residue index
	k0 = 0; // iters CA-model atoms (first-NH and last-CO included)

	// Compute Inertia Matrix (Inertia Momentum) (coords CoM must be centered 0,0,0)
	// Ec. 22 (Noguti & Go, 1983)
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		iter_frag = new pdbIter( iter_seg->get_segment() );
		num_res = iter_frag->num_fragment();
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			// first segment fragment has NH
			if(iter_frag->pos_fragment == 0)
			{
				iter->pos_atom = k0; // NH index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k1*3) * 3]; // NH position in coord
				r[1] = coord[(k1*3) * 3 + 1];
				r[2] = coord[(k1*3) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}

				// computing CoM
				if(debug)
				{
					com[0] += mta * r[0];
					com[1] += mta * r[1];
					com[2] += mta * r[2];
				}

				k0++;
			}

			// CA atom
			iter->pos_atom = k0; // CA index
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			mtot += mta; // total mass
			r[0] = coord[(k1*3+1) * 3]; // CA position in coord
			r[1] = coord[(k1*3+1) * 3 + 1];
			r[2] = coord[(k1*3+1) * 3 + 2];
			temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
			for ( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
				{
					temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
					if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
					I[m][n] += mta * temp2;
				}

			// computing CoM
			if(debug)
			{
				com[0] += mta * r[0];
				com[1] += mta * r[1];
				com[2] += mta * r[2];
			}

			k0++;

			// Last segment fragment has CO
			if(iter_frag->pos_fragment == num_res-1)
			{
				iter->pos_atom = k0; // CO index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k1*3+2) * 3]; // CO position in coord
				r[1] = coord[(k1*3+2) * 3 + 1];
				r[2] = coord[(k1*3+2) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}

				// computing CoM
				if(debug)
				{
					com[0] += mta * r[0];
					com[1] += mta * r[1];
					com[2] += mta * r[2];
				}

				k0++;
			}

			k1++;
		}
		delete iter_frag;
	}
	delete iter_seg;

	// Compute CoM
	if(debug)
	{
		com[0] /= mtot;
		com[1] /= mtot;
		com[2] /= mtot;
		printf( "Msg(kineticMCAx): Mass %8.8f Center %8.8f %8.8f %8.8f \n",mtot,com[0],com[1],com[2]);
	}

	// Compute the inverse of I
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	if(debug)
	{
		printf("TOTAL I-Inertia-Matrix:   mtot= %f\n",mtot);
		for (int k = 0; k < 3; k++ )
		{
			for (int l = 0; l < 3; l++ )
				printf("%10.4f ",I[k][l]);
			printf("\n");
		}
	}

	double *M1; // body-1 mass
	double *M3; // body-3 mass
	double **el; // e_lambda  x  (y_lambda - Y1)
	double **em; // e_mu  x  (y_mu - Y3)
	double **elMY; // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	double **emMY; // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)

	M1 = (double *) malloc( sizeof(double) * size ); // body-1 mass
	M3 = (double *) malloc( sizeof(double) * size ); // body-3 mass
	el = (double **) malloc( sizeof(double *) * size ); // e_lambda  x  (y_lambda - Y1)
	em = (double **) malloc( sizeof(double *) * size ); // e_mu  x  (y_mu - Y3)
	elMY = (double **) malloc( sizeof(double *) * size ); // M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda)
	emMY = (double **) malloc( sizeof(double *) * size ); // M3 * Y3  x  (e_mu  x  y_mu) - I3 * e_mu)
	for(int i=0; i<size; i++)
	{
		el[i] = (double *) malloc( sizeof(double) * 3);
		em[i] = (double *) malloc( sizeof(double) * 3);
		elMY[i] = (double *) malloc( sizeof(double) * 3);
		emMY[i] = (double *) malloc( sizeof(double) * 3);
	}

	// Calculate matrix I1 (Inertia Matrix 1)
	for(m=0;m<3;m++)
		for(n=0;n<3;n++)
			I1[m][n] = 0.0;

	// Pre-Computing M1, M3, el, em, elMY, emMY
	// Screen ALL residues
	k1 = 0; // residue index
	k0 = 0; // iters CA-model atoms (first-NH and last-CO included)
	Segment * seg;
	Atom * atom;
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

//		if(iter_seg->pos_segment == 719)
//			debug = true;
//		else
//			debug = false;

		if(debug)
			printf ("\nProcessing segment %d  num_res= %d  k1 %d k0 %d ***************\n"
					,iter_seg->pos_segment,num_res,k1,k0);

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			//			printf("Entering segment %d\n",iter_seg->pos_segment);

			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				// Computing MASSES (updating)
				// M1 and M2
				mb3 = mtot - mb1; // Right-half = Tot - Left_half
				if(debug)
					printf("ROT-TRANS pos_segment= %d  j= %d  mb1= %f  mb3= %f\n",iter_seg->pos_segment,j,mb1,mb3);

				// y_lambda (NH pos 0) --> NOT EXIST FOR TRASLATION
				// e_lambda --> NOT EXIST FOR TRASLATION

				// Computing CoM
				for ( m = 0; m < 3; m++ )
				{
					// [ (Previous CoM) * (Previous mass) + (NH mass)*(NH pos) ] / mb1 = body 1 CoM
					Y1[m] = my[m] / mb1;
					// [ (Remaining CoM) * (remaining mass) - (NH mass)*(NH pos) ] / mb3 = body 3 CoM
					Y3[m] = -my[m] / mb3; // Y3 = (mtot*CoM - mb1*Y1) / mb3  (given CoM=0,0,0)
				}
				//				printf("TRANS y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
				//						,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);

				// Calculate Inertia Matrix I3
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I3[m][n] = I[m][n] - I1[m][n];

				if(debug)
				{
					printf("ROT  I1-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%10.4f ",I1[k][l]);
						printf("\n");
					}
					printf("ROT  I3-Inertia-Matrix:\n");
					for (int k = 0; k < 3; k++ )
					{
						for (int l = 0; l < 3; l++ )
							printf("%9.6f ",I3[k][l]);
						printf("\n");
					}
				}
				// End Inertia matrix computation
			}

			// Adding 3 TRASLATIONS
			for(int axis=0; axis<3; axis++) // 3 translations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3; // Right-half = Tot - Left_half

					if(debug)
						printf("TRANS pos_segment= %d  j= %d  j2= %d  mb1= %f  mb3= %f  k1= %d  k0= %d\n"
								,iter_seg->pos_segment,j,j2,mb1,mb3,k1,k0);

					// Table I. --> Braun et al. (1984)
					// elMY[size] = M1 * Y1  x  (e_lambda  x  y_lambda) - I1 * e_lambda
					// elMY[size] = M1 * Y1  x  gamma_v
					for(m=0; m<3; m++)
					{
						y[m] = 0.0; // y_lambda (Phi)
						e[m] = 0.0; // e_lambda (Psi)
					}
					e[axis] = -1.0; // e_lambda (unit vector) (Psi)

					calcKine2(mb1,Y1,mb3,Y3,I1,I3,y,e,elMY[j2],emMY[j2],el[j2],em[j2]);
					j2++;
				}
				j++; // adding traslational degree of freedom
			} // 3 trans added

			// Adding 3 ROTATIONS
			for(int axis=0; axis<3; axis++) // 3 rotations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// M1 & M3 arrays
					M1[j2] = mb1;
					M3[j2] = mb3; // Right-half = Tot - Left_half

					if(debug)
						printf("ROT pos_segment= %d  j= %d  j2= %d  mb1= %f  mb3= %f  k1= %d  k0= %d\n"
								,iter_seg->pos_segment,j,j2,mb1,mb3,k1,k0);

					// Centers of Mass
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // delta_v initialization
//						y[0] = 0.0;
						y[m] = 0.0;
					}
					e[axis] = 1.0; // delta_v

					calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,y,elMY[j2],emMY[j2],el[j2],em[j2]);
					j2++;
				}
				j++; // adding rotational degree of freedom
			} // 3 rots added
		}

		// Screening residues (fragments)
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

//			if(debug)
//				printf("FRAG pos_fragment= %d (tmol= %d) %s  k1= %d  k0= %d\n"
//						,iter_frag->pos_fragment,fragtype,res->getName(),k1,k0);

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{
				// First NH
				if( iter_frag->pos_fragment == 0  )
				{
					// Iter first segment NH
					// NH position
					iter->pos_atom = k0;
					mta = ( iter->get_atom() )->getPdbocc(); // get mass
					mb1 += mta; // accumulating body 1 mass
					// Calculate Inertia Matrix (I1)
					r[0] = coord[k1*3*3];
					r[1] = coord[k1*3*3+1];
					r[2] = coord[k1*3*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
					{
						my[m] += r[m] * mta; // Accumulating body 1 momentum
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
					}
					if(debug)
					{
						printf("k1= %d  k0= %d  mb1= %f  mb2= %18.15e  my= %f %f %f\n",k1,k0,mb1,mtot-mb1,my[0],my[1],my[2]);
						printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
					}
					k0++;
				}

				// ("PHI-bond")
				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
				if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
				{
					// printf ("%d Not first not PRO -> PHI \n",j);
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						//					printf ("PHI  j= %d  j2= %d\n",j,j2);

						// Calculate Inertia Matrix (I3) from I1 (previously accumulated)
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m][n] = I[m][n] - I1[m][n];

						if(debug)
							printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I3[0][0],I3[0][1],I3[0][2],I3[1][0],I3[1][1],I3[1][2],I3[2][0],I3[2][1],I3[2][2]);

						// M1 & M3 arrays
						mb3 = mtot - mb1;
						M1[j2] = mb1;
						M3[j2] = mb3;
						// printf("PHI  mb1= %f  mb3= %f\n",mb1,mb3);

						// y_lambda (NH pos 0)
						y[0] = coord[k1*3*3];
						y[1] = coord[k1*3*3 + 1];
						y[2] = coord[k1*3*3 + 2];
						// CA pos 1
						e[0] = coord[(k1*3+1) * 3]; // CA position
						e[1] = coord[(k1*3+1) * 3 + 1];
						e[2] = coord[(k1*3+1) * 3 + 2];
						// e_lambda
						e[0] -= y[0]; // NH --> CA
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1;
							Y3[m] = -my[m] / mb3; // thanks to the pre-centering of the coordinates!
						}
//						if(debug)
//							printf("PHI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
//								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);
						if(debug)
							printf("PHI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);

						// r = Psi = e_lambda  x  y_lambda
						r[0] = e[1]*y[2] - e[2]*y[1];
						r[1] = e[2]*y[0] - e[0]*y[2];
						r[2] = e[0]*y[1] - e[1]*y[0];

						calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,r,elMY[j2],emMY[j2],el[j2],em[j2]);
						j2++;
					}
					j++; // der[] index
				}  // NOT FIRST NOT PRO

				// Current CA
				iter->pos_atom = k0;
				mta = ( iter->get_atom() )->getPdbocc(); // get mass
				mb1 += mta; // accumulating body 1 mass
				// Calculate Inertia Matrix (I1)
				r[0] = coord[(k1*3+1)*3];
				r[1] = coord[(k1*3+1)*3+1];
				r[2] = coord[(k1*3+1)*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
				{
					my[m] += r[m] * mta; // Accumulating body 1 momentum
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}
				}
				if(debug)
				{
					printf("CA_k1= %d  k0= %d  mb1= %f mta= %f mb2= %18.15e  my= %f %f %f  r= %f %f %f\n",k1,k0,mb1,mta,mtot-mb1,my[0],my[1],my[2],r[0],r[1],r[2]);
					printf("CA_I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
				}
				k0++; // CA

				// "PSI-bond" (In CA-only model, non-first and non-last)
				if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
				{
					//				printf("%d Not last residue -> PSI\n",j);

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						//					printf ("PSI  j= %d  j2= %d\n",j,j2);
						// Calculate matrix I3
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I3[m] [n] = I[m] [n] - I1[m] [n];

						if(debug)
							printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I3[0][0],I3[0][1],I3[0][2],I3[1][0],I3[1][1],I3[1][2],I3[2][0],I3[2][1],I3[2][2]);

						// M1 & M3 arrays
						mb3 = mtot - mb1;
						M1[j2] = mb1;
						M3[j2] = mb3;

						//					printf("PSI  mb1= %f  mb3= %f\n",mb1,mb3);
						// get C pos 2
						e[0] = coord[(k1*3+2) * 3];
						e[1] = coord[(k1*3+2) * 3 + 1];
						e[2] = coord[(k1*3+2) * 3 + 2];
						// get CA pos 1
						y[0] = coord[(k1*3+1) * 3]; // CA position
						y[1] = coord[(k1*3+1) * 3 + 1];
						y[2] = coord[(k1*3+1) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= y[0];
						e[1] -= y[1];
						e[2] -= y[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1;
							Y3[m] = -my[m] / mb3;
						}
//						if(debug)
//							printf("PSI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f  Y3= %f,%f,%f\n"
//								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2],Y3[0],Y3[1],Y3[2]);
						if(debug)
							printf("PSI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);

						// r = Psi = e_lambda  x  y_lambda
						r[0] = e[1]*y[2] - e[2]*y[1];
						r[1] = e[2]*y[0] - e[0]*y[2];
						r[2] = e[0]*y[1] - e[1]*y[0];

						calcKine2(mb1,Y1,mb3,Y3,I1,I3,e,r,elMY[j2],emMY[j2],el[j2],em[j2]);

						j2++;
					}
					j++; // dihedral index
				}

				// Last CO
				if( iter_frag->pos_fragment == num_res-1 )
				{
					// CO position
					iter->pos_atom = k0;
					mta = ( iter->get_atom() )->getPdbocc(); // get mass
					mb1 += mta; // accumulating body 1 mass
					// Calculate Inertia Matrix (I1)
					r[0] = coord[(k1*3+2)*3];
					r[1] = coord[(k1*3+2)*3+1];
					r[2] = coord[(k1*3+2)*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
					{
						my[m] += r[m] * mta; // Accumulating body 1 momentum
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
					}
					if(debug)
					{
//						printf("CO_k1= %d  k0= %d  mb1= %f  mb2= %18.15e  my= %f %f %f  r= %f %f %f\n",k1,k0,mb1,mtot-mb1,my[0],my[1],my[2],r[0],r[1],r[2]);
						printf("CO_k1= %d  k0= %d  mb1= %f mta= %f mb2= %18.15e  my= %f %f %f  r= %f %f %f\n",k1,k0,mb1,mta,mtot-mb1,my[0],my[1],my[2],r[0],r[1],r[2]);
						printf("CO_I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
					}
					k0++;
				}

			} // END PROTEIN FRAGMENT
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(kineticMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			// Initialize for the next residue
			k1++; // Residue index (CA index)
		} // end fragment iterator
		delete(iter_frag);
	} // end segment iterator
	delete(iter_seg);
	delete(iter);
	free(coord); // free double precision centered coordinates

	if(debug)
		printf("\nEverything was pre-computed..........\n");

//	printf("\n"); // <--- remove this

	for(i=0; i<size; i++)
	{
		if(debug)
		{
			printf("%5d Mon's elMY= ",i);
			for (int k = 0; k < 3; k++ )
				printf("%14f ",elMY[i][k]);
			printf("\n");
			printf("%5d Mon's emMY= ",i);
			for (int k = 0; k < 3; k++ )
				printf("%14f ",emMY[i][k]);
			printf("\n");
			printf("%5d Mon's el= ",i);
			for (int k = 0; k < 3; k++ )
				printf("%14f ",el[i][k]);
			printf("\n");
			printf("%5d Mon's em= ",i);
			for (int k = 0; k < 3; k++ )
				printf("%14f ",em[i][k]);
			printf("\n");
		}

		// T = elMY * I^-1    ==> (1x3) * (3x3) = (1x3)
		// T = I^-1 * emMY   ==> (3x3) * (3x1) = (3x1)
		// If it's not CHI
		for(m=0; m<3; m++)
		{
//			// WHY THIS ????
//			/////////////////////////////////
//			T[m] = 0.0;
//			for(n=0; n<3; n++)
//				T[m] += elMY[i][n] * J[m][n];
//			/////////////////////////////////

			EL[m] = el[i][m]; // stores current "el"
		}
		fm = M1[i]; // stores current "M1" (First Mass)

		j = (long) i; // Diagonal case... (it's as usual...)
		temp = 0.0;
		temp1 = 0.0;
		for(m=0; m<3; m++)
		{
			T[m] = 0.0;
			for(n=0; n<3; n++)
				T[m] += elMY[i][n] * J[m][n];
			temp += T[m] * emMY[j][m];
			temp1 += el[i][m] * em[j][m];

		}
		kine = (M1[i]*M3[j]/mtot)*temp1 + temp;

		if(!justfile) // the Kinetic Energy matrix will be dumped into RAM.
			Hmatrix[ i + j*(j+1)/2 ] = kine;
		if(file != NULL) // the Kinetic Energy matrix will be dumped into a HD disk file
			fwrite(&kine,sizeof(double),1,f_file);

		if(debug)
		{
			printf("J=  %18.15e %18.15e %18.15e  %18.15e %18.15e %18.15e  %18.15e %18.15e %18.15e\n",J[0][0],J[0][1],J[0][2],J[1][0],J[1][1],J[1][2],J[2][0],J[2][1],J[2][2]);
			printf("{elMY * I^-1}= %18.15e %18.15e %18.15e\n",T[0],T[1],T[2]);
			printf("i= %5d  j= %5ld  M1= %18.15e  M3= %18.15e  mtot= %18.15e  temp1= %18.15e  temp= %18.15e\n",i,j,M1[i],M3[j],mtot,temp1,temp);
			printf("==> K[i][i]= %18.15e\n\n",kine);
		}

		for(j=i+1; j<size; j++) // upper triangle only!
		{
			// P = T * emMY   ==> (1x3) * (3x1) = (1x1)
			temp = 0.0;
			temp1 = 0.0;
			for(m=0; m<3; m++)
			{
				//				temp += elMY[i][m] * T[m];
				temp += T[m] * emMY[j][m];
				temp1 += EL[m] * em[j][m];
				//			printf("{elMY * J * emMY}= %f,%f,%f\n",T[0],T[1],T[2]);
			}
			kine = (fm*M3[j]/mtot)*temp1 + temp;
			if(!justfile) // the Kinetic Energy matrix will be dumped into RAM.
				Hmatrix[ i + j*(j+1)/2 ] = kine;
			if(file != NULL) // the Kinetic Energy matrix will be dumped into a HD disk file
				fwrite(&kine,sizeof(double),1,f_file);

//			if(debug)
//				printf("i= %d  j= %d  M1= %f  M3= %f  mtot= %f  temp1= %f  temp= %f\n",i,j,M1[i],M3[j],mtot,temp1,temp);
		}

	}

	for(int i=0; i<3; i++)
	{
		free( I1[i] );
		free( I3[i] );
	}
	free( I1 );
	free( I3 );
	free( M1 );
	free( M3 );
	for(int i=0; i<size; i++)
	{
		free( el[i] );
		free( em[i] );
		free( elMY[i] );
		free( emMY[i] );
	}
	free( el );
	free( em );
	free( elMY );
	free( emMY );
}
