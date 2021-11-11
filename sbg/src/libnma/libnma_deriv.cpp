/************************************************************************
 *                     LIBRARY: libnma_deriv                             *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Library to compute system derivatives (K-matrix <=> Jacobian) in    *
 *   Internal Coordinate Space (ICS).                                    *
 *   Memory efficient approach based on V/W-arrays implemented            *
 *   See Go's classic papers about the so called "Naive" method.         *
 *   (It takes into account Muliple-Chains and different CG-models)      *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <libnma_deriv.h>
#include <libnma_misc.h>
#include <libnma_matrix.h>
#define DEVEL

// Computes "v" and "w" (Multi-Chain) given Phi & Psi, etc...
void calcoefM( double mb1, double mb2, double phi[3], double psi[3], double Y1[3], double I1[3][3], double I2[3][3],
		double J[3][3], double v[2][3], double w[2][3] )
{
	//     K-matrix coeficients computation
	//     It computes the 4 vectors (from Go's paper) in terms of which the derivatives
	//     of the cartesian coordinates w/r/t the dihedral coordinates are expressed as:
	//     	Dr_b/Dq_j = v[0] + w[0] x r_b  (for body 1)
	//     	Dr_c/Dq_j = v[1] + w[1] x r_c  (for body 2) */

	int m, n;
	double mtot, ta[3], tb[3], tc[3];
	double e[3], y[3], Y2[3];

	mtot = mb1 + mb2;

	// V's
	// v (body 1)
	for ( m = 0; m < 3; m++ )
		v[0][m] = (mb2 / mtot) * psi[m]; // (1-M1/M) * Psi   (ec. 40 Braun et al.)
	// (M1/M) * (Phi x Y1)
	v[0][0] += (mb1/mtot) * ( phi[1] * Y1[2] - phi[2] * Y1[1] ); // vectorial product
	v[0][1] += (mb1/mtot) * ( phi[2] * Y1[0] - phi[0] * Y1[2] );
	v[0][2] += (mb1/mtot) * ( phi[0] * Y1[1] - phi[1] * Y1[0] );
	//  printf("v1 = %f,%f,%f\n",v[0][0],v[0][1],v[0][2]);

	for ( m = 0; m < 3; m++ )
		Y2[m] =  -Y1[m] * mb1 / mb2; // (Computes Y2[] from Y1[] and y[])????
	// *****************************
	// CoM(body 2) = Y2 = (CoM*mtot - Y1*mb1)/mb2 = -Y1*mb1/mb2
	// v2 right part (ec.25) = (mb2/mtot)*Y2 = -Y1*mb1/mtot

	// v (body 2)
	for ( m = 0; m < 3; m++ )
		v[1][m] = -(mb1 / mtot) * psi[m]; // -(1-M2/M) * Psi   (ec. 40 Braun et al.)
	// -(M2/M) * (Phi x Y2)
	v[1] [0] -= (mb2/mtot) * ( phi[1] * Y2[2] - phi[2] * Y2[1] ); // vectorial product
	v[1] [1] -= (mb2/mtot) * ( phi[2] * Y2[0] - phi[0] * Y2[2] );
	v[1] [2] -= (mb2/mtot) * ( phi[0] * Y2[1] - phi[1] * Y2[0] );

	// W's
	// tb = mb1 * (Y1 x Psi)
	tb[0] = mb1 * ( Y1[1] * psi[2] - Y1[2] * psi[1] ); // vectorial product
	tb[1] = mb1 * ( Y1[2] * psi[0] - Y1[0] * psi[2] );
	tb[2] = mb1 * ( Y1[0] * psi[1] - Y1[1] * psi[0] );

	// w (body 1)
	//> (3x3 * 3x1) matrix product to obtain "tc"
	//> vector addition to obtain "ta" from "tb" & "tc"
	for ( m = 0; m < 3; m++ ) // rows
	{
		tc[m] = 0.0;
		for ( n = 0; n < 3; n++ ) // cols
			tc[m] += I2[m] [n] * phi[n]; // I2 * Phi
		ta[m] = tb[m] + tc[m]; // mb1 * (Y1 x Psi) + I2 * Phi
	}

	for ( m = 0; m < 3; m++ ) // rows
	{
		w[0][m] = 0.0;
		for ( n = 0; n < 3; n++ ) // cols
			w[0][m] -= J[m] [n] * ta[n]; // "-="
	}

	// tb = mb2 * (Y2 x Psi)
	tb[0] = mb2 * ( Y2[1] * psi[2] - Y2[2] * psi[1] ); // vectorial product
	tb[1] = mb2 * ( Y2[2] * psi[0] - Y2[0] * psi[2] );
	tb[2] = mb2 * ( Y2[0] * psi[1] - Y2[1] * psi[0] );

	// w (body 2)
	//> (3x3 * 3x1) matrix product to obtain "tc"
	//> vector addition to obtain "ta" from "tb" & "tc"
	for ( m = 0; m < 3; m++ )
	{
		tc[m] = 0.0;
		for ( n = 0; n < 3; n++ )
			tc[m] += I1[m][n] * phi[n]; // I1 * Phi
		ta[m] = tb[m] + tc[m]; // mb2 * (Y2 x Psi) + I2 * Phi
	}

	for ( m = 0; m < 3; m++ )
	{
		w[1][m] = 0.0;
		for ( n = 0; n < 3; n++ )
			w[1][m] += J[m] [n] * ta[n]; // "+="
	}
}

// Computes derivatives for an interaction pair of atoms given "v"/"w" elements
double derVW(double **v, double **w,float *ri,float *rj)
{
	double der2[3],der3[3];
	double dummy;

	// i < j --> allways!
	// i-der --> = v + (w x r) (vectorial product)
	der3[0] = v[0] [0] + w[0] [1] * ri[2] - w[0] [2] * ri[1];
	der3[1] = v[0] [1] + w[0] [2] * ri[0] - w[0] [0] * ri[2];
	der3[2] = v[0] [2] + w[0] [0] * ri[1] - w[0] [1] * ri[0];
	// j-der --> = v + (w x r) (vectorial product)
	der2[0] = v[1] [0] + w[1] [1] * rj[2] - w[1] [2] * rj[1];
	der2[1] = v[1] [1] + w[1] [2] * rj[0] - w[1] [0] * rj[2];
	der2[2] = v[1] [2] + w[1] [0] * rj[1] - w[1] [1] * rj[0];

	dummy = (rj[0]-ri[0]) * (der2[0]-der3[0]) +
	(rj[1]-ri[1]) * (der2[1]-der3[1]) +
	(rj[2]-ri[2]) * (der2[2]-der3[2]);

	return dummy;
}

// Compute derivatives of cartesian coordinates w/r/t dihedrals (Multi-Chain + Fixation)
// NEEDS: -Macromolecule (first-NH + CA-atoms + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
//        -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void dydqMCAx(Macromolecule *mol,float *coord, trd **p_der, int size, bool *fix)
{
	bool debug = false;
	double mtot,mta;
	double r[3];
	double rd[3];
	double I[3][3],J[3][3];
	double temp1,temp2;
	int m, n, num_atoms, i, k, j, j2, k2;
	Residue *res;
	int resn,num_res,num_res2,num_seg;

	pdbIter *iter = new pdbIter( mol ); //Iterador para recorrer atomos
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_seg; // Iter to screen all segments
	pdbIter *iter_frag2; // Iter to screen all fragments
	pdbIter *iter_seg2 = new pdbIter(mol); // Iter to screen all segments

	num_seg = iter->num_segment();
	num_atoms = mol->get_num_atoms();

	// Computing the PDB's Center of Mass (CoM) --> Better outside!

	// Initializing Inertia-Matrix
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			I[i][j] = 0.0;

	// Compute matrix I (Inertia Momentum) and Total Mass
	// Ec. 22 (Noguti & Go, 1983)
	i = 0;
	mtot = 0.0;
	k=0; // residue index
	k2=0; // first-NH, CA and last-C atom index
	for ( iter_seg2->pos_segment = 0; !iter_seg2->gend_segment(); iter_seg2->next_segment() )
	{
		iter_frag2 = new pdbIter( iter_seg2->get_segment() );
		num_res2 = iter_frag2->num_fragment();
		for ( iter_frag2->pos_fragment = 0; !iter_frag2->gend_fragment(); iter_frag2->next_fragment() )
		{
			if(iter_frag2->pos_fragment == 0) // first segment fragment has NH
			{
				iter->pos_atom = k2; // NH index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k*3) * 3]; // NH position in coord
				r[1] = coord[(k*3) * 3 + 1];
				r[2] = coord[(k*3) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}
				k2++;
			}

			iter->pos_atom = k2; // CA index
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			mtot += mta; // total mass
			r[0] = coord[(k*3+1) * 3]; // NH position in coord
			r[1] = coord[(k*3+1) * 3 + 1];
			r[2] = coord[(k*3+1) * 3 + 2];
			temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
			for ( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
				{
					temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
					if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
					I[m][n] += mta * temp2;
				}
			k2++;

			if(iter_frag2->pos_fragment == num_res2-1) // last segment fragment has CO
			{
				iter->pos_atom = k2; // CO index
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mtot += mta; // total mass
				r[0] = coord[(k*3+2) * 3]; // CO position in coord
				r[1] = coord[(k*3+2) * 3 + 1];
				r[2] = coord[(k*3+2) * 3 + 2];
				temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
						if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
						I[m][n] += mta * temp2;
					}
				k2++;
			}

			k++;
		}
		delete iter_frag2;
	}
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	if(debug)
	{
		printf("Inertia-total (I)       Inverse-of-it (J)  mtot= %f\n",mtot);
		for(int x=0;x<3;x++)
			printf("%8.3f %8.3f %8.3f   %8.3f %8.3f %8.3f\n",
					I[x][0],I[x][1],I[x][2],J[x][0],J[x][1],J[x][2]);
	}

	// **********************************
	// ** DERIVATIVE COMPUTATION ********
	// **********************************
	double temp;
	double my[3], e[3], y[3], Y1[3], I1[3] [3], I2[3] [3], v[2] [3], w[2] [3]; // calcoef's
	double nr2, mfa;
	double mb1, mb2;
	int k1, b, k0;
	float pos[3];

	// K-Matrix Memory Allocation: der[] --> N x size -matrix !!!! (DERIVATIVES)
	trd *der; // K-matrix (derivatives)
	if(*p_der == NULL)
	{ // Allocate memory
		if( !( der = ( trd * ) malloc( num_atoms * size * sizeof( trd ) ) ) )
		{
			printf("Msg(dydqMCAx): K-Matrix memory allocation failed! (%d Mb)\n"
					"Forcing exit!\n",(int) (num_atoms * size * sizeof( trd ) / 1000000));
			exit(1);
		}
		else
		{
			printf("Msg(dydqMCAx): K-Matrix memory allocated: %d bytes (%d Mb)\n",num_atoms * size * sizeof( trd ),(int) (num_atoms * size * sizeof( trd ) / 1000000));
			*p_der = der; // outputs K-matrix
		}
	}
	else // Use already allocated memory
		der = *p_der;

	// Intialize variables
	// Only CA-atoms derivatives will be considered! (num_atoms = number of CA-atoms)
	for( i = 0; i < num_atoms * size; i++ )
	{
		der[i].x = 0.0;
		der[i].y = 0.0;
		der[i].z = 0.0;
	}

	// inertia matrix 1 initialization
	for ( m = 0; m < 3; m++ )
		for ( n = 0; n < 3; n++ )
			I1[m] [n] = 0.0;

	j = 0; // current dihedral index
	j2 = 0; // mobile dihedral index
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */
	k1 = 0; // first-NH + CA + last-CO model index
	// N,CA,C model index
	k0 = 0; // first-NH + CA + last-CO model index
	mb1 = 0.0;

	// ********************************************************************************************
	// ************* DERIVATIVES LOOP BEGIN *******************************************************
	// ********************************************************************************************
	// printf("Msg(dydqM): Computing K-matrix for %d dihedrals and %d pseudo-atoms\n",size,num_atoms);

	Segment * seg;
	Atom * atom;
	// Screening segments
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// INTER-SEGMENT DEGREES OF FREEDOM (3 Translations + 3 Rotations)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				// mb1, mb2 and my, already updated before (between PHI and PSI)
				if(debug)
					printf("Segment %d:  mb1= %f   mb2= %f\n",iter_seg->pos_segment,mb1,mb2);

				// CoM 1
				for( m = 0; m < 3; m++ )
					Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM
				// now matrix I2
				for( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I2[m][n] = I[m][n] - I1[m][n];
			}

			// 3 TRANSLATIONS
			//  printf ("Translation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // Phi
						y[m] = 0.0; // Psi
					}
					y[axis] = -1.0; // Psi = -gamma_v
					// gamma_v --> Column unit vectors along positive x, y, z axes of S-system with respect to F-system

					// Computing some coefficients (Multi-Chain)
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					k=0; // Residue index
					k2=0; // NH and CA atom index
					for ( iter_seg2->pos_segment = 0; !iter_seg2->gend_segment(); iter_seg2->next_segment() )
					{
						iter_frag2 = new pdbIter( iter_seg2->get_segment() );
						num_res2 = iter_frag2->num_fragment();

						for ( iter_frag2->pos_fragment = 0; !iter_frag2->gend_fragment(); iter_frag2->next_fragment() )
						{
							// For NH's...
							if(iter_frag2->pos_fragment==0) // first segment fragment has NH
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3) * 3]; // CA position
								r[1] = coord[(k*3) * 3 + 1];
								r[2] = coord[(k*3) * 3 + 2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}

							// For CA's...
							if( k2 < k0 )
								b = 0;    // body 1
							else
								b = 1;	  // body 2
							r[0] = coord[(k*3+1) * 3]; // CA position
							r[1] = coord[(k*3+1) * 3 + 1];
							r[2] = coord[(k*3+1) * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
							k2++;

							// For CO's...
							if(iter_frag2->pos_fragment == num_res2-1) // last segment fragment has CO
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3+2)*3];  // CO position
								r[1] = coord[(k*3+2)*3+1];
								r[2] = coord[(k*3+2)*3+2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}

							k++;
						}
						delete iter_frag2;
					}
					j2++;
				}
				j++; // der[] index
			}

			// 3 ROTATIONS
			//  printf ("Rotation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0;
						r[m] = 0.0;
					}
					e[axis] = 1.0; // Phi

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// Computing Derivatives ( dy(b)/dTheta(lamda) ) (K-Matrix)
					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					k=0; // Residue index
					k2=0; // NH and CA atom index
					for ( iter_seg2->pos_segment = 0; !iter_seg2->gend_segment(); iter_seg2->next_segment() )
					{
						iter_frag2 = new pdbIter( iter_seg2->get_segment() );
						num_res2 = iter_frag2->num_fragment();

						for ( iter_frag2->pos_fragment = 0; !iter_frag2->gend_fragment(); iter_frag2->next_fragment() )
						{
							// For NH's...
							if(iter_frag2->pos_fragment==0) // first segment fragment has NH
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3) * 3]; // CA position
								r[1] = coord[(k*3) * 3 + 1];
								r[2] = coord[(k*3) * 3 + 2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}

							// For CA's...
							if( k2 < k0 )
								b = 0;    // body 1
							else
								b = 1;	  // body 2
							r[0] = coord[(k*3+1) * 3]; // CA position
							r[1] = coord[(k*3+1) * 3 + 1];
							r[2] = coord[(k*3+1) * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
							k2++;

							// For CO's...
							if(iter_frag2->pos_fragment == num_res2-1) // last segment fragment has CO
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3+2)*3];  // CO position
								r[1] = coord[(k*3+2)*3+1];
								r[2] = coord[(k*3+2)*3+2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}
							k++;
						}
						delete iter_frag2;
					}
					j2++;
				}
				j++; // der[] index
			}
		}

		// Screen ALL residues
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			// NH terms
			if(iter_frag->pos_fragment==0) // first segment fragment has NH
			{
				iter->pos_atom = k0; // NH
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mb1 += mta;
				mb2 = mtot - mb1;
				r[0] = coord[k1*3 * 3]; // NH atom
				r[1] = coord[k1*3 * 3 + 1];
				r[2] = coord[k1*3 * 3 + 2];
				my[0] += mta * r[0];
				my[1] += mta * r[1];
				my[2] += mta * r[2];

				// Cummulative Inertia Matrix 1
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
				if(debug)
				{
					printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mb2,my[0],my[1],my[2]);
					printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
				}
				k0++; // NH, CA and CO index (now points to the next CA)
			}

			// ********************************************************************************************
			// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
			if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
			{
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get CA pos 1
					e[0] = pos[0] = coord[(k1*3+1) * 3]; // CA position
					e[1] = pos[1] = coord[(k1*3+1) * 3 + 1];
					e[2] = pos[2] = coord[(k1*3+1) * 3 + 2];

					// get NH pos 0
					y[0] = r[0] = coord[k1*3 * 3];
					y[1] = r[1] = coord[k1*3 * 3 + 1];
					y[2] = r[2] = coord[k1*3 * 3 + 2];
					e[0] -= y[0]; // NH --> CA
					e[1] -= y[1];
					e[2] -= y[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM
					}

					// Inertia Matrix 1 (already updated)
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					if(debug)
					{
						printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I2[0][0],I2[0][1],I2[0][2],I2[1][0],I2[1][1],I2[1][2],I2[2][0],I2[2][1],I2[2][2]);
						printf("PHI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi = (e x r)
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					k=0; // CA atom index
					k2=0; // NH and CA atom index
					for ( iter_seg2->pos_segment = 0; !iter_seg2->gend_segment(); iter_seg2->next_segment() )
					{
						iter_frag2 = new pdbIter( iter_seg2->get_segment() );
						num_res2 = iter_frag2->num_fragment();

						for ( iter_frag2->pos_fragment = 0; !iter_frag2->gend_fragment(); iter_frag2->next_fragment() )
						{
							// For NH's...
							if(iter_frag2->pos_fragment==0) // first segment fragment has NH
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3) * 3]; // NH position
								r[1] = coord[(k*3) * 3 + 1];
								r[2] = coord[(k*3) * 3 + 2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}

							// For CA's...
							if( k2 < k0 )
								b = 0;    // body 1
							else
								b = 1;	  // body 2
							r[0] = coord[(k*3+1) * 3]; // CA position
							r[1] = coord[(k*3+1) * 3 + 1];
							r[2] = coord[(k*3+1) * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
							k2++;

							// For CO's...
							if(iter_frag2->pos_fragment == num_res2-1) // last segment fragment has CO
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3+2)*3];  // CO position
								r[1] = coord[(k*3+2)*3+1];
								r[2] = coord[(k*3+2)*3+2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}
							k++;
						}
						delete iter_frag2;
					}
					j2++;
				}
				j++; // der[] index
			}  // NOT FIRST NOT PRO

			// CA terms
			iter->pos_atom = k0; // CA
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			mb1 += mta;
			mb2 = mtot - mb1;
			r[0] = coord[(k1*3+1)*3];
			r[1] = coord[(k1*3+1)*3+1];
			r[2] = coord[(k1*3+1)*3+2];
			my[0] += mta * r[0];
			my[1] += mta * r[1];
			my[2] += mta * r[2];

			// Cummulative Inertia Matrix 1
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
			if(debug)
			{
				printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mb2,my[0],my[1],my[2]);
				printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
			}
			k0++; // NH, CA and CO index (now points to the next CA or CO)

			// ********************************************************************************************
			// "PSI-bond" (In CA-only model, non-first and non-last)
			if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
			{
				if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
				{
					// Gets e_lambda (CA-->) Psi-bond
					// get C pos 2
					e[0] = coord[(k1*3+2) * 3];
					e[1] = coord[(k1*3+2) * 3 + 1];
					e[2] = coord[(k1*3+2) * 3 + 2];

					// get CA pos 1
					y[0] = pos[0] = coord[(k1*3+1) * 3]; // CA position
					y[1] = pos[1] = coord[(k1*3+1) * 3 + 1];
					y[2] = pos[2] = coord[(k1*3+1) * 3 + 2];

					e[0] -= y[0]; // CA --> C
					e[1] -= y[1];
					e[2] -= y[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]/mb1; // Updating Body 1 Center of Mass: "Y1[m]"
					}

					// Inertia Matrix 1 (already updated above)
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
							I2[m][n] = I[m][n] - I1[m][n];

					if(debug)
					{
						printf("I2=  %f %f %f  %f %f %f  %f %f %f\n",I2[0][0],I2[0][1],I2[0][2],I2[1][0],I2[1][1],I2[1][2],I2[2][0],I2[2][1],I2[2][2]);
						printf("PSI j2=%d y= %f,%f,%f  e= %f,%f,%f  Y1= %f,%f,%f\n"
								,j2,y[0],y[1],y[2],e[0],e[1],e[2],Y1[0],Y1[1],Y1[2]);
					}

					y[0] = e[1] * pos[2] - e[2] * pos[1]; // Psi (e x r)
					y[1] = e[2] * pos[0] - e[0] * pos[2];
					y[2] = e[0] * pos[1] - e[1] * pos[0];

					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					k=0; // residue index
					k2=0; // NH and CA atom index
					for ( iter_seg2->pos_segment = 0; !iter_seg2->gend_segment(); iter_seg2->next_segment() )
					{
						iter_frag2 = new pdbIter( iter_seg2->get_segment() );
						num_res2 = iter_frag2->num_fragment();

						for ( iter_frag2->pos_fragment = 0; !iter_frag2->gend_fragment(); iter_frag2->next_fragment() )
						{
							// For NH's...
							if(iter_frag2->pos_fragment==0) // first segment fragment has NH
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3) * 3]; // CA position
								r[1] = coord[(k*3) * 3 + 1];
								r[2] = coord[(k*3) * 3 + 2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}

							// For CA's...
							if( k2 < k0 )
								b = 0;    // body 1
							else
								b = 1;	  // body 2
							r[0] = coord[(k*3+1) * 3]; // CA position
							r[1] = coord[(k*3+1) * 3 + 1];
							r[2] = coord[(k*3+1) * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
							k2++;

							// For CO's...
							if(iter_frag2->pos_fragment == num_res2-1) // last segment fragment has CO
							{
								if( k2 < k0 )
									b = 0;    // body 1
								else
									b = 1;	  // body 2
								r[0] = coord[(k*3+2)*3];  // CO position
								r[1] = coord[(k*3+2)*3+1];
								r[2] = coord[(k*3+2)*3+2];
								// = v + (w x r) (vectorial product)
								der[k2 * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k2 * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k2 * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
								k2++;
							}
							k++;
						}
						delete iter_frag2;
					}
					j2++;
				}
				j++; // der[] index
			}  // NOT LAST & NOT FIRST

			if(iter_frag->pos_fragment == num_res-1) // last segment fragment has CO
			{
				// CO terms
				iter->pos_atom = k0; // CO
				mta = ( iter->get_atom() )->getPdbocc(); // mass
				mb1 += mta;
				mb2 = mtot - mb1;
				r[0] = coord[(k1*3+2)*3];  // CO position
				r[1] = coord[(k1*3+2)*3+1];
				r[2] = coord[(k1*3+2)*3+2];
				my[0] += mta * r[0];
				my[1] += mta * r[1];
				my[2] += mta * r[2];

				// Cummulative Inertia Matrix 1
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
				if(debug)
				{
					printf("k1= %d  k0= %d  mb1= %f  mb2= %f  my= %f %f %f\n",k1,k0,mb1,mb2,my[0],my[1],my[2]);
					printf("I1=  %f %f %f  %f %f %f  %f %f %f\n",I1[0][0],I1[0][1],I1[0][2],I1[1][0],I1[1][1],I1[1][2],I1[2][0],I1[2][1],I1[2][2]);
				}
				k0++; // NH, CA and CO index (now points to the next NH)
			}

			k1++; // residue (or CA) index
		}
		delete iter_frag;
	}
	// end derivatives
	delete iter;
	delete iter_seg;
	delete iter_seg2;
}

// Compute derivatives of cartesian coordinates w/r/t dihedrals (Multi-Chain)
// NEEDS: -Macromolecule (N,CA,CA model) (3-atoms model),
//        -Single row (N,CA,C model) coordinates,
//        -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void dydqMCA3x(Macromolecule *mol,float *coord, trd **p_der, int size, bool *fix)
{
	bool debug = false;
	double mtot,mta;
	float r[3];
	double rd[3];
	double I[3][3],J[3][3];
	double temp1,temp2;
	int m, n, num_atoms, i, k, j, j2;
	Residue *res;
	int resn,num_res,num_seg;

	num_atoms = mol->get_num_atoms(); // Note this will be "3*num_res" (CA-model)
	num_res = mol->get_num_fragments();
	if(debug)
	{
		printf("Msg(deriv): Number of pseudo-atoms = %d\n",num_atoms);
		printf("Msg(deriv): Number of fragments = %d\n",num_res);
	}

	pdbIter *iter = new pdbIter( mol ); //Iterador para recorrer atomos
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_seg; // Iter to screen all segments
	num_seg = iter->num_segment();

	// Initializing Inertia-Matrix
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			I[i][j] = 0.0;

	// Compute matrix I (Inertia Momentum) and Total Mass
	// Ec. 22 (Noguti & Go, 1983)
	i = 0;
	mtot = 0.0;
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screens all-CA
	{
		mta = ( iter->get_atom() )->getPdbocc(); // mass
		if(mta != 0.0)
		{
			mtot += mta; // total mass
			r[0] = coord[iter->pos_atom * 3]; // atom position in coord
			r[1] = coord[iter->pos_atom * 3 + 1];
			r[2] = coord[iter->pos_atom * 3 + 2];
			temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
			for ( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
				{
					temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
					if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
					I[m][n] += mta * temp2;
				}
		}
	}
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	if(debug)
		for(int x=0;x<3;x++)
			printf("%8.3f %8.3f %8.3f   %10.6f %10.6f %10.6f\n",
					I[x][0],I[x][1],I[x][2],J[x][0],J[x][1],J[x][2]);

	// **********************************
	// ** DERIVATIVE COMPUTATION ********
	// **********************************
	double temp;
	double my[3], e[3], y[3], Y1[3], I1[3] [3], I2[3] [3], v[2] [3], w[2] [3]; // calcoef's
	double nr2, mfa;
	float mb1, mb2;
	int k1, b;
	float pos[3];

	// K-Matrix Memory Allocation: der[] --> N x size -matrix !!!! (DERIVATIVES)
	trd *der; // K-matrix (derivatives)
	if(*p_der == NULL)
	{ // Allocate memory
		if( !( der = ( trd * ) malloc( num_atoms * size * sizeof( trd ) ) ) )
		{
			printf("Msg(dydqMCA3x): K-Matrix memory allocation failed! (%.3f Mb)\n"
					"Forcing exit!\n",(float) (num_atoms * size * sizeof( trd ) / 1e6));
			exit(1);
		}
		else
		{
			printf("Msg(dydqMCA3x): K-Matrix memory allocated: %ld bytes (%ld Mb)\n",num_atoms * size * sizeof( trd ),(num_atoms * size * sizeof( trd ) / 1000000));
			*p_der = der; // outputs K-matrix
		}
	}
	else // Use already allocated memory
		der = *p_der;

	// Intialize variables
	for( i = 0; i < num_atoms * size; i++ )
	{
		der[i].x = 0.0;
		der[i].y = 0.0;
		der[i].z = 0.0;
	}

	// inertia matrix 1 initialization
	for ( m = 0; m < 3; m++ )
		for ( n = 0; n < 3; n++ )
			I1[m] [n] = 0.0;

	j = 0; // current dihedral index
	j2 = 0; // current dihedral index
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */
	k1 = 0;
	mb1 = 0.0;

	// ********************************************************************************************
	// ************* DERIVATIVES LOOP BEGIN *******************************************************
	// ********************************************************************************************

	Segment * seg;
	Atom * atom;
	// Screening segments
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments

	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		if(debug)
			printf("\nProcessing segment %d (%d): k1 %d\n",iter_seg->pos_segment,num_res,k1);

		// INTER-SEGMENT DEGREES OF FREEDOM (3 Translations + 3 Rotations)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				// CoM 1
				for ( m = 0; m < 3; m++ )
					Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM

				// now matrix I2
				for( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I2[m][n] = I[m][n] - I1[m][n];
			}

			// 3 TRANSLATIONS
			//  printf ("Translation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // Phi
						y[m] = 0.0; // Psi
					}
					y[axis] = -1.0; // Psi = -gamma_v
					// gamma_v --> Column unit vectors along positive x, y, z axes of S-system with respect to F-system

					// Computing some coefficients (Multi-Chain)
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k < k1 )
							b = 0;    // body 1
						else
							b = 1;	   // body 2

						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						// k --> atoms
						// j --> dihedrals
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++;
				}
				j++; // der[] index
			}

			// 3 ROTATIONS
			//  printf ("Rotation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{

					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0;
						r[m] = 0.0;
					}
					e[axis] = 1.0; // Phi

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// Computing Derivatives ( dy(b)/dTheta(lamda) ) (K-Matrix)
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k < k1 )
							b = 0;    // body 1
						else
							b = 1;	   // body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++;
				}
				j++; // der[] index
			}
		}

		// Screen ALL residues
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			// NH terms
			iter->pos_atom = k1; // NH
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if(mta != 0.0)
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[k1 * 3]; // NH atom
				my[1] += mta * coord[k1 * 3+1];
				my[2] += mta * coord[k1 * 3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // N,CA,C positions
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

			// ********************************************************************************************
			// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (Ca-model)
			if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
			{

				if(debug)
					printf("res= %d  IC= %d  PHI  k1= %d\n",iter_frag->pos_fragment,j,k1);

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get CA pos 1
					e[0] = pos[0] = coord[(k1+1) * 3]; // CA position
					e[1] = pos[1] = coord[(k1+1) * 3 + 1];
					e[2] = pos[2] = coord[(k1+1) * 3 + 2];
					if(debug)
						printf("PHI %d CA-pos= %f %f %f\n",j,e[0],e[1],e[2]);
					// get NH pos 0
					y[0] = r[0] = coord[k1 * 3];
					y[1] = r[1] = coord[k1 * 3 + 1];
					y[2] = r[2] = coord[k1 * 3 + 2];
					e[0] -= y[0]; // NH --> CA
					e[1] -= y[1];
					e[2] -= y[2];
					if(debug)
						printf("PHI %d NH-pos= %f %f %f\n",j,r[0],r[1],r[2]);

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM
					}
					if(debug)
						printf("j= %d  j2= %d  (%3s) PHI M1= %f  M2= %f  after Y1= %f %f %f\n",j,j2,res->getName(),mb1,mb2,Y1[0],Y1[1],Y1[2]);

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi = (e x r)
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)
					if(debug)
					{
						printf("calcoefM[res=%d dh=%d] b=0  v= %8.5f %8.5f %8.5f  w= %8.5f %8.5f %8.5f\n"
								,iter_frag->pos_fragment,j,v[0][0],v[0][1],v[0][2],w[0][0],w[0][1],w[0][2]);
						printf("calcoefM[res=%d dh=%d] b=1  v= %8.5f %8.5f %8.5f  w= %8.5f %8.5f %8.5f\n"
								,iter_frag->pos_fragment,j,v[1][0],v[1][1],v[1][2],w[1][0],w[1][1],w[1][2]);
					}

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for(int k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1 )
							b = 0;    // body 1 (NH included)
						else
							b = 1;	  // body 2
						r[0] = coord[k * 3]; // N,CA,C positions
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];

						if(debug)
							printf("der[at=%d dh=%d]= %8.5f %8.5f %8.5f\n",k,j,der[k*size+j].x,der[k*size+j].y,der[k*size+j].z);
					}
					j2++;
				}
				j++; // der[] index
			}  // NOT FIRST NOT PRO

			// CA terms
			iter->pos_atom = k1+1; // residue CA (contains the whole residue mass)
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if( mta != 0.0 )
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[(k1+1)*3]; // CA atom
				my[1] += mta * coord[(k1+1)*3+1];
				my[2] += mta * coord[(k1+1)*3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // atom positions
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

			// ********************************************************************************************
			// "PSI-bond" (In CA-only model, non-first and non-last)
			if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
			{
				if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
				{
					// Gets e_lambda (CA-->) Psi-bond
					e[0] = coord[(k1+2) * 3]; // C position
					e[1] = coord[(k1+2) * 3 + 1];
					e[2] = coord[(k1+2) * 3 + 2];

					// get CA pos 1
					y[0] = pos[0] = coord[(k1+1) * 3]; // CA position
					y[1] = pos[1] = coord[(k1+1) * 3 + 1];
					y[2] = pos[2] = coord[(k1+1) * 3 + 2];

					e[0] -= y[0]; // CA --> C
					e[1] -= y[1];
					e[2] -= y[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]/mb1; // Updating Body 1 Center of Mass: "Y1[m]"
					}

					// now matrix I2
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
							I2[m][n] = I[m][n] - I1[m][n];

					y[0] = e[1] * pos[2] - e[2] * pos[1]; // Psi (e x r)
					y[1] = e[2] * pos[0] - e[0] * pos[2];
					y[2] = e[0] * pos[1] - e[1] * pos[0];

					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// Screens ALL atoms
					for(int k = 0; k < num_atoms; k++ )
					{
						if( k <= k1+1 ) // k1+1 --> CA index
							b = 0; // body 1
						else
							b = 1;	// body 2
						r[0] = coord[k*3]; // N,CA,C positions
						r[1] = coord[k*3 + 1];
						r[2] = coord[k*3 + 2];
						// vectorial product
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++;
				}
				j++; // der[] index
			}  // NOT LAST & NOT FIRST

			// CO terms
			iter->pos_atom = k1+2; // residue CA (contains the whole residue mass)
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if( mta != 0.0 )
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[(k1+2)*3]; // CA atom
				my[1] += mta * coord[(k1+2)*3+1];
				my[2] += mta * coord[(k1+2)*3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // atom positions
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

			k1 += 3; // current atom index (first atom of the residue)
		}
		delete iter_frag;
	}// end derivatives
	// *********************************************************
	// ************* DERIVATIVES LOOP END **********************
	// *********************************************************
	iter_seg->clean_virtual();
	delete iter_seg;
	delete iter;
}

// Compute derivatives of cartesian coordinates w/r/t dihedrals
// (Multi-Chain & FULL-ATOM & Protein/RNA/DNA/SMOL)
// NEEDS: Macromolecule, Single row PDB coordinates,
// -Properties structure array (props[]), -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: type = 0 --> phi,psi,... type = 1 & 2 --> phi,chi,psi...
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
void dydqMFAx(Macromolecule *molr,float *coord, tri *props, trd **p_der, int type, int model, int size, bool *fix, bool *addrot)
{
	bool debug=false;
	bool debug2=false;
	bool debug3=false;
	double mtot,mta;
	double r[3];
	double rd[3];
	double I[3][3],J[3][3];
	double temp1,temp2;
	int m, n, num_atoms, i, k, j;
	Residue *res;
	int resn,num_res,num_seg,res_index=0;
	int CBindex;
	bool has_oxt = false;

	pdbIter *iter = new pdbIter( molr ); //Iterador para recorrer atomos
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_seg; // Iter to screen all segments
	num_atoms = molr->num_atoms();

	// Initializing Inertia-Matrix
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			I[i][j] = 0.0;

	// Compute matrix I (Inertia Momentum) and Total Mass
	// Ec. 22 (Noguti & Go, 1983)
	i = 0;
	mtot = 0.0;
	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screens all-atoms
	{
		mta = ( iter->get_atom() )->getPdbocc(); // mass
		mtot += mta; // total mass
		r[0] = coord[i * 3];
		r[1] = coord[i * 3 + 1];
		r[2] = coord[i * 3 + 2];
		i++;
		temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
		for ( m = 0; m < 3; m++ )
			for ( n = 0; n < 3; n++ )
			{
				temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
				if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
				I[m][n] += mta * temp2;
			}
	}
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	// ************************************************
	// ************ DERIVATIVE COMPUTATION ************
	// ************************************************
	double mpr,temp, my[3], e[3], y[3], Y1[3], I1[3] [3], I2[3] [3], v[2] [3], w[2] [3];
	double nr2, mfa;
	double mb1, mb2;
	int k1, k2, b;
	float pos[3];
	TMOL fragtype;

	// K-Matrix Memory Allocation: der[] --> N x size -matrix !!!! (DERIVATIVES)
	trd *der; // K-matrix (derivatives)
	if(*p_der == NULL)
	{ // Allocate memory
		if( !( der = ( trd * ) malloc( num_atoms * size * sizeof( trd ) ) ) )
		{
			printf("Msg(dydqM): K-Matrix memory allocation failed! (%d Mb)\nForcing exit!\n",(int) (num_atoms * size * sizeof( trd ) / 1000000));
			exit(1);
		}
		else
		{
			if(debug)
				printf("MEM-USAGE: atoms=%d by size=%d K-Matrix size (Mb) = %.3f\n",num_atoms,size,num_atoms * size * sizeof( trd )/1e6);
			*p_der = der; // outputs K-matrix
		}
	}
	else // Use already allocated memory
		der = *p_der;

	// Intialize variables
	for ( i = 0; i < num_atoms * size; i++ )
	{
		der[i].x = 0.0;
		der[i].y = 0.0;
		der[i].z = 0.0;
	}

	for ( m = 0; m < 3; m++ )
		for ( n = 0; n < 3; n++ )
			I1[m] [n] = 0.0;

	j = 0; // current dihedral index (screens all dihedrals)
	int j2 = 0; // active dihedral index (screens all mobile dihedrals)
	mpr=0; /* mpr = mass of previous residues */
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */
	k1 = 0;
	k2 = props[0].nat - 1; /* 1st and last indices of atoms of residue i */ // Residue 0 ?
	int indexbase;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

	// ********************************************************************************************
	// ************* DERIVATIVES LOOP BEGIN *******************************************************
	// ********************************************************************************************

	iter_seg = new pdbIter( molr, true, true, true, true ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	Segment *seg;
	Atom *atom;

	j = 0; // must be reseted
	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// M1 and M2
			// (mb1 already updated)
			mb2 = mtot - mb1; // Right-half = Tot - Left_half

			// CoM 1
			for( m = 0; m < 3; m++ )
				Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM

			// Inertia Matrix 1
			for(m=0;m<3;m++)
				for(n=0;n<3;n++)
					I1[m][n] = 0.0;
			for( iter->pos_atom = 0; iter->pos_atom < k1; iter->next_atom() )
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
			// now matrix I2
			for( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
					I2[m][n] = I[m][n] - I1[m][n];
			// End Inertia matrix computation

			// 3 TRANSLATIONS
			//  printf ("Translation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // Phi = 0 (translation)
						y[m] = 0.0; // Psi
					}
					y[axis] = -1.0; // Psi = -gamma_v
					// gamma_v --> Column unit vectors along positive x, y, z axes of S-system with respect to F-system

					if(debug3)
						printf("TRA  seg= %d  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",iter_seg->pos_segment,k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

					// Computing some coefficients (Multi-Chain)
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k < k1 )
							b = 0;	// body 1
						else
							b = 1;	// body 2

						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						// k --> atoms
						// j --> dihedrals
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++;
				}
				j++; // der[] index
			}

			// 3 ROTATIONS
			//  printf ("Rotation \n");
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{
						// e_lambda
						for ( m = 0; m < 3; m++ )
						{
							e[m] = 0.0; // Phi
							r[m] = 0.0; // Psi (in rotation is the origin: 0,0,0)
						}
						e[axis] = 1.0; // Phi

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi = 0 ????
						y[1] = e[2] * r[0] - e[0] * r[2]; // Put Psi=0 directly whenever!
						y[2] = e[0] * r[1] - e[1] * r[0];

						if(debug3)
							printf("ROT  seg= %d  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",iter_seg->pos_segment,k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// Computing Derivatives ( dy(b)/dTheta(lamda) ) (K-Matrix)
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k < k1 )
								b = 0;	// body 1
							else
								b = 1;	// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++;
					}
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // der[] index
				}
		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// Take into account: N-atom
				//
				// Updating M1 and M2 masses...
				iter->pos_atom = k1; // N-atom
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// updating previous masses*positions
				y[0] = r[0] = coord[k1 * 3];
				y[1] = r[1] = coord[k1 * 3 + 1];
				y[2] = r[2] = coord[k1 * 3 + 2];
				for ( m = 0; m < 3; m++ )
					my[m] +=  mta * y[m];

				// PHI --> NOT FIRST, NOT PRO
				if ( iter_frag->pos_fragment != 0 && resn != PRO )
				{
					//  printf (" Not first not PRO -> PHI \n");
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get CA pos 1
						e[0] = coord[(k1+1) * 3];
						e[1] = coord[(k1+1) * 3 + 1];
						e[2] = coord[(k1+1) * 3 + 2];

						// get NH pos 0
						e[0] -= y[0]; // NH --> CA
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
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
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						if(debug3)
							printf("PHI  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// K-matrix --> der[]
						// Computing Derivatives ( dy(b)/dTheta(lamda) )
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k <= k1 )
								b = 0;    // body 1 (NH of residue included)
							else
								b = 1;				// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index
				}  // NOT FIRST NOT PRO

				if(model==1)
					// 3BB2R
					CBindex = 3;
				else
				{	// Full-Atom
					CBindex = 4;

					// Checking if last residue has OXT (Ot)
					if(iter_frag->pos_fragment == num_res-1)
					{
						iter->pos_atom = k2;
						if( strcmp(iter->get_atom()->getName()," OXT ")==0 )
						{
							k2--; // this avoids including OXT in CHI's body 2.
							has_oxt = true;
						}
					}
				}

				//  LATERAL CHAIN-->CHI
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if(type == 2) // phi,chi,psi,...
					if( props[res_index].nan==3 ||
							(props[res_index].nan==2 &&
									(iter_frag->pos_fragment==0 ||
											(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
						{
							// M1 and M2
							mb2 = 0.0;
							// screening CB & -R atoms (CB: "k1+3")
							for( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
								mb2 += ( iter->get_atom() )->getPdbocc(); // atom mass
							// Now CB & -R is the whole 2nd body, the remaining atoms will be the 1st one!
							mpr = mb1; // body 1 mass buffer
							mb1 = mtot - mb2;

							// e_lambda
							// get CB pos 4
							e[0] = pos[0] = coord[(k1+CBindex) * 3 + 0];
							e[1] = pos[1] = coord[(k1+CBindex) * 3 + 1];
							e[2] = pos[2] = coord[(k1+CBindex) * 3 + 2];
							// get CA pos 1
							r[0] = y[0] = pos[0] = coord[(k1+1) * 3 + 0];
							r[1] = y[1] = pos[1] = coord[(k1+1) * 3 + 1];
							r[2] = y[2] = pos[2] = coord[(k1+1) * 3 + 2];
							e[0] -= y[0]; // CA-->CB unit vector
							e[1] -= y[1];
							e[2] -= y[2];
							temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
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

							y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
							y[1] = e[2] * r[0] - e[0] * r[2];
							y[2] = e[0] * r[1] - e[1] * r[0];

							// this should be expressed different!
							// screening [k1+4,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc(); // atom mass
								// Substracting to the CoM (total) the body 2 (CB-R) contribution
								Y1[0] -= mta * coord[iter->pos_atom*3];
								Y1[1] -= mta * coord[iter->pos_atom*3+1];
								Y1[2] -= mta * coord[iter->pos_atom*3+2];
							}
							for(m=0;m<3;m++)
								Y1[m] /= mb1;

							// Calculate matrix I2
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I2[m][n] = 0;

							// screening [k1+4,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
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
										I2[m] [n] += mta * mfa;
									}
								}
							}
							// Now matrix I1
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I1[m][n] = I[m][n]-I2[m][n];

							calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

							for( k = 0; k < num_atoms; k++ )
							{
								if ( k <k1 || k>k2 || k<k1+CBindex )
									b = 0; // body 1 + res but CB-R
								else
									b = 1; 							 // body 2 next residues // ( CB & -R )
								r[0] = coord[k * 3];
								r[1] = coord[k * 3 + 1];
								r[2] = coord[k * 3 + 2];
								// = v + (w x r) (vectorial product)
								der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
								der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
								der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
							}
							j2++; // selected dihedral index
							mb1 = mpr; // restoring mb1's value!
						}
						j++; // der[] index
					}

				// Take into account: Not C and O atoms
				for ( iter->pos_atom = k1+1; iter->pos_atom <= k2; iter->next_atom() ) // Not N-atom
					if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
					{
						// Masses
						mta = iter->get_atom()->getPdbocc(); // Non-(C/O) atoms mass
						mb1 += mta;
						// Masses * positions
						my[0] += mta * coord[iter->pos_atom*3];
						my[1] += mta * coord[iter->pos_atom*3+1];
						my[2] += mta * coord[iter->pos_atom*3+2];
					}
				// Current body 2 mass (M2)
				mb2 = mtot - mb1;

				// NOT LAST RESIDUE--> PSI
				// "PSI-bond" (in Full-Atom, every residue has PSI)
				if ( iter_frag->pos_fragment != num_res - 1 || model==2 ) // Full-Atom allways has PSI
				{
					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// get C pos +2
						e[0] = coord[(k1+2) * 3 + 0];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// get CA pos +1
						y[0] = pos[0] = coord[(k1+1) * 3 + 0];
						y[1] = pos[1] = coord[(k1+1) * 3 + 1];
						y[2] = pos[2] = coord[(k1+1) * 3 + 2];
						e[0] -= y[0]; // CA --> C
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// Updating Body 1 Center of Mass: "Y1[m]"
						for(m=0;m<3;m++) // 3D coords
							Y1[m] = my[m]/mb1;

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
						// now matrix I2
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						if(debug3)
							printf("PSI  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

						y[0] = e[1] * pos[2] - e[2] * pos[1]; // Psi
						y[1] = e[2] * pos[0] - e[0] * pos[2];
						y[2] = e[0] * pos[1] - e[1] * pos[0];

						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// Screens ALL atoms
						for ( k = 0; k < num_atoms; k++ )
						{
							// this could be expresed clearly!!!
							if (  k <= k2 && !(k == k1+2 || (model == 2 && k == k1+3 ) ) ) // note that OXT is not in body 1
								b = 0; // body 1 + all except C/O of current residue
							else
								b = 1;	// body 2 next residues (C/O-atoms belong to 2nd body)
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// vectorial product
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index
				}

				// Updating "mb1" and "my" for the next residue!
				// C/O-atoms
				iter->pos_atom = k1+2;
				mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
				mb1 += mta;
				my[0] += mta * coord[iter->pos_atom*3];
				my[1] += mta * coord[iter->pos_atom*3+1];
				my[2] += mta * coord[iter->pos_atom*3+2];

				if(model == 2) // then O-atom must be taken into account!
				{			   // (if model==1 --> already accounted for...
					iter->pos_atom = k1+3;
					mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
					mb1 += mta;
					my[0] += mta * coord[iter->pos_atom*3];
					my[1] += mta * coord[iter->pos_atom*3+1];
					my[2] += mta * coord[iter->pos_atom*3+2];

					// Taking into account whether last residue has OXT (Ot) (Full-atom only)
					if(iter_frag->pos_fragment == num_res-1)
					{
						if(has_oxt)
						{
							k2++;
							iter->pos_atom = k2; // now it's "k2+1" (above: k2--)
							mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
							mb1 += mta;
							my[0] += mta * coord[iter->pos_atom*3];
							my[1] += mta * coord[iter->pos_atom*3+1];
							my[2] += mta * coord[iter->pos_atom*3+2];
						}
						has_oxt = false; // initialice for the next segment
					}
				}
			}
			// if RNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)
				// ALPHA (bond between P and O5*)
				for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
				{	// M1 += P + O1P + O2P
					// M1 and M2
					mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
					mb1 += mta; // Left_half = previous + actual (NH)

					// Updating previous masses*positions
					for ( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
				}
				mb2 = mtot - mb1; // Right-half = Tot - Left_half

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = r[0] = coord[k1*3];
					y[1] = r[1] = coord[k1*3 + 1];
					y[2] = r[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+2 )
							b = 0;	// body 1 (O2P included)
						else
							b = 1;	// body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				// BETA (bond between O5* and C5*)
				iter->pos_atom = k1+3; // O5* index
				// M1 += O5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = r[0] = coord[(k1+3)*3];
					y[1] = r[1] = coord[(k1+3)*3 + 1];
					y[2] = r[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+3 )
							b = 0;	// body 1 (O5* included)
						else
							b = 1;	// body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				// GAMMA (bond between C5* and C4*)
				iter->pos_atom = k1+4; // C5* index
				// M1 += C5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = r[0] = coord[(k1+4)*3];
					y[1] = r[1] = coord[(k1+4)*3 + 1];
					y[2] = r[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+4 )
							b = 0;    // body 1 (C5* included)
						else
							b = 1;	  // body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				//  LATERAL CHAIN --> CHI
				if(type == 2)
				{
					//printf (" Lateral Chain --> CHI\n");

					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// Screening Nucleotide-Base atoms
						// WARNING when switching to DNA (12 --- ���)
						Y1[0] = 0.0;
						Y1[1] = 0.0;
						Y1[2] = 0.0;
						// The CoM of the whole protein is = 0.0 !!!
						// We'll substract to the CoM (total) the contribution
						// of the 2nd body (CB-R). (See below)

						// M1 and M2
						mb2 = 0.0;
						for( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc();
							mb2 += mta; // atom mass
							// Substracting to the CoM (total) the body 2 (CB-R) contribution
							Y1[0] -= mta * coord[iter->pos_atom*3];
							Y1[1] -= mta * coord[iter->pos_atom*3+1];
							Y1[2] -= mta * coord[iter->pos_atom*3+2];
						}
						mpr = mb1; // body 1 mass buffer
						mb1 = mtot - mb2;
						// Computing CoM
						for(m=0;m<3;m++)
							Y1[m] /= mb1;

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
						r[0] = y[0] = coord[(k1+11)*3];
						r[1] = y[1] = coord[(k1+11)*3 + 1];
						r[2] = y[2] = coord[(k1+11)*3 + 2];
						e[0] -= y[0]; // C1* --> N1/N9 unit vector
						e[1] -= y[1];
						e[2] -= y[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Calculate matrix I2
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I2[m][n] = 0;
						// Screening Nucleotide-Base atoms
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
									I2[m] [n] += mta * mfa;
								}
							}
						}
						// Now matrix I1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = I[m][n]-I2[m][n];

						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						for( k = 0; k < num_atoms; k++ )
						{
							if ( !(k >= k1+12 && k <= k2) ) // if not body 2
								b = 0;	// body 1 + (all except Nucleotide-Base)
							else
								b = 1;	// body 2 (Nucleotide-Base)
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
						mb1 = mpr; // restoring previous mb1's value!
					}
					j++; // der[] index
				}

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{

					// EPSILON (bond between C3* and O3*)
					// M1 and M2
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
						if( iter->pos_atom != k1+8 ) // if not O3*
						{	// M1 += Sugar + Base
							mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
							mb1 += mta; // Left_half = previous + actual (NH)

							// Updating previous masses*positions
							for( m = 0; m < 3; m++ )
							{
								r[m] = coord[iter->pos_atom*3 + m];
								my[m] +=  mta * r[m];
							}
						}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = r[0] = coord[(k1+7)*3];
						y[1] = r[1] = coord[(k1+7)*3 + 1];
						y[2] = r[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0;
						for( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
							if( iter->pos_atom != k1+8 ) // if not O3*
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
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// K-matrix --> der[]
						// Computing Derivatives ( dy(b)/dTheta(lamda) )
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k <= k2 && k != k1+8 ) // not O3*
								b = 0;	// body 1 (O2P included)
							else
								b = 1;	// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index

					// ZETA (bond between O3* and next-P)
					// M1 and M2
					iter->pos_atom = k1+8; // O3* index
					mta = ( iter->get_atom() )->getPdbocc(); // N-atom mass
					// M1 += O3*
					mb1 += mta; // Left_half = previous + current (O3*)
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					// Updating previous masses*positions
					for( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = r[0] = coord[(k1+8)*3];
						y[1] = r[1] = coord[(k1+8)*3 + 1];
						y[2] = r[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0;
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
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
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// K-matrix --> der[]
						// Computing Derivatives ( dy(b)/dTheta(lamda) )
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k <= k2 )
								b = 0;	// body 1 (O3* included)
							else
								b = 1;	// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index
				}
				else
				{
					// Computing M1/M2 and "my" due to EPSILON and ZETA bodies
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
					{	// M1 += Sugar + Base
						mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
						mb1 += mta; // Left_half = previous + actual (NH)
						// Updating previous masses*positions
						for( m = 0; m < 3; m++ )
						{
							r[m] = coord[iter->pos_atom*3 + m];
							my[m] +=  mta * r[m];
						}
					}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
				}
			} // END RNA
			// if DNA fragment
			else if( resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// ********************************************************************************************
				//				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)

				// ALPHA (bond between P and O5*)
				for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
				{	// M1 += P + O1P + O2P
					// M1 and M2
					mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
					mb1 += mta; // Left_half = previous + actual (NH)

					// Updating previous masses*positions
					for ( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
				}
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = r[0] = coord[k1*3];
					y[1] = r[1] = coord[k1*3 + 1];
					y[2] = r[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+2 )
							b = 0;	// body 1 (O2P included)
						else
							b = 1;	// body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				// BETA (bond between O5* and C5*)
				iter->pos_atom = k1+3; // O5* index
				// M1 += O5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = r[0] = coord[(k1+3)*3];
					y[1] = r[1] = coord[(k1+3)*3 + 1];
					y[2] = r[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+3 )
							b = 0;	// body 1 (O5* included)
						else
							b = 1;	// body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				// GAMMA (bond between C5* and C4*)
				iter->pos_atom = k1+4; // C5* index
				// M1 += C5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = r[0] = coord[(k1+4)*3];
					y[1] = r[1] = coord[(k1+4)*3 + 1];
					y[2] = r[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Inertia Matrix 1
					for(m=0;m<3;m++)
						for(n=0;n<3;n++)
							I1[m][n] = 0.0;
					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
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
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// K-matrix --> der[]
					// Computing Derivatives ( dy(b)/dTheta(lamda) )
					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
					{
						if( k <= k1+4 )
							b = 0;    // body 1 (C5* included)
						else
							b = 1;	  // body 2
						r[0] = coord[k * 3];
						r[1] = coord[k * 3 + 1];
						r[2] = coord[k * 3 + 2];
						// = v + (w x r) (vectorial product)
						der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
						der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
						der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
					}
					j2++; // selected dihedral index
				}
				j++; // der[] index

				//  LATERAL CHAIN --> CHI
				if(type == 2)
				{
					//printf (" Lateral Chain --> CHI\n");

					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// Screening Nucleotide-Base atoms
						// WARNING when switching to DNA (12 --- ���)
						Y1[0] = 0.0;
						Y1[1] = 0.0;
						Y1[2] = 0.0;
						// The CoM of the whole protein is = 0.0 !!!
						// We'll substract to the CoM (total) the contribution
						// of the 2nd body (CB-R). (See below)

						// M1 and M2
						mb2 = 0.0;
						for( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() ) // DNA
						{
							mta = ( iter->get_atom() )->getPdbocc();
							mb2 += mta; // atom mass
							// Substracting to the CoM (total) the body 2 (CB-R) contribution
							Y1[0] -= mta * coord[iter->pos_atom*3];
							Y1[1] -= mta * coord[iter->pos_atom*3+1];
							Y1[2] -= mta * coord[iter->pos_atom*3+2];
						}
						mpr = mb1; // body 1 mass buffer
						mb1 = mtot - mb2;
						// Computing CoM
						for(m=0;m<3;m++)
							Y1[m] /= mb1;

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
						r[0] = y[0] = coord[(k1+10)*3];
						r[1] = y[1] = coord[(k1+10)*3 + 1];
						r[2] = y[2] = coord[(k1+10)*3 + 2];
						e[0] -= y[0]; // C1* --> N1/N9 unit vector
						e[1] -= y[1];
						e[2] -= y[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Calculate matrix I2
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I2[m][n] = 0;
						// Screening Nucleotide-Base atoms
						for ( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() )
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
									I2[m] [n] += mta * mfa;
								}
							}
						}
						// Now matrix I1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = I[m][n]-I2[m][n];

						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						for( k = 0; k < num_atoms; k++ )
						{
							if ( !(k >= k1+11 && k <= k2) ) // if not body 2
								b = 0;	// body 1 + (all except Nucleotide-Base)
							else
								b = 1;	// body 2 (Nucleotide-Base)
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
						mb1 = mpr; // restoring previous mb1's value!
					}
					j++; // der[] index
				}

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{

					// EPSILON (bond between C3* and O3*)
					// M1 and M2
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
						if( iter->pos_atom != k1+8 ) // if not O3*
						{	// M1 += Sugar + Base
							mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
							mb1 += mta; // Left_half = previous + actual (NH)

							// Updating previous masses*positions
							for( m = 0; m < 3; m++ )
							{
								r[m] = coord[iter->pos_atom*3 + m];
								my[m] +=  mta * r[m];
							}
						}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = r[0] = coord[(k1+7)*3];
						y[1] = r[1] = coord[(k1+7)*3 + 1];
						y[2] = r[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0;
						for( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
							if( iter->pos_atom != k1+8 ) // if not O3*
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
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// K-matrix --> der[]
						// Computing Derivatives ( dy(b)/dTheta(lamda) )
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k <= k2 && k != k1+8 ) // not O3*
								b = 0;	// body 1 (O2P included)
							else
								b = 1;	// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index

					// ZETA (bond between O3* and next-P)
					// M1 and M2
					iter->pos_atom = k1+8; // O3* index
					mta = ( iter->get_atom() )->getPdbocc(); // N-atom mass
					// M1 += O3*
					mb1 += mta; // Left_half = previous + current (O3*)
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					// Updating previous masses*positions
					for( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = r[0] = coord[(k1+8)*3];
						y[1] = r[1] = coord[(k1+8)*3 + 1];
						y[2] = r[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I1[m][n] = 0.0;
						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
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
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// K-matrix --> der[]
						// Computing Derivatives ( dy(b)/dTheta(lamda) )
						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
						{
							if( k <= k2 )
								b = 0;	// body 1 (O3* included)
							else
								b = 1;	// body 2
							r[0] = coord[k * 3];
							r[1] = coord[k * 3 + 1];
							r[2] = coord[k * 3 + 2];
							// = v + (w x r) (vectorial product)
							der[k * size + j2].x = v[b] [0] + w[b] [1] * r[2] - w[b] [2] * r[1];
							der[k * size + j2].y = v[b] [1] + w[b] [2] * r[0] - w[b] [0] * r[2];
							der[k * size + j2].z = v[b] [2] + w[b] [0] * r[1] - w[b] [1] * r[0];
						}
						j2++; // selected dihedral index
					}
					j++; // der[] index
				}
				else
				{
					// Computing M1/M2 and "my" due to EPSILON and ZETA bodies
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
					{	// M1 += Sugar + Base
						mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
						mb1 += mta; // Left_half = previous + actual (NH)
						// Updating previous masses*positions
						for( m = 0; m < 3; m++ )
						{
							r[m] = coord[iter->pos_atom*3 + m];
							my[m] +=  mta * r[m];
						}
					}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
				}
			} // END DNA
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype == tmol_smol )
			{
				if(debug3)
					printf("SMOL  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

				for( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
				{
					mta = iter->get_atom()->getPdbocc(); // Non-(C/O) atoms mass
					mb1 += mta;
					my[0] += mta * coord[iter->pos_atom*3];
					my[1] += mta * coord[iter->pos_atom*3+1];
					my[2] += mta * coord[iter->pos_atom*3+2];
				}
			}
			else
			{
			    // NOT FOUND MOL-TYPE
				printf("Msg(dydqMFAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			if( !( (iter_frag->pos_fragment == num_res-1) && (iter_seg->pos_segment == num_seg - 1) ) )
			{
				k1 += props[res_index].nat;
				k2 += props[res_index + 1].nat;
			}
			res_index++;
		}
		delete iter_frag;
	}
	// end derivatives
	iter_seg->clean_virtual();
	delete iter_seg;
	delete iter;
}

// Compute V/W-arrays - Memory Efficient (see ec. 40 Braun et al.) (Multi-Chain + Fixation)
// NEEDS: -Macromolecule (N,CA,CA model) (3-atoms model),
//        -Single row (N,CA,C)-PDB coordinates,
//        -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void vwMCA3x(Macromolecule *mol,float *coord, double ****p_V, double ****p_W, int ***p_body1, int size, bool *fix, int model)
{
	bool debug = false;
	double mtot,mta;
	float r[3];
	double rd[3];
	double I[3][3],J[3][3];
	double temp1,temp2;
	int m, n, num_atoms, i, k, j, j2;
	Residue *res;
	int resn,num_res,num_seg;

	num_atoms = mol->get_num_atoms(); // Note this will be "3*num_res" (CA-model)
	num_res = mol->get_num_fragments();

	pdbIter *iter = new pdbIter( mol ); //Iterador para recorrer atomos
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_seg; // Iter to screen all segments

	num_seg = iter->num_segment();

	// V and W vector arrays Memory Allocation: der[] --> N x size -matrix !!!! (DERIVATIVES)
	double ***V;
	double ***W;
	int **body1;

	// V-array memory allocation
	if(*p_V == NULL)
	{ // Allocate memory
		if( !( V = ( double *** ) malloc( size * sizeof( double ** ) ) ) )
		{
			printf("Msg(vwMCA3x): V-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_V = V; // outputs V-matrix
		// Allocating V-arrays
		for ( i = 0; i < size; i++ )
		{
			V[i] = (double **) malloc( sizeof(double *) * 2);
			V[i][0] = (double *) malloc( sizeof(double) * 3 );
			V[i][1] = (double *) malloc( sizeof(double) * 3 );
		}
	}
	else // Use already allocated memory
		V = *p_V;
	// Initializing V-arrays
	for ( i = 0; i < size; i++ )
	{
		V[i][0][0] = 0.0;
		V[i][0][1] = 0.0;
		V[i][0][2] = 0.0;
		V[i][1][0] = 0.0;
		V[i][1][1] = 0.0;
		V[i][1][2] = 0.0;
	}

	// W-array memory allocation
	if(*p_W == NULL)
	{ // Allocate memory
		if( !( W = ( double *** ) malloc( size * sizeof( double ** ) ) ) )
		{
			printf("Msg(vwMCA3x): W-array memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_W = W; // outputs V-array
		// Allocating W-arrays
		for ( i = 0; i < size; i++ )
		{
			W[i] = (double **) malloc( sizeof(double *) * 2);
			W[i][0] = (double *) malloc( sizeof(double) * 3 );
			W[i][1] = (double *) malloc( sizeof(double) * 3 );
		}
	}
	else // Use already allocated memory
		W = *p_W;
	// Initializing W-arrays
	for ( i = 0; i < size; i++ )
	{
		W[i][0][0] = 0.0;
		W[i][0][1] = 0.0;
		W[i][0][2] = 0.0;
		W[i][1][0] = 0.0;
		W[i][1][1] = 0.0;
		W[i][1][2] = 0.0;
	}

	// body1 array memory allocation
	if(*p_body1 == NULL)
	{ // Allocate memory
		if( !( body1 = ( int ** ) malloc( size * sizeof( int *) ) ) )
		{
			printf("Msg(vwMCA3x): body1-array memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_body1 = body1; // outputs body1-matrix
		// Allocating V/W-arrays
		for ( i = 0; i < size; i++ )
			if( !(body1[i] = (int *) malloc( sizeof(int) * 1)) )
			{
				printf("Msg(vwMCA3x): body1-array memory allocation failed!\nForcing exit!\n");
				exit(1);
			}
	}
	else // Use already allocated memory
		body1 = *p_body1;
	// Initializing body1-arrays
	for ( i = 0; i < size; i++ )
		body1[i][0] = 0;

	// Initializing Inertia-Matrix
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			I[i][j] = 0.0;

	// Compute matrix I (Inertia Momentum) and Total Mass
	// Ec. 22 (Noguti & Go, 1983)
	i = 0;
	mtot = 0.0;
	for ( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screens all atoms
	{
		mta = ( iter->get_atom() )->getPdbocc(); // mass
		if(mta != 0.0) // Only "real" atoms (to be valid for CA and NCAC models)
		{
			mtot += mta; // total mass
			r[0] = coord[iter->pos_atom * 3]; // atom position in coord
			r[1] = coord[iter->pos_atom * 3 + 1];
			r[2] = coord[iter->pos_atom * 3 + 2];
			temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
			for ( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
				{
					temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
					if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
					I[m][n] += mta * temp2;
				}
		}
	}
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	// **********************************
	// ** V / W COMPUTATION *************
	// **********************************
	double temp;
	double my[3], e[3], y[3], Y1[3], I1[3] [3], I2[3] [3], v[2] [3], w[2] [3]; // calcoef's
	double nr2, mfa;
	float mb1, mb2;
	int k1, b, kCA;
	float pos[3];

	// inertia matrix 1 initialization
	for ( m = 0; m < 3; m++ )
		for ( n = 0; n < 3; n++ )
			I1[m] [n] = 0.0;

	j = 0; // current dihedral index
	j2 = 0; // current dihedral index
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */
	k1 = 0;
	mb1 = 0.0;
	kCA = 0; // CA-model index

	// ********************************************************************************************
	// ************* DERIVATIVES LOOP BEGIN *******************************************************
	// ********************************************************************************************

	Segment * seg;
	Atom * atom;
	// Screening segments
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments

	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// INTER-SEGMENT DEGREES OF FREEDOM (3 Translations + 3 Rotations)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				// CoM 1
				for ( m = 0; m < 3; m++ )
					Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM

				// now matrix I2
				for( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
						I2[m][n] = I[m][n] - I1[m][n];
			}

			// 3 TRANSLATIONS
			//  printf ("Translation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{

					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // Phi
						y[m] = 0.0; // Psi
					}
					y[axis] = -1.0; // Psi = -gamma_v
					// gamma_v --> Column unit vectors along positive x, y, z axes of S-system with respect to F-system

					// Computing some coefficients (Multi-Chain)
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					// Storing which body an atom belongs to.
					if(model == 0) // If CA-model things are different
						body1[j2][0] = kCA; // atom_index < kCA --> body 1
					else
						body1[j2][0] = k1; // atom_index < k1 --> body 1

					j2++;
				}
				j++; // der[] index
			}

			// 3 ROTATIONS
			//  printf ("Rotation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0;
						r[m] = 0.0;
					}
					e[axis] = 1.0; // Phi

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					// Storing which body an atom belongs to.
					if(model == 0) // If CA-model things are different
						body1[j2][0] = kCA; // atom_index < kCA --> body 1
					else
						body1[j2][0] = k1; // atom_index < k1 --> body 1

					j2++;
				}
				j++; // der[] index
			}
		}

		// Screen ALL residues
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			// NH terms
			iter->pos_atom = k1; // NH
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if(mta != 0.0) // Only "real" atoms (to be valid for CA and NCAC models)
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[k1 * 3]; // NH atom
				my[1] += mta * coord[k1 * 3+1];
				my[2] += mta * coord[k1 * 3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // N,CA,C positions
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

				if(model == 0)
					kCA++; // Permits CA-model compatibility
			}

			// ********************************************************************************************
			// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (Ca-model)
			if( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
			{
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// mb1, mb2 and my, already updated before

					// e_lambda
					// get CA pos 1
					e[0] = pos[0] = coord[(k1+1) * 3]; // CA position
					e[1] = pos[1] = coord[(k1+1) * 3 + 1];
					e[2] = pos[2] = coord[(k1+1) * 3 + 2];

					// get NH pos 0
					y[0] = r[0] = coord[k1 * 3];
					y[1] = r[1] = coord[k1 * 3 + 1];
					y[2] = r[2] = coord[k1 * 3 + 2];
					e[0] -= y[0]; // NH --> CA
					e[1] -= y[1];
					e[2] -= y[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi = (e x r)
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					// Storing which body an atom belongs to.
					if(model == 0) // If CA-model things are different
						body1[j2][0] = kCA; // atom_index < kCA --> body 1
					else
						body1[j2][0] = k1; // atom_index < k1 --> body 1

					j2++;
				}
				j++; // der[] index
			}  // NOT FIRST NOT PRO

			// CA terms
			iter->pos_atom = k1+1; // residue CA (contains the whole residue mass)
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if( mta != 0.0 ) // Only "real" atoms (to be valid for CA and NCAC models)
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[(k1+1)*3]; // CA atom
				my[1] += mta * coord[(k1+1)*3+1];
				my[2] += mta * coord[(k1+1)*3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // atom positions
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

				if(model == 0)
					kCA++; // Permits CA-model compatibility
			}

			// ********************************************************************************************
			// "PSI-bond" (In CA-only model, non-first and non-last)
			if( !(iter_frag->pos_fragment == num_res-1 || iter_frag->pos_fragment == 0) )
			{
				if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
				{
					// mb1, mb2 and my, already updated above

					// Gets e_lambda (CA-->) Psi-bond
					e[0] = coord[(k1+2) * 3]; // C position
					e[1] = coord[(k1+2) * 3 + 1];
					e[2] = coord[(k1+2) * 3 + 2];

					// get CA pos 1
					y[0] = pos[0] = coord[(k1+1) * 3]; // CA position
					y[1] = pos[1] = coord[(k1+1) * 3 + 1];
					y[2] = pos[2] = coord[(k1+1) * 3 + 2];

					e[0] -= y[0]; // CA --> C
					e[1] -= y[1];
					e[2] -= y[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m]/mb1; // Updating Body 1 Center of Mass: "Y1[m]"
					}

					// now matrix I2
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
							I2[m][n] = I[m][n] - I1[m][n];

					y[0] = e[1] * pos[2] - e[2] * pos[1]; // Psi (e x r)
					y[1] = e[2] * pos[0] - e[0] * pos[2];
					y[2] = e[0] * pos[1] - e[1] * pos[0];

					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					// Storing which body an atom belongs to.
					if(model == 0) // If CA-model things are different
						body1[j2][0] = kCA; // atom_index < kCA --> body 1
					else
						body1[j2][0] = k1; // atom_index < k1 --> body 1

					j2++;
				}
				j++; // der[] index
			}  // NOT LAST & NOT FIRST

			// CO terms
			iter->pos_atom = k1+2; // residue CA (contains the whole residue mass)
			mta = ( iter->get_atom() )->getPdbocc(); // mass
			if( mta != 0.0 ) // Only "real" atoms (to be valid for CA and NCAC models)
			{
				// Masses
				mb1 += mta;
				mb2 = mtot - mb1;
				// Momentum
				my[0] += mta * coord[(k1+2)*3]; // CA atom
				my[1] += mta * coord[(k1+2)*3+1];
				my[2] += mta * coord[(k1+2)*3+2];
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3]; // atom positions
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

				if(model == 0)
					kCA++; // Permits CA-model compatibility
			}

			k1 += 3; // current atom index (first atom of the residue)
		}
		delete iter_frag;
	}

	// end derivatives
	iter_seg->clean_virtual();
	delete iter_seg;
	delete iter;
}


// Compute V/W-arrays - Memory Efficient (see ec. 40 Braun et al.) (Multi-Chain + Fixation)
// (Multi-Chain & FULL-ATOM & Protein/RNA/DNA)
// NEEDS: Macromolecule, Single row PDB coordinates,
// -Properties structure array (props[]), -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: type = 0 --> phi,psi,... type = 1 & 2 --> phi,chi,psi...
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void vwMFAx(Macromolecule *molr,float *coord, tri *props, double ****p_V, double ****p_W, int ***p_body1, int type, int model, int size, bool *fix, bool *addrot)
{
	bool debug=false;
	bool debug2=false;
	double mtot,mta;
	double r[3];
	double rd[3];
	double I[3][3],J[3][3];
	double temp1,temp2;
	int m, n, num_atoms, i, k, j;
	Residue *res;
	int resn,num_res,num_seg,res_index=0;
	int CBindex;
	int lim;
	bool has_oxt = false;

	// Needed below for proper indexing of C5 and HA CG-models...
	if(model==1) // C5
		CBindex = 3;
	else // HA
		CBindex = 4;

	num_atoms = molr->get_num_atoms();
	num_res = molr->get_num_fragments();

	pdbIter *iter = new pdbIter( molr ); //Iterador para recorrer atomos
	pdbIter *iter_frag; // Iter to screen all fragments
	pdbIter *iter_seg; // Iter to screen all segments

	// V and W vector arrays Memory Allocation: der[] --> N x size -matrix !!!! (DERIVATIVES)
	double ***V;
	double ***W;
	int **body1; // integer (num_atoms x 4) matrix to store body boundary atom index ( atom_index < body0[i] will belong to body 1)

	if(*p_V == NULL)
	{ // Allocate memory
		if( !( V = ( double *** ) malloc( size * sizeof( double ** ) ) ) )
		{
			printf("Msg(vwMFAx): V-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_V = V; // outputs V-matrix
		// Allocating V-arrays
		for ( i = 0; i < size; i++ )
		{
			V[i] = (double **) malloc( sizeof(double *) * 2);
			V[i][0] = (double *) malloc( sizeof(double) * 3 );
			V[i][1] = (double *) malloc( sizeof(double) * 3 );
		}
	}
	else // Use already allocated memory
		V = *p_V;
	// Initializing V-arrays
	for ( i = 0; i < size; i++ )
	{
		V[i][0][0] = 0.0;
		V[i][0][1] = 0.0;
		V[i][0][2] = 0.0;
		V[i][1][0] = 0.0;
		V[i][1][1] = 0.0;
		V[i][1][2] = 0.0;
	}

	if(*p_W == NULL)
	{ // Allocate memory
		if( !( W = ( double *** ) malloc( size * sizeof( double ** ) ) ) )
		{
			printf("Msg(vwMFAx): W-array memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_W = W; // outputs V-array
		// Allocating W-arrays
		for ( i = 0; i < size; i++ )
		{
			W[i] = (double **) malloc( sizeof(double *) * 2);
			W[i][0] = (double *) malloc( sizeof(double) * 3 );
			W[i][1] = (double *) malloc( sizeof(double) * 3 );
		}
	}
	else // Use already allocated memory
		W = *p_W;
	// Initializing W-arrays
	for ( i = 0; i < size; i++ )
	{
		W[i][0][0] = 0.0;
		W[i][0][1] = 0.0;
		W[i][0][2] = 0.0;
		W[i][1][0] = 0.0;
		W[i][1][1] = 0.0;
		W[i][1][2] = 0.0;
	}

	if(*p_body1 == NULL)
	{ // Allocate memory
		if( !( body1 = ( int ** ) malloc( size * sizeof( int * ) ) ) )
		{
			printf("Msg(vwMFAx): body1-array memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		else
			*p_body1 = body1; // outputs body1-matrix
		// Allocating V/W-arrays
		for ( i = 0; i < size; i++ )
		{
//			if( !(body1[i] = (int *) malloc( sizeof(int) * num_atoms * 3)) )
			if( !(body1[i] = (int *) malloc( sizeof(int) * 3)) )
			{
				printf("Msg(vwMFAx): body1-array memory allocation failed!\nForcing exit!\n");
				exit(1);
			}
		}
	}
	else // Use already allocated memory
		body1 = *p_body1;

	// Initializing Inertia-Matrix
	for ( i = 0; i < 3; i++ )
		for ( j = 0; j < 3; j++ )
			I[i][j] = 0.0;

	// Compute matrix I (Inertia Momentum) and Total Mass
	// Ec. 22 (Noguti & Go, 1983)
	i = 0;
	mtot = 0.0;
	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screens all-atoms
	{
		mta = ( iter->get_atom() )->getPdbocc(); // mass
		mtot += mta; // total mass
		r[0] = coord[i * 3];
		r[1] = coord[i * 3 + 1];
		r[2] = coord[i * 3 + 2];
		i++;
		temp1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2]; // Sum,k( y^2(ka) ) // ¿Diagonal?
		for ( m = 0; m < 3; m++ )
			for ( n = 0; n < 3; n++ )
			{
				temp2 = -r[m] * r[n]; // -y(ia)*y(ja)
				if ( m == n ) temp2 += temp1; // delta_Dirac(ij)* Sum,k( y^2(ka) )
				I[m][n] += mta * temp2;
			}
	}
	inverse( I, J ); // J = Inverse of I (3x3 matrix)

	// *************************************************
	// ** now the actual DERIVATIVE COMPUTATION ********
	// (... ranging over the residues is more practical)
	// *************************************************
	double mpr,temp, my[3], e[3], y[3], Y1[3], I1[3][3], I2[3][3], Ix[3][3], v[2] [3], w[2] [3];
	double nr2, mfa;
	double mb1, mb2;
	int k1, k2, b;
	float pos[3];
	TMOL fragtype;

	for ( m = 0; m < 3; m++ )
		for ( n = 0; n < 3; n++ )
			I1[m] [n] = 0.0;

	j = 0; // current dihedral index (screens all dihedrals)
	int j2 = 0; // active dihedral index (screens all mobile dihedrals)
	mpr=0; /* mpr = mass of previous residues */
	my[0] = my[1] = my[2] = 0.0; /* Sum(mass*coord) of previous residues */
	k1 = 0;
	k2 = props[0].nat - 1; /* 1st and last indices of atoms of residue i */ // Residue 0 ?
	int indexbase;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

	// ********************************************************************************************
	// ************* DERIVATIVES LOOP BEGIN *******************************************************
	// ********************************************************************************************

	iter_seg = new pdbIter( molr, true, true, true, true ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	Segment *seg;
	Atom *atom;

	// Screening segments
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// M1 and M2
			// (mb1 already updated)
			mb2 = mtot - mb1; // Right-half = Tot - Left_half

			// CoM 1
			for( m = 0; m < 3; m++ )
				Y1[m] = my[m] / mb1; // Previous CoM = body 1 CoM

//			// Inertia Matrix 1
//			for(m=0;m<3;m++)
//				for(n=0;n<3;n++)
//					I1[m][n] = 0.0;
//			for( iter->pos_atom = 0; iter->pos_atom < k1; iter->next_atom() )
//			{
//				mta = ( iter->get_atom() )->getPdbocc(); // mass
//				r[0] = coord[iter->pos_atom*3];
//				r[1] = coord[iter->pos_atom*3+1];
//				r[2] = coord[iter->pos_atom*3+2];
//				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//				for ( m = 0; m < 3; m++ )
//				{
//					for ( n = 0; n < 3; n++ )
//					{
//						mfa = -r[m] * r[n];
//						if ( m == n ) mfa += nr2;
//						I1[m][n] +=  mta * mfa;
//					}
//				}
//			}
			// now matrix I2
			for( m = 0; m < 3; m++ )
				for ( n = 0; n < 3; n++ )
					I2[m][n] = I[m][n] - I1[m][n];
			// End Inertia matrix computation

			// 3 TRANSLATIONS
			//  printf ("Translation \n");
			for(int axis=0; axis<3; axis++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's fixed
				{
					// e_lambda
					for ( m = 0; m < 3; m++ )
					{
						e[m] = 0.0; // Phi = 0 (translation)
						y[m] = 0.0; // Psi
					}
					y[axis] = -1.0; // Psi = -gamma_v
					// gamma_v --> Column unit vectors along positive x, y, z axes of S-system with respect to F-system

					// Computing some coefficients (Multi-Chain)
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1; // k < k1 --> body1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//					{
//						if( k < k1 )
//							body1[j2][k] = true;    // body 1
//						else
//							body1[j2][k] = false;    // body 2
//					}

					j2++;
				}
				j++; // der[] index
			}

			// 3 ROTATIONS
			//  printf ("Rotation \n");
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{
						// e_lambda
						for ( m = 0; m < 3; m++ )
						{
							e[m] = 0.0; // Phi
							r[m] = 0.0; // Psi (in rotation is the origin: 0,0,0)
						}
						e[axis] = 1.0; // Phi

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi = 0 ????
						y[1] = e[2] * r[0] - e[0] * r[2]; // Put Psi=0 directly whenever!
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k1; // k < k1 --> body1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = -1;   // No "!=" condition
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k < k1 )
//								body1[j2][k] = true;    // body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++;
					}
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // der[] index
				}

		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				// ********************************************************************************************
				// PHI --> NOT FIRST, NOT PRO
				if ( iter_frag->pos_fragment != 0 && resn != PRO )
				{
					//  printf (" Not first not PRO -> PHI \n");

					// M1 and M2
					// Computing masses... (updating)
					iter->pos_atom = k1; // N-atom
					mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
					mb1 += mta; // Left_half = previous + actual (NH)
					mb2 = mtot - mb1; // Right-half = Tot - Left_half

					// updating previous masses*positions
					y[0] = r[0] = coord[k1 * 3];
					y[1] = r[1] = coord[k1 * 3 + 1];
					y[2] = r[2] = coord[k1 * 3 + 2];
					for ( m = 0; m < 3; m++ )
						my[m] +=  mta * y[m];

					// Body 1 Inertia matrix accumulation (I1)
//					r[0] = coord[iter->pos_atom*3];
//					r[1] = coord[iter->pos_atom*3+1];
//					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get CA pos 1
						e[0] = coord[(k1+1) * 3];
						e[1] = coord[(k1+1) * 3 + 1];
						e[2] = coord[(k1+1) * 3 + 2];
						// get NH pos 0
						//					y[0] = r[0] = coord[k1 * 3];
						//					y[1] = r[1] = coord[k1 * 3 + 1];
						//					y[2] = r[2] = coord[k1 * 3 + 2];
						e[0] -= y[0]; // NH --> CA
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}
						// printf("PHI M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0;
//						for ( iter->pos_atom = 0; iter->pos_atom <= k1; iter->next_atom() )
//						{
//							mta = ( iter->get_atom() )->getPdbocc(); // mass
//							r[0] = coord[iter->pos_atom*3];
//							r[1] = coord[iter->pos_atom*3+1];
//							r[2] = coord[iter->pos_atom*3+2];
//							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//							for ( m = 0; m < 3; m++ )
//							{
//								for ( n = 0; n < 3; n++ )
//								{
//									mfa = -r[m] * r[n];
//									if ( m == n ) mfa += nr2;
//									I1[m][n] +=  mta * mfa;
//								}
//							}
//						}

//						// Adding NH contribution to I1
//						iter->pos_atom = k1;
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k1+1; // k <= k1 --> body1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = -1;   // No "!=" condition
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k <= k1 )
//								body1[j2][k] = true;    // body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}
					j++; // der[] index
				}  // NOT FIRST NOT PRO

				// Checking if last residue has OXT (Ot)
				if(model==2)
					if(iter_frag->pos_fragment == num_res-1)
					{
						iter->pos_atom = k2;
						if( strcmp(iter->get_atom()->getName()," OXT ")==0 )
						{
							k2--;
							has_oxt = true;
						}
					}

				// ********************************************************************************************
				//  LATERAL CHAIN-->CHI
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if(type == 2) // phi,chi,psi,...
					if( props[res_index].nan==3 ||
							(props[res_index].nan==2 &&
									(iter_frag->pos_fragment==0 ||
											(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						//printf (" Lateral Chain --> CHI\n");

						if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
						{
							// M1 and M2
							mb2 = 0.0;
							// screening CB & -R atoms (CB: "k1+3")
							for( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
								mb2 += ( iter->get_atom() )->getPdbocc(); // atom mass
							// Now CB & -R is the whole 2nd body, the remaining atoms will be the 1st one!
							mpr = mb1; // body 1 mass buffer
							mb1 = mtot - mb2;

							// e_lambda
							// get CB pos 4
							e[0] = pos[0] = coord[(k1+CBindex) * 3 + 0];
							e[1] = pos[1] = coord[(k1+CBindex) * 3 + 1];
							e[2] = pos[2] = coord[(k1+CBindex) * 3 + 2];
							// get CA pos 1
							r[0] = y[0] = pos[0] = coord[(k1+1) * 3 + 0];
							r[1] = y[1] = pos[1] = coord[(k1+1) * 3 + 1];
							r[2] = y[2] = pos[2] = coord[(k1+1) * 3 + 2];
							e[0] -= y[0]; // CA-->CB unit vector
							e[1] -= y[1];
							e[2] -= y[2];
							temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
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

							y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
							y[1] = e[2] * r[0] - e[0] * r[2];
							y[2] = e[0] * r[1] - e[1] * r[0];

							// this should be expressed different!
							// screening [k1+4,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
							{
								mta = ( iter->get_atom() )->getPdbocc(); // atom mass
								// Substracting to the CoM (total) the body 2 (CB-R) contribution
								Y1[0] -= mta * coord[iter->pos_atom*3];
								Y1[1] -= mta * coord[iter->pos_atom*3+1];
								Y1[2] -= mta * coord[iter->pos_atom*3+2];
							}
							for(m=0;m<3;m++)
								Y1[m] /= mb1;

							// Calculate matrix I2
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									I2[m][n] = 0.0;
							// screening [k1+4,k2]   CB & -R
							for ( iter->pos_atom = k1+CBindex; iter->pos_atom <= k2; iter->next_atom() )
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
										I2[m] [n] += mta * mfa;
									}
								}
							}
							// Now matrix I1
							for(m=0;m<3;m++)
								for(n=0;n<3;n++)
									Ix[m][n] = I[m][n]-I2[m][n];

							calcoefM( mb1, mb2, e, y, Y1, Ix, I2, J, v, w ); // (ec. 40 Braun et al.)

							// loading "v" and "w" into V/W-arrays
							V[j2][0][0] = v[0][0];
							V[j2][0][1] = v[0][1];
							V[j2][0][2] = v[0][2];
							V[j2][1][0] = v[1][0];
							V[j2][1][1] = v[1][1];
							V[j2][1][2] = v[1][2];
							W[j2][0][0] = w[0][0];
							W[j2][0][1] = w[0][1];
							W[j2][0][2] = w[0][2];
							W[j2][1][0] = w[1][0];
							W[j2][1][1] = w[1][1];
							W[j2][1][2] = w[1][2];

							body1[j2][0] = k1+CBindex; // k < k1+CBindex --> body1
							body1[j2][1] = k2; // k > k2 --> body1
							body1[j2][2] = -1;   // No "!=" condition
//							for( k = 0; k < num_atoms; k++ )
//								if ( k < k1+CBindex || k > k2 )
//									body1[j2][k] = true;    // body 1
//								else
//									body1[j2][k] = false;    // body 2

							j2++; // selected dihedral index
							mb1 = mpr; // restoring mb1's value!
						}
						j++; // der[] index
					}

				// ********************************************************************************************
				// NOT LAST RESIDUE--> PSI
//				int lim = k1+1; // not includes N (it was previously included in PHI)
				if(iter_frag->pos_fragment==0 || resn == PRO) // if no-PHI (First or PRO)
					lim = k1; // includes N
				else
					lim = k1+1; // not includes N (it was previously included in PHI)

				for ( iter->pos_atom = lim; iter->pos_atom <= k2; iter->next_atom() )
					if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
					{
						mta = iter->get_atom()->getPdbocc(); // Non-(C/O) atoms mass
						mb1 += mta;
						my[0] += mta * coord[iter->pos_atom*3];
						my[1] += mta * coord[iter->pos_atom*3+1];
						my[2] += mta * coord[iter->pos_atom*3+2];

						// Body 1 Inertia matrix accumulation (I1)
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
					}

				// "PSI-bond" (in Full-Atom, every residue has PSI)
				if ( iter_frag->pos_fragment != num_res - 1 || model==2 ) // Full-Atom always has PSI
				{
					// Current body 2 mass (M2)
					mb2 = mtot - mb1;

					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// get C pos +2
						e[0] = coord[(k1+2) * 3 + 0];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// get CA pos +1
						y[0] = pos[0] = coord[(k1+1) * 3 + 0];
						y[1] = pos[1] = coord[(k1+1) * 3 + 1];
						y[2] = pos[2] = coord[(k1+1) * 3 + 2];
						e[0] -= y[0]; // CA --> C
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// Updating Body 1 Center of Mass: "Y1[m]"
						for(m=0;m<3;m++) // 3D coords
							Y1[m] = my[m]/mb1;

						// printf("PSI M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0; // Initialization
						// Screens ALL body 1
//						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )

//						for ( iter->pos_atom = k1+1; iter->pos_atom <= k2; iter->next_atom() )
//							if( !(iter->pos_atom == k1+2 || (model == 2 && iter->pos_atom == k1+3) ) ) // NOT (C or O)
//							{
//								mta = ( iter->get_atom() )->getPdbocc(); // mass
//								// Do this with iterators: **********************
//								r[0] = coord[iter->pos_atom*3];
//								r[1] = coord[iter->pos_atom*3+1];
//								r[2] = coord[iter->pos_atom*3+2];
//
//								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//								for( m = 0; m < 3; m++ )
//									for( n = 0; n < 3; n++ )
//									{
//										mfa = -r[m] * r[n];
//										if ( m == n ) mfa += nr2;
//										I1[m][n] +=  mta * mfa;
//									}
//							}

						// now matrix I2
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						y[0] = e[1] * pos[2] - e[2] * pos[1]; // Psi
						y[1] = e[2] * pos[0] - e[0] * pos[2];
						y[2] = e[0] * pos[1] - e[1] * pos[0];

						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						if(model == 2) // HA
						{
							body1[j2][0] = k2+1; // k <= k2 --> body 1
							body1[j2][1] = -1;   // No ">" condition
							body1[j2][2] = k1+3; // k != k1+3 --> body 1
						}
						else
						{
							body1[j2][0] = k2+1; // k <= k2 --> body 1
							body1[j2][1] = -1;   // No ">" condition
							body1[j2][2] = -1;   // No "!=" condition
						}

//						if(model == 2) // HA
//							for ( k = 0; k < num_atoms; k++ ) // Screens ALL atoms
//								if (  k <= k2 && k != k1+3 ) // note that OXT is not in body 1
//									body1[j2][k] = true;    // body 1
//								else
//									body1[j2][k] = false;    // body 2
//						else // C5
//							for ( k = 0; k < num_atoms; k++ ) // Screens ALL atoms
//								if (  k <= k2  ) // In C5-model there is no OXT!
//									body1[j2][k] = true;    // body 1
//								else
//									body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}

					iter->pos_atom = k1+2;
					mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
					mb1 += mta;
					my[0] += mta * coord[iter->pos_atom*3];
					my[1] += mta * coord[iter->pos_atom*3+1];
					my[2] += mta * coord[iter->pos_atom*3+2];
					// Body 1 Inertia matrix accumulation (I1)
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}

					if(model == 2) // then O-atom must be taken into account!
					{			   // (if model==1 --> already accounted for...
						iter->pos_atom = k1+3;
						mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
						mb1 += mta;
						my[0] += mta * coord[iter->pos_atom*3];
						my[1] += mta * coord[iter->pos_atom*3+1];
						my[2] += mta * coord[iter->pos_atom*3+2];
						// Body 1 Inertia matrix accumulation (I1)
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for( m = 0; m < 3; m++ )
							for( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
					}

					// Taking into account whether last residue has OXT (Ot)
					if(iter_frag->pos_fragment == num_res-1)
					{
						if(has_oxt)
						{
							k2++;
							iter->pos_atom = k2; // now it's "k2+1" (above: k2--)
							mta = iter->get_atom()->getPdbocc(); // Non-(N/C/O) atoms mass
							mb1 += mta;
							my[0] += mta * coord[iter->pos_atom*3];
							my[1] += mta * coord[iter->pos_atom*3+1];
							my[2] += mta * coord[iter->pos_atom*3+2];
							// Body 1 Inertia matrix accumulation (I1)
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for( m = 0; m < 3; m++ )
								for( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
						}
						has_oxt = false; // initialice for the next segment
					}

					j++; // der[] index
				}

				// Updating last C-atom if 3BB2R...
				if(model==1 && iter_frag->pos_fragment == num_res - 1)
				{
					iter->pos_atom = k1+2; // C-atom index
					mta = iter->get_atom()->getPdbocc(); // C-atom mass
					mb1 += mta;
					my[0] += mta * coord[iter->pos_atom*3];
					my[1] += mta * coord[iter->pos_atom*3+1];
					my[2] += mta * coord[iter->pos_atom*3+2];
					// Body 1 Inertia matrix accumulation (I1)
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for( m = 0; m < 3; m++ )
						for( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
				}

			}
			// if RNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				// ********************************************************************************************
				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)

				// ALPHA (bond between P and O5*)
				for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
				{	// M1 += P + O1P + O2P
					// M1 and M2
					mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
					mb1 += mta; // Left_half = previous + actual (NH)

					// Updating previous masses*positions
					for ( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}

					// Inertia Matrix 1
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
				}
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = r[0] = coord[k1*3];
					y[1] = r[1] = coord[k1*3 + 1];
					y[2] = r[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+3; // k <= k1+2 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+2 )
//							body1[j2][k] = true;    // body 1  (O2P included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				// BETA (bond between O5* and C5*)
				iter->pos_atom = k1+3; // O5* index
				// M1 += O5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3];
				r[1] = coord[iter->pos_atom*3+1];
				r[2] = coord[iter->pos_atom*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = r[0] = coord[(k1+3)*3];
					y[1] = r[1] = coord[(k1+3)*3 + 1];
					y[2] = r[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+4; // k <= k1+3 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+3 )
//							body1[j2][k] = true;	// body 1 (O5* included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				// GAMMA (bond between C5* and C4*)
				iter->pos_atom = k1+4; // C5* index
				// M1 += C5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3];
				r[1] = coord[iter->pos_atom*3+1];
				r[2] = coord[iter->pos_atom*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = r[0] = coord[(k1+4)*3];
					y[1] = r[1] = coord[(k1+4)*3 + 1];
					y[2] = r[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+4; // k <= k1+3 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "1=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+4 )
//							body1[j2][k] = true;    // body 1 (C5* included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				//  LATERAL CHAIN --> CHI
				if(type == 2)
				{
					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// Screening Nucleotide-Base atoms
						// WARNING when switching to DNA (12 --- ���)
						Y1[0] = 0.0;
						Y1[1] = 0.0;
						Y1[2] = 0.0;
						// The CoM of the whole protein is = 0.0 !!!
						// We'll substract to the CoM (total) the contribution
						// of the 2nd body (CB-R). (See below)

						// M1 and M2
						mb2 = 0.0;
						for( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc();
							mb2 += mta; // atom mass
							// Substracting to the CoM (total) the body 2 (CB-R) contribution
							Y1[0] -= mta * coord[iter->pos_atom*3];
							Y1[1] -= mta * coord[iter->pos_atom*3+1];
							Y1[2] -= mta * coord[iter->pos_atom*3+2];
						}
						mpr = mb1; // body 1 mass buffer
						mb1 = mtot - mb2;
						// Computing CoM
						for(m=0;m<3;m++)
							Y1[m] /= mb1;

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
						r[0] = y[0] = coord[(k1+11)*3];
						r[1] = y[1] = coord[(k1+11)*3 + 1];
						r[2] = y[2] = coord[(k1+11)*3 + 2];
						e[0] -= y[0]; // C1* --> N1/N9 unit vector
						e[1] -= y[1];
						e[2] -= y[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Calculate matrix I2
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								Ix[m][n] = 0;
						// Screening Nucleotide-Base atoms
						for ( iter->pos_atom = k1+12; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc();
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for ( m = 0; m < 3; m++ )
								for ( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									Ix[m] [n] += mta * mfa;
								}
						}
						// Now matrix I1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								Ix[m][n] = I[m][n]-I2[m][n];

						calcoefM( mb1, mb2, e, y, Y1, Ix, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k1+13; // k <= k1+12 --> body 1
						body1[j2][1] = k2;   // k > k2 --> body 1
						body1[j2][2] = -1;   // No "1=" condition
//						for( k = 0; k < num_atoms; k++ )
//							if ( k < k1+12 || k > k2 ) // if body 1
//								body1[j2][k] = true;	// body 1 + (all except Nucleotide-Base)
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
						mb1 = mpr; // restoring previous mb1's value!
					}
					j++; // der[] index
				}

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{

					// EPSILON (bond between C3* and O3*)
					// M1 and M2
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
						if( iter->pos_atom != k1+8 ) // if not O3*
						{	// M1 += Sugar + Base
							mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
							mb1 += mta; // Left_half = previous + actual (NH)

							// Updating previous masses*positions
							for( m = 0; m < 3; m++ )
							{
								r[m] = coord[iter->pos_atom*3 + m];
								my[m] +=  mta * r[m];
							}
							// Inertia Matrix 1
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for ( m = 0; m < 3; m++ )
								for ( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
						}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = r[0] = coord[(k1+7)*3];
						y[1] = r[1] = coord[(k1+7)*3 + 1];
						y[2] = r[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

//						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0;
//						for( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
//							if( iter->pos_atom != k1+8 ) // if not O3*
//							{
//								mta = ( iter->get_atom() )->getPdbocc(); // mass
//								r[0] = coord[iter->pos_atom*3];
//								r[1] = coord[iter->pos_atom*3+1];
//								r[2] = coord[iter->pos_atom*3+2];
//								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//								for ( m = 0; m < 3; m++ )
//								{
//									for ( n = 0; n < 3; n++ )
//									{
//										mfa = -r[m] * r[n];
//										if ( m == n ) mfa += nr2;
//										I1[m][n] +=  mta * mfa;
//									}
//								}
//							}
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k2+1; // k <= k2 --> body 1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = k1+8; // k != k1+8 --> body 1
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k <= k2 && k != k1+8 ) // not O3*
//								body1[j2][k] = true;	// body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}
					j++; // der[] index

					// ZETA (bond between O3* and next-P)
					// M1 and M2
					iter->pos_atom = k1+8; // O3* index
					mta = ( iter->get_atom() )->getPdbocc(); // N-atom mass
					// M1 += O3*
					mb1 += mta; // Left_half = previous + current (O3*)
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					// Updating previous masses*positions
					for( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
					// Inertia Matrix 1
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = r[0] = coord[(k1+8)*3];
						y[1] = r[1] = coord[(k1+8)*3 + 1];
						y[2] = r[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

//						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0;
//						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
//						{
//							mta = ( iter->get_atom() )->getPdbocc(); // mass
//							r[0] = coord[iter->pos_atom*3];
//							r[1] = coord[iter->pos_atom*3+1];
//							r[2] = coord[iter->pos_atom*3+2];
//							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//							for ( m = 0; m < 3; m++ )
//							{
//								for ( n = 0; n < 3; n++ )
//								{
//									mfa = -r[m] * r[n];
//									if ( m == n ) mfa += nr2;
//									I1[m][n] +=  mta * mfa;
//								}
//							}
//						}
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k2+1; // k <= k2 --> body 1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = -1;   // No "!=" condition
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k <= k2 )
//								body1[j2][k] = true;	// body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}
					j++; // der[] index
				}
				else
				{
					// Computing M1/M2 and "my" due to EPSILON and ZETA bodies
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
					{	// M1 += Sugar + Base
						mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
						mb1 += mta; // Left_half = previous + actual (NH)
						// Updating previous masses*positions
						for( m = 0; m < 3; m++ )
						{
							r[m] = coord[iter->pos_atom*3 + m];
							my[m] +=  mta * r[m];
						}
						// Inertia Matrix 1
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
					}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
				}
			} // END RNA
			// if DNA fragment
			else if( resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// ********************************************************************************************
				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)

				// ALPHA (bond between P and O5*)
				for ( iter->pos_atom = k1; iter->pos_atom <= k1+2; iter->next_atom() )
				{	// M1 += P + O1P + O2P
					// M1 and M2
					mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
					mb1 += mta; // Left_half = previous + actual (NH)

					// Updating previous masses*positions
					for ( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
					// Inertia Matrix 1
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}
				}
				mb2 = mtot - mb1; // Right-half = Tot - Left_half

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					y[0] = r[0] = coord[k1*3];
					y[1] = r[1] = coord[k1*3 + 1];
					y[2] = r[2] = coord[k1*3 + 2];
					e[0] -= y[0]; // P --> O5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+2; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+3; // k <= k1+2 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+2 )
//							body1[j2][k] = true;    // body 1  (O2P included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				// BETA (bond between O5* and C5*)
				iter->pos_atom = k1+3; // O5* index
				// M1 += O5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3];
				r[1] = coord[iter->pos_atom*3+1];
				r[2] = coord[iter->pos_atom*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					y[0] = r[0] = coord[(k1+3)*3];
					y[1] = r[1] = coord[(k1+3)*3 + 1];
					y[2] = r[2] = coord[(k1+3)*3 + 2];
					e[0] -= y[0]; // O5* --> C5*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+3; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+4; // k <= k1+3 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+3 )
//							body1[j2][k] = true;	// body 1 (O5* included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				// GAMMA (bond between C5* and C4*)
				iter->pos_atom = k1+4; // C5* index
				// M1 += C5*
				// M1 and M2
				mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
				mb1 += mta; // Left_half = previous + actual (NH)
				mb2 = mtot - mb1; // Right-half = Tot - Left_half
				// Updating previous masses*positions
				for( m = 0; m < 3; m++ )
				{
					r[m] = coord[iter->pos_atom*3 + m];
					my[m] +=  mta * r[m];
				}
				// Inertia Matrix 1
				r[0] = coord[iter->pos_atom*3];
				r[1] = coord[iter->pos_atom*3+1];
				r[2] = coord[iter->pos_atom*3+2];
				nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
				for ( m = 0; m < 3; m++ )
					for ( n = 0; n < 3; n++ )
					{
						mfa = -r[m] * r[n];
						if ( m == n ) mfa += nr2;
						I1[m][n] +=  mta * mfa;
					}

				if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
				{
					// e_lambda
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					y[0] = r[0] = coord[(k1+4)*3];
					y[1] = r[1] = coord[(k1+4)*3 + 1];
					y[2] = r[2] = coord[(k1+4)*3 + 2];
					e[0] -= y[0]; // C5* --> C4*
					e[1] -= y[1];
					e[2] -= y[2];
					temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
					for ( m = 0; m < 3; m++ )
					{
						e[m] /= temp; // Unit vector normalization
						Y1[m] = my[m] / mb1; // updating body 1 CoM
					}

					y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
					y[1] = e[2] * r[0] - e[0] * r[2];
					y[2] = e[0] * r[1] - e[1] * r[0];

//					// Inertia Matrix 1
//					for(m=0;m<3;m++)
//						for(n=0;n<3;n++)
//							I1[m][n] = 0.0;
//					for ( iter->pos_atom = 0; iter->pos_atom <= k1+4; iter->next_atom() )
//					{
//						mta = ( iter->get_atom() )->getPdbocc(); // mass
//						r[0] = coord[iter->pos_atom*3];
//						r[1] = coord[iter->pos_atom*3+1];
//						r[2] = coord[iter->pos_atom*3+2];
//						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//						for ( m = 0; m < 3; m++ )
//						{
//							for ( n = 0; n < 3; n++ )
//							{
//								mfa = -r[m] * r[n];
//								if ( m == n ) mfa += nr2;
//								I1[m][n] +=  mta * mfa;
//							}
//						}
//					}
					// now matrix I2
					for( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
							I2[m] [n] = I[m] [n] - I1[m] [n];

					// Computing some coefficients
					calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

					// loading "v" and "w" into V/W-arrays
					V[j2][0][0] = v[0][0];
					V[j2][0][1] = v[0][1];
					V[j2][0][2] = v[0][2];
					V[j2][1][0] = v[1][0];
					V[j2][1][1] = v[1][1];
					V[j2][1][2] = v[1][2];
					W[j2][0][0] = w[0][0];
					W[j2][0][1] = w[0][1];
					W[j2][0][2] = w[0][2];
					W[j2][1][0] = w[1][0];
					W[j2][1][1] = w[1][1];
					W[j2][1][2] = w[1][2];

					body1[j2][0] = k1+5; // k <= k1+4 --> body 1
					body1[j2][1] = -1;   // No ">" condition
					body1[j2][2] = -1;   // No "!=" condition
//					for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//						if( k <= k1+4 )
//							body1[j2][k] = true;    // body 1 (C5* included)
//						else
//							body1[j2][k] = false;    // body 2

					j2++; // selected dihedral index
				}
				j++; // der[] index

				//  LATERAL CHAIN --> CHI
				if(type == 2)
				{
					if( fix == NULL || fix[ j ] ) // Check whether the variable it's fixed
					{
						// Screening Nucleotide-Base atoms
						// WARNING when switching to DNA (12 --- ���)
						Y1[0] = 0.0;
						Y1[1] = 0.0;
						Y1[2] = 0.0;
						// The CoM of the whole protein is = 0.0 !!!
						// We'll substract to the CoM (total) the contribution
						// of the 2nd body (CB-R). (See below)

						// M1 and M2
						mb2 = 0.0;
						for( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() ) // DNA
						{
							mta = ( iter->get_atom() )->getPdbocc();
							mb2 += mta; // atom mass
							// Substracting to the CoM (total) the body 2 (CB-R) contribution
							Y1[0] -= mta * coord[iter->pos_atom*3];
							Y1[1] -= mta * coord[iter->pos_atom*3+1];
							Y1[2] -= mta * coord[iter->pos_atom*3+2];
						}
						mpr = mb1; // body 1 mass buffer
						mb1 = mtot - mb2;
						// Computing CoM
						for(m=0;m<3;m++)
							Y1[m] /= mb1;

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
						r[0] = y[0] = coord[(k1+10)*3];
						r[1] = y[1] = coord[(k1+10)*3 + 1];
						r[2] = y[2] = coord[(k1+10)*3 + 2];
						e[0] -= y[0]; // C1* --> N1/N9 unit vector
						e[1] -= y[1];
						e[2] -= y[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

						// Calculate matrix I2
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								I2[m][n] = 0;
						// Screening Nucleotide-Base atoms
						for ( iter->pos_atom = k1+11; iter->pos_atom <= k2; iter->next_atom() )
						{
							mta = ( iter->get_atom() )->getPdbocc();
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for ( m = 0; m < 3; m++ )
								for ( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I2[m] [n] += mta * mfa;
								}
						}
						// Now matrix I1
						for(m=0;m<3;m++)
							for(n=0;n<3;n++)
								Ix[m][n] = I[m][n]-I2[m][n];

						calcoefM( mb1, mb2, e, y, Y1, Ix, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k1+12; // k < k1+11 --> body 1
						body1[j2][1] = k2;    // k > k2 --> body 1
						body1[j2][2] = -1;    // No "!=" condition
//						for( k = 0; k < num_atoms; k++ )
//							if ( k < k1+11 || k > k2 ) // if body 1 (DNA)
//								body1[j2][k] = true;	// body 1 + (all except Nucleotide-Base)
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
						mb1 = mpr; // restoring previous mb1's value!
					}
					j++; // der[] index
				}

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{

					// EPSILON (bond between C3* and O3*)
					// M1 and M2
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
						if( iter->pos_atom != k1+8 ) // if not O3*
						{	// M1 += Sugar + Base
							mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
							mb1 += mta; // Left_half = previous + actual (NH)
							// Updating previous masses*positions
							for( m = 0; m < 3; m++ )
							{
								r[m] = coord[iter->pos_atom*3 + m];
								my[m] +=  mta * r[m];
							}
							// Inertia Matrix 1
							r[0] = coord[iter->pos_atom*3];
							r[1] = coord[iter->pos_atom*3+1];
							r[2] = coord[iter->pos_atom*3+2];
							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
							for ( m = 0; m < 3; m++ )
								for ( n = 0; n < 3; n++ )
								{
									mfa = -r[m] * r[n];
									if ( m == n ) mfa += nr2;
									I1[m][n] +=  mta * mfa;
								}
						}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						y[0] = r[0] = coord[(k1+7)*3];
						y[1] = r[1] = coord[(k1+7)*3 + 1];
						y[2] = r[2] = coord[(k1+7)*3 + 2];
						e[0] -= y[0]; // C3* --> O3*
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

//						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0;
//						for( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
//							if( iter->pos_atom != k1+8 ) // if not O3*
//							{
//								mta = ( iter->get_atom() )->getPdbocc(); // mass
//								r[0] = coord[iter->pos_atom*3];
//								r[1] = coord[iter->pos_atom*3+1];
//								r[2] = coord[iter->pos_atom*3+2];
//								nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//								for ( m = 0; m < 3; m++ )
//								{
//									for ( n = 0; n < 3; n++ )
//									{
//										mfa = -r[m] * r[n];
//										if ( m == n ) mfa += nr2;
//										I1[m][n] +=  mta * mfa;
//									}
//								}
//							}
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k2+1; // k <= k2 --> body 1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = k1+8; // k != k1+8 --> body 1
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k <= k2 && k != k1+8 ) // not O3*
//								body1[j2][k] = true;	// body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}
					j++; // der[] index

					// ZETA (bond between O3* and next-P)
					// M1 and M2
					iter->pos_atom = k1+8; // O3* index
					mta = ( iter->get_atom() )->getPdbocc(); // N-atom mass
					// M1 += O3*
					mb1 += mta; // Left_half = previous + current (O3*)
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
					// Updating previous masses*positions
					for( m = 0; m < 3; m++ )
					{
						r[m] = coord[iter->pos_atom*3 + m];
						my[m] +=  mta * r[m];
					}
					// Inertia Matrix 1
					r[0] = coord[iter->pos_atom*3];
					r[1] = coord[iter->pos_atom*3+1];
					r[2] = coord[iter->pos_atom*3+2];
					nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
					for ( m = 0; m < 3; m++ )
						for ( n = 0; n < 3; n++ )
						{
							mfa = -r[m] * r[n];
							if ( m == n ) mfa += nr2;
							I1[m][n] +=  mta * mfa;
						}

					if( fix == NULL || fix[ j ] ) // Check whether current variable it's fixed
					{
						// e_lambda
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						y[0] = r[0] = coord[(k1+8)*3];
						y[1] = r[1] = coord[(k1+8)*3 + 1];
						y[2] = r[2] = coord[(k1+8)*3 + 2];
						e[0] -= y[0]; // O3* --> next-P
						e[1] -= y[1];
						e[2] -= y[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
						{
							e[m] /= temp; // Unit vector normalization
							Y1[m] = my[m] / mb1; // updating body 1 CoM
						}

						y[0] = e[1] * r[2] - e[2] * r[1]; // Psi
						y[1] = e[2] * r[0] - e[0] * r[2];
						y[2] = e[0] * r[1] - e[1] * r[0];

//						// Inertia Matrix 1
//						for(m=0;m<3;m++)
//							for(n=0;n<3;n++)
//								I1[m][n] = 0.0;
//						for ( iter->pos_atom = 0; iter->pos_atom <= k2; iter->next_atom() )
//						{
//							mta = ( iter->get_atom() )->getPdbocc(); // mass
//							r[0] = coord[iter->pos_atom*3];
//							r[1] = coord[iter->pos_atom*3+1];
//							r[2] = coord[iter->pos_atom*3+2];
//							nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
//							for ( m = 0; m < 3; m++ )
//							{
//								for ( n = 0; n < 3; n++ )
//								{
//									mfa = -r[m] * r[n];
//									if ( m == n ) mfa += nr2;
//									I1[m][n] +=  mta * mfa;
//								}
//							}
//						}
						// now matrix I2
						for( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
								I2[m] [n] = I[m] [n] - I1[m] [n];

						// Computing some coefficients
						calcoefM( mb1, mb2, e, y, Y1, I1, I2, J, v, w ); // (ec. 40 Braun et al.)

						// loading "v" and "w" into V/W-arrays
						V[j2][0][0] = v[0][0];
						V[j2][0][1] = v[0][1];
						V[j2][0][2] = v[0][2];
						V[j2][1][0] = v[1][0];
						V[j2][1][1] = v[1][1];
						V[j2][1][2] = v[1][2];
						W[j2][0][0] = w[0][0];
						W[j2][0][1] = w[0][1];
						W[j2][0][2] = w[0][2];
						W[j2][1][0] = w[1][0];
						W[j2][1][1] = w[1][1];
						W[j2][1][2] = w[1][2];

						body1[j2][0] = k2+1; // k <= k2 --> body 1
						body1[j2][1] = -1;   // No ">" condition
						body1[j2][2] = -1;   // No "!=" condition
//						for( k = 0; k < num_atoms; k++ ) // screening ALL atoms
//							if( k <= k2 )
//								body1[j2][k] = true;	// body 1
//							else
//								body1[j2][k] = false;    // body 2

						j2++; // selected dihedral index
					}
					j++; // der[] index
				}
				else
				{
					// Computing M1/M2 and "my" due to EPSILON and ZETA bodies
					for( iter->pos_atom = k1+5; iter->pos_atom <= k2; iter->next_atom() )
					{	// M1 += Sugar + Base
						mta= ( iter->get_atom() )->getPdbocc(); // N-atom mass
						mb1 += mta; // Left_half = previous + actual (NH)
						// Updating previous masses*positions
						for( m = 0; m < 3; m++ )
						{
							r[m] = coord[iter->pos_atom*3 + m];
							my[m] +=  mta * r[m];
						}
						// Inertia Matrix 1
						r[0] = coord[iter->pos_atom*3];
						r[1] = coord[iter->pos_atom*3+1];
						r[2] = coord[iter->pos_atom*3+2];
						nr2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
						for ( m = 0; m < 3; m++ )
							for ( n = 0; n < 3; n++ )
							{
								mfa = -r[m] * r[n];
								if ( m == n ) mfa += nr2;
								I1[m][n] +=  mta * mfa;
							}
					}
					mb2 = mtot - mb1; // Right-half = Tot - Left_half
				}
			} // END DNA
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype == tmol_smol )
			{
				if(debug2)
					printf("SMOL  k1= %d  k2= %d  mb1= %f  mb2= %f  my= %f %f %f  Y1= %f %f %f\n",k1,k2,mb1,mb2,my[0],my[1],my[2],Y1[0],Y1[1],Y1[2]);

				for( iter->pos_atom = k1; iter->pos_atom <= k2; iter->next_atom() )
				{
					mta = iter->get_atom()->getPdbocc(); // Non-(C/O) atoms mass
					mb1 += mta;
					my[0] += mta * coord[iter->pos_atom*3];
					my[1] += mta * coord[iter->pos_atom*3+1];
					my[2] += mta * coord[iter->pos_atom*3+2];
				}
			}
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(vwMFAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			if( !( (iter_frag->pos_fragment == num_res-1) && (iter_seg->pos_segment == num_seg - 1) ) )
			{
				k1 += props[res_index].nat;
				k2 += props[res_index + 1].nat;
			}
			res_index++;
		}
		delete iter_frag;
	}
	// end derivatives
	iter_seg->clean_virtual();
	delete iter_seg;
	delete iter;
}

// Translate from IC-modes to Mass-weighted/unweighted Cartesian modes using K-matrix
// IC-modes ("hess_matix") --> CCS-modes (*mass = NULL)
//                         --> Mass-Weighted CCS-modes, if *mass = array with SQRT(mass_i)
// Warning: "evec" memory should be already allocated!
void di2cart(floating *hess_matrix,trd *der,floating *evec,int size,int num_atoms,int modes_saved, double *mass)
{
	for(int i=0; i< num_atoms*3*modes_saved; i++)
		evec[i] = 0.0; // initialization

	// Transforming d(dihedral) into Cartesian coordinates (dx,dy,dz) with K-matrix
	// Eigenvalues not translated to Cartesian. How we deal with this????
	// Both, CCS and DAS eigenvalues are the same--> See Go's paper about NMA theorem CCS & DAS.
	int der_index,evec_index,hess_index;

	if(mass==NULL) // Cartesian coordinates
		for(int n=0; n<modes_saved; n++) // screening only modes to save
		{
			evec_index = n*3*num_atoms;
			hess_index = n*size;
			for(int k=0; k<num_atoms; k++) // screen pseudo-atoms
			{
				der_index = k*size;
				for(int l=0; l<size; l++) // screen dihedrals
				{   // "Each atom accumulates the perturbation caused by the Dihedral NM"
					evec[evec_index] += der[der_index + l].x * hess_matrix[hess_index + l]; // x cartesian component
					evec[evec_index + 1] += der[der_index + l].y * hess_matrix[hess_index + l]; // y cartesian component
					evec[evec_index + 2] += der[der_index + l].z * hess_matrix[hess_index + l]; // z cartesian component
				}
				evec_index += 3;
			}
		}
	else // Mass-Weighted cartesian coordinates
		for(int n=0; n<modes_saved; n++) // screening only modes to save
		{
			evec_index = n*3*num_atoms;
			hess_index = n*size;
			for(int k=0; k<num_atoms; k++) // screen pseudo-atoms
			{
				der_index = k*size;
				for(int l=0; l<size; l++) // screen dihedrals
				{   // "Each atom accumulates the perturbation caused by the Dihedral NM"
					evec[evec_index]     += mass[k] * der[der_index + l].x * hess_matrix[hess_index + l]; // x cartesian component
					evec[evec_index + 1] += mass[k] * der[der_index + l].y * hess_matrix[hess_index + l]; // y cartesian component
					evec[evec_index + 2] += mass[k] * der[der_index + l].z * hess_matrix[hess_index + l]; // z cartesian component
				}
				evec_index += 3;
			}
		}
}


// Translate from IC-modes to Mass-weighted/unweighted Cartesian modes using V/W-arrays
// IC-modes ("hess_matix") --> CCS-modes
// Warning: "evec" memory should be already allocated!
void di2cartVW(floating *hess_matrix,float *coord,double ***V,double ***W,int **body1,floating *evec,int size,int num_atoms,int modes_saved, int model)
{
	float rk[3];
	double der[3];

	if(evec != NULL )
		for(int i=0; i< num_atoms*3*modes_saved; i++)
			evec[i] = 0.0; // initialization

	// Is not neccessary to translate the eigenvalues from ICS to CCS,
	// see Go's theorem paper elsewere...
	int evec_index,hess_index,atom_index,k,i,n;

	// fprintf(stderr,"Begin di2cartVW\n");

	if( model == 1 || model == 2 ) // C5/HA models
		for(k=0,atom_index=0; k<num_atoms; k++,atom_index+=3) // screen pseudo-atoms
		{
			rk[0] = coord[k*3]; // k-atom position
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			// Compute "derivative" for current atom and IC
			for(i=0; i<size; i++) // screen ICs
			{
				if( (k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) // if body 1 atom
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					der[0] = V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1];
					der[1] = V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2];
					der[2] = V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0];
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					der[0] = V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1];
					der[1] = V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2];
					der[2] = V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0];
				}
				// this adds each IC contribution to every mode...
				for(n = 0, evec_index = atom_index, hess_index = 0; n<modes_saved; n++, evec_index += 3*num_atoms, hess_index += size) // screening only the modes to save
				{
					// Each atom accumulates the perturbation caused by each IC
					evec[evec_index]     += der[0] * hess_matrix[hess_index + i]; // x cartesian component
					evec[evec_index + 1] += der[1] * hess_matrix[hess_index + i]; // y cartesian component
					evec[evec_index + 2] += der[2] * hess_matrix[hess_index + i]; // z cartesian component
				}
			}
		}
	else // CA/CA3 models
	{
		// fprintf(stderr,"CA/CA3 model\n");
		for(k=0,atom_index=0; k<num_atoms; k++,atom_index+=3) // screen pseudo-atoms
		{
			rk[0] = coord[k*3]; // k-atom position
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			// Compute "derivative" for current atom and IC
			for(i=0; i<size; i++) // screen ICs
			{
				if(k < body1[i][0]) // if body 1 atom
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					der[0] = V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1];
					der[1] = V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2];
					der[2] = V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0];
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					der[0] = V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1];
					der[1] = V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2];
					der[2] = V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0];
				}
				// this adds each IC contribution to every mode...
				for(n = 0, evec_index = atom_index, hess_index = 0; n<modes_saved; n++, evec_index += 3*num_atoms, hess_index += size) // screening only the modes to save
				{
					// Each atom accumulates the perturbation caused by each IC
					evec[evec_index]     += der[0] * hess_matrix[hess_index + i]; // x cartesian component
					evec[evec_index + 1] += der[1] * hess_matrix[hess_index + i]; // y cartesian component
					evec[evec_index + 2] += der[2] * hess_matrix[hess_index + i]; // z cartesian component
				}
			}
		}
	}
}

// Multi-thread routine to convert CC modes into IC modes in parallel
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void *di2cartVW_thread(void *threadarg)
{
	float rk[3];
	double der[3];
	// Is not neccessary to translate the eigenvalues from ICS to CCS,
	// see Go's theorem paper elsewere...
	int evec_index,hess_index,atom_index,k,i,n;

	double *evec = ((di2cartVW_data *) threadarg)->CCevec;
	double *hess_matrix = ((di2cartVW_data *) threadarg)->ICevec;
	float *coord = ((di2cartVW_data *) threadarg)->coord;
	int nmodes = ((di2cartVW_data *) threadarg)->nmodes;
	int num_atoms = ((di2cartVW_data *) threadarg)->natoms;
	int model = ((di2cartVW_data *) threadarg)->model;
	int size = ((di2cartVW_data *) threadarg)->size;
	double ***V = ((di2cartVW_data *) threadarg)->V;
	double ***W = ((di2cartVW_data *) threadarg)->W;
	int **body1 = ((di2cartVW_data *) threadarg)->body1;
	int first = ((di2cartVW_data *) threadarg)->first;
	int last = ((di2cartVW_data *) threadarg)->last;

	if( model == 1 || model == 2 ) // C5/HA models
		for(k=first,atom_index=first*3; k<=last; k++,atom_index+=3) // screen pseudo-atoms
		{
			rk[0] = coord[k*3]; // k-atom position
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			// Compute "derivative" for current atom and IC
			for(i=0; i<size; i++) // screen ICs
			{
				if( (k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) // if body 1 atom
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					der[0] = V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1];
					der[1] = V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2];
					der[2] = V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0];
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					der[0] = V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1];
					der[1] = V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2];
					der[2] = V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0];
				}
				// this adds each IC contribution to every mode...
				for(n = 0, evec_index = atom_index, hess_index = 0; n<nmodes; n++, evec_index += 3*num_atoms, hess_index += size) // screening only the modes to save
				{
					// Each atom accumulates the perturbation caused by each IC
					evec[evec_index]     += der[0] * hess_matrix[hess_index + i]; // x cartesian component
					evec[evec_index + 1] += der[1] * hess_matrix[hess_index + i]; // y cartesian component
					evec[evec_index + 2] += der[2] * hess_matrix[hess_index + i]; // z cartesian component
				}
			}
		}
	else // CA/CA3 models
		for(k=first,atom_index=first*3; k<=last; k++,atom_index+=3) // screen pseudo-atoms
		{
			rk[0] = coord[k*3]; // k-atom position
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			// Compute "derivative" for current atom and IC
			for(i=0; i<size; i++) // screen ICs
			{
				if(k < body1[i][0]) // if body 1 atom
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					der[0] = V[i][0][0] + W[i][0][1] * rk[2] - W[i][0][2] * rk[1];
					der[1] = V[i][0][1] + W[i][0][2] * rk[0] - W[i][0][0] * rk[2];
					der[2] = V[i][0][2] + W[i][0][0] * rk[1] - W[i][0][1] * rk[0];
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					der[0] = V[i][1][0] + W[i][1][1] * rk[2] - W[i][1][2] * rk[1];
					der[1] = V[i][1][1] + W[i][1][2] * rk[0] - W[i][1][0] * rk[2];
					der[2] = V[i][1][2] + W[i][1][0] * rk[1] - W[i][1][1] * rk[0];
				}
				// this adds each IC contribution to every mode...
				for(n = 0, evec_index = atom_index, hess_index = 0; n<nmodes; n++, evec_index += 3*num_atoms, hess_index += size) // screening only the modes to save
				{
					// Each atom accumulates the perturbation caused by each IC
					evec[evec_index]     += der[0] * hess_matrix[hess_index + i]; // x cartesian component
					evec[evec_index + 1] += der[1] * hess_matrix[hess_index + i]; // y cartesian component
					evec[evec_index + 2] += der[2] * hess_matrix[hess_index + i]; // z cartesian component
				}
			}
		}
	pthread_exit(NULL);
}

// Parallel routine to convert CC modes into IC modes
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void di2cartVW_par(floating *hess_matrix,float *coord,double ***V,double ***W,int **body1,floating *evec,int size,int num_atoms,int modes_saved, int model, int nthreads)
{
	bool debug = false;
	di2cartVW_data *threads_data;
	pthread_t *threads;
	float rk[3];
	double der[3];
	// Is not neccessary to translate the eigenvalues from ICS to CCS,
	// see Go's theorem paper elsewere...
	int evec_index,hess_index,atom_index,k,i,n,rc,current,init=0;
	void *status;

	if(evec != NULL )
		for(int i=0; i< num_atoms*3*modes_saved; i++)
			evec[i] = 0.0; // initialization

	// Allocating threads data
	if( !(threads_data = (di2cartVW_data *) malloc(sizeof(di2cartVW_data) * nthreads)) )
	{
		fprintf(stderr,"Msg(move_vwMFAx_par): I'm sorry, thread memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(threads = (pthread_t *) malloc(sizeof(pthread_t) * nthreads)) )
	{
		fprintf(stderr,"Msg(move_vwMFAx_par): I'm sorry, thread allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Creating threads input data
	current = num_atoms / nthreads; // balancing number of atoms to be moved equally (watch out odd numbers!)
	if(debug)
		fprintf(stderr,"Creating thread data:  num_atoms= %d  current= %d atoms/thread\n",num_atoms, current);


	for(i = 0; i < nthreads; i++)
	{
		threads_data[i].CCevec = evec;
		threads_data[i].ICevec = hess_matrix;
		threads_data[i].coord = coord;
		threads_data[i].nmodes = modes_saved;
		threads_data[i].natoms = num_atoms;
		threads_data[i].model = model;
		threads_data[i].size = size;
		threads_data[i].V = V;
		threads_data[i].W = W;
		threads_data[i].body1 = body1;
		threads_data[i].first = init;
		threads_data[i].last = init + current - 1;
		init += current;
		if(debug)
			fprintf(stderr,"Thread %d:  first= %d  last= %d\n",i,threads_data[i].first,threads_data[i].last);
	}
	threads_data[nthreads-1].last = num_atoms - 1; // this should ensure that all atoms will be moved!
	if(debug)
		fprintf(stderr,"Thread %d:  first= %d  last= %d\n",nthreads-1,threads_data[nthreads-1].first,threads_data[nthreads-1].last);

	// Creating threads
	for(i = 0; i < nthreads; i++)
	{
		if(debug)
			fprintf(stderr,"Creating thread (move_vwMCAx_thread): %d\n", i);
		rc = pthread_create(&threads[i], NULL, di2cartVW_thread, (void *) &threads_data[i]);
		if (rc)
		{
			fprintf(stderr,"ERROR; return code from pthread_create() is %d\n", rc);
			exit(-1);
		}
	}

	// This pauses the main thread until all sub-threads finish
	for(i = 0; i < nthreads; i++)
	{
		if(debug)
			fprintf(stderr,"Joining thread: %d\n", i);
		rc = pthread_join(threads[i], &status);
		if (rc)
		{
			fprintf(stderr,"ERROR; return code from pthread_join() is %d\n", rc);
			exit(-1);
		}
	}

	free(threads_data);
	free(threads);
}


// Translates from Cartesian normal modes to Mass-Weighted ones
// (mass <-> mass square root)
void cart2wcart(floating *evec,int num_atoms,int modes_saved,double *mass)
{
	int evec_index,atom_index,k,n;
	for(k=0,atom_index=0; k<num_atoms; k++,atom_index+=3)
	{
		for(n=0,evec_index=atom_index; n<modes_saved; n++,evec_index+=3*num_atoms)
		{
			evec[evec_index]     *= mass[k];
			evec[evec_index + 1] *= mass[k];
			evec[evec_index + 2] *= mass[k];
		}
	}
}

