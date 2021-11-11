/************************************************************************
*                     LIBRARY: libnma_move                              *
*************************************************************************
* Program is part of the ADP package URL: http://sbg.cib.csic.es        *
* (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
* Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
*************************************************************************
*                                                                       *
*   NMA motion in Internal Coordinate Space (ICS) library using         *
*   different Coarse-Graining models.                                   *
*   (It takes into account Muliple-Chains and different CG-models)      *
*                                                                       *
*************************************************************************
* This library is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation; either version 2 of the License, or     *
* (at your option) any later version.                                   *
************************************************************************/
#include <pthread.h>
#include <libnma_move.h>
#include <libnma_matrix.h>
#define DEVEL

//// GOOD REFERENCE OK
//// Moves atoms given an Internal Coordinates Normal Mode and a Step
//// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
//// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
//// If fix=NULL, no fixed dihedrals (equal to previous-version)
//void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix)
//{
//	bool debug = false;
//	int j=0;
//	int j2=0;
//	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
//	Residue *res;
//	Atom *atom;
//	pdbIter *iter,*iter_res_atom,*iter_seg;
//	double **Ri; // current rotation (Ri-matrix);
//	double **Si; // accummulated rotation (Si-matrix)
//	double **dummy; // dummy rotation (dummy-matrix)
//	double **temp;
//	int num_res = 0;
//	int index_res=0;
//	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
//	int indexbase;
//	int resn;
//
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	int num_seg = iter_seg->num_segment();
//
//	// Centering Macromolecule (Temporal)
//	// Computing the PDB's Center of Mass (CoM)
//	double mtot,mta;
//	double com[3];
//	pdbIter *iter_com;
//	Tcoor pos;
//
//	if(num_seg > 1)
//	{
//		iter_com = new pdbIter(mol);
//		com[0] = com[1] = com[2] = 0.0;
//		mtot = 0.0;
//		for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
//		{
//			atom = ( iter_com->get_atom() );
//			mta = atom->getPdbocc(); // Load mass...
//			mtot += mta;
//			atom->getPosition(pos);
//			/* Sum(mass*coord) before putting the CoM at 0 */
//			com[0] += mta * pos[0];
//			com[1] += mta * pos[1];
//			com[2] += mta * pos[2];
//		}
//		com[0] /= mtot;
//		com[1] /= mtot;
//		com[2] /= mtot;
//		if(debug)
//			printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );
//		delete iter_com;
//	}
//
//	// Ri-matrix, Si-matrix, & dummy-matrix initialization
//	Ri = (double **) malloc( sizeof(double *) * 3 );
//	Si = (double **) malloc( sizeof(double *) * 3 );
//	dummy = (double **) malloc( sizeof(double *) * 3 );
//	// 1st Si --> I-matrix
//	// 1st Ti --> (0,0,0) vector
//	for(int i=0; i<3; i++)
//	{
//		Ri[i] = (double *) malloc( sizeof(double) * 3 );
//		Si[i] = (double *) malloc( sizeof(double) * 3 );
//		dummy[i] = (double *) malloc( sizeof(double) * 3 );
//		Ti[i] = 0.0;
//		for(int j=0; j<3; j++)
//			if(i==j)
//			{
//				Si[i][j] = 1.0;
//				dummy[i][j] = 1.0;
//			}
//			else
//			{
//				Si[i][j] = 0.0;
//				dummy[i][j] = 0.0;
//			}
//	}
//
//	if(debug)
//		printf("Msg(move_dihedralM): Moving dihedrals:\n");
//
//	Segment * seg;
//	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter = new pdbIter( seg );
//		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)
//
//		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			T[0] = Ti[0];
//			T[1] = Ti[1];
//			T[2] = Ti[2];
//
//			// 3 Translations
//			for(int axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					Ti2[axis] = uu[j2] * step;
//					if(debug)
//						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//					j2++;
//				}
//				j++;
//			}
//
//			// 3 Rotations
//			for(int axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					e1[0] = 0;
//					e1[1] = 0;
//					e1[2] = 0;
//					e1[axis] = 1.0; // e1 ==> Rotation axis
//
//					// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
//					rotmat(uu[j2] * step, e1, Ri);
//
//					if(debug)
//						printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//					// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//					mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//					temp = Si; // swaping...
//					Si = dummy;
//					dummy = temp;
//
//					// This avoids pre-centering...
//					T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
//					T[1] -= com[1];
//					T[2] -= com[2];
//					// Rotating the "rotation origin position" -->
//					multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)
//
//					T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
//					T[1] = Ti3[1] + com[1];
//					T[2] = Ti3[2] + com[2];
//
//					j2++;
//				}
//				j++;
//			}
//
//			Ti[0] = T[0] + Ti2[0];
//			Ti[1] = T[1] + Ti2[1];
//			Ti[2] = T[2] + Ti2[2];
//		}
//
//		if(debug)
//			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);
//
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			res = ( Residue * ) iter->get_fragment();
//			iter_res_atom = new pdbIter( res ); // iter residue atoms
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType(); // gets current fragment moltype
//
//			if(fragtype == tmol_protein)
//			{
//				// We should move now "non-first segment" first residue
//				// (first segment first residue does not move!)
//				if(iter_seg->pos_segment != 0 && iter->pos_fragment==0) // non-first segment && first residue
//				{
//					// rotating current residue atoms (due to backwards internal coords)
//					for ( iter_res_atom->pos_atom = 0;
//						!iter_res_atom->gend_atom();
//						iter_res_atom->next_atom()   )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//				}
//
//				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
//				if( !(iter->pos_fragment == 0 || iter->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					// rotating N and CA atoms (due to previous residues)
//					for ( iter_res_atom->pos_atom = 0; // N and CA
//						iter_res_atom->pos_atom < 2;
//						iter_res_atom->next_atom()   )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PHI
//						// Computing axis
//						iter_res_atom->pos_atom = 0; // N
//						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
//						e1[1] = ca[1] - nh[1];
//						e1[2] = ca[2] - nh[2];
//
//						// Computing PHI rotation matrix given "e1" and the current "phi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - ca[1];
//						T[2] = Ti[2] - ca[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += ca[0]; // translating it back to the Qi position
//						Ti[1] += ca[1];
//						Ti[2] += ca[2];
//
//						j2++;
//					}
//
//					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
//					iter_res_atom->pos_atom = 2; // CO moves only!
//					atom = (Atom *) iter_res_atom->get_atom();
//					atom->getPosition( tv ) ; // position before rotation
//					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//					atom->setPosition( tvp ) ; // position after rotation
//
//					j++;
//				}
//
//				// If it's PRO we must move now !!!
//				// (Watch out this: When switching to Full-atom model)
//				if( iter->pos_fragment == num_res-1 || (resn == PRO && iter->pos_fragment != 0 ) )
//				{
//					// rotating current residue atoms (due to backwards dihedrals)
//					for( iter_res_atom->pos_atom = 0;
//						!iter_res_atom->gend_atom();
//						iter_res_atom->next_atom() )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//				}
//
//				// PSI
//				// First and Last residue never have PSI in CA-model
//				if( iter->pos_fragment != num_res-1 && iter->pos_fragment != 0 )
//				{
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PSI
//						// Computing axis
//						iter_res_atom->pos_atom = 2; // C
//						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
//						e1[1] = co[1] - ca[1];
//						e1[2] = co[2] - ca[2];
//
//						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - co[1];
//						T[2] = Ti[2] - co[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += co[0]; // translating it back to the Qi position
//						Ti[1] += co[1];
//						Ti[2] += co[2];
//
//						j2++;
//					}
//					// During PSI-rotation, no current residue atoms are moved !!!
//					j++;
//				}
//			}
//			// END PROTEIN.....
//			else
//			{
//				printf("Msg(move_dihedralMCAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
//				exit(1);
//			}
//			delete iter_res_atom;
//			index_res++; // counts residues
//		}
//		delete iter;
//	}
//	delete iter_seg;
//
//	// Freeing memory
//	for(int i=0; i<3; i++)
//	{
//		free( Ri[i] );
//		free( Si[i] );
//		free( dummy[i] );
//	}
//	free( Ri );
//	free( Si );
//	free( dummy );
//}

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix)
{
	bool debug = false;
	int j=0;
	int j2=0;
	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
	Residue *res;
	Atom *atom;
	pdbIter *iter,*iter_res_atom,*iter_seg;
	double **Ri; // current rotation (Ri-matrix);
	double **Si; // accummulated rotation (Si-matrix)
	double **dummy; // dummy rotation (dummy-matrix)
	double **temp;
	int num_res = 0;
	int index_res=0;
	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
	int indexbase;
	int resn;

	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	int num_seg = iter_seg->num_segment();

	// Centering Macromolecule (Temporal)
	// Computing the PDB's Center of Mass (CoM)
	double mtot,mta;
	double com[3];
	com[0] = com[1] = com[2] = 0.0;
	pdbIter *iter_com;
	Tcoor pos;

	if(num_seg > 1)
	{
		iter_com = new pdbIter(mol);
		mtot = 0.0;
		for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
		{
			atom = ( iter_com->get_atom() );
			mta = atom->getPdbocc(); // Load mass...
			mtot += mta;
			atom->getPosition(pos);
			/* Sum(mass*coord) before putting the CoM at 0 */
			com[0] += mta * pos[0];
			com[1] += mta * pos[1];
			com[2] += mta * pos[2];
		}
		com[0] /= mtot;
		com[1] /= mtot;
		com[2] /= mtot;
		if(debug)
			printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );
		delete iter_com;
	}

	// Ri-matrix, Si-matrix, & dummy-matrix initialization
	Ri = (double **) malloc( sizeof(double *) * 3 );
	Si = (double **) malloc( sizeof(double *) * 3 );
	dummy = (double **) malloc( sizeof(double *) * 3 );
	// 1st Si --> I-matrix
	// 1st Ti --> (0,0,0) vector
	for(int i=0; i<3; i++)
	{
		Ri[i] = (double *) malloc( sizeof(double) * 3 );
		Si[i] = (double *) malloc( sizeof(double) * 3 );
		dummy[i] = (double *) malloc( sizeof(double) * 3 );
		Ti[i] = 0.0;
		for(int j=0; j<3; j++)
			if(i==j)
			{
				Si[i][j] = 1.0;
				dummy[i][j] = 1.0;
			}
			else
			{
				Si[i][j] = 0.0;
				dummy[i][j] = 0.0;
			}
	}

	if(debug)
		printf("Msg(move_dihedralM): Moving dihedrals:\n");

	Segment * seg;
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter = new pdbIter( seg );
		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			T[0] = Ti[0];
			T[1] = Ti[1];
			T[2] = Ti[2];

			// 3 Translations
			for(int axis=0; axis<3; axis++)
			{
				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
				{
					Ti2[axis] = uu[j2] * step;
					if(debug)
						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
					j2++;
				}
				j++;
			}

			// 3 Rotations
			for(int axis=0; axis<3; axis++)
			{
				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
				{
					e1[0] = 0;
					e1[1] = 0;
					e1[2] = 0;
					e1[axis] = 1.0; // e1 ==> Rotation axis

					// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
					rotmat(uu[j2] * step, e1, Ri);

					if(debug)
						printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

					// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
					mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
					temp = Si; // swaping...
					Si = dummy;
					dummy = temp;

					// This avoids pre-centering...
					T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
					T[1] -= com[1];
					T[2] -= com[2];
					// Rotating the "rotation origin position" -->
					multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)

					T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
					T[1] = Ti3[1] + com[1];
					T[2] = Ti3[2] + com[2];

					j2++;
				}
				j++;
			}

			Ti[0] = T[0] + Ti2[0];
			Ti[1] = T[1] + Ti2[1];
			Ti[2] = T[2] + Ti2[2];
		}

		if(debug)
			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);

		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		{
			res = ( Residue * ) iter->get_fragment();
			iter_res_atom = new pdbIter( res ); // iter residue atoms
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType(); // gets current fragment moltype

			if(fragtype == tmol_protein)
			{
				// We should move now "non-first segment" first residue
				// (first segment first residue does not move!)
				if(iter->pos_fragment==0 && iter_seg->pos_segment != 0) // non-first segment && first residue
				{
					// rotating current residue atoms (due to backwards internal coords)
					for ( iter_res_atom->pos_atom = 0;
						!iter_res_atom->gend_atom();
						iter_res_atom->next_atom()   )
					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
					}
				}

				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
				if( !(iter->pos_fragment == 0 || iter->pos_fragment == num_res-1 || resn == PRO) )
				{
					// rotating N and CA atoms (due to previous residues)
					for ( iter_res_atom->pos_atom = 0; // N and CA
						iter_res_atom->pos_atom < 2;
						iter_res_atom->next_atom()   )
					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
					}

					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PHI
						// Computing axis
						iter_res_atom->pos_atom = 0; // N
						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
						iter_res_atom->pos_atom = 1; // CA
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing PHI rotation matrix given "e1" and the current "phi" angle
						rotmat(uu[j2] * step, e1, Ri);

						if(debug)
							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];

						j2++;
					}

					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
					iter_res_atom->pos_atom = 2; // CO moves only!
//					atom = (Atom *) iter_res_atom->get_atom();
//					atom->getPosition( tv ) ; // position before rotation
//					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//					atom->setPosition( tvp ) ; // position after rotation
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.

					j++;
				}

				// If it's PRO we must move now !!!
				// (Watch out this: When switching to Full-atom model)
				if( iter->pos_fragment == num_res-1 || (resn == PRO && iter->pos_fragment != 0 ) )
				{
					// rotating current residue atoms (due to backwards dihedrals)
					for( iter_res_atom->pos_atom = 0;
						!iter_res_atom->gend_atom();
						iter_res_atom->next_atom() )
					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
					}
				}

				// PSI
				// First and Last residue never have PSI in CA-model
				if( iter->pos_fragment != num_res-1 && iter->pos_fragment != 0 )
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PSI
						// Computing axis
						iter_res_atom->pos_atom = 2; // C
						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
						iter_res_atom->pos_atom = 1; // CA
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
						e1[1] = co[1] - ca[1];
						e1[2] = co[2] - ca[2];

						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
						rotmat(uu[j2] * step, e1, Ri);

						if(debug)
							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - co[1];
						T[2] = Ti[2] - co[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += co[0]; // translating it back to the Qi position
						Ti[1] += co[1];
						Ti[2] += co[2];

						j2++;
					}
					// During PSI-rotation, no current residue atoms are moved !!!
					j++;
				}
			}
			// END PROTEIN.....
			else
			{
				printf("Msg(move_dihedralMCAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
				exit(1);
			}
			delete iter_res_atom;
			index_res++; // counts residues
		}
		delete iter;
	}
	delete iter_seg;

	// Freeing memory
	for(int i=0; i<3; i++)
	{
		free( Ri[i] );
		free( Si[i] );
		free( dummy[i] );
	}
	free( Ri );
	free( Si );
	free( dummy );
}

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
// "iters" array of residue-level (pdbIter *)'s
void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix, pdbIter **iters)
{
	bool debug = false;
	int j=0;
	int j2=0;
	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
	Residue *res;
	Atom *atom;
	pdbIter *iter,*iter_res_atom,*iter_seg;
	double **Ri; // current rotation (Ri-matrix);
	double **Si; // accummulated rotation (Si-matrix)
	double **dummy; // dummy rotation (dummy-matrix)
	double **temp;
	int num_res = 0;
	int index_res=0;
	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
	int indexbase;
	int resn;

	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	int num_seg = iter_seg->num_segment();

	// Centering Macromolecule (Temporal)
	// Computing the PDB's Center of Mass (CoM)
	double mtot,mta;
	double com[3];
	pdbIter *iter_com;
	Tcoor pos;

	if(num_seg > 1)
	{
		iter_com = new pdbIter(mol);
		com[0] = com[1] = com[2] = 0.0;
		mtot = 0.0;
		for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
		{
			atom = ( iter_com->get_atom() );
			mta = atom->getPdbocc(); // Load mass...
			mtot += mta;
			atom->getPosition(pos);
			/* Sum(mass*coord) before putting the CoM at 0 */
			com[0] += mta * pos[0];
			com[1] += mta * pos[1];
			com[2] += mta * pos[2];
		}
		com[0] /= mtot;
		com[1] /= mtot;
		com[2] /= mtot;
		if(debug)
			printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );
		delete iter_com;
	}

	// Ri-matrix, Si-matrix, & dummy-matrix initialization
	Ri = (double **) malloc( sizeof(double *) * 3 );
	Si = (double **) malloc( sizeof(double *) * 3 );
	dummy = (double **) malloc( sizeof(double *) * 3 );
	// 1st Si --> I-matrix
	// 1st Ti --> (0,0,0) vector
	for(int i=0; i<3; i++)
	{
		Ri[i] = (double *) malloc( sizeof(double) * 3 );
		Si[i] = (double *) malloc( sizeof(double) * 3 );
		dummy[i] = (double *) malloc( sizeof(double) * 3 );
		Ti[i] = 0.0;
		for(int j=0; j<3; j++)
			if(i==j)
			{
				Si[i][j] = 1.0;
				dummy[i][j] = 1.0;
			}
			else
			{
				Si[i][j] = 0.0;
				dummy[i][j] = 0.0;
			}
	}

	if(debug)
		printf("Msg(move_dihedralM): Moving dihedrals:\n");

	Segment * seg;
	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter = new pdbIter( seg );
//		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)
		num_res=524;

		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			T[0] = Ti[0];
			T[1] = Ti[1];
			T[2] = Ti[2];

			// 3 Translations
			for(int axis=0; axis<3; axis++)
			{
				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
				{
					Ti2[axis] = uu[j2] * step;
					if(debug)
						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
					j2++;
				}
				j++;
			}

			// 3 Rotations
			for(int axis=0; axis<3; axis++)
			{
				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
				{
					e1[0] = 0;
					e1[1] = 0;
					e1[2] = 0;
					e1[axis] = 1.0; // e1 ==> Rotation axis

					// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
					rotmat(uu[j2] * step, e1, Ri);

					if(debug)
						printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

					// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
					mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
					temp = Si; // swaping...
					Si = dummy;
					dummy = temp;

					// This avoids pre-centering...
					T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
					T[1] -= com[1];
					T[2] -= com[2];
					// Rotating the "rotation origin position" -->
					multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)

					T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
					T[1] = Ti3[1] + com[1];
					T[2] = Ti3[2] + com[2];

					j2++;
				}
				j++;
			}

			Ti[0] = T[0] + Ti2[0];
			Ti[1] = T[1] + Ti2[1];
			Ti[2] = T[2] + Ti2[2];
		}

		if(debug)
			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);

		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		{
			res = ( Residue * ) iter->get_fragment();
//			iter_res_atom = new pdbIter( res ); // iter residue atoms
			iter_res_atom = iters[ iter->pos_fragment ];
//			resn = resnum_from_resname( res->getName() );
			resn = 9999;
//			fragtype = res->getMolType(); // gets current fragment moltype

//			if(fragtype == tmol_protein)
			if(true)
			{
				// We should move now "non-first segment" first residue
				// (first segment first residue does not move!)
				if(iter->pos_fragment==0 && iter_seg->pos_segment != 0) // non-first segment && first residue
				{
					// rotating current residue atoms (due to backwards internal coords)
					iter_res_atom->pos_atom = 0;
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
//					iter_res_atom->pos_atom = 1;
					iter_res_atom->next_atom();
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
//					iter_res_atom->pos_atom = 2;
					iter_res_atom->next_atom();
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
				}

				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
				if( !(iter->pos_fragment == 0 || iter->pos_fragment == num_res-1 || resn == PRO) )
				{
					// rotating N and CA atoms (due to previous residues)
					iter_res_atom->pos_atom = 0;
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
//					iter_res_atom->pos_atom = 1;
					iter_res_atom->next_atom();
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.

					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PHI
						// Computing axis
						iter_res_atom->pos_atom = 0; // N
						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
//						iter_res_atom->pos_atom = 1; // CA
						iter_res_atom->next_atom();
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing PHI rotation matrix given "e1" and the current "phi" angle
						rotmat(uu[j2] * step, e1, Ri);

						if(debug)
							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];

						j2++;
					}

					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
					iter_res_atom->pos_atom = 2; // CO moves only!
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.

					j++;
				}

				// If it's PRO we must move now !!!
				// (Watch out this: When switching to Full-atom model)
				if( iter->pos_fragment == num_res-1 || (resn == PRO && iter->pos_fragment != 0 ) )
				{
					// rotating current residue atoms (due to backwards dihedrals)
					iter_res_atom->pos_atom = 0;
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
//					iter_res_atom->pos_atom = 1;
					iter_res_atom->next_atom();
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
//					iter_res_atom->pos_atom = 2;
					iter_res_atom->next_atom();
					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti,tvp); // direct rot.trans.
				}

				// PSI
				// First and Last residue never have PSI in CA-model
				if( iter->pos_fragment != num_res-1 && iter->pos_fragment != 0 )
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PSI
						// Computing axis
						iter_res_atom->pos_atom = 2; // C
						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
//						iter_res_atom->pos_atom = 1; // CA
						iter_res_atom->back_atom();
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
						e1[1] = co[1] - ca[1];
						e1[2] = co[2] - ca[2];

						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
						rotmat(uu[j2] * step, e1, Ri);

						if(debug)
							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - co[1];
						T[2] = Ti[2] - co[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += co[0]; // translating it back to the Qi position
						Ti[1] += co[1];
						Ti[2] += co[2];

						j2++;
					}
					// During PSI-rotation, no current residue atoms are moved !!!
					j++;
				}
			}
			// END PROTEIN.....
			else
			{
				printf("Msg(move_dihedralMCAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
				exit(1);
			}
//			delete iter_res_atom;
			index_res++; // counts residues
		}
		delete iter;
	}
	delete iter_seg;

	// Freeing memory
	for(int i=0; i<3; i++)
	{
		free( Ri[i] );
		free( Si[i] );
		free( dummy[i] );
	}
	free( Ri );
	free( Si );
	free( dummy );
}

//// Moves atoms given an Internal Coordinates Normal Mode and a Step
//// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
//// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
//// If fix=NULL, no fixed dihedrals (equal to previous-version)
//void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix)
//{
//	bool debug = false;
//	int i=0,j=0,j2=0;
//	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
//	Residue *res;
//	Atom *atom;
//	pdbIter *iter,*iter_res_atom,*iter_seg;
//	double **Ri; // current rotation (Ri-matrix);
//	double **Si; // accummulated rotation (Si-matrix)
//	double **dummy; // dummy rotation (dummy-matrix)
//	double **temp;
//	int num_res = 0;
//	int index_res=0;
//	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
//	int indexbase;
//	int resn;
//	int axis;
//
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	int num_seg = iter_seg->num_segment();
//
//	// Centering Macromolecule (Temporal)
//	// Computing the PDB's Center of Mass (CoM)
//	double mtot,mta;
//	double com[3];
//	pdbIter *iter_com;
//	Tcoor pos;
//
//	if(num_seg > 1)
//	{
//		iter_com = new pdbIter(mol);
//		com[0] = com[1] = com[2] = 0.0;
//		mtot = 0.0;
//		for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
//		{
//			atom = ( iter_com->get_atom() );
//			mta = atom->getPdbocc(); // Load mass...
//			mtot += mta;
//			atom->getPosition(pos);
//			/* Sum(mass*coord) before putting the CoM at 0 */
//			com[0] += mta * pos[0];
//			com[1] += mta * pos[1];
//			com[2] += mta * pos[2];
//		}
//		com[0] /= mtot;
//		com[1] /= mtot;
//		com[2] /= mtot;
//		if(debug)
//			printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );
//		delete iter_com;
//	}
//
//	// Ri-matrix, Si-matrix, & dummy-matrix initialization
//	Ri = (double **) malloc( sizeof(double *) * 3 );
//	Si = (double **) malloc( sizeof(double *) * 3 );
//	dummy = (double **) malloc( sizeof(double *) * 3 );
//	// 1st Si --> I-matrix
//	// 1st Ti --> (0,0,0) vector
//	for(i=0; i<3; i++)
//	{
//		Ri[i] = (double *) malloc( sizeof(double) * 3 );
//		Si[i] = (double *) malloc( sizeof(double) * 3 );
//		dummy[i] = (double *) malloc( sizeof(double) * 3 );
//		Ti[i] = 0.0;
//		for(j=0; j<3; j++)
//			if(i==j)
//			{
//				Si[i][j] = 1.0;
//				dummy[i][j] = 1.0;
//			}
//			else
//			{
//				Si[i][j] = 0.0;
//				dummy[i][j] = 0.0;
//			}
//	}
//
//	if(debug)
//		printf("Msg(move_dihedralM): Moving dihedrals:\n");
//
//	Segment * seg;
//	j = 0;
//	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter = new pdbIter( seg );
//		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)
//
//		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			T[0] = Ti[0];
//			T[1] = Ti[1];
//			T[2] = Ti[2];
//
//			// 3 Translations
//			for(axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					Ti2[axis] = uu[j2] * step;
//					if(debug)
//						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//					j2++;
//				}
//				j++;
//			}
//
//			// 3 Rotations
//			for(axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					e1[0] = 0;
//					e1[1] = 0;
//					e1[2] = 0;
//					e1[axis] = 1.0; // e1 ==> Rotation axis
//
//					// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
//					rotmat(uu[j2] * step, e1, Ri);
//
//					if(debug)
//						printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//					// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//					mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//					temp = Si; // swaping...
//					Si = dummy;
//					dummy = temp;
//
//					// This avoids pre-centering...
//					T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
//					T[1] -= com[1];
//					T[2] -= com[2];
//					// Rotating the "rotation origin position" -->
//					multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)
//
//					T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
//					T[1] = Ti3[1] + com[1];
//					T[2] = Ti3[2] + com[2];
//
//					j2++;
//				}
//				j++;
//			}
//
//			Ti[0] = T[0] + Ti2[0];
//			Ti[1] = T[1] + Ti2[1];
//			Ti[2] = T[2] + Ti2[2];
//		}
//
//		if(debug)
//			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);
//
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			res = ( Residue * ) iter->get_fragment();
//			iter_res_atom = new pdbIter( res ); // iter residue atoms
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType(); // gets current fragment moltype
//
//			if(fragtype == tmol_protein)
//			{
//				// We should move now "non-first segment" first residue
//				// (first segment first residue does not move!)
//				if(iter->pos_fragment==0 && iter_seg->pos_segment != 0 ) // non-first segment && first residue
//					// rotating current residue atoms (due to backwards internal coords)
//					for( iter_res_atom->pos_atom = 0; !iter_res_atom->gend_atom(); iter_res_atom->next_atom() )
//						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti); // direct rot.trans.
//
//				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
//				if( !(iter->pos_fragment == 0 || iter->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					// rotating N and CA atoms (due to previous residues)
//					for(iter_res_atom->pos_atom = 0; iter_res_atom->pos_atom < 2; iter_res_atom->next_atom() )
//						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti); // direct rot.trans.
//
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PHI
//						// Computing axis
//						iter_res_atom->pos_atom = 0; // N
//						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
//						e1[1] = ca[1] - nh[1];
//						e1[2] = ca[2] - nh[2];
//
//						// Computing PHI rotation matrix given "e1" and the current "phi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - ca[1];
//						T[2] = Ti[2] - ca[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += ca[0]; // translating it back to the Qi position
//						Ti[1] += ca[1];
//						Ti[2] += ca[2];
//
//						j2++;
//					}
//
//					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
//					iter_res_atom->pos_atom = 2; // CO moves only!
//					( iter_res_atom->get_atom() )->rotate_trans(Si,Ti); // direct rot.trans.
//					j++;
//				}
//
//				// If it's PRO we must move now !!!
//				// (Watch out this: When switching to Full-atom model)
//				if( iter->pos_fragment == num_res-1 || (resn == PRO && iter->pos_fragment != 0 ) )
//					// rotating current residue atoms (due to backwards dihedrals)
//					for(iter_res_atom->pos_atom = 0; !iter_res_atom->gend_atom(); iter_res_atom->next_atom() )
//						( iter_res_atom->get_atom() )->rotate_trans(Si,Ti); // direct rot.trans.
//
//				// PSI
//				// First and Last residue never have PSI in CA-model
//				if( iter->pos_fragment != num_res-1 && iter->pos_fragment != 0 )
//				{
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PSI
//						// Computing axis
//						iter_res_atom->pos_atom = 2; // C
//						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
//						e1[1] = co[1] - ca[1];
//						e1[2] = co[2] - ca[2];
//
//						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - co[1];
//						T[2] = Ti[2] - co[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += co[0]; // translating it back to the Qi position
//						Ti[1] += co[1];
//						Ti[2] += co[2];
//
//						j2++;
//					}
//					// During PSI-rotation, no current residue atoms are moved !!!
//					j++;
//				}
//			}
//			// END PROTEIN.....
//			else
//			{
//				printf("Msg(move_dihedralMCAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
//				exit(1);
//			}
//			delete iter_res_atom;
//			index_res++; // counts residues
//		}
//		delete iter;
//	}
//	delete iter_seg;
//
//	// Freeing memory
//	for(i=0; i<3; i++)
//	{
//		free( Ri[i] );
//		free( Si[i] );
//		free( dummy[i] );
//	}
//	free( Ri );
//	free( Si );
//	free( dummy );
//}

//// Moves atoms given an Internal Coordinates Normal Mode and a Step
//// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
//// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
//// If fix=NULL, no fixed dihedrals (equal to previous-version)
//void move_dihedralMCAx_old(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix)
//{
//	bool debug = false;
//	int j=0;
//	int j2=0;
//	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
//	Residue *res;
//	Atom *atom;
//	pdbIter *iter,*iter_res_atom,*iter_seg;
//	double **Ri; // current rotation (Ri-matrix);
//	double **Si; // accummulated rotation (Si-matrix)
//	double **dummy; // dummy rotation (dummy-matrix)
//	double **temp;
//	int num_res = 0;
//	int index_res=0;
//	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
//	int indexbase;
//	int resn;
//
//	// Centering Macromolecule (Temporal)
//	// Computing the PDB's Center of Mass (CoM)
//	double mtot,mta;
//	double com[3];
//	pdbIter *iter_com = new pdbIter(mol);
//	Tcoor pos;
//	com[0] = com[1] = com[2] = 0.0;
//	mtot = 0.0;
//	for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
//	{
//		atom = ( iter_com->get_atom() );
//		mta = atom->getPdbocc(); // Load mass...
//		mtot += mta;
//		atom->getPosition(pos);
//		/* Sum(mass*coord) before putting the CoM at 0 */
//		com[0] += mta * pos[0];
//		com[1] += mta * pos[1];
//		com[2] += mta * pos[2];
//	}
//	com[0] /= mtot;
//	com[1] /= mtot;
//	com[2] /= mtot;
//	if(debug)
//		printf( "Msg(move_dihedralM): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );
//
//	// Ri-matrix, Si-matrix, & dummy-matrix initialization
//	Ri = (double **) malloc( sizeof(double *) * 3 );
//	Si = (double **) malloc( sizeof(double *) * 3 );
//	dummy = (double **) malloc( sizeof(double *) * 3 );
//	// 1st Si --> I-matrix
//	// 1st Ti --> (0,0,0) vector
//	for(int i=0; i<3; i++)
//	{
//		Ri[i] = (double *) malloc( sizeof(double) * 3 );
//		Si[i] = (double *) malloc( sizeof(double) * 3 );
//		dummy[i] = (double *) malloc( sizeof(double) * 3 );
//		Ti[i] = 0.0;
//		Ti2[i] = 0.0;
//		for(int j=0; j<3; j++)
//			if(i==j)
//			{
//				Si[i][j] = 1.0;
//				dummy[i][j] = 1.0;
//			}
//			else
//			{
//				Si[i][j] = 0.0;
//				dummy[i][j] = 0.0;
//			}
//	}
//
//	if(debug)
//		printf("Msg(move_dihedralM): Moving dihedrals:\n");
//
//	Segment * seg;
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	int num_seg = iter_seg->num_segment();
//	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter = new pdbIter( seg );
//		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)
//
//		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			T[0] = Ti[0];
//			T[1] = Ti[1];
//			T[2] = Ti[2];
//
//			// 3 Translations
//			for(int axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					Ti2[axis] = uu[j2] * step;
//					if(debug)
//						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//					j2++;
//				}
//				j++;
//			}
//
//			// 3 Rotations
//			for(int axis=0; axis<3; axis++)
//			{
//				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//				{
//					e1[0] = 0;
//					e1[1] = 0;
//					e1[2] = 0;
//					e1[axis] = 1.0; // e1 ==> Rotation axis
//
//					// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
//					rotmat(uu[j2] * step, e1, Ri);
//
//					if(debug)
//						printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//					// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//					mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//					temp = Si; // swaping...
//					Si = dummy;
//					dummy = temp;
//
//					// This avoids pre-centering...
//					T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
//					T[1] -= com[1];
//					T[2] -= com[2];
//					// Rotating the "rotation origin position" -->
//					multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)
//
//					T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
//					T[1] = Ti3[1] + com[1];
//					T[2] = Ti3[2] + com[2];
//
//					j2++;
//				}
//				j++;
//			}
//
//			Ti[0] = T[0] + Ti2[0];
//			Ti[1] = T[1] + Ti2[1];
//			Ti[2] = T[2] + Ti2[2];
//		}
//
//		if(debug)
//			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);
//
//		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
//		{
//			res = ( Residue * ) iter->get_fragment();
//			iter_res_atom = new pdbIter( res ); // iter residue atoms
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType(); // gets current fragment moltype
//
//			if(fragtype == tmol_protein)
//			{
//				// We should move now "non-first segment" first residue
//				// (first segment first residue does not move!)
//				if(iter_seg->pos_segment != 0 && iter->pos_fragment==0) // non-first segment && first residue
//				{
//					// rotating current residue atoms (due to backwards internal coords)
//					for ( iter_res_atom->pos_atom = 0;
//						!iter_res_atom->gend_atom();
//						iter_res_atom->next_atom()   )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//				}
//
//				// PHI --> NOT FIRST, NOT PRO, NOT LAST too! (CA-model)
//				if( !(iter->pos_fragment == 0 || iter->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					// rotating N and CA atoms (due to previous residues)
//					for ( iter_res_atom->pos_atom = 0; // N and CA
//						iter_res_atom->pos_atom < 2;
//						iter_res_atom->next_atom()   )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PHI
//						// Computing axis
//						iter_res_atom->pos_atom = 0; // N
//						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
//						e1[1] = ca[1] - nh[1];
//						e1[2] = ca[2] - nh[2];
//
//						// Computing PHI rotation matrix given "e1" and the current "phi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - ca[1];
//						T[2] = Ti[2] - ca[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += ca[0]; // translating it back to the Qi position
//						Ti[1] += ca[1];
//						Ti[2] += ca[2];
//
//						j2++;
//					}
//
//					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
//					iter_res_atom->pos_atom = 2; // CO moves only!
//					atom = (Atom *) iter_res_atom->get_atom();
//					atom->getPosition( tv ) ; // position before rotation
//					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//					atom->setPosition( tvp ) ; // position after rotation
//
//					j++;
//				}
//
//				// If it's PRO we must move now !!!
//				// (Watch out this: When switching to Full-atom model)
//				if( iter->pos_fragment == num_res-1 ||
//						( strcmp(res->getName(), "PRO") == 0  && iter->pos_fragment != 0 ) )
//				{
//					// rotating current residue atoms (due to backwards dihedrals)
//					for( iter_res_atom->pos_atom = 0;
//						!iter_res_atom->gend_atom();
//						iter_res_atom->next_atom() )
//					{
//						atom = (Atom *) iter_res_atom->get_atom();
//						atom->getPosition( tv ) ; // position before rotation
//						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
//						atom->setPosition( tvp ) ; // position after rotation
//					}
//				}
//
//				// PSI
//				// First and Last residue never have PSI in CA-model
//				if( iter->pos_fragment != num_res-1 && iter->pos_fragment != 0 )
//				{
//					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
//					{
//						// Rotate PSI
//						// Computing axis
//						iter_res_atom->pos_atom = 2; // C
//						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
//						iter_res_atom->pos_atom = 1; // CA
//						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
//						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
//						e1[1] = co[1] - ca[1];
//						e1[2] = co[2] - ca[2];
//
//						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
//						rotmat(uu[j2] * step, e1, Ri);
//
//						if(debug)
//							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
//
//						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
//						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
//						temp = Si; // swaping...
//						Si = dummy;
//						dummy = temp;
//						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
//						T[1] = Ti[1] - co[1];
//						T[2] = Ti[2] - co[2];
//						// Rotating the "rotation origin position" -->
//						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
//						Ti[0] += co[0]; // translating it back to the Qi position
//						Ti[1] += co[1];
//						Ti[2] += co[2];
//
//						j2++;
//					}
//					// During PSI-rotation, no current residue atoms are moved !!!
//					j++;
//				}
//			}
//			// END PROTEIN.....
//			else
//			{
//				printf("Msg(move_dihedralMCAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
//				exit(1);
//			}
//			delete iter_res_atom;
//			index_res++; // counts residues
//		}
//		delete iter;
//	}
//	delete iter_seg;
//
//	// Freeing memory
//	for(int i=0; i<3; i++)
//	{
//		free( Ri[i] );
//		free( Si[i] );
//		free( dummy[i] );
//	}
//	free( Ri );
//	free( Si );
//	free( dummy );
//	delete iter_com;
//}

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain/Prot/DNA/RNA/SMOL): 3BB2R & Full-Atom models
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
void move_dihedralMFAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix, bool *addrot)
{
	bool debug = false;
	int j=0;
	int j2=0;
	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
	Residue *res;
	Atom *atom;
	pdbIter *iter,*iter_res_atom,*iter_seg;
	double **Ri; // current rotation (Ri-matrix);
	double **Si; // accummulated rotation (Si-matrix)
	double **dummy; // dummy rotation (dummy-matrix)
	double **temp;
	int num_res = 0;
	int index_res=0;
	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
	int indexbase;
	int resn;

	// Centering Macromolecule (Temporal)
	// Computing the PDB's Center of Mass (CoM)
	double mtot,mta;
	double com[3];
	pdbIter *iter_com = new pdbIter(mol);
	Tcoor pos;
	com[0] = com[1] = com[2] = 0.0;
	mtot = 0.0;
	for( iter_com->pos_atom = 0; !iter_com->gend_atom(); iter_com->next_atom() )  // screens all-atoms
	{
		atom = ( iter_com->get_atom() );
		mta = atom->getPdbocc(); // Load mass...
		mtot += mta;
		atom->getPosition(pos);
		// Sum(mass*coord) before putting the CoM at 0
		com[0] += mta * pos[0];
		com[1] += mta * pos[1];
		com[2] += mta * pos[2];
	}
	com[0] /= mtot;
	com[1] /= mtot;
	com[2] /= mtot;
	if(debug)
		printf( "Msg(move_dihedralMFAx): Mass %8.8f Center %8.8f %8.8f %8.8f (shifting it to 0,0,0)\n", mtot, com[0], com[1], com[2] );

	// Ri-matrix, Si-matrix, & dummy-matrix initialization
	Ri = (double **) malloc( sizeof(double *) * 3 );
	Si = (double **) malloc( sizeof(double *) * 3 );
	dummy = (double **) malloc( sizeof(double *) * 3 );
	// 1st Si --> I-matrix
	// 1st Ti --> (0,0,0) vector
	for(int i=0; i<3; i++)
	{
		Ri[i] = (double *) malloc( sizeof(double) * 3 );
		Si[i] = (double *) malloc( sizeof(double) * 3 );
		dummy[i] = (double *) malloc( sizeof(double) * 3 );
		Ti[i] = 0.0;
		Ti2[i] = 0.0;
		for(int j=0; j<3; j++)
			if(i==j)
			{
				Si[i][j] = 1.0;
				dummy[i][j] = 1.0;
			}
			else
			{
				Si[i][j] = 0.0;
				dummy[i][j] = 0.0;
			}
	}

	if(debug)
		printf("Msg(move_dihedralM): Moving dihedrals:\n");

	Segment * seg;
	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	int num_seg;
	//,seg_atoms;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

	num_seg = iter_seg->num_segment();

	for ( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() ) // Screening segments
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter = new pdbIter( seg );
		num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)

		// Fixes bug?
		fragtype = seg->getMolType(); // gets current fragment moltype

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		// Moving the inter-segment degrees of freedom (3 Trans and 3 Rots)
		if(iter_seg->pos_segment != 0) // non-first segment
		{
			T[0] = Ti[0];
			T[1] = Ti[1];
			T[2] = Ti[2];

//			// Checking single atom segments... SMOL's
//			seg_atoms = seg->num_atoms();

			// 3 Translations
			for(int axis=0; axis<3; axis++)
			{
				if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
				{
					Ti2[axis] = uu[j2] * step;
					if(debug)
						printf("%4d Translations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
					j2++;
				}
				j++;
			}

			// 3 Rotations
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{
						e1[0] = 0;
						e1[1] = 0;
						e1[2] = 0;
						e1[axis] = 1.0; // e1 ==> Rotation axis

						// Computing Rot[axis] rotation matrix given "e1" and the current R[axis] angle
						rotmat(uu[j2] * step, e1, Ri);

						if(debug)
							printf("%4d Rotations (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						// This avoids macromolecule pre-centering...
						T[0] -= com[0]; // translating to the origin the Translational part of M due to rotation
						T[1] -= com[1];
						T[2] -= com[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti3); // Ri x T = Ti (Applies a rotation to a position vector)
						T[0] = Ti3[0] + com[0]; // translating it back to the Qi position
						T[1] = Ti3[1] + com[1];
						T[2] = Ti3[2] + com[2];
						j2++;
					}
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++;
				}

			// moving after
			Ti[0] = T[0] + Ti2[0];
			Ti[1] = T[1] + Ti2[1];
			Ti[2] = T[2] + Ti2[2];
		}

		if(debug)
			printf("\nProcessing segment %d (%d)\n",iter_seg->pos_segment,num_res);

		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
		{
			res = ( Residue * ) iter->get_fragment();
			iter_res_atom = new pdbIter( res ); // iter residue atoms
//			fragtype = res->getMolType(); // gets current fragment moltype
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			if(fragtype == tmol_protein)
			{
				// We should move now "non-first segment" first residue
				// (first segment first residue does not move!)
				if(iter_seg->pos_segment != 0 && iter->pos_fragment==0) // non-first segment && first residue
				{
					// rotating current residue atoms (due to backwards internal coords)
					for ( iter_res_atom->pos_atom = 0;
					!iter_res_atom->gend_atom();
					iter_res_atom->next_atom()   )
						if( !(iter_res_atom->pos_atom == 3 && model == 2) ) // if not O-atom in Full-atom
						{
							atom = (Atom *) iter_res_atom->get_atom();
							atom->getPosition( tv ) ; // position before rotation
							rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
							atom->setPosition( tvp ) ; // position after rotation
						}
				}

				// PHI (No first, No PRO)
				if( iter->pos_fragment != 0 && (strcmp(res->getName(), "PRO") != 0) )
				{
					// rotating N and CA atoms (due to previous residues)
					for ( iter_res_atom->pos_atom = 0; // N and CA
					iter_res_atom->pos_atom < 2;
					iter_res_atom->next_atom()   )
					{
						atom = (Atom *) iter_res_atom->get_atom();
						atom->getPosition( tv ) ; // position before rotation
						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
						atom->setPosition( tvp ) ; // position after rotation
					}

					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PHI
						// Computing axis
						iter_res_atom->pos_atom = 0; // N
						( iter_res_atom->get_atom() )->getPosition( nh ); // N position
						iter_res_atom->pos_atom = 1; // CA
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = ca[0] - nh[0]; // e1 ==> vector N-->CA (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing PHI rotation matrix given "e1" and the current "phi" angle
						rotmat(uu[j2] * step, e1, Ri);
						// d_rotmat(uu[j2] * step / (float) maxrot, e1, Ri);

						if(debug)
							printf("%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];


						j2++;
					}

					// rotating current-PHi affected residue atoms (due to backwards dihedrals and current-Phi)
					for(iter_res_atom->pos_atom = 2;!iter_res_atom->gend_atom();iter_res_atom->next_atom())
						if( !(iter_res_atom->pos_atom == 3 && model == 2) ) // if not O-atom in Full-atom
						{
							atom = (Atom *) iter_res_atom->get_atom();
							atom->getPosition( tv ) ; // position before rotation
							rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
							atom->setPosition( tvp ) ; // position after rotation
						}

					j++;
				}

				// If it's PRO we must move now !!!
				// (Watch out this: When switching to Full-atom model)
				if( strcmp(res->getName(), "PRO") == 0 && iter->pos_fragment != 0 )
				{
					// rotating current residue atoms (due to backwards dihedrals)
					for( iter_res_atom->pos_atom = 0;
					!iter_res_atom->gend_atom();
					iter_res_atom->next_atom() )
						if( !(iter_res_atom->pos_atom == 3 && model==2) ) // excepting O-atom in Full-atom model
						{
							atom = (Atom *) iter_res_atom->get_atom();
							atom->getPosition( tv ) ; // position before rotation
							rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
							atom->setPosition( tvp ) ; // position after rotation
						}
				}

				// CHI (3BB2R)
				// Note that "resn" is available, you can check numericaly Aminoacid identities!!! (faster)
				if( type == 2 )
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
											(iter->pos_fragment==0 ||
													(iter->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
						{
							// Rotate CHI
							// Computing axis
							if(model==2)
								iter_res_atom->pos_atom = 4; // CB full-atom
							else
								iter_res_atom->pos_atom = 3; // CB 3bb2r
							( iter_res_atom->get_atom() )->getPosition( co ); // CB position (R) (further needed: Qi)
							iter_res_atom->pos_atom = 1; // CA
							( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
							e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->CB (P-->Q)
							e1[1] = co[1] - ca[1];
							e1[2] = co[2] - ca[2];

							// Computing CHI rotation matrix given "e1" (O-->Y) and the current "chi" angle
							rotmat(uu[j2] * step, e1, Ri);
							// d_rotmat(uu[j2] * step / (float) maxrot, e1, Ri);

							if(debug)
								printf("%4d CHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

							// NOT-Updating acummulated rotation matrix (Si). (Rotational part of Mi)
							// We will rotate CHI after PHI rotation-traslation is applied!
							// ( We are moveing just one atom: O )

							// translating to the origin the Translational part of M due to rotation
							T[0] = - co[0]; // T = -CB-position
							T[1] = - co[1];
							T[2] = - co[2];
							// Rotating the "rotation origin position" -->
							// Ri x T = Ti (Applies a rotation to a position vector)
							multvec3(Ri, T, nh);
							nh[0] += co[0];
							nh[1] += co[1];
							nh[2] += co[2];

							// Move CHI
							if(model==2)
							{
								// With Full-Atom, only atoms after CB move!!!
								for ( iter_res_atom->pos_atom = 5;
								!iter_res_atom->gend_atom();
								iter_res_atom->next_atom()   )
								{
									atom = (Atom *) iter_res_atom->get_atom();
									atom->getPosition( tv ) ; // position before rotation
									rmot(Ri, nh, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
									atom->setPosition( tvp ) ; // position after rotation
								}
							}
							else
							{
								// With 3BB2R, only "O" moves!!!
								iter_res_atom->pos_atom = 4;
								atom = (Atom *) iter_res_atom->get_atom();
								atom->getPosition( tv ) ; // position before rotation
								rmot(Ri, nh, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
								atom->setPosition( tvp ) ; // position after rotation
							}

							j2++;
						}
						j++;
					}

				// PSI
				if( iter->pos_fragment != num_res-1 || model == 2 ) // Segment Last residue never has PSI
				{													// (excepting in Full-Atom!)
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Rotate PSI
						// Computing axis
						iter_res_atom->pos_atom = 2; // C
						( iter_res_atom->get_atom() )->getPosition( co ); // C position (further needed: Qi)
						iter_res_atom->pos_atom = 1; // CA
						( iter_res_atom->get_atom() )->getPosition( ca ); // CA position
						e1[0] = co[0] - ca[0]; // e1 ==> vector CA-->C (P-->Q)
						e1[1] = co[1] - ca[1];
						e1[2] = co[2] - ca[2];

						// Computing PSI rotation matrix given "e1" (O-->Y) and the current "psi" angle
						rotmat(uu[j2] * step, e1, Ri);
						// d_rotmat(uu[j2] * step / (float) maxrot, e1, Ri);

						if(debug)
							printf("%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - co[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - co[1];
						T[2] = Ti[2] - co[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += co[0]; // translating it back to the Qi position
						Ti[1] += co[1];
						Ti[2] += co[2];

						j2++;
					}
					// During PSI-rotation, no current residue atoms are moved !!!

					// Excepting in Full-Atom model!
					// In Full-Atom, the O-atom moves
					if(model == 2)
					{
						// rotating O-atom
						iter_res_atom->pos_atom = 3;
						atom = (Atom *) iter_res_atom->get_atom();
						atom->getPosition( tv ) ; // position before rotation
						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
						atom->setPosition( tvp ) ; // position after rotation
					}
					j++;
				}
			}
			// END PROTEIN.....
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				// rotating P (due to backwards internal coords: Epsilon or 3Trans + 3Rot)
				if( !(iter_seg->pos_segment == 0 && iter->pos_fragment == 0) ) // non-first segment && first residue
				{
					if(debug)
						printf("Rotating P (due to backwards internal coords: Epsilon or 3Trans + 3Rot)\n");
					iter_res_atom->pos_atom=0;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}

				// if Not-first fragment
				// previous residue ETA (the O3* position is stored in "nh", see below)
				// (the last fragment doesn't have EPSILON and ETA)
				if( iter->pos_fragment != 0 )
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 0; // Final current-P
						( iter_res_atom->get_atom() )->getPosition( ca );
						// Initial O3* already computed below... (stored in "ca")
						e1[0] = ca[0] - nh[0]; // e1 ==> vector O5* --> C5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing previous-ETA rotation matrix given "e1" and the current "eta" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d ETA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];

					}
					j++;
				}

				// First segment and first fragment's O1P,O2P,O5* don't move!
				// if Not-First segment and Not-First fragment
				if( !(iter_seg->pos_segment == 0 && iter->pos_fragment == 0) ) // non-first segment && first residue
				{
					// rotating P,O1P,O2P,O5* (due to backwards internal coords)
					if(debug)
						printf("Rotating O1P,O2P,O5* (due to backwards internal coords: Epsilon + Eta or 3Trans + 3Rot)\n");
					for( iter_res_atom->pos_atom=1; iter_res_atom->pos_atom<=3; iter_res_atom->next_atom() )
					{
						atom = (Atom *) iter_res_atom->get_atom();
						atom->getPosition( tv ) ; // position before rotation
						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
						atom->setPosition( tvp ) ; // position after rotation
					}
				}

				// ALPHA
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 3; // Final O5*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 0; // Initial P
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector P --> O5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing ALPHA rotation matrix given "e1" and the current "alpha" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d ALPHA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// BETA
				{
					// Moving C5* atom (+4) before updating BETA's Ti and Si
					if(debug)
						printf("Rotating C5* atom (+4) before updating BETA's Ti and Si\n");
					iter_res_atom->pos_atom = 4;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation

					// we should get Final and Initial positions first! (before moving)
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 4; // Final C5*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 3; // Initial O5*
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector O5* --> C5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing BETA rotation matrix given "e1" and the current "beta" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d BETA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// GAMMA
				{
					// Moving C4* atom (+5) before updating GAMMA's Ti and Si
					if(debug)
						printf("Rotating C4* atom (+5) before updating GAMMA's Ti and Si\n");
					iter_res_atom->pos_atom = 5;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation

					// we should get Final and Initial positions first! (before moving)
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 5; // Final C4*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 4; // Initial C5*
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector C5* --> C4* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing GAMMA rotation matrix given "e1" and the current "gamma" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d GAMMA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// Moving all remaining atoms of the nucleotide
				// (before updating "LOCAL" CHI's T and Ri)
				if(debug)
					printf("Rotating all remaining atoms of the nucleotide\n");
				for( iter_res_atom->pos_atom=6; !iter_res_atom->gend_atom(); iter_res_atom->next_atom() )
				{
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}

				// Moving CHI (with Ri and nh, without Si and Ti)
				if( type == 2 )
				{
					// we should get Final and Initial positions after moving!
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
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

						// Computing axis
						iter_res_atom->pos_atom = indexbase; // Final N1/N9
						( iter_res_atom->get_atom() )->getPosition( co ); // Final N1/N9
						iter_res_atom->pos_atom = 11; // Initial C1*
						( iter_res_atom->get_atom() )->getPosition( ca );
						e1[0] = co[0] - ca[0]; // e1 ==> vector C1* --> N1/N9 (P-->Q)
						e1[1] = co[1] - ca[1];
						e1[2] = co[2] - ca[2];

						// Computing CHI rotation matrix given "e1" (O-->Y) and the current "chi" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d CHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// NOT-Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						// We will rotate CHI after PHI rotation-traslation is applied!
						// ( We are moveing just one atom: O )

						// translating to the origin the Translational part of M due to rotation
						T[0] = -co[0]; // T = -CB-position
						T[1] = -co[1];
						T[2] = -co[2];
						// Rotating the "rotation origin position" -->
						// Ri x T = Ti (Applies a rotation to a position vector)
						multvec3(Ri, T, nh);
						nh[0] += co[0];
						nh[1] += co[1];
						nh[2] += co[2];

						// With Full-Atom, only atoms after N1/N9 move due to CHI
						if(debug)
							printf("Rotating only atoms after N1/N9 move due to CHI\n");
						for ( iter_res_atom->pos_atom = 12;
						!iter_res_atom->gend_atom();
						iter_res_atom->next_atom()   )
							if(iter_res_atom->pos_atom != indexbase)
							{
								atom = (Atom *) iter_res_atom->get_atom();
								atom->getPosition( tv ) ; // position before rotation
								rmot(Ri, nh, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
								atom->setPosition( tvp ) ; // position after rotation
							}
					}
					j++;
				}

				// if Not-Last
				if( iter->pos_fragment != num_res-1 )
				{
					// EPSILON
					{
						if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
						{
							// Computing axis
							iter_res_atom->pos_atom = 8; // Final O3*
							( iter_res_atom->get_atom() )->getPosition( ca );
							iter_res_atom->pos_atom = 7; // Initial C3*
							( iter_res_atom->get_atom() )->getPosition( nh );
							e1[0] = ca[0] - nh[0]; // e1 ==> vector C3* --> O3* (P-->Q)
							e1[1] = ca[1] - nh[1];
							e1[2] = ca[2] - nh[2];

							// Computing EPSILON rotation matrix given "e1" and the current "epsilon" angle
							rotmat(uu[j2] * step, e1, Ri);
							if(debug)
								printf("%4d EPSILON (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
							j2++;

							// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
							mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
							temp = Si; // swaping...
							Si = dummy;
							dummy = temp;
							T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
							T[1] = Ti[1] - ca[1];
							T[2] = Ti[2] - ca[2];
							// Rotating the "rotation origin position" -->
							multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
							Ti[0] += ca[0]; // translating it back to the Qi position
							Ti[1] += ca[1];
							Ti[2] += ca[2];
						}
						j++;
					}

					// ETA
					{
						if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
						{
							// Computing axis
							iter_res_atom->pos_atom = 8; // Initial O3*
							( iter_res_atom->get_atom() )->getPosition( nh );
							// the Final position will be computed before ALPHA (see above)
							// (because next-P position is not available yet (in "res") !!!!)
						}
					}
				}
			} // END RNA
			else if( resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// rotating P (due to backwards internal coords: Epsilon or 3Trans + 3Rot)
				if( !(iter_seg->pos_segment == 0 && iter->pos_fragment == 0) ) // non-first segment && first residue
				{
					if(debug)
						printf("Rotating P (due to backwards internal coords: Epsilon or 3Trans + 3Rot)\n");
					iter_res_atom->pos_atom=0;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}

				// if Not-first fragment
				// previous residue ETA (the O3* position is stored in "nh", see below)
				// (the last fragment doesn't have EPSILON and ETA)
				if( iter->pos_fragment != 0 )
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (otherwise --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 0; // Final current-P
						( iter_res_atom->get_atom() )->getPosition( ca );
						// Initial O3* already computed below... (stored in "ca")
						e1[0] = ca[0] - nh[0]; // e1 ==> vector O5* --> C5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing previous-ETA rotation matrix given "e1" and the current "eta" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d ETA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// First segment and first fragment's O1P,O2P,O5* don't move!
				// if Not-First segment and Not-First fragment
				if( !(iter_seg->pos_segment == 0 && iter->pos_fragment == 0) ) // non-first segment && first residue
				{
					// rotating P,O1P,O2P,O5* (due to backwards internal coords)
					if(debug)
						printf("Rotating O1P,O2P,O5* (due to backwards internal coords: Epsilon + Eta or 3Trans + 3Rot)\n");
					for( iter_res_atom->pos_atom=1; iter_res_atom->pos_atom<=3; iter_res_atom->next_atom() )
					{
						atom = (Atom *) iter_res_atom->get_atom();
						atom->getPosition( tv ) ; // position before rotation
						rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
						atom->setPosition( tvp ) ; // position after rotation
					}
				}

				// ALPHA
				{
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 3; // Final O5*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 0; // Initial P
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector P --> O5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing ALPHA rotation matrix given "e1" and the current "alpha" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d ALPHA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// BETA
				{
					// Moving C5* atom (+4) before updating BETA's Ti and Si
					if(debug)
						printf("Rotating C5* atom (+4) before updating BETA's Ti and Si\n");
					iter_res_atom->pos_atom = 4;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation

					// we should get Final and Initial positions first! (before moving)
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 4; // Final C5*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 3; // Initial O5*
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector O5* --> C5* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing BETA rotation matrix given "e1" and the current "beta" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d BETA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// GAMMA
				{
					// Moving C4* atom (+5) before updating GAMMA's Ti and Si
					if(debug)
						printf("Rotating C4* atom (+5) before updating GAMMA's Ti and Si\n");
					iter_res_atom->pos_atom = 5;
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation

					// we should get Final and Initial positions first! (before moving)
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
						// Computing axis
						iter_res_atom->pos_atom = 5; // Final C4*
						( iter_res_atom->get_atom() )->getPosition( ca );
						iter_res_atom->pos_atom = 4; // Initial C5*
						( iter_res_atom->get_atom() )->getPosition( nh );
						e1[0] = ca[0] - nh[0]; // e1 ==> vector C5* --> C4* (P-->Q)
						e1[1] = ca[1] - nh[1];
						e1[2] = ca[2] - nh[2];

						// Computing GAMMA rotation matrix given "e1" and the current "gamma" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d GAMMA (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
						temp = Si; // swaping...
						Si = dummy;
						dummy = temp;
						T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
						T[1] = Ti[1] - ca[1];
						T[2] = Ti[2] - ca[2];
						// Rotating the "rotation origin position" -->
						multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
						Ti[0] += ca[0]; // translating it back to the Qi position
						Ti[1] += ca[1];
						Ti[2] += ca[2];
					}
					j++;
				}

				// Moving all remaining atoms of the nucleotide
				// (before updating "LOCAL" CHI's T and Ri)
				if(debug)
					printf("Rotating all remaining atoms of the nucleotide\n");
				for( iter_res_atom->pos_atom=6; !iter_res_atom->gend_atom(); iter_res_atom->next_atom() )
				{
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}

				// Moving CHI (with Ri and nh, without Si and Ti)
				if( type == 2 )
				{
					// we should get Final and Initial positions after moving!
					if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
					{
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

						// Computing axis
						iter_res_atom->pos_atom = indexbase; // Final N1/N9
						( iter_res_atom->get_atom() )->getPosition( co ); // Final N1/N9
						iter_res_atom->pos_atom = 10; // Initial C1* (DNA)
						( iter_res_atom->get_atom() )->getPosition( ca );
						e1[0] = co[0] - ca[0]; // e1 ==> vector C1* --> N1/N9 (P-->Q)
						e1[1] = co[1] - ca[1];
						e1[2] = co[2] - ca[2];

						// Computing CHI rotation matrix given "e1" (O-->Y) and the current "chi" angle
						rotmat(uu[j2] * step, e1, Ri);
						if(debug)
							printf("%4d CHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
						j2++;

						// NOT-Updating acummulated rotation matrix (Si). (Rotational part of Mi)
						// We will rotate CHI after PHI rotation-traslation is applied!
						// ( We are moveing just one atom: O )

						// translating to the origin the Translational part of M due to rotation
						T[0] = -co[0]; // T = -CB-position
						T[1] = -co[1];
						T[2] = -co[2];
						// Rotating the "rotation origin position" -->
						// Ri x T = Ti (Applies a rotation to a position vector)
						multvec3(Ri, T, nh);
						nh[0] += co[0];
						nh[1] += co[1];
						nh[2] += co[2];

						// With Full-Atom, only atoms after N1/N9 move due to CHI
						if(debug)
							printf("Rotating only atoms after N1/N9 move due to CHI\n");
						for( iter_res_atom->pos_atom = 11; !iter_res_atom->gend_atom(); iter_res_atom->next_atom() ) // DNA
							if(iter_res_atom->pos_atom != indexbase)
							{
								atom = (Atom *) iter_res_atom->get_atom();
								atom->getPosition( tv ) ; // position before rotation
								rmot(Ri, nh, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
								atom->setPosition( tvp ) ; // position after rotation
							}
					}
					j++;
				}

				// if Not-Last
				if( iter->pos_fragment != num_res-1 )
				{
					// EPSILON
					{
						if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
						{
							// Computing axis
							iter_res_atom->pos_atom = 8; // Final O3*
							( iter_res_atom->get_atom() )->getPosition( ca );
							iter_res_atom->pos_atom = 7; // Initial C3*
							( iter_res_atom->get_atom() )->getPosition( nh );
							e1[0] = ca[0] - nh[0]; // e1 ==> vector C3* --> O3* (P-->Q)
							e1[1] = ca[1] - nh[1];
							e1[2] = ca[2] - nh[2];

							// Computing EPSILON rotation matrix given "e1" and the current "epsilon" angle
							rotmat(uu[j2] * step, e1, Ri);
							if(debug)
								printf("%4d EPSILON (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);
							j2++;

							// Updating acummulated rotation matrix (Si). (Rotational part of Mi)
							mult3(Ri,Si,dummy); // Si = Ri x S(i-1)
							temp = Si; // swaping...
							Si = dummy;
							dummy = temp;
							T[0] = Ti[0] - ca[0]; // translating to the origin the Translational part of M due to rotation
							T[1] = Ti[1] - ca[1];
							T[2] = Ti[2] - ca[2];
							// Rotating the "rotation origin position" -->
							multvec3(Ri, T, Ti); // Ri x T = Ti (Applies a rotation to a position vector)
							Ti[0] += ca[0]; // translating it back to the Qi position
							Ti[1] += ca[1];
							Ti[2] += ca[2];
						}
						j++;
					}

					// ZETA
					{
						if(fix == NULL || fix[j]) // "fix == NULL" must go 1st!!! (if not --> segmentation fault!)
						{
							// Computing axis
							iter_res_atom->pos_atom = 8; // Initial O3*
							( iter_res_atom->get_atom() )->getPosition( nh );
							// the Final position will be computed before ALPHA (see above)
							// (because next-P position is not available yet (in "res") !!!!)
						}
					}
				}
			} // END DNA
			// SMOL - Small MOLecules (ligands) (we should move here)
			else if( fragtype == tmol_smol )
			{
				for ( iter_res_atom->pos_atom = 0; // N and CA
				!iter_res_atom->gend_atom();
				iter_res_atom->next_atom()   )
				{
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}
			}
			else
			{
				printf("Msg(move_dihedralMFAx): Sorry, unknown MolType %d\nForcing exit\n",fragtype);
				exit(1);
			}

			delete iter_res_atom;
			index_res++; // counts residues
		}
		delete iter;
	}
	iter_seg->clean_virtual();
	delete iter_seg;

	// Freeing memory
	for(int i=0; i<3; i++)
	{
		free( Ri[i] );
		free( Si[i] );
		free( dummy[i] );
	}
	free( Ri );
	free( Si );
	free( dummy );
	delete iter_com;
}

// It applies a CCS Normal Mode with "step" amplitude to an input structure "mol".
// "uu"-array has 3N elements, where N= number of atoms. (Valid for all atomic models)
void move_cart(Macromolecule *mol, double *uu, double step)
{
	int l;
	Tcoor pos;
	pdbIter *iter_mol;
	iter_mol = new pdbIter(mol);
	// (from cartesian Normal Modes)
	for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom() )
	{
		( iter_mol->get_atom() )->getPosition(pos); // reads "mol" atom positions
		for(l=0; l< 3; l++) // x, y, z
			// updating coord.
			pos[l] += uu[iter_mol->pos_atom * 3 + l] * step;
		( iter_mol->get_atom() )->setPosition(pos); // writes "mol"
	}
	delete iter_mol;
}

// It applies many consecutive CCS Normal Modes with their amplitudes stored
// in "amp"-array to an input structure "mol".
// "uu"-array has 3N elements, where N= number of atoms. (Valid for all atomic models)
void move_cart(Macromolecule *mol, double *uu, double *amp, int nmodes)
{
	int l,k,size;
	Tcoor pos;
	double delta[3];
	pdbIter *iter_mol;
	iter_mol = new pdbIter(mol);
	size = iter_mol->num_atom() * 3; // number of DoFs
	// (from cartesian Normal Modes)
	for(iter_mol->pos_atom=0; !iter_mol->gend_atom(); iter_mol->next_atom() )
	{
		delta[0] = 0.0;
		delta[1] = 0.0;
		delta[2] = 0.0;
		( iter_mol->get_atom() )->getPosition(pos); // reads "mol" atom positions
		for(k=0; k<nmodes; k++)
		{
			delta[0] += uu[size*k + iter_mol->pos_atom * 3 ] * amp[k];
			delta[1] += uu[size*k + iter_mol->pos_atom * 3 + 1] * amp[k];
			delta[2] += uu[size*k + iter_mol->pos_atom * 3 + 2] * amp[k];
		}
		pos[0] += delta[0];
		pos[1] += delta[1];
		pos[2] += delta[2];
		( iter_mol->get_atom() )->setPosition(pos); // writes "mol"
	}
	delete iter_mol;
}

//// Moves atoms given an Internal Coordinates Normal Mode and a Step (CA atomic model)
//// Using V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
//void move_vwMCAx(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1 )
//{
//	double rk[3],temp;
//	int l;
//	Atom *atom;
//	Tcoor pos;
//	pdbIter *iter_atoms = new pdbIter( mol ); // atom iterator
//
//	for(iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
//	{ // screen atoms
//		atom = iter_atoms->get_atom();
//		atom->getPosition(pos); // get position
//		rk[0] = pos[0];
//		rk[1] = pos[1];
//		rk[2] = pos[2];
//		for(l=0; l< size; l++) // screen Dihedrals
//		{
//			temp = evec[index + l] * step;
//			if( iter_atoms->pos_atom < body1[l][0] )
//			{
//				// i-der for k-atom --> = v + (w x r) (vectorial product)
//				pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
//				pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
//				pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
//			}
//			else // if body 2 atom
//			{
//				// i-der for k-atom  --> = v + (w x r) (vectorial product)
//				pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
//				pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
//				pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
//			}
//		}
//		atom->setPosition(pos); // set new position
//	}
//
//	delete iter_atoms;
//}
//
//// Multi-thread routine to move atoms given an Internal Coordinates Normal Mode and a Step (CA atomic model)
//// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
//void *move_vwMCAx_thread(void *threadarg)
//{
//	double rk[3],temp;
//	int l;
//	Atom *atom;
//	Tcoor pos;
//	pdbIter *iter_atoms;
//
//	// Loading data from "thread argument"
////	Macromolecule *mol = ((moveVW_data *) threadarg)->mol;
//	iter_atoms = ((moveVW_data *) threadarg)->iter;
//	double step = ((moveVW_data *) threadarg)->step;
//	double *evec = ((moveVW_data *) threadarg)->evec;
//	int index = ((moveVW_data *) threadarg)->index;
//	int size = ((moveVW_data *) threadarg)->size;
//	double ***V = ((moveVW_data *) threadarg)->V;
//	double ***W = ((moveVW_data *) threadarg)->W;
//	int **body1 = ((moveVW_data *) threadarg)->body1;
//	int first = ((moveVW_data *) threadarg)->first;
//	int last = ((moveVW_data *) threadarg)->last;
//
////	pdbIter *iter_atoms = new pdbIter( mol ); // atom iterator
//	for(iter_atoms->pos_atom = first; iter_atoms->pos_atom <= last; iter_atoms->next_atom())
//	{ // screen atoms
//		atom = iter_atoms->get_atom();
//		atom->getPosition(pos); // get position
//		rk[0] = pos[0];
//		rk[1] = pos[1];
//		rk[2] = pos[2];
//		for(l=0; l< size; l++) // screen Dihedrals
//		{
//			temp = evec[index + l] * step;
//			if( iter_atoms->pos_atom < body1[l][0] )
//			{
//				// i-der for k-atom --> = v + (w x r) (vectorial product)
//				pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
//				pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
//				pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
//			}
//			else // if body 2 atom
//			{
//				// i-der for k-atom  --> = v + (w x r) (vectorial product)
//				pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
//				pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
//				pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
//			}
//		}
//		atom->setPosition(pos); // set new position
//	}
////	delete iter_atoms;
//	pthread_exit(NULL);
//}

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Using V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void move_VW(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1, int model )
{
	double rk[3],temp;
	int l;
	Atom *atom;
	Tcoor pos;
	pdbIter *iter_atoms = new pdbIter( mol, true, true, true, true ); // atom iterator

	if( model == 1 || model == 2 ) // C5/HA models
		for(iter_atoms->pos_atom=0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
		{ // screen atoms
			atom = iter_atoms->get_atom();
			atom->getPosition(pos); // get position
			rk[0] = pos[0];
			rk[1] = pos[1];
			rk[2] = pos[2];
			for(l=0; l< size; l++) // screen Dihedrals
			{
				temp = evec[index + l] * step;
				if( (iter_atoms->pos_atom < body1[l][0] || (body1[l][1]>=0 && iter_atoms->pos_atom > body1[l][1])) && (iter_atoms->pos_atom != body1[l][2]) )
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
					pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
					pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
					pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
					pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
				}
			}
			atom->setPosition(pos); // set new position
		}
	else // CA/CA3 models
		for(iter_atoms->pos_atom = 0; !iter_atoms->gend_atom(); iter_atoms->next_atom())
		{ // screen atoms
			atom = iter_atoms->get_atom();
			atom->getPosition(pos); // get position
			rk[0] = pos[0];
			rk[1] = pos[1];
			rk[2] = pos[2];
			for(l=0; l< size; l++) // screen Dihedrals
			{
				temp = evec[index + l] * step;
				if( iter_atoms->pos_atom < body1[l][0] )
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
					pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
					pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
					pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
					pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
				}
			}
			atom->setPosition(pos); // set new position
		}

	delete iter_atoms;
}

// Multi-thread routine to move atoms given an Internal Coordinates Normal Mode and a Step
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void *move_VWthread(void *threadarg)
{
	double rk[3],temp;
	int l,model;
	Atom *atom;
	Tcoor pos;
	pdbIter *iter_atoms;

	// Loading data from "thread argument"
	model = ((moveVW_data *) threadarg)->model;
	iter_atoms = ((moveVW_data *) threadarg)->iter;
	double step = ((moveVW_data *) threadarg)->step;
	double *evec = ((moveVW_data *) threadarg)->evec;
	int index = ((moveVW_data *) threadarg)->index;
	int size = ((moveVW_data *) threadarg)->size;
	double ***V = ((moveVW_data *) threadarg)->V;
	double ***W = ((moveVW_data *) threadarg)->W;
	int **body1 = ((moveVW_data *) threadarg)->body1;
	int first = ((moveVW_data *) threadarg)->first;
	int last = ((moveVW_data *) threadarg)->last;

	if( model == 1 || model == 2 ) // C5/HA models
		for(iter_atoms->pos_atom = first; iter_atoms->pos_atom <= last; iter_atoms->next_atom())
		{ // screen atoms
			atom = iter_atoms->get_atom();
			atom->getPosition(pos); // get position
			rk[0] = pos[0];
			rk[1] = pos[1];
			rk[2] = pos[2];
			for(l=0; l< size; l++) // screen Dihedrals
			{
				temp = evec[index + l] * step;
				if( (iter_atoms->pos_atom < body1[l][0] || (body1[l][1]>=0 && iter_atoms->pos_atom > body1[l][1])) && (iter_atoms->pos_atom != body1[l][2]) )
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
					pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
					pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
					pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
					pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
				}
			}
			atom->setPosition(pos); // set new position
		}
	else
		for(iter_atoms->pos_atom = first; iter_atoms->pos_atom <= last; iter_atoms->next_atom())
		{ // screen atoms
			atom = iter_atoms->get_atom();
			atom->getPosition(pos); // get position
			rk[0] = pos[0];
			rk[1] = pos[1];
			rk[2] = pos[2];
			for(l=0; l< size; l++) // screen Dihedrals
			{
				temp = evec[index + l] * step;
				if( iter_atoms->pos_atom < body1[l][0] )
				{
					// i-der for k-atom --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][0][0] + W[l][0][1] * rk[2] - W[l][0][2] * rk[1]);
					pos[1] += temp * (V[l][0][1] + W[l][0][2] * rk[0] - W[l][0][0] * rk[2]);
					pos[2] += temp * (V[l][0][2] + W[l][0][0] * rk[1] - W[l][0][1] * rk[0]);
				}
				else // if body 2 atom
				{
					// i-der for k-atom  --> = v + (w x r) (vectorial product)
					pos[0] += temp * (V[l][1][0] + W[l][1][1] * rk[2] - W[l][1][2] * rk[1]);
					pos[1] += temp * (V[l][1][1] + W[l][1][2] * rk[0] - W[l][1][0] * rk[2]);
					pos[2] += temp * (V[l][1][2] + W[l][1][0] * rk[1] - W[l][1][1] * rk[0]);
				}
			}
			atom->setPosition(pos); // set new position
		}
	pthread_exit(NULL);
}

// Moves atoms given an Internal Coordinates Normal Mode and a Step (C5 & HA atomic models)
// Using V/W-arrays and multi-threaded parallel schedule (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void move_VWpar(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1, int nthreads, int model)
{
	bool debug = false;
	moveVW_data *threads_data;
	pthread_t *threads;
	int i,rc,num_atoms,current,init=0;
	void *status;

	num_atoms = mol->get_num_atoms();

	// Allocating threads data
	if( !(threads_data = (moveVW_data *) malloc(sizeof(moveVW_data) * nthreads)) )
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
		threads_data[i].model = model;
		threads_data[i].iter = (pdbIter *) new pdbIter(mol,true,true,true,true); // pdbIter seems not to be "thread-safe"...
		threads_data[i].step = step;
		threads_data[i].evec = evec;
		threads_data[i].index = index;
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
		rc = pthread_create(&threads[i], NULL, move_VWthread, (void *)	&threads_data[i]);
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

	// Deleting iterators
	for(i = 0; i < nthreads; i++)
	{
		threads_data[i].iter->clean_virtual();
		delete threads_data[i].iter;
	}
	free(threads_data);
	free(threads);
}


