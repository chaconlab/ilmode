/*************************************************************************
 *                     LIBRARY: libnma_hessian                           *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Library to compute Hessian Matrices using several methods:          *
 *   	-In Internal Coordinate Space (ICS):                             *
 *   		-naive K-matrix (O^4),                                       *
 * 			-naive V/W-arrays (memmory efficient) (O^4),                 *
 * 			-fast and memmory efficient Go's method (O^2).               *
 *       -In Cartesian Coordinate Space:                                 *
 *           -Elastic Network Model hessian (also mass-weighted)         *
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
#include <libnma_deriv.h>
#include <libnma_misc.h>  // this should be removed when seg_props() is passed as an argument
#include <libnma_hessian.h>
#include <libnma_io.h>
#include <unistd.h>

#define DEVEL

// Computing Tij element "on the fly"
// (ec.21) Noguti & Go 1983 pp.3685-90
// (T 6x6 matrix must be already allocated!)
void calcTij( double **T, int *index, int num, twid *decint, float *coord )
{
	int k,l,ki,li,in,i,j,m,n;
	double r_alpha[3];
	double r_beta[3];
	double v[6];
	double d;

	for(i=0; i<6; i++)
		for(j=0; j<6; j++)
			T[i][j] = 0.0; // T initialization

	for(i=1; i < num; i++) // screens inter-unit contacts (k vs. l)
	{
		in = index[i];
		k = decint[ in ].k;
		l = decint[ in ].l;
		ki = 3*k;
		li = 3*l;

		// k-atom position retrieval (alpha)
		r_alpha[0] = coord[ki];
		r_alpha[1] = coord[ki + 1];
		r_alpha[2] = coord[ki + 2];

		// l-atom position retrieval (beta)
		r_beta[0] = coord[li];
		r_beta[1] = coord[li + 1];
		r_beta[2] = coord[li + 2];

		// D-Matrix
		v[0] = r_alpha[1] * r_beta[2] - r_alpha[2] * r_beta[1];
		v[1] = r_alpha[2] * r_beta[0] - r_alpha[0] * r_beta[2];
		v[2] = r_alpha[0] * r_beta[1] - r_alpha[1] * r_beta[0];
		v[3] = r_alpha[0] - r_beta[0];
		v[4] = r_alpha[1] - r_beta[1];
		v[5] = r_alpha[2] - r_beta[2];
		d = decint[in].C / pow(decint[in].d,2);

		for(n=0;n<6;n++) // Diadic product
			for(m=0;m<6;m++)
				T[n][m] += d * v[n] * v[m]; // Storing S-matrix
	}
}

// Computing Tij element "on the fly" from "sorted ipas" array.
// (ec.21) Noguti & Go 1983 pp.3685-90
// (T 6x6 matrix must be already allocated!)
void calcTij_new(double **T, twid **sipas, int *index, float *coord, int nipa)
{
	int ki,li,m,n;
	double r_alpha[3],r_beta[3],v[6];
	double d;

	for(m=0; m<6; m++)
		for(n=0; n<6; n++)
			T[m][n] = 0.0; // T initialization

	// This screens inter-unit contacts (k vs. l) for the units pair (i,j)
//	do {	// while "sipas"'s belong to the same interacting pair of units: i,j
do
//	for( ; (*index) == 0 || ( (*index) < nipa && (sipas[*index]->i == sipas[(*index)-1]->i) && (sipas[*index]->j == sipas[(*index)-1]->j) ); (*index)++ )
	{
		ki = 3*sipas[*index]->k;
		li = 3*sipas[*index]->l;

		// k-atom position retrieval (alpha)
		r_alpha[0] = coord[ki];
		r_alpha[1] = coord[ki + 1];
		r_alpha[2] = coord[ki + 2];

		// l-atom position retrieval (beta)
		r_beta[0] = coord[li];
		r_beta[1] = coord[li + 1];
		r_beta[2] = coord[li + 2];

		// D-Matrix
		v[0] = r_alpha[1] * r_beta[2] - r_alpha[2] * r_beta[1];
		v[1] = r_alpha[2] * r_beta[0] - r_alpha[0] * r_beta[2];
		v[2] = r_alpha[0] * r_beta[1] - r_alpha[1] * r_beta[0];
		v[3] = r_alpha[0] - r_beta[0];
		v[4] = r_alpha[1] - r_beta[1];
		v[5] = r_alpha[2] - r_beta[2];
		d = sipas[*index]->C / pow(sipas[*index]->d,2);

		for(n=0;n<6;n++) // Diadic product
			for(m=0;m<6;m++)
				T[n][m] += d * v[n] * v[m]; // Storing S-matrix

		(*index) ++; // next "sipa"
	}
	while( (*index) < nipa && (sipas[*index]->i == sipas[(*index)-1]->i) && (sipas[*index]->j == sipas[(*index)-1]->j) );
//	while( (sipas[*index]->i == sipas[(*index)-1]->i) && (sipas[*index]->j == sipas[(*index)-1]->j) );
	// while the interacting pair of atoms belong to the same pair of units
	// (note that *index should never be equal to nipa, otherwise overflow-bug)
}


// Cartesian Coordinate Space (CCS) Hessian Building
// (Mass-weighted CCS if mass!=NULL; i.e. including kinetic energy matrix into the hessian)
// Triangular matrix storage and Hessian memory allocation (if *p_hess_matrix==NULL)
void hessianC(double **p_hess_matrix,double *mass,float *coord,double *dist_matrix,double *cont_matrix, int size)
{
	double *hess_matrix;
	double r[3], rsqu, dummy;
	int i, j, m, n, index, num_atoms;
	num_atoms = size/3;

	// Hessian Matrix memory allocation
	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hessianC): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
		printf("Msg(hessianC): Hessian Matrix memory allocated (rank=%d, %.3f Mb)\n",size,(sizeof(double)*size*(size+1)/2)/1e6);
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian matrix initialization
	for( i = 0; i < size*(size+1)/2; i++ )
		hess_matrix[i] = 0.0;

	// Off-diagonal terms
	for ( i = 0; i < num_atoms; ++i )
		for ( j = i+1; j < num_atoms; ++j )
		{
			r[0] = coord[i * 3] - coord[j * 3];
			r[1] = coord[i * 3 + 1] - coord[j * 3 + 1];
			r[2] = coord[i * 3 + 2] - coord[j * 3 + 2];
			index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2;
			rsqu = dist_matrix[index] * dist_matrix[index];
			for ( m = 0; m < 3; ++m )
				for ( n = 0; n < 3; ++n )
				{
					// if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
					// i + j*(j+1)/2 --> for 0<=i<=j<size (begining with 0)
					hess_matrix[3*i+m + (3*j+n)*(3*j+n+1)/2] = -cont_matrix[index] * r[m] * r[n] / rsqu;
				}
		}

	// On-diagonal terms
	for ( i = 0; i < num_atoms; ++i )
		for ( j = i+1; j < num_atoms; ++j )
		{
			r[0] = coord[i * 3] - coord[j * 3];
			r[1] = coord[i * 3 + 1] - coord[j * 3 + 1];
			r[2] = coord[i * 3 + 2] - coord[j * 3 + 2];

			index = num_atoms * i + j - 1 - i * ( i + 3 ) / 2;
			rsqu = dist_matrix[index] * dist_matrix[index];

			for ( m = 0; m < 3; ++m )
				for ( n = m; n < 3; ++n ) // <-- diagonal sub-matrices are symmetric,
				{						  //     (watch out triangular packed storage!)
					dummy = cont_matrix[index] * r[m] * r[n] / rsqu;
					hess_matrix[3*i+m + (3*i+n)*(3*i+n+1)/2] += dummy; // i-term
					hess_matrix[3*j+m + (3*j+n)*(3*j+n+1)/2] += dummy; // j-term
				}
		}

	// Including KINETIC ENERGY MATRIX into HESSIAN
	// (i.e. transforming the generalized eigenvalue problem into the standard one)
	// (in cartesian space, simply... mass-weighting)
	// multiply hess_matrix by mass matrix, i.e., form F = M^(-1/2).H.M^(-1/2)

	if( mass != NULL )
	{
		// Off-diagonal
		for(i=0;i<num_atoms;++i)
			for(j=i+1;j<num_atoms;++j)
				for(m=0;m<3;++m)
					for(n=0;n<3;++n)
						// i + j*(j+1)/2 --> for 0<=i<=j<size (begining with 0)
						hess_matrix[3*i+m + (3*j+n)*(3*j+n+1)/2] /= sqrt(mass[i]*mass[j]);

		// On-diagonal
		for(i=0;i<num_atoms;++i)
			for(m=0;m<3;++m)
				for(n=m;n<3;++n)
					// i + j*(j+1)/2 --> for 0<=i<=j<size (begining with 0)
					hess_matrix[3*i+m + (3*i+n)*(3*i+n+1)/2] /= mass[i];
	}

}

// Dihedral Angle Space (DAS) Hessian Building (Multi-chain/Protein/DNA/RNA/SMOL)
// Adds torsional springs to a previously built hessian matrix (triangular packing)
// (See second term of ec.5 from: Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471)
void hessianMDHx(Macromolecule *mol,tri *props,floating *hess,double omega,int size,bool *fix,bool *addrot)
{
	pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	pdbIter *iter_frag;
	Segment *seg;
	int index_frag=0,j=0;
	long int j2=0,i=0; // Needed for huge systems!
	double min,aux;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

//	double avg=0.0,sig=0.0;

	// Ma's ec.5 condition? Ambiguous. We understand this:
	min = hess[ 0 ];
	for(i = 0; i < size; i++)
	{
		aux = hess[ i + i*(i+1)/2 ];
//		avg += aux;
		if( min > aux && aux > 0.0000001) // neccessary when segments are not linked
			min = aux;
	}
//	avg/=size;

//	for(i = 0; i < size; i++)
//		sig += pow( hess[ i + i*(i+1)/2 ] - avg, 2 );
//	sig = sqrt(sig/size);

	// Applying "omega" torsional-stiffness
	omega *= 3*min;
//	fprintf(stderr,"Msg(hessianMDHx): min= %f  omega= %f  avg= %f  sig= %f\n",min,omega,avg,sig);

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Taking into account TRASLATIONS (there are always 3 translations between segments)
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					j2++; // No torsional contribution is added for Rotational/Traslational DoFs
				j++; // adding traslational degree of freedom
			}

			// Taking into account 3 ROTATIONS
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
						j2++; // No rotational contribution is added for Rotational/Traslational DoFs
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // adding rotational DoF
				}
		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			for(i = 0; i < props[index_frag].nan; i++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// Only diagonal terms!
					hess[ j2 + j2*(j2+1)/2 ] += (floating) omega; // dV/dTheta,j * dV/dTheta,j
					j2++;
				}
				j++;
			}
			index_frag++;
		}
		delete iter_frag;
	}
	iter_seg->clean_virtual();
	delete iter_seg;

	// Checking
	if(size != j2)
	{
		printf("Msg(hessianMDHx): Internal coordinates mismatch! (size= %d  j= %d  j2= %d)\n",size,j,j2);
		exit(1);
	}
}

// Dihedral Angle Space (DAS) Hessian Building (Multi-chain/Protein/DNA/RNA/SMOL)
// Adds torsional springs to a previously built hessian matrix from file (triangular packing)
// (See second term of ec.5 from: Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471)
void hessianMDHxHD(Macromolecule *mol,tri *props,char *file,double omega,bool *fix,bool *addrot)
{
	pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	pdbIter *iter_frag;
	Segment *seg;
	int index_frag=0,j=0;
	long int j2=0,i=0; // Needed for huge systems!
	double min,aux,hess;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;
	int size,offset;
	FILE *f_file;

	if( !(f_file = fopen(file, "r+") ) ) // Open "read and write"
	{
		printf("Msg(hessianMDHxHD): I'm sorry, unable to open %s file\nForcing exit\n",file);
		exit(1);
	}
	fread(&aux,sizeof(double),1,f_file); // getting size
	size = (int) aux; // casting double into integer

	// TO-DO LIST:
	// UN POCO MAS "COHERENTE"... "Se busca el minimo en todas y solo se aplica a los DH" ¿?
	//
	// Ma's ec.5 condition? Ambiguous. We understand this:
	// min = hess[ 0 ];
	fread(&hess,sizeof(double),1,f_file); // getting Hessian element
	min = hess; // Hessian's diagonal elements are always positive
	fseek( f_file, -1*sizeof(double), SEEK_CUR ); // one element back (because we've just read one)
	for(i = 0; i < size; i++)
	{
		fread(&hess,sizeof(double),1,f_file); // getting Hessian element
//		fprintf(stderr,"i= %ld  --> hess_diag= %f\n",i,(float)hess);
		aux = hess;
		if( min > aux && aux > 0.0000001) // necessary when segments are not linked
			min = aux;
		fseek( f_file, (size-i-1)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
	}

	// Applying "omega" torsional-stiffness
	omega *= 3*min;
//	fprintf(stderr,"min= %f --> omega= %f\n",min,omega);

	fseek( f_file, 1*sizeof(double), SEEK_SET ); // positioning at the file start

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			// Taking into account TRASLATIONS (there are always 3 translations between segments)
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
//					fread(&hess,sizeof(double),1,f_file); // getting Hessian element
//					fprintf(stderr,"j2= %ld  --> hess_diag= %f\n",j2,(float)hess);
//					fseek( f_file, (size-j2-1)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
					fseek( f_file, (size-j2)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
					j2++; // No torsional contribution is added for Rotational/Traslational DoFs
				}
				j++; // adding traslational degree of freedom
			}

			// Taking into account 3 ROTATIONS
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{
//						fread(&hess,sizeof(double),1,f_file); // getting Hessian element
//						fprintf(stderr,"j2= %ld  --> hess_diag= %f\n",j2,(float)hess);
//						fseek( f_file, (size-j2-1)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
						fseek( f_file, (size-j2)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
						j2++; // No rotational contribution is added for Rotational/Traslational DoFs
					}
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // adding rotational DoF
				}
		}

		// Screen ALL fragments
		for( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			for(i = 0; i < props[index_frag].nan; i++)
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// Only diagonal terms!
					fread(&hess,sizeof(double),1,f_file); // getting Hessian element
//					fprintf(stderr,"j2= %ld  --> hess_diag= %f  ",j2,(float)hess);
					hess += omega;
//					fprintf(stderr," --> %f\n",(float)hess);
					fseek( f_file, -1*sizeof(double), SEEK_CUR ); // one element back (because we've just read one)
					fwrite(&hess,sizeof(double),1,f_file); // setting Hessian element
					fseek( f_file, (size-j2-1)*sizeof(double), SEEK_CUR ); // positioning next diagonal element
					j2++;
				}
				j++;
			}
			index_frag++;
		}
		delete iter_frag;
	}
	iter_seg->clean_virtual();
	delete iter_seg;

	fclose(f_file);

	// Checking
	if(size != j2)
	{
		printf("Msg(hessianMDHx): Internal coordinates mismatch! (size= %d  j= %d  j2= %d)\n",size,j,j2);
		exit(1);
	}
}


// *Pending: Do it like CA-one!
// Fast Hessian Matrix computation for:  Multi-Chain & FULL-ATOM & Protein/RNA/DNA/SMOL
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing) (Efficient memory allocation)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
// Noguti & Go (1983) pp. 3685-90
void hessianMFAx(Macromolecule *mol, twid *decint,int nipa,int size,float *coord,floating **p_hess_matrix,tri *props, int *unat, int type, int model, bool *fix, bool *addrot)
{
	bool debug = false;
	bool debug_mem = false;
	FILE *out=stderr;

	double prod,prod1;
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	int prin, fin,i, j, j2, m, k1, k2, buff;
	double r_alpha[3],r_beta[3],r[3],e[3];
	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
	double c,d,temp;
	double v[6];
	int resn,num_res,num_atoms,num_seg;
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	int index_res = 0;
	floating *hess_matrix;
	Htimer ht_timer; // timer
	int CBindex;
	TMOL fragtype;

	// **************************************
	// Storing==> (ea, ea x ra) (== (eb, eb x rb) )
	// **************************************
	double **erx;
	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	for(int i=0;i<6;i++)
		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]

	bool *ischi; // CHIs are different...
	bool *istra; // Traslations are different...
	bool *isrot; // Rotations are different...

	ischi = (bool *) malloc( sizeof(bool) * size ); // CHIs are different...
	istra = (bool *) malloc( sizeof(bool) * size ); // Traslations are different...
	isrot = (bool *) malloc( sizeof(bool) * size ); // Rotations are different...

	int *undh; // returns the first unit-index on the left side of the dihedral
	int un_index = 0;
	undh = (int *) malloc( sizeof(int) * size );

	int *dhup; // returns the closest dihedral-index towards the top Uab matrix part
	int *dhright; // returns the closest dihedral-index towards the right Uab matrix part
	int *cc; // c(a)=d, array in Noguti & Go (1983) pp. 3685-90
	int dhup_index = 0;
	int dhright_index = size-1;
	dhup = (int *) malloc( sizeof(int) * size );
	dhright = (int *) malloc( sizeof(int) * size );
	cc = (int *) malloc( sizeof(int) * size );
	bool last_chi = false;
	bool last_psi = false;
	bool any_rotrans = false;
	int indexbase;

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	j = 0; // current dihedral index
	j2 = 0; // mobile dihedral index
	k1 = 0;
	k2 = props[0].nat - 1; /* 1st and last indices of atoms of residue i */ // Residue 0 ?
	int n_inter; // number of mobile intersegment variables
	int available_inter,seg_atoms;
	int num_units=0;
	int n_fixed=0; // counter of fixed internal coordinates
	int fix_offset=0; // unit-offset (joins units)
	int n_bb=0;
	int n_chi=0;
	int n_bb2=0;
	int natom=0;
	int seg_atoms_old;
	int seg_atoms_current = 0; // should be initialized
	bool may_have_rot = false;

	// Screen ALL residues
	Segment * seg;
	Atom * atom;
	// Screening segments
	iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(debug)
			fprintf(out,"\nProcessing segment %d (%d): k1 %d k2 %d\n",iter_seg->pos_segment,num_res,k1,k2);

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms_current;
		seg_atoms_current = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			n_inter=0; // number of mobile intersegment variables (should be knowk "a priori")
			// Counting TRASLATIONS
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j+axis] ) // Check whether current variable it's mobile
					n_inter++;
			}
			// Counting ROTATIONS
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=3; axis<6; axis++)
				{
					if( fix == NULL || fix[j+axis] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
						n_inter++;
				}

			any_rotrans = false;
			// Adding 3 TRASLATIONS (there are always 3 translations between segments)
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					any_rotrans = true;
					// ea (Phi)  (0)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					// ea x ra (Psi)  (-gamma_v)
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;
					erx[3+axis][j2] = -1.0;

					ischi[j2] = false;
					istra[j2] = true;
					isrot[j2] = false;
					undh[j2] = un_index;
					dhup[j2] = dhup_index-1;
					dhright[j2] = dhup_index+n_inter;

					if(debug)
						fprintf(out,"TRANS j= %d  j2= %d  k1= %d  k2=%d  dhup[j]= %d  dhright[j]= %d  erx= %f %f %f %f %f %f\n",j,j2,k1,k2,dhup[j2],dhright[j2],erx[0][j2],erx[1][j2],erx[2][j2],erx[3][j2],erx[4][j2],erx[5][j2]);

					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}

				j++; // adding translational degree of freedom
			} // 3 trans added

			// Adding 3 ROTATIONS
			if( (seg_atoms_current != 1 && may_have_rot) || (addrot != NULL && addrot[iter_seg->pos_segment]))
				for(int axis=0; axis<3; axis++)
				{
					if( fix == NULL || fix[j] || addrot[iter_seg->pos_segment] ) // Check whether current variable it's fixed
					{
						any_rotrans = true;
						// ea (Phi)  (delta_v)
						erx[0][j2] = 0;
						erx[1][j2] = 0;
						erx[2][j2] = 0;
						erx[axis][j2] = 1.0;
						// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
						erx[3][j2] = 0;
						erx[4][j2] = 0;
						erx[5][j2] = 0;

						ischi[j2] = false;
						istra[j2] = false;
						isrot[j2] = true;
						undh[j2] = un_index;
						dhup[j2] = dhup_index-1;
						dhright[j2] = dhup_index+n_inter;

						if(debug)
							fprintf(out,"ROT j= %d  j2= %d  k1= %d  k2=%d  dhup[j]= %d  dhright[j]= %d  erx= %f %f %f %f %f %f\n",j,j2,k1,k2,dhup[j2],dhright[j2],erx[0][j2],erx[1][j2],erx[2][j2],erx[3][j2],erx[4][j2],erx[5][j2]);

						j2++;
					}
					else
					{
						fix_offset++;
						n_fixed++;
					}

// MON: I think this is not used... check it!!! (otherwise "j" and "j2" would unbalance...) It seems to be always true...
					if(seg_atoms_current != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
						j++; // adding rotational DoF
				} // 3 rots added

			// Whether any of the 6D inter-segment variables are mobile
			if(any_rotrans)
			{
				// Defining rigid units
				n_bb += n_chi + 1; // +number of chis of preceding unit
				n_chi = 0; // parent unit breaks!

				un_index++; // count units (due to segment ending)
				dhup[j2-1] = dhup_index; // the last Rot it's different!!! (needed)
				dhright[j2-n_inter] = dhup_index; // the last Rot it's different!!! (needed)
				dhup_index += n_inter; // jumps "n_inter" internal coordinates (xT + xR)
			}
		}

		// Screen ALL fragments
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			if(debug)
				fprintf(out,"CURRENT FRAGMENT (tmol= %d): %s  k1= %d  k2= %d  index_res= %d\n",fragtype,res->getName(),k1,k2,index_res);

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{
				if(model==1)
					CBindex = 3;
				else
					CBindex = 4;

				unat[k1] = n_bb; // bb unit index

				// NOT FIRST, NOT PRO -> PHI
				// ("PHI-bond")
				if ( iter_frag->pos_fragment != 0 && resn != PRO )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						//			  printf ("%d Not first not PRO -> PHI \n",j);
						iter->pos_atom = k1; // first residue atomic index (NH pseudo-atom)
						// ("k2" is updated below (for the 1st not PRO -> PHI)

						// y_lambda (NH pos 0)
						e[0] = -coord[k1 * 3];
						e[1] = -coord[k1 * 3 + 1];
						e[2] = -coord[k1 * 3 + 2];
						// CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// e_lambda
						e[0] += r[0]; // NH --> CA
						e[1] += r[1];
						e[2] += r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
						// e x ta
						// v[1] [0] = e[1] * ta[2] - e[2] * ta[1]; // vectorial product (from calcoef)
						// v[1] [1] = e[2] * ta[0] - e[0] * ta[2];
						// v[1] [2] = e[0] * ta[1] - e[1] * ta[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit

						//					cc[j2] = n_bb-1;
						n_chi = 0; // parent unit breaks!

						if(debug)
							fprintf(out,"PHI j= %d  j2= %d  k1= %d  k2=%d  dhup[j2]= %d  dhright[j2]= %d  erx= %f %f %f %f %f %f\n",j,j2,k1,k2,dhup[j2],dhright[j2],erx[0][j2],erx[1][j2],erx[2][j2],erx[3][j2],erx[4][j2],erx[5][j2]);

						j2++;
					}
					else
					{
						n_fixed++;
						if(debug)
							fprintf(out,"PHI j= %d --> No (n_fixed= %3d) j= %d\n",j,n_fixed,j);
					}

					j++; // der[] index
				}  // NOT FIRST NOT PRO

				unat[k1+1] = n_bb; // CA-atom
				// if needed, will be overwritten further (due to ALA and PRO)
				if(resn == ALA)
					unat[k1+CBindex] = n_bb;
				if(resn == PRO)
				{
					unat[k1+CBindex] = n_bb;
					if(model==1) // 3BB2R
						unat[k1+4] = n_bb; // R- (O-atom)
					else
					{	// Full-Atom
						unat[k1+CBindex+1] = n_bb;
						unat[k1+CBindex+2] = n_bb;
					}
				}
				for(int n=CBindex; n<props[index_res].nat; n++)
					unat[k1+n] = n_bb;

				//  LATERAL CHAIN-->CHI
				// 3 dihedrals (normal residue) or 2 dihedrals (ending residue)
				if( type == 2 )
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frag->pos_fragment==0 ||
											(iter_frag->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
						{
							// CA pos 1
							r[0] = coord[(k1+1) * 3];
							r[1] = coord[(k1+1) * 3 + 1];
							r[2] = coord[(k1+1) * 3 + 2];
							// CB pos 3
							e[0] = coord[(k1+CBindex) * 3];
							e[1] = coord[(k1+CBindex) * 3 + 1];
							e[2] = coord[(k1+CBindex) * 3 + 2];
							// e_lambda ==> CA-->CB unit vector
							e[0] -= r[0];
							e[1] -= r[1];
							e[2] -= r[2];

							temp = sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
							for(m=0;m<3;m++)
								e[m] /= temp; // unit vector normalization
							// ea
							erx[0][j2] = e[0];
							erx[1][j2] = e[1];
							erx[2][j2] = e[2];
							// ea x ra
							erx[3][j2] = e[1] * r[2] - e[2] * r[1];
							erx[4][j2] = e[2] * r[0] - e[0] * r[2];
							erx[5][j2] = e[0] * r[1] - e[1] * r[0];

							ischi[j2] = true; // it is chi

							undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
							un_index++; // count units
							istra[j2] = false;
							isrot[j2] = false;
							dhup[j2] = dhup_index;
							dhright[j2] = dhup_index;
							dhup_index++; // counts "up" dihedrals (1st body)

							// Defining rigid units
							if(n_bb == n_bb2) // if both belong to the same unit
								n_chi++;
							else
							{
								n_chi = 1;
								n_bb2 = n_bb;
							}

							for(int n=CBindex; n<props[index_res].nat; n++)
								unat[k1+n] = n_bb + n_chi;

							if(debug)
								fprintf(out,"CHI j= %d  j2= %d  k1= %d  k2=%d  dhup[j2]= %d  dhright[j2]= %d  erx= %f %f %f %f %f %f\n",j,j2,k1,k2,dhup[j2],dhright[j2],erx[0][j2],erx[1][j2],erx[2][j2],erx[3][j2],erx[4][j2],erx[5][j2]);

							// ESTO SOBRA, CHECK PLEASE!
							if(fix != NULL)
								last_chi = true; // last chi is mobile
							j2++;
						}
						else
						{
							n_fixed++;
							if(debug)
								fprintf(out,"CHI j= %d --> No (n_fixed= %3d)\n",j,n_fixed);
						}

						j++; // der[] index
					}

				// NOT LAST RESIDUE--> PSI
				// ("PSI-bond")
				if ( iter_frag->pos_fragment != num_res - 1 || model==2 ) // Full-Atom allways has PSI
				{
					//			 fprintf(out,"%d Not last residue -> PSI\n",j);
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					{
						// get CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
						e[0] = coord[(k1+2) * 3];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= r[0];
						e[1] -= r[1];
						e[2] -= r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit
						cc[j2] = n_bb-1;
						n_chi = 0; // parent unit brokes!

						if(debug)
							fprintf(out,"PSI j= %d  j2= %d  k1= %d  k2=%d  dhup[j2]= %d  dhright[j2]= %d  erx= %f %f %f %f %f %f\n",j,j2,k1,k2,dhup[j2],dhright[j2],erx[0][j2],erx[1][j2],erx[2][j2],erx[3][j2],erx[4][j2],erx[5][j2]);

						j2++;
					}
					else
					{
						n_fixed++;
						if(debug)
							fprintf(out,"PSI j= %d--> No (n_fixed= %3d)\n",j,n_fixed);
					}

					j++; // der[] dihedral index
				}  // NOT LAST

				if(model==1) // 3BB2R
					unat[k1+2] = n_bb; // just C
				else // Full-Atom
					for(int n=2; n<=3; n++)
						unat[k1+n] = n_bb; // C and O

				// Check OXT
				iter->pos_atom = k2; // OXT position
				if( strcmp(iter->get_atom()->getName()," OXT ")==0 )
					unat[k2] = unat[k1+2]; // OXT is in the same unit as C

			} // END PROTEIN FRAGMENT
			// RNA FRAGMENT
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA )
			{
				// ********************************************************************************************
				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)

				// Updating "unat" for "alpha" atoms
				unat[k1] = n_bb; // P bb unit index
				unat[k1+1] = n_bb; // O1P bb unit index
				unat[k1+2] = n_bb; // O2P bb unit index

				// ALPHA (bond between P and O5*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					r[0] = coord[k1*3];
					r[1] = coord[k1*3 + 1];
					r[2] = coord[k1*3 + 2];
					e[0] -= r[0]; // P --> O5*
					e[1] -= r[1];
					e[2] -= r[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END ALPHA

				// Updating "unat" for "beta" atoms
				unat[k1+3] = n_bb; // O5*

				// BETA (bond between O5* and C5*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					r[0] = coord[(k1+3)*3];
					r[1] = coord[(k1+3)*3 + 1];
					r[2] = coord[(k1+3)*3 + 2];
					e[0] -= r[0]; // O5* --> C5*
					e[1] -= r[1];
					e[2] -= r[2];
					// fprintf(out,"BETA M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END BETA

				// Updating "unat" for "gamma" atoms
				unat[k1+4] = n_bb; // C5*

				// GAMMA (bond between C5* and C4*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					r[0] = coord[(k1+4)*3];
					r[1] = coord[(k1+4)*3 + 1];
					r[2] = coord[(k1+4)*3 + 2];
					e[0] -= r[0]; // C5* --> C4*
					e[1] -= r[1];
					e[2] -= r[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END GAMMA

				// Updating "unat" for "chi" atoms
				for(int n=5; n<props[index_res].nat; n++)
					unat[k1+n] = n_bb;

				//  LATERAL CHAIN-->CHI
				if( type == 2 )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
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
							fprintf(out,"Unknown indexbase (N1/N9 index)\n");
							exit(3);
						}
						// get N1/N9 position
						e[0] = coord[(k1+indexbase)*3];
						e[1] = coord[(k1+indexbase)*3 + 1];
						e[2] = coord[(k1+indexbase)*3 + 2];
						// get C1* position (+11)
						r[0] = coord[(k1+11)*3];
						r[1] = coord[(k1+11)*3 + 1];
						r[2] = coord[(k1+11)*3 + 2];
						e[0] -= r[0]; // C1* --> N1/N9 unit vector
						e[1] -= r[1];
						e[2] -= r[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						ischi[j2] = true; // it is chi
						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						if(n_bb == n_bb2) // if both belong to the same unit
							n_chi++;
						else
						{
							n_chi = 1;
							n_bb2 = n_bb;
						}

						// ESTO SOBRA? CHECK PLEASE!
						if(fix != NULL)
							last_chi = true; // last chi is mobile

						// Updating "unat" for "chi" atoms
						for(int n=12; n<props[index_res].nat; n++)
							unat[k1+n] = n_bb + n_chi;

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
				}

				// the O3*-unat will be overwritten below

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{
					// EPSILON (bond between C3* and O3*)
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						r[0] = coord[(k1+7)*3];
						r[1] = coord[(k1+7)*3 + 1];
						r[2] = coord[(k1+7)*3 + 2];
						e[0] -= r[0]; // C3* --> O3*
						e[1] -= r[1];
						e[2] -= r[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization
						// fprintf(out,"EPSILON M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit
						n_chi = 0; // parent unit breaks!

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
					// END EPSILON

					// Updating "unat" for "epsilon" atoms
					unat[k1+8] = n_bb; // O3*

					// ZETA (bond between O3* and next-P)
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						r[0] = coord[(k1+8)*3];
						r[1] = coord[(k1+8)*3 + 1];
						r[2] = coord[(k1+8)*3 + 2];
						e[0] -= r[0]; // O3* --> next-P
						e[1] -= r[1];
						e[2] -= r[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit
						n_chi = 0; // parent unit breaks!

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
					// END ZETA
				}
			} // END RNA
			// DNA FRAGMENT
			else if( resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				// ********************************************************************************************
				// STANDARD NUCLEOTIDE BACKBONE --> 6-IC's (NOT LAST)

				// Updating "unat" for "alpha" atoms
				unat[k1] = n_bb; // P bb unit index
				unat[k1+1] = n_bb; // O1P bb unit index
				unat[k1+2] = n_bb; // O2P bb unit index

				// ALPHA (bond between P and O5*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get O5* position (+3)
					e[0] = coord[(k1+3)*3];
					e[1] = coord[(k1+3)*3 + 1];
					e[2] = coord[(k1+3)*3 + 2];
					// get P position (+0)
					r[0] = coord[k1*3];
					r[1] = coord[k1*3 + 1];
					r[2] = coord[k1*3 + 2];
					e[0] -= r[0]; // P --> O5*
					e[1] -= r[1];
					e[2] -= r[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END ALPHA

				// Updating "unat" for "beta" atoms
				unat[k1+3] = n_bb; // O5*

				// BETA (bond between O5* and C5*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get C5* position (+4)
					e[0] = coord[(k1+4)*3];
					e[1] = coord[(k1+4)*3 + 1];
					e[2] = coord[(k1+4)*3 + 2];
					// get O5* position (+3)
					r[0] = coord[(k1+3)*3];
					r[1] = coord[(k1+3)*3 + 1];
					r[2] = coord[(k1+3)*3 + 2];
					e[0] -= r[0]; // O5* --> C5*
					e[1] -= r[1];
					e[2] -= r[2];
					// fprintf(out,"BETA M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END BETA

				// Updating "unat" for "gamma" atoms
				unat[k1+4] = n_bb; // C5*

				// GAMMA (bond between C5* and C4*)
				if( fix == NULL || fix[j] ) // Check whether current variable is mobile
				{
					// get C4* position (+5)
					e[0] = coord[(k1+5)*3];
					e[1] = coord[(k1+5)*3 + 1];
					e[2] = coord[(k1+5)*3 + 2];
					// get C5* position (+4)
					r[0] = coord[(k1+4)*3];
					r[1] = coord[(k1+4)*3 + 1];
					r[2] = coord[(k1+4)*3 + 2];
					e[0] -= r[0]; // C5* --> C4*
					e[1] -= r[1];
					e[2] -= r[2];

					temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
					for ( m = 0; m < 3; m++ )
						e[m] /= temp; // Unit vector normalization

					// ea
					erx[0][j2] = e[0];
					erx[1][j2] = e[1];
					erx[2][j2] = e[2];
					// ea x ra
					erx[3][j2] = e[1] * r[2] - e[2] * r[1];
					erx[4][j2] = e[2] * r[0] - e[0] * r[2];
					erx[5][j2] = e[0] * r[1] - e[1] * r[0];

					undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
					un_index++; // count units
					ischi[j2] = false; // it's not chi
					istra[j2] = false;
					isrot[j2] = false;
					dhup[j2] = dhup_index;
					dhright[j2] = dhup_index;
					dhup_index++; // counts "up" dihedrals (1st body)

					// Defining rigid units
					n_bb += n_chi + 1; // +number of chis of preceeding unit
					n_chi = 0; // parent unit breaks!

					j2++;
				}
				else
					n_fixed++;

				j++; // der[] index
				// END GAMMA

				// Updating "unat" for "chi" atoms
				for(int n=5; n<props[index_res].nat; n++)
					unat[k1+n] = n_bb;

				//  LATERAL CHAIN-->CHI
				if( type == 2 )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
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
							fprintf(out,"Unknown indexbase (N1/N9 index)\n");
							exit(3);
						}
						// get N1/N9 position
						e[0] = coord[(k1+indexbase)*3];
						e[1] = coord[(k1+indexbase)*3 + 1];
						e[2] = coord[(k1+indexbase)*3 + 2];
						// get C1* position (+10) // DNA
						r[0] = coord[(k1+10)*3];
						r[1] = coord[(k1+10)*3 + 1];
						r[2] = coord[(k1+10)*3 + 2];
						e[0] -= r[0]; // C1* --> N1/N9 unit vector
						e[1] -= r[1];
						e[2] -= r[2];
						temp=sqrt(e[0]*e[0]+e[1]*e[1]+e[2]*e[2]);
						for(m=0;m<3;m++)
							e[m] /= temp; // unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						ischi[j2] = true; // it is chi
						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						if(n_bb == n_bb2) // if both belong to the same unit
							n_chi++;
						else
						{
							n_chi = 1;
							n_bb2 = n_bb;
						}

						// ESTO SOBRA? CHECK PLEASE!
						if(fix != NULL)
							last_chi = true; // last chi is mobile

						// Updating "unat" for "chi" atoms
						for(int n=11; n<props[index_res].nat; n++)
							unat[k1+n] = n_bb + n_chi;

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
				}

				// the O3*-unat will be overwritten below

				// NOT LAST NUCLEOTIDE (4-IC's)
				if ( iter_frag->pos_fragment != num_res-1 )
				{
					// EPSILON (bond between C3* and O3*)
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// get O3* position (+8)
						e[0] = coord[(k1+8)*3];
						e[1] = coord[(k1+8)*3 + 1];
						e[2] = coord[(k1+8)*3 + 2];
						// get C3* position (+7)
						r[0] = coord[(k1+7)*3];
						r[1] = coord[(k1+7)*3 + 1];
						r[2] = coord[(k1+7)*3 + 2];
						e[0] -= r[0]; // C3* --> O3*
						e[1] -= r[1];
						e[2] -= r[2];
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization
						// fprintf(out,"EPSILON M1= %f  M2= %f  after Y1= %f %f %f\n",mb1,mb2,Y1[0],Y1[1],Y1[2]);

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit
						n_chi = 0; // parent unit breaks!

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
					// END EPSILON

					// Updating "unat" for "epsilon" atoms
					unat[k1+8] = n_bb; // O3*

					// ZETA (bond between O3* and next-P)
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// get next-P position (k2+1)
						e[0] = coord[(k2+1)*3];
						e[1] = coord[(k2+1)*3 + 1];
						e[2] = coord[(k2+1)*3 + 2];
						// get O3* position (+8)
						r[0] = coord[(k1+8)*3];
						r[1] = coord[(k1+8)*3 + 1];
						r[2] = coord[(k1+8)*3 + 2];
						e[0] -= r[0]; // O3* --> next-P
						e[1] -= r[1];
						e[2] -= r[2];
						temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
						un_index++; // count units
						ischi[j2] = false; // it's not chi
						istra[j2] = false;
						isrot[j2] = false;
						dhup[j2] = dhup_index;
						dhright[j2] = dhup_index;
						dhup_index++; // counts "up" dihedrals (1st body)

						// Defining rigid units
						n_bb += n_chi + 1; // +number of chis of preceeding unit
						n_chi = 0; // parent unit breaks!

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
					// END ZETA
				}
			} // END DNA
			else if( fragtype == tmol_smol )
			{
				for(int n=0; n<props[index_res].nat; n++)
					unat[k1+n] = n_bb;
			}
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else
			{
			    // NOT FOUND MOL-TYPE
				fprintf(out,"Msg(hessianMFAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			// initialize for the next residue
			if( !( (iter_frag->pos_fragment == num_res-1) && (iter_seg->pos_segment == num_seg-1) ) )
			{	// if not the last fragment of the last segment
				k1 += props[index_res].nat;
				k2 += props[index_res + 1].nat;
			}
			index_res++; // counts residues (due to "props")
		}
		delete(iter_frag);
	}
	if(debug)
		fprintf(out,"Msg(hessianMFAx): loop finished!\n");
	iter_seg->clean_virtual();
	delete(iter_seg);

	// Inter-unit variables (rot-trans) dont increase units number!
	num_units = un_index + 1;
	if(debug_mem)
		fprintf(out,"Msg(hessianMFAx): Number of units: %d (%d)\n",num_units,un_index);

	// Tij (triangular) related pre-computations...
	int unipa_size=0;
	int unipa_index;
	int **unipa; // Stores "ipa"s indices of interacting atoms belonging to the same pair of units
	int *p_int;

	unipa = (int **) malloc( sizeof(int *) * num_units*(num_units+1)/2 );
	unipa_size += sizeof(int *) * num_units*(num_units+1)/2;
	for(int i=0; i<num_units*(num_units+1)/2; i++)
		unipa[i] = NULL; // initialization

	if(debug_mem)
		fprintf(out,"Msg(hessianMFAx): Triangular unipa_size= %.2f Mb (nipa=%d)\n",(float)unipa_size/1e6,nipa);

	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
	{
		k = decint[index].k; // k
		l = decint[index].l; // l
		i = unat[k]; // "k" atom unit index
		j = unat[l]; // "l" atom unit index
		if(i != j) // the same unit is rigid, so inner contacts must not be taken into account!
		{		   // (as well as the non contacting pairs!)
			if( i > j ) // The following must be allways true:  (k < l) && (i < j)
			{
				buff = j;
				j = i;
				i = buff;
			}
			unipa_index = i + j*(j+1)/2; // translates from "squared" to "triangular" matrix elements

			// If not allocated yet... then it does it!
			if(unipa[unipa_index] == NULL)
			{
				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
				{
					fprintf(out,"Msg(hessianMFAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
				unipa[unipa_index][0] = 2;
				unipa_size +=  2 * sizeof(int);
			}
			// If already allocated, then increases it size one int.
			else
			{
				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms in (i,j) units
				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
				{
					fprintf(out,"Msg(hessianMFAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
				unipa_size +=  sizeof(int);
			}

			unipa[unipa_index][ unipa[unipa_index][0]-1 ] = index; // storing ipa index for k,l atoms interacting pair
		}
		if(debug)
			fprintf(out,"index= %d  k= %d  l= %d  C= %f  d= %f\n",index,k,l,decint[index].C,decint[index].d);
	}

	if(debug)
		fprintf(out,"Msg(hessianMFAx): Final UNIPA size = %d bytes (%.2f Mb)\n",unipa_size,(float)unipa_size/1e6);

	// Hessian memory allocation
	if(debug_mem)
		fprintf(out,"Msg(hessianMFAx): Hessian-Matrix memory allocation %.3f Mb (%d)\n",(float) (sizeof(floating) * size*(size+1)/2 )/1e6, size*(size+1)/2);
	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
		{
			fprintf(out,"Msg(hessianMFAx): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			//			getchar();
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// ******************************
	// Uab - Matrix Computation
	// (ec.23) Noguti & Go 1983 pp.3685-90
	// ******************************
	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)

	int Usize=0;
	int uab_index,uab1,uab2,uab3,uab4,uab5,uab6,uab_void;
	double ***Uab;
	if( Uab = (double ***) malloc( sizeof(double **) * (1 + size*(size+1)/2) ) ) // Uab[size]
		for(int i=0; i<size*(size+1)/2; i++)
			Uab[i] = NULL; // initialization
	else
	{
		fprintf(out,"Msg(hessianMFAx): Unable to allocate Uab-matrix!\n");
		exit(1);
	}
	// void element (when there are just Traslational/Rotational DoFs)
	uab_void = size*(size+1)/2; // void element

	// VOID Uab element definition
	uab_index = size*(size+1)/2; // first Uab element index ("top-right" element)
	if( Uab[uab_index] = (double **) malloc( sizeof(double *) * 6 ) )
		for(int m=0;m<6;m++)
		{
			if( Uab[uab_index][m] = (double *) malloc( sizeof(double) * 6 ) )
				for(int n=0; n<6; n++)
					Uab[uab_index][m][n] = 0.0;
		}

	if(debug_mem)
	{
		Usize += sizeof(double **) * size*(size+1)/2;
		fprintf(out,"Msg(hessianMFAx): MEM-USAGE: Uab size = %f Mb\n",(float)Usize/1e6);
	}
	// WARNING: REMIND TO FREE "Uab" MEMORY!!!

	// Allocating first Uab element (6x6 sub-matrix)
	// First Uab element: 0,n
	uab_index = (size-1)*size/2; // first Uab element index ("top-right" element)
	if( Uab[uab_index] = (double **) malloc( sizeof(double *) * 6 ) )
		for(int m=0;m<6;m++)
		{
			if( Uab[uab_index][m] = (double *) malloc( sizeof(double) * 6 ) )
				for(int n=0; n<6; n++)
					Uab[uab_index][m][n] = 0.0;
			else
			{
				fprintf(out,"Msg(hessianMFAx): Unable to allocate an Uab-matrix element!\n");
				exit(1);
			}
		}
	else
	{
		fprintf(out,"Msg(hessianMFAx): Unable to allocate an Uab-matrix element!\n");
		exit(1);
	}
	Usize += sizeof(double *)*6 + sizeof(double)*36;

	// T memory allocation (6x6 matrix)
	double **T;
	T = (double **) malloc( sizeof(double *) * 6 );
	for(int i=0; i<6; i++)
		T[i] = (double *) malloc( sizeof(double) * 6 );

	// First Uab element: 0,n   (Unipa's "top-right" element = Tij)
	unipa_index = (undh[size-1]+1)*(undh[size-1]+2)/2; // i + j*(j+1)/2 --> translates from "squared" to "triangular" matrix elements
	if( unipa[unipa_index] != NULL )
	{
		if(debug)
			fprintf(out,"Computing Tij.  unipa_index= %d  i= %d  j= %d\n",unipa_index,0,undh[size-1]+1);
		calcTij( Uab[uab_index], unipa[unipa_index], unipa[unipa_index][0], decint, coord );
	}

	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
	for(int b=size-1, c= size+5; b >= 0; b--,c--) // b --> dihedrals (column)
	{											  // c --> deleteable col. index
		for(int a=0; a <= b+1; a++)	 // a --> dihedrals (row)
		{ 							 // it fills upper triangular part (including diagonal)
			// IMPORTANT:
			//	Delaying 1 "a" in Rab computation.
			//  (because we need U(a+1,b) element to compute Rab)
			//	(a-1 = a, a = a+1,... and so on...) (thus: a <= b+1, and when a=0 no Rab computation!)
			// NOTE: if we want to add more chis, we should delay "a" accordingly during Rab calculation!

			// Allocating Uab elements as needed
			if( a <= b ) // if it's in Uab range
			{
				uab_index = a + b*(b+1)/2; // translates from "squared" to "triangular" matrix elements

				// "in situ" Uab memory allocation
				if( Uab[uab_index] == NULL )
					if( Uab[uab_index] = (double **) malloc( sizeof(double *) * 6 ) )
					{
						for(int m=0;m<6;m++)
						{
							if( Uab[uab_index][m] = (double *) malloc( sizeof(double) * 6 ) )
								for(int n=0; n<6; n++)
									Uab[uab_index][m][n] = 0.0;
							else
							{ fprintf(out,"Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
						}
						Usize += sizeof(double *)*6 + sizeof(double)*36;
					}
					else
					{ fprintf(out,"Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }

				// unipa_index = i + j*(j+1)/2;
				// (translates from "squared" to "triangular" matrix elements)
				unipa_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2; // unipa[ a ][ undh[b]+1 ]

				// Uab element computation (NEW & robust: no valgrind's "invalids")
				if(b<size-1 && dhright[b+1] < size) // if valid "uab1"
					uab1 = a + dhright[b+1] * ( dhright[b+1] + 1 ) / 2; // Uab[ a ][ dhright[b+1] ]
				else
					uab1 = uab_void; // void element

				if(a>0 && dhup[a-1] >= 0)
					uab2 = dhup[a-1] + b * (b+1) / 2; // Uab[ dhup[a-1] ][ b ]
				else
					uab2 = uab_void; // void element

				if(a>0 && dhup[a-1] >=0 && b<size-1 && dhright[b+1] < size)
					uab4 = dhup[a-1] + dhright[b] * (dhright[b]+1) / 2; // Uab[ dhup[a-1] ][ dhright[b] ]
				else
					uab4 = uab_void; // void element

				if(a>0 && dhup[a-1]>=0 && b<size-1 && dhright[b+1]<size)
					uab3 = dhup[a-1] + dhright[b+1] * (dhright[b+1]+1) / 2; // Uab[ dhup[a-1] ][ dhright[b+1] ]
				else
					uab3 = uab_void; // void element

				// Computing Uab element
				if( unipa[unipa_index] != NULL ) // if "a" and "b" 's units interact
				{								 // (whether they have ipas...) --> then, "T" must be computed
					// Computing Tij element "on the fly"
					// (ec.21) Noguti & Go 1983 pp.3685-90
					if(debug)
						fprintf(out,"Computing Tij.  unipa_index= %d  i= %d  j= %d\n",unipa_index,undh[a],undh[b]+1);
					calcTij( T, unipa[unipa_index], unipa[unipa_index][0], decint, coord );

					// (ec.23) Noguti & Go 1983 pp.3685-90
					if( a > 0 && b < size-1 ) // In the middle... (most of them)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uab[uab_index][n][m] = Uab[uab1][n][m]
								                                    + Uab[uab2][n][m]
								                                                   - Uab[uab3][n][m]
								                                                                  + T[n][m];
					}
					else if( a == 0 && b < size-1 ) // In the top row (execepting: 0,n  already computed above)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uab[uab_index][n][m] = Uab[uab1][n][m]
								                                    + T[n][m];
					}
					else if( a > 0 && b == size-1 ) // In the right column (execepting: 0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
							{
								Uab[uab_index][n][m] = Uab[uab2][n][m]
								                                    + T[n][m];
							}
					}
				}
				else // if "a" and "b" 's units don't interact
				{	 // (no ipas...) --> then, "T" computation is not necessary!

					if( a > 0 && b < size-1 ) // In the middle... (most of them)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uab[uab_index][n][m] = Uab[uab1][n][m]
								                                    + Uab[uab2][n][m]
								                                                   - Uab[uab3][n][m];
					}
					else if( a == 0 && b < size-1 ) // In the top row (execepting: 0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uab[uab_index][n][m] = Uab[uab1][n][m];
					}
					else if( a > 0 && b == size-1 ) // In the right column (execepting: 0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uab[uab_index][n][m] = Uab[uab2][n][m];
					}
				}
			} // at this point, the Uab element is fully computed!

			// ******************************
			// Rab - Matrix Computation
			// Noguti & Go (1983) pp. 3685-90
			// ******************************
			// Rab Computation ( the same Single- or Multi- Chain)

			// IMPORTANT:
			//	Delaying 1 "a" in Rab computation.
			//  (because we need U(a+1,b) element to compute Rab)
			//	(a-1 = a, a = a+1,... and so on...) (thus: a <= b+1, and when a=0 no Rab computation!)
			// NOTE: if we want to add more chis, we should delay "a" accordingly during Rab calculation!

			if( a > 0 )
			{
				// Uab indices
				uab_index = a + b*(b+1)/2; // translates from "squared" to "triangular" matrix elements
				// unipa_index = i + j*(j+1)/2; // translates from "squared" to "triangular" matrix elements

				// Table I. Noguti & Go 1983 pp.3685-90
				//                         		(DELAYED)  		 (REAL)			(PAPER's)
				uab1 = a-1 + b*(b+1)/2; 	// Uab[a-1][b] 		Uab[a][b]		Uab[a-1][b]
				uab2 = a + b*(b+1)/2; 		// Uab[a][b] 		Uab[a+1][b]		U[c(a)][b]
				uab3 = a-1 + (b+1)*(b+2)/2; // Uab[a-1][b+1] 	Uab[a][b+1]		U[a][c(b)]
				uab4 = a + (b+1)*(b+2)/2; 	// Uab[a][b+1] 		Uab[a+1][b+1]	U[c(a)][c(b)]
				uab5 = a-1 + a*(a+1)/2; 	// Uab[a-1][a] 		Uab[a][a+1]		U[a][c(a)]
				uab6 = a + a*(a+1)/2; 		// Uab[a][a] 		Uab[a+1][a+1]	U[c(a)][c(a)]

				// backbone vs. backbone
				for(int n=0;n<6;n++)
					for(int m=0;m<6;m++)
						U[n][m] = Uab[uab1][n][m]; // Uab[a-1][b] --> Uab[a][b]

				// Table I. Noguti & Go 1983 pp.3685-90

				// chi vs. backbone ("a" not last row)
				if( (ischi[a-1] && !ischi[b]) && a-1 < size-1 )
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							U[n][m] -= Uab[uab2][n][m]; // -= U[c(a)][b] // -= Uab[a][b] --> Uab[a+1][b]

				// backbone vs. chi ("b" not last column)
				else if( (!ischi[a-1] && ischi[b]) && b < size-1 )
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							U[n][m] -= Uab[uab3][n][m]; // -= U[a][c(b)] // -= Uab[a-1][b+1] --> Uab[a][b+1]

				// chi vs. chi (a != b, a+1 in range, and b+1 in range)
				else if( (ischi[a-1] && ischi[b]) && a-1 != b && a < size && b+1 < size ) // && a+1 != b+1
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
						{
							U[n][m] -= Uab[uab2][n][m]; // -= U[c(a)][b]    //  -= Uab[a][b] --> Uab[a+1][b]
							U[n][m] -= Uab[uab3][n][m]; // -= U[a][c(b)]    //  -= Uab[a-1][b+1] --> Uab[a][b+1]
							U[n][m] += Uab[uab4][n][m]; // -= U[c(a)][c(b)] //  += Uab[a][b+1] --> Uab[a+1][b+1]
						}

				// chi vs. chi (a != b, a+1 in range, and b+1 is last column)
				else if( (ischi[a-1] && ischi[b]) && a-1 != b && a < size && b+1 == size) // && a+1 != b+1
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							U[n][m] -= Uab[uab2][n][m]; // -= U[c(a)][b]   //  -= Uab[a][b] --> Uab[a+1][b]

				// chi vs. chi (a != b, a+1 is outside, and b+1 is in range)
				else if( (ischi[a-1] && ischi[b]) && a-1 != b && a == size && b+1 < size) // && a+1 != b+1
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							U[n][m] -= Uab[uab3][n][m]; // -= U[a][c(b)]   //  -= Uab[a-1][b+1] --> Uab[a][b+1]

				// chi vs. chi (a == b, a+1 and b are in range )
				else if( ischi[a-1] && a-1 == b && a < size ) // && a+1 == b+1
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
						{
							U[n][m] -= Uab[uab5][n][m]; // -= U[a][c(a)]       // -= Uab[a-1][a] --> Uab[a][a+1]
							U[n][m] += Uab[uab6][m][n]; // += U[c(a)][c(a)]^t  // += Uab[a][a]^t --> Uab[a+1][a+1]^t
							U[n][m] -= Uab[uab5][m][n]; // -= U[b][c(a)]^t     // -= Uab[a-1][a] --> Uab[a][a+1]^t
						}

				// (ec.16) Noguti & Go 1983 pp.3685-90
				for(int n=0; n<6; n++)
				{
					R[n] = 0.0;
					for(int m=0; m<6; m++)
						R[n] += erx[m][a-1] * U[m][n]; // era x Rab
				}

				temp = 0.0; // dummy buffer
				for(int n=0; n<6; n++)
					temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

				hess_matrix[a-1 + b*(b+1)/2] = temp; // Rab' x erb = Hessian
			}
		}

		// Deleting unnecessary colums! (it saves much memory!)
		if( c < size ) // if col. is deleteable
		{
			// Deleting "c" column from Uab (up-diagonal part)
			for(int n=0; n <= c; n++)
			{
				uab_index = n + c*(c+1)/2; // translates from "squared" to "triangular" matrix elements

				// Deleting Uab element
				for(int m=0; m<6; m++)
					free( Uab[uab_index][m] );
				free( Uab[uab_index] );
				Usize -= sizeof(double *)*6 + sizeof(double)*36;
			}
		}
	}
	// Deleting 6 last Uab columns-rows
	if(size>6)
		for(int a=0; a<6; a++)
		{
			for(int b=a; b<6; b++)
			{
				uab_index = a + b*(b+1)/2; // translates from "squared" to "triangular" matrix elements
				for(int m=0; m<6; m++)
					free( Uab[uab_index][m] );
				free( Uab[uab_index] );
				Usize -= sizeof(double *)*6 + sizeof(double)*36;
			}
		}
	free( Uab );
	for(int i=0; i<num_units*(num_units+1)/2; i++)
		free( unipa[i] );
	free( unipa );

	for(int i=0;i<6;i++)
	{
		free( erx[i] ); // erx[6][size]
		free( T[i] );
	}
	free( erx );
	free( T );
	free( ischi );
	free( istra );
	free( isrot );
	free( undh );
	free( dhup );
	free( dhright );
}

//// Fast Hessian Matrix computation // (Multi-Chain & CA-model)
//// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
////        -Single row (N,CA,C)-PDB coordinates,
//// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
//// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
//// Noguti & Go (1983) pp. 3685-90
//// (3/2/2010)
//void hessianMCAx_ok(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,floating **p_hess_matrix,tri *props, int *unat, bool *fix)
//{
//	bool debug = false;
//	double prod,prod1;
//	int k,l,ind2,ind3,ls,ks;
//	double *dummy;
//	int prin, fin,i, j, j2, m, k1, k2, buff;
//	double r_alpha[3],r_beta[3],r[3],e[3];
//	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
//	double c,d,temp;
//	double v[6];
//	int resn,num_res,num_atoms,num_seg;
//	Residue *res;
//	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
//	pdbIter *iter_frag; // Iterator to screen fragments
//	pdbIter *iter_seg; // Iterator to screen segments
//	num_atoms = mol->get_num_atoms();
//	num_res = mol->get_num_fragments();
//	int index_res = 0;
//	floating *hess_matrix;
//	Htimer ht_timer; // timer
//	TMOL fragtype;
//
//	// ********************************************
//	// Storing==> (ea, ea x ra) (== (eb, eb x rb) )
//	// ********************************************
//	double **erx;
//	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
//	for(int i=0;i<6;i++)
//		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
//
//	int *undh; // returns the first unit-index on the left side of the dihedral
//	int un_index = 0;
//	undh = (int *) malloc( sizeof(int) * size );
//
//	bool any_rotrans = false;
//
//	// ****************************************************************************
//	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
//	// ****************************************************************************
//	j = 0; // current dihedral index
//	j2 = 0; // mobile dihedral index
//	k1 = 0;
//	int k0 = 0; // first-NH + CA + last-CO model index
//
//	int n_inter; // number of mobile intersegment variables
//	int num_units=0;
//	int n_fixed=0; // counter of fixed internal coordinates
//	int fix_offset=0; // unit-offset (joins units)
//	int natom=0;
//	bool break_unit = false;
//	bool already_broken = false;
//
//	// ---------------------------------------------------
//	// 1. BUILDING "erx" array and other auxiliars: "undh"
//	// ---------------------------------------------------
//
//	Segment * seg;
//	Atom * atom;
//	// Screening segments
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	num_seg = iter_seg->num_segment();
//	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter_frag = new pdbIter( seg );
//		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
//
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			n_inter=0; // number of mobile intersegment variables
//			for(int n=0;n<6;n++)
//				if(fix == NULL || fix[j+n]) // if mobile
//					n_inter++;
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
//			{
//				any_rotrans = true;
//				// Defining rigid units
//			}
//			else
//				any_rotrans = false;
//
//			// ********************
//			// Adding 3 TRASLATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (0)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					// ea x ra (Psi)  (-gamma_v)
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//					erx[3+axis][j2] = -1.0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 trans added
//
//			// ********************
//			// Adding 3 ROTATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (delta_v)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					erx[axis][j2] = 1.0;
//					// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 rots added
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if(any_rotrans)
//			{
//				un_index++; // count units (due to segment ending)
//			}
//		}
//
//		// Screen ALL fragments
//		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
//		{
//			res = ( Residue * ) iter_frag->get_fragment();
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();
//
//			// PROTEIN FRAGMENT
//			if( fragtype == tmol_protein )
//			{
//
//				if(iter_frag->pos_fragment == 0) // has NH
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//				// NOT FIRST, NOT LAST, NOT PRO -> PHI
//				// ("PHI-bond")
//				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
//					{
//						//			  printf ("%d Not first not PRO -> PHI \n",j);
//						// ("k2" is updated below (for the 1st not PRO -> PHI)
//
//						// y_lambda (NH pos 0)
//						e[0] = -coord[k1 * 3];
//						e[1] = -coord[k1 * 3 + 1];
//						e[2] = -coord[k1 * 3 + 2];
//						// CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// e_lambda
//						e[0] += r[0]; // NH --> CA
//						e[1] += r[1];
//						e[2] += r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] index
//				}  // NOT FIRST NOT PRO
//
//				if(break_unit)
//				{
//					un_index++;
//					break_unit = false; // by default not breaking unit
//				}
//
//				unat[k0] = un_index;
//				k0++;
//
//				already_broken = false;
//
//				// NOT LAST, NOT FIRST RESIDUE--> PSI
//				// ("PSI-bond")
//				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
//				{
//					//			 printf("%d Not last residue -> PSI\n",j);
//					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//					{
//						// get CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
//						e[0] = coord[(k1+2) * 3];
//						e[1] = coord[(k1+2) * 3 + 1];
//						e[2] = coord[(k1+2) * 3 + 2];
//						// e_lambda ==> CA --> C (unit vector)
//						e[0] -= r[0];
//						e[1] -= r[1];
//						e[2] -= r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] dihedral index
//				}  // NOT LAST
//
//				if(iter_frag->pos_fragment == num_res-1) // has CO
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//			} // END PROTEIN FRAGMENT
//			else // NOT FOUND MOL-TYPE
//			{
//				printf("Msg(hessianMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
//				exit(2);
//			}
//
//			k1+=3; // residue index // first residue atom index
//			index_res++; // counts residues (due to "props")
//		}
//		delete(iter_frag);
//	}
//	delete(iter_seg);
//
//	// ---------------------------------------------------------------------------
//	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
//	// ---------------------------------------------------------------------------
//	// Tij related pre-computations (to get direct Tij computation)
//
//	num_units = un_index + 1; // Inter-unit variables (rot-trans) dont increase units number
//
//	if(debug)
//		printf("Msg(hessianMCAx): Number of units: %d (%d)\n",num_units,un_index);
//
//	int unipa_size=0;
//	int unipa_index;
//	int *p_int;
//	int **unipa;
//
//	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
//	// (triangular packing storage)
//	unipa = (int **) malloc( sizeof(int *) * num_units*(num_units+1)/2 );
//	unipa_size += sizeof(int *) * num_units*(num_units+1)/2;
//
//	// "unipa" initialization
//	for(int i=0; i<num_units*(num_units+1)/2; i++)
//		unipa[i] = NULL;
//
//	if(debug)
//		printf("Msg(hessianMCAx): Triangular unipa_size= %d (nipa=%d)\n",unipa_size,nipa);
//
//	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
//	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
//	{
//		k = decint[index].k; // k
//		l = decint[index].l; // l
//		i = unat[k]; // "k" atom unit index
//		j = unat[l]; // "l" atom unit index
//		if(i != j) // the same unit is rigid, so inner contacts must not be taken into account!
//		{		   // (as well as the non contacting pairs!)
//			// The following must be allways true:  (k < l) && (i < j)
//			if( i > j )
//			{
//				buff = j;
//				j = i;
//				i = buff;
//			}
//			// translates from "squared" to "triangular" matrix elements
//			unipa_index = i + j*(j+1)/2;
//
//			// Allocating memory for the new element (memmory efficient)
//			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
//			{
//				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa[unipa_index][0] = 2;
//				unipa_size +=  2 * sizeof(int);
//			}
//			else // If already allocated, then increases it size one int.
//			{
//				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
//				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa_size +=  sizeof(int);
//			}
//
//			// Storing ipa index for k,l interacting pair of atoms
//			unipa[unipa_index][ unipa[unipa_index][0]-1 ] = index;
//		}
//	}
//
//	if(debug)
//		printf("Msg(hessianMCAx): Final UNIPA size = %d bytes (%.3f Mb)\n",unipa_size,(float)unipa_size/1e6);
//
//	// -----------------------------------------------------------------------------------
//	// 3. HESSIAN BUILDING from Uab elements:
//	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
//	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//	// -----------------------------------------------------------------------------------
//
//	// Hessian memory allocation
//	if(*p_hess_matrix == NULL)
//	{
//		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
//		{
//			printf("Msg(hessianMCAx): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
//			exit(1);
//		}
//		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
//	}
//	else
//		hess_matrix = *p_hess_matrix;
//
//	// Hessian initialization
//	for(i=0; i<size*(size+1)/2; i++)
//		hess_matrix[i]=0.0;
//
//	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
//	int Usize=0;
//	int uij_index,uij1,uij2,uij3,uij4,uij_void;
//	double ***Uij;
//	// +1 due to void element (when there are just Traslational/Rotational DoFs)
//	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
//		for(int i=0; i<(num_units*(num_units+1)/2); i++)
//			Uij[i] = NULL; // initialization
//	else
//	{
//		printf("Msg(hessianMCAx): Unable to allocate Uij-matrix!\n");
//		exit(1);
//	}
//
//	Usize += sizeof(double **) * num_units*(num_units+1)/2;
//
//	// T memory allocation (6x6 matrix)
//	double **T;
//	T = (double **) malloc( sizeof(double *) * 6 );
//	for(int i=0; i<6; i++)
//		T[i] = (double *) malloc( sizeof(double) * 6 );
//
//	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
//	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
//	for(int b=size-1, c= size+5; b >= 0; b--,c--) // b --> dihedrals (column)
//		//	for(int b=size-1, c= size+2; b >= 0; b--,c--) // b --> dihedrals (column)
//	{											  // c --> deleteable col. index
//		for(int a=0; a <= b; a++)	 // a --> dihedrals (row)
//		{ 							 // it fills upper triangular part (including diagonal)
//			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
//			uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2;
//
//			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguos positions
//			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
//			// the Uij for the 4 ICs combination should be computed only once!
//			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
//			{
//				// "in situ" Uij memory allocation
//				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
//				{
//					for(int m=0;m<6;m++)
//					{
//						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
//							for(int n=0; n<6; n++)
//								Uij[uij_index][m][n] = 0.0;
//						else
//						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
//					}
//					Usize += sizeof(double *)*6 + sizeof(double)*36;
//				}
//				else
//				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
//
//				// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//				// Computing valid indices
//
//				// U(a,b+1)
//				if(undh[b]+1 < num_units-1 ) // if valid "uab1"
//					uij1 = undh[a] + (undh[b]+2)*(undh[b]+3)/2; // i= undh[a]  j= (undh[b]+1)+1
//
//				// U(a-1,b)
//				if( undh[a] > 0 )
//					uij2 = undh[a]-1 + (undh[b]+1)*(undh[b]+2) / 2; // Uab[ dhup[a-1] ][ dhright[b] ]
//
//				// U(a-1,b+1)
//				if( undh[a] > 0 && undh[b]+2 < num_units )
//					uij3 = undh[a]-1 + (undh[b]+2)*(undh[b]+3) / 2; // Uab[ dhup[a-1] ][ dhright[b+1] ]
//
//				// Computing Uij element
//				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
//				{								 // (whether they have ipas...) --> then, "T" must be computed
//					// Computing Tij element "on the fly"
//					// (ec.21) Noguti & Go 1983 pp.3685-90
//					calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );
//
//					// (ec.23) Noguti & Go 1983 pp.3685-90
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + Uij[uij2][n][m]
//								                                                   - Uij[uij3][n][m]
//								                                                                  + T[n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + T[n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij2][n][m]
//								                                    + T[n][m];
//					}
//					else // a == 0 && b == size-1  (0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = T[n][m];
//					}
//				}
//				else // if "a" and "b" 's units don't interact
//				{	 // (no ipas...) --> then, "T" computation is not necessary!
//
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + Uij[uij2][n][m]
//								                                                   - Uij[uij3][n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij2][n][m];
//					}
//					// if no interaction, then no Tij addition to Uab needed!
//				}
//				// at this point, the Uij element is fully computed!
//			}   // only if not computed yet!
//
//			// ****************************************************
//			// Rab - Matrix Computation
//			// Noguti & Go (1983) pp. 3685-90
//			// ****************************************************
//			// Rab Computation ( the same Single- or Multi- Chain)
//			// Table I. Noguti & Go 1983 pp.3685-90
//			// backbone vs. backbone --> Rab = Uab
//
//			// (ec.16) Noguti & Go 1983 pp.3685-90
//			for(int n=0; n<6; n++)
//			{
//				R[n] = 0.0;
//				for(int m=0; m<6; m++)
//					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
//			}
//			temp = 0.0;
//			for(int n=0; n<6; n++)
//				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian
//			hess_matrix[a + b*(b+1)/2] = temp; // Rab' x erb = Hessian
//		}
//
//		// Deleting unnecessary colums! (it saves much memory!)
//		if( c < size ) // if col. is deleteable
//		{
//			// Deleting "c" column from Uab (up-diagonal part)
//			for(int n=0; n <= c; n++)
//			{
//				uij_index = undh[n] + (undh[c]+1)*(undh[c]+2)/2; // translates from "squared" to "triangular" matrix elements
//
//				if(Uij[uij_index]!=NULL)
//				{
//					// Deleting Uab element
//					for(int m=0; m<6; m++)
//						free( Uij[uij_index][m] );
//					free( Uij[uij_index] );
//					Uij[uij_index]=NULL;
//					Usize -= sizeof(double *)*6 + sizeof(double)*36;
//				}
//			}
//		}
//	}
//	// Deleting 6 last Uab columns-rows
//	if(size>6)
//		//	if(size>3)
//		for(int a=0; a<6; a++)
//			//	for(int a=0; a<3; a++)
//		{
//			for(int b=a; b<6; b++)
//				//		for(int b=a; b<3; b++)
//			{
//				uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2; // translates from "squared" to "triangular" matrix elements
//
//				if(Uij[uij_index]!=NULL)
//				{
//					// Deleting Uab element
//					for(int m=0; m<6; m++)
//						free( Uij[uij_index][m] );
//					free( Uij[uij_index] );
//					Uij[uij_index]=NULL;
//					Usize -= sizeof(double *)*6 + sizeof(double)*36;
//				}
//			}
//		}
//	free( Uij );
//	for(int i=0; i<num_units*(num_units+1)/2; i++)
//		free( unipa[i] );
//	free( unipa );
//
//	for(int i=0;i<6;i++)
//	{
//		free( erx[i] ); // erx[6][size]
//		free( T[i] );
//	}
//	free( erx );
//	free( T );
//	free( undh );
//}


//// BACKUP of hessianMCAx() working ok (28/3/2012)
//// Fast Hessian Matrix computation (Multi-Chain & CA-model) (Huge-systems)
//// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
////        -Single row (N,CA,C)-PDB coordinates,
//// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
//// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
//// Noguti & Go (1983) pp. 3685-90
//// (3/2/2010)
//void hessianMCAx_backup(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,floating **p_hess_matrix,tri *props, int *unat, bool *fix)
//{
//	bool debug = false;
//	double prod,prod1;
//	int k,l,ind2,ind3,ls,ks;
//	double *dummy;
//	int prin, fin, j2, m, k1, k2, buff;
//	long int i,j; // Needed for huge systems!
//	long int sizel = size; // Needed for huge systems!
//	double r_alpha[3],r_beta[3],r[3],e[3];
//	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
//	double c,d,temp;
//	double v[6];
//	int resn,num_res,num_atoms,num_seg;
//	Residue *res;
//	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
//	pdbIter *iter_frag; // Iterator to screen fragments
//	pdbIter *iter_seg; // Iterator to screen segments
//	num_atoms = mol->get_num_atoms();
//	num_res = mol->get_num_fragments();
//	int index_res = 0;
//	floating *hess_matrix;
//	Htimer ht_timer; // timer
//	TMOL fragtype;
//
//	// ********************************************
//	// Storing==> (ea, ea x ra) (== (eb, eb x rb) )
//	// ********************************************
//	double **erx;
//	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
//	for(int i=0;i<6;i++)
//		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
//
//	int *undh; // returns the first unit-index on the left side of the dihedral
//	int un_index = 0;
//	undh = (int *) malloc( sizeof(int) * size );
//
//	bool any_rotrans = false;
//
//	// ****************************************************************************
//	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
//	// ****************************************************************************
//	j = 0; // current dihedral index
//	j2 = 0; // mobile dihedral index
//	k1 = 0;
//	int k0 = 0; // first-NH + CA + last-CO model index
//
//	int n_inter; // number of mobile intersegment variables
//	long int num_units=0; // Needed for huge systems!
//	int n_fixed=0; // counter of fixed internal coordinates
//	int fix_offset=0; // unit-offset (joins units)
//	int natom=0;
//	bool break_unit = false;
//	bool already_broken = false;
//
//	// ---------------------------------------------------
//	// 1. BUILDING "erx" array and other auxiliars: "undh"
//	// ---------------------------------------------------
//
//	Segment * seg;
//	Atom * atom;
//	// Screening segments
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	num_seg = iter_seg->num_segment();
//	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter_frag = new pdbIter( seg );
//		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
//
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			n_inter=0; // number of mobile intersegment variables
//			for(int n=0;n<6;n++)
//				if(fix == NULL || fix[j+n]) // if mobile
//					n_inter++;
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
//			{
//				any_rotrans = true;
//				// Defining rigid units
//			}
//			else
//				any_rotrans = false;
//
//			// ********************
//			// Adding 3 TRASLATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (0)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					// ea x ra (Psi)  (-gamma_v)
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//					erx[3+axis][j2] = -1.0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 trans added
//
//			// ********************
//			// Adding 3 ROTATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (delta_v)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					erx[axis][j2] = 1.0;
//					// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 rots added
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if(any_rotrans)
//			{
//				un_index++; // count units (due to segment ending)
//			}
//		}
//
//		// Screen ALL fragments
//		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
//		{
//			res = ( Residue * ) iter_frag->get_fragment();
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();
//
//			// PROTEIN FRAGMENT
//			if( fragtype == tmol_protein )
//			{
//
//				if(iter_frag->pos_fragment == 0) // has NH
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//				// NOT FIRST, NOT LAST, NOT PRO -> PHI
//				// ("PHI-bond")
//				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
//					{
//						//			  printf ("%d Not first not PRO -> PHI \n",j);
//						// ("k2" is updated below (for the 1st not PRO -> PHI)
//
//						// y_lambda (NH pos 0)
//						e[0] = -coord[k1 * 3];
//						e[1] = -coord[k1 * 3 + 1];
//						e[2] = -coord[k1 * 3 + 2];
//						// CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// e_lambda
//						e[0] += r[0]; // NH --> CA
//						e[1] += r[1];
//						e[2] += r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] index
//				}  // NOT FIRST NOT PRO
//
//				if(break_unit)
//				{
//					un_index++;
//					break_unit = false; // by default not breaking unit
//				}
//
//				unat[k0] = un_index;
//				k0++;
//
//				already_broken = false;
//
//				// NOT LAST, NOT FIRST RESIDUE--> PSI
//				// ("PSI-bond")
//				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
//				{
//					//			 printf("%d Not last residue -> PSI\n",j);
//					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//					{
//						// get CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
//						e[0] = coord[(k1+2) * 3];
//						e[1] = coord[(k1+2) * 3 + 1];
//						e[2] = coord[(k1+2) * 3 + 2];
//						// e_lambda ==> CA --> C (unit vector)
//						e[0] -= r[0];
//						e[1] -= r[1];
//						e[2] -= r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] dihedral index
//				}  // NOT LAST
//
//				if(iter_frag->pos_fragment == num_res-1) // has CO
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//			} // END PROTEIN FRAGMENT
//			else // NOT FOUND MOL-TYPE
//			{
//				printf("Msg(hessianMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
//				exit(2);
//			}
//
//			k1+=3; // residue index // first residue atom index
//			index_res++; // counts residues (due to "props")
//		}
//		delete(iter_frag);
//	}
//	delete(iter_seg);
//
//	// ---------------------------------------------------------------------------
//	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
//	// ---------------------------------------------------------------------------
//	// Tij related pre-computations (to get direct Tij computation)
//
//	num_units = un_index + 1; // Inter-unit variables (rot-trans) dont increase units number
//
//	if(debug)
//		printf("Msg(hessianMCAx): Number of units: %ld (%d)\n",num_units,un_index);
//
//	long int unipa_size=0; // Needed for huge systems!
//	long int unipa_index;
//	int *p_int;
//	int **unipa;
//
//	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
//	// (triangular packing storage)
//	unipa = (int **) malloc( sizeof(int *) * num_units*(num_units+1)/2 );
//	unipa_size += sizeof(int *) * num_units*(num_units+1)/2;
//
//	// "unipa" initialization
//	for(i=0; i<num_units*(num_units+1)/2; i++)
//		unipa[i] = NULL;
//
//	if(debug)
//		printf("Msg(hessianMCAx): Triangular unipa_size= %ld (nipa=%d)\n",unipa_size,nipa);
//
//	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
//	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
//	{
//		k = decint[index].k; // k
//		l = decint[index].l; // l
//		i = unat[k]; // "k" atom unit index
//		j = unat[l]; // "l" atom unit index
//		if(i != j) // the same unit is rigid, so inner contacts must not be taken into account!
//		{		   // (as well as the non contacting pairs!)
//			// The following must be allways true:  (k < l) && (i < j)
//			if( i > j )
//			{
//				buff = j;
//				j = i;
//				i = buff;
//			}
//			// translates from "squared" to "triangular" matrix elements
//			unipa_index = i + j*(j+1)/2;
//
//			// Allocating memory for the new element (memmory efficient)
//			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
//			{
//				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa[unipa_index][0] = 2;
//				unipa_size +=  2 * sizeof(int);
//			}
//			else // If already allocated, then increases it size one int.
//			{
//				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
//				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa_size +=  sizeof(int);
//			}
//
//			// Storing ipa index for k,l interacting pair of atoms
//			unipa[unipa_index][ unipa[unipa_index][0]-1 ] = index;
//		}
//	}
//
//	if(debug)
//		printf("Msg(hessianMCAx): Final UNIPA size = %ld bytes (%.3f Mb)\n",unipa_size,(float)unipa_size/1e6);
//
//	// -----------------------------------------------------------------------------------
//	// 3. HESSIAN BUILDING from Uab elements:
//	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
//	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//	// -----------------------------------------------------------------------------------
//
//	// Hessian memory allocation
//	if(*p_hess_matrix == NULL)
//	{
////fprintf(stderr,"BEFORE 1) malloc\n");
//		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * sizel*(sizel+1)/2 )) ) // Triangular
//		{
//			printf("Msg(hessianMCAx): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
//			exit(1);
//		}
//		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
////fprintf(stderr,"AFTER 1) malloc\n");
//	}
//	else
//		hess_matrix = *p_hess_matrix;
//
////fprintf(stderr,"BEFORE 2) initialization\n");
//	// Hessian initialization
//	for(i=0; i<sizel*(sizel+1)/2; i++)
//		hess_matrix[i]=0.0;
////fprintf(stderr,"AFTER 2) initialization\n");
//
////fprintf(stderr,"BEFORE 3) Uab\n");
////fprintf(stderr,"BEFORE 4) Uab allocation\n");
//	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
//	long int Usize=0; // Needed for huge systems!
//	long int uij_index,uij1,uij2,uij3,uij4,uij_void; // Needed for huge systems!
//	double ***Uij;
//	// +1 due to void element (when there are just Traslational/Rotational DoFs)
//	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
//		for(i=0; i<(num_units*(num_units+1)/2); i++)
//			Uij[i] = NULL; // initialization
//	else
//	{
//		printf("Msg(hessianMCAx): Unable to allocate Uij-matrix!\n");
//		exit(1);
//	}
////fprintf(stderr,"AFTER 4) Uab allocation\n");
//
//	Usize += sizeof(double **) * num_units*(num_units+1)/2;
//
//	// T memory allocation (6x6 matrix)
//	double **T;
//	T = (double **) malloc( sizeof(double *) * 6 );
//	for(int i=0; i<6; i++)
//		T[i] = (double *) malloc( sizeof(double) * 6 );
//
//	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
//	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
////fprintf(stderr,"BEFORE 5) Uab main-loop\n");
//	for(long int b=size-1, c= size+5; b >= 0; b--,c--) // b --> dihedrals (column)
//		//	for(int b=size-1, c= size+2; b >= 0; b--,c--) // b --> dihedrals (column)
//	{											  // c --> deleteable col. index
//		for(long int a=0; a <= b; a++)	 // a --> dihedrals (row)
//		{ 							 // it fills upper triangular part (including diagonal)
//			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
//			uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2;
//			// fprintf(stderr,"a= %d  b= %d  undh[a]= %d  undh[b]+1= %d  uij_index= %d\n",a,b,undh[a],undh[b]+1,uij_index);
////fprintf(stderr,"%ld ",uij_index);
//			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguos positions
//			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
//			// the Uij for the 4 ICs combination should be computed only once!
//			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
//			{
//				// "in situ" Uij memory allocation
////fprintf(stderr,"BEFORE 6) Uij not-computed\n");
//
////fprintf(stderr,"BEFORE 7) Uij element allocation\n");
//				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
//				{
//					for(int m=0;m<6;m++)
//					{
//						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
//							for(int n=0; n<6; n++)
//								Uij[uij_index][m][n] = 0.0;
//						else
//						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
//					}
//					Usize += sizeof(double *)*6 + sizeof(double)*36;
//				}
//				else
//				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
////fprintf(stderr,"AFTER 7) Uij element allocation\n");
//
//				// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//				// Computing valid indices
//
//				// U(a,b+1)
//				if(undh[b]+1 < num_units-1 ) // if valid "uab1"
//					uij1 = undh[a] + (undh[b]+2)*(undh[b]+3)/2; // i= undh[a]  j= (undh[b]+1)+1
//
//				// U(a-1,b)
//				if( undh[a] > 0 )
//					uij2 = undh[a]-1 + (undh[b]+1)*(undh[b]+2) / 2; // Uab[ dhup[a-1] ][ dhright[b] ]
//
//				// U(a-1,b+1)
//				if( undh[a] > 0 && undh[b]+2 < num_units )
//					uij3 = undh[a]-1 + (undh[b]+2)*(undh[b]+3) / 2; // Uab[ dhup[a-1] ][ dhright[b+1] ]
//
//				// fprintf(stderr,"\tUij: uij1(row=%d col=%d) uij2(row=%d col=%d) uij3(row=%d col=%d) \n",undh[a],undh[b]+2,undh[a]-1,undh[b]+1,undh[a]-1,undh[b]+2);
//
//				// Computing Uij element
////fprintf(stderr,"BEFORE 8) Uij element computation\n");
//				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
//				{								 // (whether they have ipas...) --> then, "T" must be computed
//					// Computing Tij element "on the fly"
//					// (ec.21) Noguti & Go 1983 pp.3685-90
//					calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );
//					// fprintf(stderr,"\tTij computed\n");
//
//					// (ec.23) Noguti & Go 1983 pp.3685-90
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + Uij[uij2][n][m]
//								                                                   - Uij[uij3][n][m]
//								                                                                  + T[n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + T[n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij2][n][m]
//								                                    + T[n][m];
//					}
//					else // a == 0 && b == size-1  (0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = T[n][m];
//					}
//				}
//				else // if "a" and "b" 's units don't interact
//				{	 // (no ipas...) --> then, "T" computation is not necessary!
//					// fprintf(stderr,"\tTij not-computed\n");
//
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m]
//								                                    + Uij[uij2][n][m]
//								                                                   - Uij[uij3][n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij1][n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uij_index][n][m] = Uij[uij2][n][m];
//					}
//					// if no interaction, then no Tij addition to Uab needed!
//				}
////fprintf(stderr,"AFTER 8) Uij element computation\n");
//
//				// at this point, the Uij element is fully computed!
//			}   // only if not computed yet!
////fprintf(stderr,"AFTER 6) Uij not-computed\n");
//
//			// ****************************************************
//			// Rab - Matrix Computation
//			// Noguti & Go (1983) pp. 3685-90
//			// ****************************************************
//			// Rab Computation ( the same Single- or Multi- Chain)
//			// Table I. Noguti & Go 1983 pp.3685-90
//			// backbone vs. backbone --> Rab = Uab
//
//			// (ec.16) Noguti & Go 1983 pp.3685-90
//			// fprintf(stderr,"\tRij= ");
////fprintf(stderr,"BEFORE 9) R from Uij\n");
//			for(int n=0; n<6; n++)
//			{
//				R[n] = 0.0;
//				for(int m=0; m<6; m++)
//					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
//				// fprintf(stderr,"%f ",R[n]);
//			}
////fprintf(stderr,"AFTER 9) R from Uij\n");
////fprintf(stderr,"BEFORE 10) Hess from R\n");
//			temp = 0.0;
//			for(int n=0; n<6; n++)
//				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian
////			indexl = a+b*(b+1)/2; // Needed for huge systems!
////fprintf(stderr,"BEFORE 11) Hess writting (a= %ld  b= %ld  index= %ld\n",a,b,a+b*(b+1)/2);
//			hess_matrix[a+b*(b+1)/2] = temp; // Rab' x erb = Hessian
////fprintf(stderr,"AFTER 11) Hess writting\n");
//			// fprintf(stderr," Hij= %f\n",hess_matrix[a + b*(b+1)/2]);
////fprintf(stderr,"AFTER 10) Hess from R\n");
//		}
////fprintf(stderr,"AFTER 5) Uab main-loop\n");
//
//		// Deleting unnecessary colums! (it saves much memory!)
//		if( c < size ) // if col. is deleteable
//		{
//			// Deleting "c" column from Uab (up-diagonal part)
//			for(int n=0; n <= c; n++)
//			{
//				uij_index = undh[n] + (undh[c]+1)*(undh[c]+2)/2; // translates from "squared" to "triangular" matrix elements
//
//				if(Uij[uij_index]!=NULL)
//				{
//					// Deleting Uab element
//					for(int m=0; m<6; m++)
//						free( Uij[uij_index][m] );
//					free( Uij[uij_index] );
//					Uij[uij_index]=NULL;
//					Usize -= sizeof(double *)*6 + sizeof(double)*36;
//				}
//			}
//		}
////fprintf(stderr,"AFTER 3) Uab\n");
//	}
//	// Deleting 6 last Uab columns-rows
//	if(size>6)
//		//	if(size>3)
//		for(int a=0; a<6; a++)
//			//	for(int a=0; a<3; a++)
//		{
//			for(int b=a; b<6; b++)
//				//		for(int b=a; b<3; b++)
//			{
//				uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2; // translates from "squared" to "triangular" matrix elements
//
//				if(Uij[uij_index]!=NULL)
//				{
//					// Deleting Uab element
//					for(int m=0; m<6; m++)
//						free( Uij[uij_index][m] );
//					free( Uij[uij_index] );
//					Uij[uij_index]=NULL;
//					Usize -= sizeof(double *)*6 + sizeof(double)*36;
//				}
//			}
//		}
//	free( Uij );
//	for(i=0; i<num_units*(num_units+1)/2; i++)
//		free( unipa[i] );
//	free( unipa );
//
//	for(i=0;i<6;i++)
//	{
//		free( erx[i] ); // erx[6][size]
//		free( T[i] );
//	}
//	free( erx );
//	free( T );
//	free( undh );
//}

// Fast Hessian Matrix computation (Multi-Chain & CA-model) (Huge-systems)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// Noguti & Go (1983) pp. 3685-90
// (3/2/2010)
void hessianMCAx(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,floating **p_hess_matrix,tri *props, int *unat, bool *fix)
{
	// TESTING variables
	long int allocated_mem = 0;

	bool debug = false;
	double prod,prod1;
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	int prin, fin, j2, m, k1, k2, buff;
	long int i,j; // Needed for huge systems!
	long int sizel = size; // Needed for huge systems!
	double r_alpha[3],r_beta[3],r[3],e[3];
	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
	double c,d,temp;
	double v[6];
	int resn,num_res,num_atoms,num_seg;
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();
	int index_res = 0;
	floating *hess_matrix;
	Htimer ht_timer; // timer
	TMOL fragtype;

	if(debug)
		printf("debug> WARNING! Using working: hessianMCAx()\n");

	// ********************************************
	// Storing==> (ea, ea x ra) (== (eb, eb x rb) )
	// ********************************************
	double **erx;
	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	allocated_mem += sizeof(double *) * 6;
	for(int i=0;i<6;i++)
	{
		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
		allocated_mem += sizeof(double *) * size;
	}

	int *undh; // returns the first unit-index on the left side of the dihedral
	int un_index = 0;
	undh = (int *) malloc( sizeof(int) * size );
	allocated_mem += sizeof(int) * size;
	bool any_rotrans = false;

	if(debug)
		fprintf(stderr,"debug> Memory allocated by hessianMCAx() before (ea, ea x ra) and (eb, eb x rb) building: %ld\n",allocated_mem);

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	j = 0; // current dihedral index
	j2 = 0; // mobile dihedral index
	k1 = 0;
	int k0 = 0; // first-NH + CA + last-CO model index

	int n_inter; // number of mobile intersegment variables
	long int num_units=0; // Needed for huge systems!
	int n_fixed=0; // counter of fixed internal coordinates
	int fix_offset=0; // unit-offset (joins units)
	int natom=0;
	bool break_unit = false;
	bool already_broken = false;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary: "undh"
	// ---------------------------------------------------

	Segment * seg;
	Atom * atom;
	// Screening segments
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			n_inter=0; // number of mobile intersegment variables
			for(int n=0;n<6;n++)
				if(fix == NULL || fix[j+n]) // if mobile
					n_inter++;

			// Check whether any of the 6D inter-segment variables are fixed
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
			{
				any_rotrans = true;
				// Defining rigid units
			}
			else
				any_rotrans = false;

			// ********************
			// Adding 3 TRASLATIONS
			// ********************
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// ea (Phi)  (0)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					// ea x ra (Psi)  (-gamma_v)
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;
					erx[3+axis][j2] = -1.0;

					undh[j2] = un_index;

					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}

				j++; // adding traslational degree of freedom
			} // 3 trans added

			// ********************
			// Adding 3 ROTATIONS
			// ********************
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// ea (Phi)  (delta_v)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					erx[axis][j2] = 1.0;
					// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;

					undh[j2] = un_index;

					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}

				j++; // adding traslational degree of freedom
			} // 3 rots added

			// Check whether any of the 6D inter-segment variables are fixed
			if(any_rotrans)
			{
				un_index++; // count units (due to segment ending)
			}
		}

		// Screen ALL fragments
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{

				if(iter_frag->pos_fragment == 0) // has NH
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

				// NOT FIRST, NOT LAST, NOT PRO -> PHI
				// ("PHI-bond")
				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						//			  printf ("%d Not first not PRO -> PHI \n",j);
						// ("k2" is updated below (for the 1st not PRO -> PHI)

						// y_lambda (NH pos 0)
						e[0] = -coord[k1 * 3];
						e[1] = -coord[k1 * 3 + 1];
						e[2] = -coord[k1 * 3 + 2];
						// CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// e_lambda
						e[0] += r[0]; // NH --> CA
						e[1] += r[1];
						e[2] += r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
				}  // NOT FIRST NOT PRO

				if(break_unit)
				{
					un_index++;
					break_unit = false; // by default not breaking unit
				}

				unat[k0] = un_index;
				k0++;

				already_broken = false;

				// NOT LAST, NOT FIRST RESIDUE--> PSI
				// ("PSI-bond")
				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
				{
					//			 printf("%d Not last residue -> PSI\n",j);
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					{
						// get CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
						e[0] = coord[(k1+2) * 3];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= r[0];
						e[1] -= r[1];
						e[2] -= r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] dihedral index
				}  // NOT LAST

				if(iter_frag->pos_fragment == num_res-1) // has CO
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

			} // END PROTEIN FRAGMENT
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(hessianMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			k1+=3; // residue index // first residue atom index
			index_res++; // counts residues (due to "props")
		}
		delete(iter_frag);
	}
	delete(iter_seg);

	if(debug)
		fprintf(stderr,"debug> Memory allocated by hessianMCAx() after (ea, ea x ra) and (eb, eb x rb) building: %ld\n",allocated_mem);

	// ---------------------------------------------------------------------------
	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
	// ---------------------------------------------------------------------------
	// Tij related pre-computations (to get direct Tij computation)

	num_units = un_index + 1; // Inter-unit variables (rot-trans) dont increase units number

	if(debug)
		printf("Msg(hessianMCAx): Number of units: %ld (%d)\n",num_units,un_index);

	long int unipa_size=0; // Needed for huge systems!
	long int unipa_index;
	int *p_int;
	int **unipa;

	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
	// (triangular packing storage)
	unipa = (int **) malloc( sizeof(int *) * num_units*(num_units+1)/2 );
	allocated_mem += sizeof(int *) * num_units*(num_units+1)/2;

	unipa_size += sizeof(int *) * num_units*(num_units+1)/2;

	// "unipa" initialization
	for(i=0; i<num_units*(num_units+1)/2; i++)
		unipa[i] = NULL;

	if(debug)
		printf("Msg(hessianMCAx): Triangular unipa_size= %ld (nipa=%d)\n",unipa_size,nipa);

	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
	{
		k = decint[index].k; // k
		l = decint[index].l; // l
		i = unat[k]; // "k" atom unit index
		j = unat[l]; // "l" atom unit index
		if(i != j) // the same unit is rigid, so inner contacts must not be taken into account!
		{		   // (as well as the non contacting pairs!)
			// The following must be always true:  (k < l) && (i < j)
			if( i > j )
			{
				buff = j;
				j = i;
				i = buff;
			}
			// translates from "squared" to "triangular" matrix elements
			unipa_index = i + j*(j+1)/2;

			// Allocating memory for the new element (memmory efficient)
			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
			{
				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
				{
					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
				unipa[unipa_index][0] = 2;
				unipa_size +=  2 * sizeof(int);
				allocated_mem += sizeof(int) * 2;
			}
			else // If already allocated, then increases it size one int.
			{
				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
				{
					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
				unipa_size +=  sizeof(int);
				allocated_mem += sizeof(int);
			}

			// Storing ipa index for k,l interacting pair of atoms
			unipa[unipa_index][ unipa[unipa_index][0]-1 ] = index;
		}
	}

	if(debug)
	{
		printf("Msg(hessianMCAx): Final UNIPA size = %ld bytes (%.3f Mb)\n",unipa_size,(float)unipa_size/1e6);
		fprintf(stderr,"debug> Memory allocated by hessianMCAx() after UNIPA building: %ld\n",allocated_mem);
	}

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------

	// Hessian memory allocation
	if(*p_hess_matrix == NULL)
	{
//fprintf(stderr,"BEFORE 1) malloc\n");
		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * sizel*(sizel+1)/2 )) ) // Triangular
		{
			printf("Msg(hessianMCAx): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
		allocated_mem += sizeof(floating) * sizel*(sizel+1)/2;
		if(debug)
			fprintf(stderr,"debug> Memory allocated by hessianMCAx() HESSIAN MATRIX: %ld\n",sizeof(floating) * sizel*(sizel+1)/2 );
//fprintf(stderr,"AFTER 1) malloc\n");
	}
	else
		hess_matrix = *p_hess_matrix;

//fprintf(stderr,"BEFORE 2) initialization\n");
	// Hessian initialization
	for(i=0; i<sizel*(sizel+1)/2; i++)
		hess_matrix[i]=0.0;
//fprintf(stderr,"AFTER 2) initialization\n");

//fprintf(stderr,"BEFORE 3) Uab\n");
//fprintf(stderr,"BEFORE 4) Uab allocation\n");
	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
	long int Usize=0; // Needed for huge systems!
	long int uij_index,uij1,uij2,uij3,uij4,uij_void; // Needed for huge systems!
	double ***Uij;
	// +1 due to void element (when there are just Traslational/Rotational DoFs)
	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
	{
		allocated_mem += sizeof(double **) * (num_units*(num_units+1)/2);
		for(i=0; i<(num_units*(num_units+1)/2); i++)
			Uij[i] = NULL; // initialization
	}
	else
	{
		printf("Msg(hessianMCAx): Unable to allocate Uij-matrix!\n");
		exit(1);
	}
//fprintf(stderr,"AFTER 4) Uab allocation\n");

	Usize += sizeof(double **) * num_units*(num_units+1)/2;

	if(debug)
		fprintf(stderr,"debug> Memory allocated by hessianMCAx() after Uij initialization: %ld\n",allocated_mem);

	// T memory allocation (6x6 matrix)
	double **T;
	T = (double **) malloc( sizeof(double *) * 6 );
	allocated_mem += sizeof(double *) * 6;
	for(int i=0; i<6; i++)
	{
		allocated_mem += sizeof(double) * 6;
		T[i] = (double *) malloc( sizeof(double) * 6 );
	}

	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
//fprintf(stderr,"BEFORE 5) Uab main-loop\n");
	for(long int b=size-1, c= size+5; b >= 0; b--,c--) // b --> dihedrals (column)
		//	for(int b=size-1, c= size+2; b >= 0; b--,c--) // b --> dihedrals (column)
	{											  // c --> deleteable col. index
		for(long int a=0; a <= b; a++)	 // a --> dihedrals (row)
		{ 							 // it fills upper triangular part (including diagonal)
			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
			uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2;
			// fprintf(stderr,"a= %d  b= %d  undh[a]= %d  undh[b]+1= %d  uij_index= %d\n",a,b,undh[a],undh[b]+1,uij_index);
//fprintf(stderr,"%ld ",uij_index);
			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
			// the Uij for the 4 ICs combination should be computed only once!
			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
			{
				// "in situ" Uij memory allocation
//fprintf(stderr,"BEFORE 6) Uij not-computed\n");

//fprintf(stderr,"BEFORE 7) Uij element allocation\n");
				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
				{
					allocated_mem += sizeof(double *) * 6;
					for(int m=0;m<6;m++)
					{
						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
						{
							allocated_mem += sizeof(double) * 6;
							for(int n=0; n<6; n++)
								Uij[uij_index][m][n] = 0.0;
						}
						else
						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
					}
					Usize += sizeof(double *)*6 + sizeof(double)*36;
				}
				else
				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
//fprintf(stderr,"AFTER 7) Uij element allocation\n");

				// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
				// Computing valid indices

				// U(a,b+1)
				if(undh[b]+1 < num_units-1 ) // if valid "uab1"
					uij1 = undh[a] + (undh[b]+2)*(undh[b]+3)/2; // i= undh[a]  j= (undh[b]+1)+1

				// U(a-1,b)
				if( undh[a] > 0 )
					uij2 = undh[a]-1 + (undh[b]+1)*(undh[b]+2) / 2; // Uab[ dhup[a-1] ][ dhright[b] ]

				// U(a-1,b+1)
				if( undh[a] > 0 && undh[b]+2 < num_units )
					uij3 = undh[a]-1 + (undh[b]+2)*(undh[b]+3) / 2; // Uab[ dhup[a-1] ][ dhright[b+1] ]

				// fprintf(stderr,"\tUij: uij1(row=%d col=%d) uij2(row=%d col=%d) uij3(row=%d col=%d) \n",undh[a],undh[b]+2,undh[a]-1,undh[b]+1,undh[a]-1,undh[b]+2);

				// Computing Uij element
//fprintf(stderr,"BEFORE 8) Uij element computation\n");
				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
				{								 // (whether they have ipas...) --> then, "T" must be computed
					// Computing Tij element "on the fly"
					// (ec.21) Noguti & Go 1983 pp.3685-90
					calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );
					// fprintf(stderr,"\tTij computed\n");

					// (ec.23) Noguti & Go 1983 pp.3685-90
					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij1][n][m]
								                                    + Uij[uij2][n][m]
								                                                   - Uij[uij3][n][m]
								                                                                  + T[n][m];
					}
					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij1][n][m]
								                                    + T[n][m];
					}
					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij2][n][m]
								                                    + T[n][m];
					}
					else // a == 0 && b == size-1  (0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = T[n][m];
					}
				}
				else // if "a" and "b" 's units don't interact
				{	 // (no ipas...) --> then, "T" computation is not necessary!
					// fprintf(stderr,"\tTij not-computed\n");

					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij1][n][m]
								                                    + Uij[uij2][n][m]
								                                                   - Uij[uij3][n][m];
					}
					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij1][n][m];
					}
					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uij_index][n][m] = Uij[uij2][n][m];
					}
					// if no interaction, then no Tij addition to Uab needed!
				}
//fprintf(stderr,"AFTER 8) Uij element computation\n");

				// at this point, the Uij element is fully computed!
			}   // only if not computed yet!
//fprintf(stderr,"AFTER 6) Uij not-computed\n");

			// ****************************************************
			// Rab - Matrix Computation
			// Noguti & Go (1983) pp. 3685-90
			// ****************************************************
			// Rab Computation ( the same Single- or Multi- Chain)
			// Table I. Noguti & Go 1983 pp.3685-90
			// backbone vs. backbone --> Rab = Uab

			// (ec.16) Noguti & Go 1983 pp.3685-90
			// fprintf(stderr,"\tRij= ");
//fprintf(stderr,"BEFORE 9) R from Uij\n");
			for(int n=0; n<6; n++)
			{
				R[n] = 0.0;
				for(int m=0; m<6; m++)
					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
				// fprintf(stderr,"%f ",R[n]);
			}
//fprintf(stderr,"AFTER 9) R from Uij\n");
//fprintf(stderr,"BEFORE 10) Hess from R\n");
			temp = 0.0;
			for(int n=0; n<6; n++)
				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian
//			indexl = a+b*(b+1)/2; // Needed for huge systems!
//fprintf(stderr,"BEFORE 11) Hess writting (a= %ld  b= %ld  index= %ld\n",a,b,a+b*(b+1)/2);
			hess_matrix[a+b*(b+1)/2] = temp; // Rab' x erb = Hessian
//fprintf(stderr,"AFTER 11) Hess writting\n");
			// fprintf(stderr," Hij= %f\n",hess_matrix[a + b*(b+1)/2]);
//fprintf(stderr,"AFTER 10) Hess from R\n");
		}
//fprintf(stderr,"AFTER 5) Uab main-loop\n");

		// Deleting unnecessary colums! (it saves much memory!)
		if( c < size ) // if col. is deleteable
		{
			// Deleting "c" column from Uab (up-diagonal part)
			for(int n=0; n <= c; n++)
			{
				uij_index = undh[n] + (undh[c]+1)*(undh[c]+2)/2; // translates from "squared" to "triangular" matrix elements

				if(Uij[uij_index]!=NULL)
				{
					// Deleting Uab element
					for(int m=0; m<6; m++)
						free( Uij[uij_index][m] );
					free( Uij[uij_index] );
					Uij[uij_index]=NULL;
					Usize -= sizeof(double *)*6 + sizeof(double)*36;
				}
			}
		}
//fprintf(stderr,"AFTER 3) Uab\n");
	}

	// Deleting 6 last Uab columns-rows
	if(size>6)
		//	if(size>3)
		for(int a=0; a<6; a++)
			//	for(int a=0; a<3; a++)
		{
			for(int b=a; b<6; b++)
				//		for(int b=a; b<3; b++)
			{
				uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2; // translates from "squared" to "triangular" matrix elements

				if(Uij[uij_index]!=NULL)
				{
					// Deleting Uab element
					for(int m=0; m<6; m++)
						free( Uij[uij_index][m] );
					free( Uij[uij_index] );
					Uij[uij_index]=NULL;
					Usize -= sizeof(double *)*6 + sizeof(double)*36;
				}
			}
		}

	free( Uij );
	for(i=0; i<num_units*(num_units+1)/2; i++)
		free( unipa[i] );
	free( unipa );

	for(i=0;i<6;i++)
	{
		free( erx[i] ); // erx[6][size]
		free( T[i] );
	}
	free( erx );
	free( T );
	free( undh );
}

// Fast Hessian Matrix computation: Multi-Chain & fix & CA-model & Hard-Disk output. (4/4/2012)
// Minimal memory requirements. (Efficient "Uab" computation, i.e. not using "unipa")
// Input: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// [Ref.] Noguti & Go (1983) pp. 3685-90
void hessianMCAxHD(Macromolecule *mol,twid *decint,int nipa,long int size,float *coord,float *coordCA,
		floating **p_hess_matrix,tri *props, int *unat, bool *fix, char *file, bool justfile)
{
	bool debug = false;
	int j2,m,k1;
	long int i,j;

	double r[3],e[3],R[6];
	double temp;
	int resn,num_res,num_atoms,num_seg;
	int num_res_old; // previous segment number of residues
	int index_res = 0;
	floating *hess_matrix;
	TMOL fragtype;
	twid **sipas; // "Sorted ipas" array
	int isipas = 0; // "Sorted ipas" array index
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	num_atoms = mol->get_num_atoms();
// MON: "num_res" is loaded while iterating each segment, initialization may be avoided...
//	num_res = mol->get_num_fragments();
	num_res = 0;

	if(debug)
		printf("debug> WARNING! Using: hessianMCAxHD()\n");

	if(justfile && file == NULL) // Stupid options...
	{
		printf("Msg(hessianMCAxHD): Sorry, you should introduce some Hessian matrix file name!\nForcing exit!\n");
		exit(2);
	}

	// ********************************************
	// Creating --> (ea, ea x ra)  <-- [equal to (eb, eb x rb)]
	// ********************************************
	double **erx;
	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	for(int i=0;i<6;i++)
		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]

	bool any_rotrans = false;
	int *undh; // returns the first unit-index on the left side of the dihedral
	int un_index = 0;
	undh = (int *) malloc( sizeof(int) * size );

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	j = 0; // current dihedral index
	j2 = 0; // mobile dihedral index
	k1 = 0;
	int k0 = 0; // first-NH + CA + last-CO model index

	int n_inter; // number of mobile intersegment variables
//	int n_rotrans; // number of rotational variables
	int num_units=0;
	int n_fixed=0; // counter of fixed internal coordinates
	int fix_offset=0; // unit-offset (joins units)
	int natom=0;
	bool break_unit = false;
	bool already_broken = false;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary ones: "undh"
	// ---------------------------------------------------
	Segment * seg;
	Atom * atom;
	// Screening segments
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
//		num_res_old = num_res;
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
		fragtype = seg->getMolType(); // This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			n_inter=0; // number of mobile intersegment variables
			for(int n=0; n<6; n++)
				if(fix == NULL || fix[j+n]) // if mobile
					n_inter++;
//			if(num_res_old == 1)
//				n_rotrans = 3;
//			else
//				n_rotrans = 6;
//			for(int n=0; n < n_rotrans; n++)
//				if(fix == NULL || fix[j+n]) // if mobile
//					n_inter++;

//			// Check whether any of the 6D inter-segment variables are fixed (to define rigid units)
//			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
//				any_rotrans = true;
//			else
//				any_rotrans = false;

			any_rotrans = n_inter > 0; // "true" if any inter-segment variable is mobile (to properly define rigid units)

			// ********************
			//  Adding 3 TRASLATIONS
			// ********************
			for(int axis=0; axis<3; axis++) // 3 translations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// ea (Phi)  (0)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					// ea x ra (Psi)  (-gamma_v)
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;
					erx[3+axis][j2] = -1.0;

					undh[j2] = un_index;
					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}
				j++; // adding traslational degree of freedom
			} // 3 trans added

			// ********************
			//  Adding 3 ROTATIONS
			// ********************
//			if(num_res_old != 1)
				for(int axis=0; axis<3; axis++) // 3 rotations
				{
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					{
						// ea (Phi)  (delta_v)
						erx[0][j2] = 0;
						erx[1][j2] = 0;
						erx[2][j2] = 0;
						erx[axis][j2] = 1.0;
						// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
						erx[3][j2] = 0;
						erx[4][j2] = 0;
						erx[5][j2] = 0;

						undh[j2] = un_index;
						j2++;
					}
					else
					{
						fix_offset++;
						n_fixed++;
					}
					j++; // adding traslational degree of freedom
				} // 3 rots added

			// Check whether any of the 6D inter-segment variables are fixed
			if(any_rotrans)
				un_index++; // count units (due to segment ending)
		}

		// Screen ALL fragments
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{
				if(iter_frag->pos_fragment == 0) // has NH
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

				// NOT FIRST, NOT LAST, NOT PRO -> PHI
				// ("PHI-bond")
				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// printf ("%d Not first not PRO -> PHI \n",j);

						// y_lambda (NH pos 0)
						e[0] = -coord[k1 * 3];
						e[1] = -coord[k1 * 3 + 1];
						e[2] = -coord[k1 * 3 + 2];
						// CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// e_lambda
						e[0] += r[0]; // NH --> CA
						e[1] += r[1];
						e[2] += r[2];

						// Unit vector normalization
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );	// Normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp;

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
				}  // NOT FIRST NOT PRO

				if(break_unit)
				{
					un_index++;
					break_unit = false; // by default not breaking unit
				}

				unat[k0] = un_index;
				k0++;
				already_broken = false;

				// NOT LAST, NOT FIRST RESIDUE--> PSI
				// ("PSI-bond")
				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					{
						// get CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
						e[0] = coord[(k1+2) * 3];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= r[0];
						e[1] -= r[1];
						e[2] -= r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] dihedral index
				}  // NOT LAST

				if(iter_frag->pos_fragment == num_res-1) // has CO
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

			} // END PROTEIN FRAGMENT
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(hessianMCAxHD): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			k1+=3; // residue index // first residue atom index
			index_res++; // counts residues (due to "props")
		}
		delete(iter_frag);
	}
	num_units = un_index + 1; // Inter-unit variables (rot-trans) dont increase units number
	delete(iter_seg);

	if(debug)
		printf("Msg(hessianMCAxHD): Number of units: %d (%d)\n",num_units,un_index);

	// ---------------------------------------------------------------------------
	// 2. Sorting "ipas" instead of building "unipa" (mandatory for huge systems!)
	// ---------------------------------------------------------------------------
	int x_nipa = nipa;
	if(debug)
		printf("Msg(hessianMCAxHD): Sorting IPAS...");
	sipas = sort_ipas(decint,&x_nipa,unat,num_units); // sipas --> "Sorted ipas" array
	if(debug)
		printf(" OK!\n");

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------
	FILE *f_file = NULL;
	if(file != NULL) // the Hessian will be dumped into a HD disk file
	{
		if( !(f_file = fopen(file, "w")) ) // if file creation error...
		{
			printf("Msg(hessianMCAxHD): Sorry, Hessian-Matrix file creation failed!\nForcing exit!\n");
			exit(1);
		}
		double dsize = (double) size;
		fwrite(&dsize,sizeof(double),1,f_file); // first number is the matrix rank
	}

	if(!justfile) // the Hessian will be dumped into RAM.
	{
		// Hessian memory allocation
		if(*p_hess_matrix == NULL)
		{
			if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
			{
				printf("Msg(hessianMCAxHD): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
				exit(1);
			}
			*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
		}
		else
			hess_matrix = *p_hess_matrix;

		// Hessian initialization
		if(debug)
			printf("Msg(hessianMCAxHD): hess_matrix size= %ld bytes (size= %ld)\n",(long int) sizeof(floating) * size*(size+1)/2, size);
		for(i=0; i< size*(size+1)/2; i++) // "size" is "long int"
			hess_matrix[i]=0.0;
	}

	// Uabmn minimal-memory allocation (full-square + 2 cols. + Tij "on the fly" computation)
	long int Usize=0;

	// Now indexing is a bit different
	int uij1x,uij1y,uij2x,uij2y,uij3x,uij3y,uijx,uijy;
	long int col, col_buf;

	// Now we will allocate just a (2 cols. x num_units rows.) matrix
	// with 6x6 elements on each position.
	double ****Uij;
	double ***Uitemp;
	if( Uij = (double ****) malloc( sizeof(double ***) * 2 ) ) // Just 2 columns!
		for(int i=0; i<2; i++)
		{
			if( Uij[i] = (double ***) malloc( sizeof(double **) * num_units) ) // rows
				for(int j=0; j<num_units; j++)
					if( Uij[i][j] = (double **) malloc( sizeof(double *) * 6 ) )
					{
						for(int m=0;m<6;m++)
						{
							if( Uij[i][j][m] = (double *) malloc( sizeof(double) * 6 ) )
								for(int n=0; n<6; n++)
									Uij[i][j][m][n] = 0.0; // initialization
							else
							{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
						}
						Usize += sizeof(double *)*6 + sizeof(double)*36;
					}
					else
					{ printf("Msg(hessianM): Unable to allocate an Uij-matrix element!\n"); exit(1); }
			else
			{ printf("Msg(hessianM): Unable to allocate an Uij-matrix column!\n"); exit(1); }
		}
	else
	{
		printf("Msg(hessianMCAxHD): Unable to allocate Uij-matrix (2 x num_units) !\n");
		exit(1);
	}
	Usize += ( sizeof(double **) * num_units * 2 ) + ( sizeof(double ***) * 2 );
	if(debug)
		printf("Msg(hessianMCAxHD): Final Uij size = %ld bytes (%.3f MB)\n",Usize,(float)Usize/1e6);

	// T memory allocation (6x6 matrix)
	double **T;
	T = (double **) malloc( sizeof(double *) * 6 );
	for(int i=0; i<6; i++)
		T[i] = (double *) malloc( sizeof(double) * 6 );

	// Uab Computation Main-Loop (Computing Tij elements "on the fly") <--- (EASILY PARALELIZABLE)
	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
	col_buf = undh[size-1]+1; // chieck this! It is senseless... Initialize col=0.
	for(long int b=size-1; b >= 0; b--) // b --> dihedrals (column)
	{
		col_buf = col; // Uij column buffer
		col = undh[b]+1; // Current Uij column

		if( col != col_buf && b < size-1) // If we changed Uij column (i.e. unit)
		{
			// Displacing Uij one position to the right...
			// First, last element initialization
			for(int j=0; j<num_units; j++)
				for(int m=0;m<6;m++)
					for(int n=0; n<6; n++)
						Uij[1][j][m][n] = 0.0; // initialization

			// Second, displacing 1 position to the right... (Column swapping)
			Uitemp = Uij[1]; // buffer
			Uij[1] = Uij[0]; // current element displacement
			Uij[0] = Uitemp;
		}

		for(long int a=0; a <= b; a++)	 // a --> dihedrals (row)
		{ 							     // it fills upper triangular part (including diagonal)
			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
			//
			//   a0   a1   a2
			// O----O----O----O
			// u0   u1   u2   u3
			//
			// U(a,b)
			uijx = undh[a];
			uijy = 0;

			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
			// the Uij for the 4 ICs combination should be computed only once!
			if( Uij[uijy][uijx][0][0] == 0.0 ) // If "Uij" was not computed yet (compute Uij only once!)
			{
				// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
				// Computing valid indices

				// U(a,b+1)
				if(undh[b]+1 < num_units-1 ) // if valid "uij1"
				{
					uij1x = undh[a]; // i= undh[a] // (row)
					uij1y = 1 ; // j= (undh[b]+1)+1 // (col.)
				}

				// U(a-1,b)
				if( undh[a] > 0 ) // if valid "uij2"
				{
					uij2x = undh[a]-1;
					uij2y = 0;
				}

				// U(a-1,b+1)
				if( undh[a] > 0 && undh[b]+2 < num_units ) // if valid "uij3"
				{
					uij3x = undh[a]-1;
					uij3y = 1;
				}

//				fprintf(stderr,"Current sipas[%d]->i= %d  sipas[%d]->j= %d  undh[a]= %d  undh[b]+1= %d\n",isipas,sipas[isipas]->i,isipas,sipas[isipas]->j,undh[a],undh[b]+1 );

				// Computing Uij element
				if( sipas[isipas]->i == undh[a] && sipas[isipas]->j == undh[b]+1 ) // if "a"'s and "b"'s units interact
				{	// If "i and "j" units have ipas... then, "T" must be computed
					// Computing Tij element "on the fly" from "sorted ipas" array.
					// (ec.21) Noguti & Go 1983 pp.3685-90
					// calcTij_new( T, sipas, &isipas, coordCA);
					calcTij_new( T, sipas, &isipas, coordCA, x_nipa);

//					fprintf(stderr,"Current sipas[%d]->i= %d  sipas[%d]->j= %d  undh[a]= %d  undh[b]+1= %d\n",isipas,sipas[isipas]->i,isipas,sipas[isipas]->j,undh[a],undh[b]+1 );

//					fprintf(stderr,"\nT-matrix: %d vs. %d\n",sipas[isipas]->i,sipas[isipas]->j);
//					for(int n=0;n<6;n++)
//					{
//						for(int m=0;m<6;m++)
//							fprintf(stderr," %f",T[n][m]);
//						fprintf(stderr,"\n");
//					}

					// (ec.23) Noguti & Go 1983 pp.3685-90
					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
					{
//						fprintf(stderr,"middle (using: uij1 uij2 uij3)\n");
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
													   + Uij[uij2y][uij2x][n][m]
													   - Uij[uij3y][uij3x][n][m]
													   + T[n][m];
					}
					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row, but (0,n) (computed below)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
													   + T[n][m];
					}
					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column, but (0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m]
													   + T[n][m];
					}
					else // a == 0 && b == size-1  --> (0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = T[n][m];
					}
				}
				else // if "a" and "b" 's units don't interact
				{	 // (no ipas...) --> then, "T" computation is not necessary!
//					fprintf(stderr,"NO interaccionan\n");
					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
													   + Uij[uij2y][uij2x][n][m]
													   - Uij[uij3y][uij3x][n][m];
					}
					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row, but (0,n) (computed above)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m];
					}
					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column, but (0,n)
					{
						for(int n=0;n<6;n++)
							for(int m=0;m<6;m++)
								Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m];
					}
					// if there is no interaction, then not Tij addition to Uab needed!
				}
			} // At this point, the Uij element is fully computed (if needed).
			// fprintf(stderr,"a=%ld, b=%ld, uijy=%d, uijx=%d, Uij[uijy][uijx][0][0]= %f\n",a,b,uijy,uijx,Uij[uijy][uijx][0][0]);
//			fprintf(stderr,"a=%ld, b=%ld, uijy=%d, uijx=%d, Uij[uijy][uijx][0][0]= %f, T[0][0]= %f\n",a,b,uijy,uijx,Uij[uijy][uijx][0][0],T[0][0]);

			// *********************************************************
			// Rab Computation ( the same Single- or Multi- Chain)
			// Table I. Noguti & Go 1983 pp.3685-90
			// (In CA-model, only "backbone vs. backbone" --> Rab = Uab)
			// *********************************************************

			// (ec.16) Noguti & Go 1983 pp.3685-90
//			fprintf(stderr,"\nR-vector: uijy= %d  uijx= %d\n",uijy,uijx);
			for(int n=0; n<6; n++)
			{
				R[n] = 0.0;
				for(int m=0; m<6; m++)
				{
					R[n] += erx[m][a] * Uij[uijy][uijx][m][n]; // era x Rab
//					fprintf(stderr," %f*%f",erx[m][a],Uij[uijy][uijx][m][n]);
				}
//				fprintf(stderr," %f",R[n]);
			}
//			fprintf(stderr,"\n");

			temp = 0.0;
			for(int n=0; n<6; n++)
				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

			if(!justfile) // the Hessian will be dumped into RAM.
				hess_matrix[a + b*(b+1)/2] = temp; // Rab' x erb = Hessian

//			fprintf(stderr,"a=%ld b=%ld hess= %f\n",a,b,hess_matrix[a + b*(b+1)/2]);
//
//			if(a>=58 && a<=61 && b==943)
//			{
//				fprintf(stderr,"This is the full current Uij:\n");
//				for(int ii=0; ii<6; ii++)
//				{
//					for(int jj=0; jj<6; jj++)
//						fprintf(stderr,"%f ",Uij[uijy][uijx][ii][jj]);
//					fprintf(stderr,"\n");
//				}
//
//				fprintf(stderr,"This is the full current erx[a]:\n");
//				for(int ii=0; ii<6; ii++)
//					fprintf(stderr,"%f ",erx[ii][a]);
//				fprintf(stderr,"\n");
//
//				fprintf(stderr,"This is the full current erx[b]:\n");
//				for(int ii=0; ii<6; ii++)
//					fprintf(stderr,"%f ",erx[ii][b]);
//				fprintf(stderr,"\n");
//			}

			if(file != NULL) // the Hessian will be dumped into a HD disk file
				fwrite(&temp,sizeof(double),1,f_file);
//					printf("%5d  %21.18e\n",a+1,temp);
		}
	} // loop end

	// Free "Uij"
	for(int i=0; i<2; i++)
	{
		for(int j=0; j<num_units; j++)
		{
			for(int m=0;m<6;m++)
				free(Uij[i][j][m]);
			free(Uij[i][j]);
			Usize -= sizeof(double *)*6 + sizeof(double)*36;
		}
		free(Uij[i]);
		Usize -= sizeof(double **)*num_units;
	}
	free(Uij);
	Usize -= sizeof(double ***)*2;

	if(debug)
		printf("Msg(hessianMCAxHD): Final \"after free\" Uij size = %ld bytes (should be =0)\n",Usize);

	for(int i=0;i<6;i++)
	{
		free( erx[i] ); // erx[6][size]
		free( T[i] );
	}
	free( erx );
	free( T );
	free( undh );

	if(file != NULL) // the Hessian will be dumped into a HD disk file
		fclose(f_file); // closing file
}

// Fast Hessian Matrix computation: Multi-Chain & fix & CA-model & Hard-Disk output. (4/4/2012)
// Minimal memory requirements. (Efficient "Uab" computation, i.e. not using "unipa")
// Input: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "nthreads" --> Number of threads to be used.
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// [Ref.] Noguti & Go (1983) pp. 3685-90
void hessianMCAxHD_par(Macromolecule *mol,twid *decint,int nipa,long int size,float *coord,float *coordCA,
		floating **p_hess_matrix,tri *props, int *unat, bool *fix, int nthreads, char *file, bool justfile)
{
	bool debug = false;
	int j2,m,k1;
	long int i,j;

	double r[3],e[3],R[6];
	double temp;
	int resn,num_res,num_atoms,num_seg;
	int index_res = 0;
	floating *hess_matrix;
	TMOL fragtype;
	twid **sipas; // "Sorted ipas" array
	int isipas = 0; // "Sorted ipas" array index
	Residue *res;
	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
	pdbIter *iter_frag; // Iterator to screen fragments
	pdbIter *iter_seg; // Iterator to screen segments
	num_atoms = mol->get_num_atoms();
	num_res = mol->get_num_fragments();

	if(debug)
		printf("debug> WARNING! Using: hessianMCAxHD()\n");

	if(justfile && file == NULL) // Stupid options...
	{
		printf("Msg(hessianMCAxHD): Sorry, you should introduce some Hessian matrix file name!\nForcing exit!\n");
		exit(2);
	}

	// ********************************************
	// Creating --> (ea, ea x ra)  <-- [equal to (eb, eb x rb)]
	// ********************************************
	double **erx;
	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	for(int i=0;i<6;i++)
		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]

	bool any_rotrans = false;
	int *undh; // returns the first unit-index on the left side of the dihedral
	int un_index = 0;
	undh = (int *) malloc( sizeof(int) * size );

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	j = 0; // current dihedral index
	j2 = 0; // mobile dihedral index
	k1 = 0;
	int k0 = 0; // first-NH + CA + last-CO model index

	int n_inter; // number of mobile intersegment variables
	int num_units=0;
	int n_fixed=0; // counter of fixed internal coordinates
	int fix_offset=0; // unit-offset (joins units)
	int natom=0;
	bool break_unit = false;
	bool already_broken = false;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary ones: "undh"
	// ---------------------------------------------------
	Segment * seg;
	Atom * atom;
	// Screening segments
	iter_seg = new pdbIter( mol ); // Iterator to screen segments
	num_seg = iter_seg->num_segment();
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = ( Segment * ) iter_seg->get_segment();
		iter_frag = new pdbIter( seg );
		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
		fragtype = seg->getMolType(); // This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?

		if(iter_seg->pos_segment != 0) // non-first segment
		{
			n_inter=0; // number of mobile intersegment variables
			for(int n=0;n<6;n++)
				if(fix == NULL || fix[j+n]) // if mobile
					n_inter++;

			// Check whether any of the 6D inter-segment variables are fixed (to define rigid units)
			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
				any_rotrans = true;
			else
				any_rotrans = false;

			// ********************
			//  Adding 3 TRASLATIONS
			// ********************
			for(int axis=0; axis<3; axis++) // 3 traslations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// ea (Phi)  (0)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					// ea x ra (Psi)  (-gamma_v)
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;
					erx[3+axis][j2] = -1.0;

					undh[j2] = un_index;
					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}
				j++; // adding traslational degree of freedom
			} // 3 trans added

			// ********************
			//  Adding 3 ROTATIONS
			// ********************
			for(int axis=0; axis<3; axis++) // 3 rotations
			{
				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
				{
					// ea (Phi)  (delta_v)
					erx[0][j2] = 0;
					erx[1][j2] = 0;
					erx[2][j2] = 0;
					erx[axis][j2] = 1.0;
					// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
					erx[3][j2] = 0;
					erx[4][j2] = 0;
					erx[5][j2] = 0;

					undh[j2] = un_index;
					j2++;
				}
				else
				{
					fix_offset++;
					n_fixed++;
				}
				j++; // adding traslational degree of freedom
			} // 3 rots added

			// Check whether any of the 6D inter-segment variables are fixed
			if(any_rotrans)
				un_index++; // count units (due to segment ending)
		}

		// Screen ALL fragments
		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
		{
			res = ( Residue * ) iter_frag->get_fragment();
			resn = resnum_from_resname( res->getName() );

			// PROTEIN FRAGMENT
			if( fragtype == tmol_protein )
			{
				if(iter_frag->pos_fragment == 0) // has NH
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

				// NOT FIRST, NOT LAST, NOT PRO -> PHI
				// ("PHI-bond")
				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
					{
						// printf ("%d Not first not PRO -> PHI \n",j);

						// y_lambda (NH pos 0)
						e[0] = -coord[k1 * 3];
						e[1] = -coord[k1 * 3 + 1];
						e[2] = -coord[k1 * 3 + 2];
						// CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// e_lambda
						e[0] += r[0]; // NH --> CA
						e[1] += r[1];
						e[2] += r[2];

						// Unit vector normalization
						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );	// Normalization factor
						for ( m = 0; m < 3; m++ )
							e[m] /= temp;

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] index
				}  // NOT FIRST NOT PRO

				if(break_unit)
				{
					un_index++;
					break_unit = false; // by default not breaking unit
				}

				unat[k0] = un_index;
				k0++;
				already_broken = false;

				// NOT LAST, NOT FIRST RESIDUE--> PSI
				// ("PSI-bond")
				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
				{
					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
					{
						// get CA pos 1
						r[0] = coord[(k1+1) * 3];
						r[1] = coord[(k1+1) * 3 + 1];
						r[2] = coord[(k1+1) * 3 + 2];
						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
						e[0] = coord[(k1+2) * 3];
						e[1] = coord[(k1+2) * 3 + 1];
						e[2] = coord[(k1+2) * 3 + 2];
						// e_lambda ==> CA --> C (unit vector)
						e[0] -= r[0];
						e[1] -= r[1];
						e[2] -= r[2];

						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
						for ( m = 0; m < 3; m++ )
							e[m] /= temp; // Unit vector normalization

						// ea
						erx[0][j2] = e[0];
						erx[1][j2] = e[1];
						erx[2][j2] = e[2];
						// ea x ra
						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
						erx[5][j2] = e[0] * r[1] - e[1] * r[0];

						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)

						if(!already_broken)
						{
							break_unit = true;
							already_broken = true;
						}

						j2++;
					}
					else
						n_fixed++;

					j++; // der[] dihedral index
				}  // NOT LAST

				if(iter_frag->pos_fragment == num_res-1) // has CO
				{
					unat[k0] = un_index; // bb unit index
					k0++;
				}

			} // END PROTEIN FRAGMENT
			else // NOT FOUND MOL-TYPE
			{
				printf("Msg(hessianMCAxHD): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
				exit(2);
			}

			k1+=3; // residue index // first residue atom index
			index_res++; // counts residues (due to "props")
		}
		delete(iter_frag);
	}
	num_units = un_index + 1; // Inter-unit variables (rot-trans) don't increase units number
	delete(iter_seg);

	if(debug)
		printf("Msg(hessianMCAxHD): Number of units: %d (%d)\n",num_units,un_index);

	// ---------------------------------------------------------------------------
	// 2. Sorting "ipas" instead of building "unipa" (mandatory for huge systems!)
	// ---------------------------------------------------------------------------
	int x_nipa = nipa; // Because sort_ipas() removes intra-unit contacts (units are rigid!)
	if(debug)
		printf("Msg(hessianMCAxHD): Sorting IPAS...");
	// Sorting "ipas", removeing intra-unit contacts, and setting units ownership by using "unat".
	sipas = sort_ipas(decint,&x_nipa,unat,num_units); // sipas --> "Sorted ipas" array
	if(debug)
		printf(" OK!\n");

	// 2.1. In parallel, interacting pair of atoms list ("sipas") should be split into chunks corresponding to interactions belonging exclusively to a single
	//      Hessian (a,b) box.
	//      Sorting should be maintained, so it will be screened like during Hessian computation.
	//      To avoid memory overflows in the following sections, keeping note about how many "ipas" belong to each box is mandatory! This is done in "nasipas" array.
	//
	bool debug2 = true; // show some debug-info related to this section
	int icbox = 100; // Hessian parallel computation box size in internal coordinates (a,b-dimensions: rows,cols)
	int iasipas,k,l,nbox,maxbox;
	nbox = (long int) ceil((double)size/icbox); // number of boxes per box-matrix dimension ("ceiled")
	maxbox = nbox*(nbox+1)/2; // maximum number of boxes (or jobs)

	twid ***asipas; // Array with pointers to an array of "sorted ipas" belonging only to a single Hessian box.
	int *nasipas; // Array storing the number of "ipas" within each "asipas" element
	if( !(asipas = (twid ***) malloc( sizeof(twid **) * maxbox ) ) ) // Triangular-Packed-Storage matrix
	{
		fprintf(stderr,"Sorry, asipas malloc failed!\n");
		exit(1);
	}
	if( !(nasipas = (int *) malloc( sizeof(int) * maxbox ) ) )
	{
		fprintf(stderr,"Sorry, nasipas malloc failed!\n");
		exit(1);
	}
	for(int n = 0; n < maxbox; n++)
	{	// Initialization
		asipas[n] = NULL;
		nasipas[n] = 0;
	}

	for(long int b=size-1; b >= 0; b--) // b --> dihedrals (column)
		for(long int a=0; a <= b; a++)	// a --> dihedrals (row)
			if( sipas[isipas]->i == undh[a] && sipas[isipas]->j == undh[b]+1 ) // if "a"'s and "b"'s units interact
			{
				// a + b*(b+1)/2
				k = a / icbox; // [0,icbox-1] --> 1st row, ...
				l = nbox - ((size-1)-b)/icbox - 1; // (boxes are taken into account from most right column)

				iasipas = k + l*(l+1)/2; // asipas index

				// This is efficient given "ipas" were already sorted into "sipas"
				do {// while "sipas"'s "ipas" belong to the same interacting pair of units: i,j
					nasipas[iasipas]++;
					if(asipas[iasipas] == NULL) // Not allocated yet, then allocating first element
					{
						if( !(asipas[iasipas] = (twid **) malloc( sizeof(twid *) ) ) ) // allocating first block of elements
						{
							fprintf(stderr,"Sorry, first asipas[iasipas] element allocation failed!\n");
							exit(1);
						}
					}
					else // Realloc...
						if(!(asipas[iasipas] = ( twid **) realloc(asipas[iasipas] , nasipas[iasipas] * sizeof( twid *)))) // resizes contact-list structure
						{
							fprintf(stderr,"Sorry, realloc failed!\n");
							exit(1);
						}

					asipas[iasipas][nasipas[iasipas]-1] = sipas[isipas]; // "sipas" are referenced to current box "asipas"
					isipas++; // next "sipa"
				}
//				while( (sipas[isipas]->i == sipas[isipas-1]->i) && (sipas[isipas]->j == sipas[isipas-1]->j) );
				while( isipas < x_nipa && (sipas[isipas]->i == sipas[isipas-1]->i) && (sipas[isipas]->j == sipas[isipas-1]->j) );
				// while the interacting pair of atoms belong to the same pair of units
			}
	isipas = 0; // reset isipas
	// "asipas" computation finished...

//	// Showing "asipas"...
//	if(debug2)
//	{
//		printf("\nx_nipa= %d (number of nipas after intra-unit redundance)\n",x_nipa);
//		for(int y = nbox-1; y >= 0; y--) // l (cols.)
//			for(int x = 0; x <= y; x++) // k (rows)
//			{
//				iasipas = x + y*(y+1)/2; // asipas index
//				printf("Showing contacts for: k= %d  l= %d\n",x,y);
//				for(int n=0; n<nasipas[iasipas]; n++)
//					printf("\t(k,l)=(%4d,%4d) --> %4d %4d %4d %4d %8.3f %8.3f\n",x,y,asipas[iasipas][n]->i,asipas[iasipas][n]->j,asipas[iasipas][n]->k,
//							asipas[iasipas][n]->l,asipas[iasipas][n]->d,asipas[iasipas][n]->C );
//			}
//	}
//	exit(0);

	// 2.2.	In parallel, a "readiness" matrix should be created ("isready"). It will tell whether a Hessian-Box(k,l) is ready to be computed, i.e. if needed "U-segments"
	//		were either computed or not-needed (near the borders).
	//
	bool *isready; // This bool array tells if box(k,l) is ready to be computed.
	bool *isdone; // This bool array tells if box(k,l) was already computed.
	double ****Ut,****Ur; // Ut = U-top, Ur = U-right. Both matrices will store corresponding U interaction elements fore every Hessian-Box.

	// Memory allocation
	if( !(isready = (bool *) malloc( sizeof(bool) * maxbox )) ) // Triangular-Packed-Storage matrix
	{
		printf("Msg(hessianMCAxHD): I'm sorry, memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(isdone = (bool *) malloc( sizeof(bool) * maxbox )) ) // Triangular-Packed-Storage matrix
	{
		printf("Msg(hessianMCAxHD): I'm sorry, memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(Ut = (double ****) malloc( sizeof(double ***) * maxbox )) ) // Triangular-Packed-Storage matrix
	{
		printf("Msg(hessianMCAxHD): I'm sorry, memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(Ur = (double ****) malloc( sizeof(double ***) * maxbox )) ) // Triangular-Packed-Storage matrix
	{
		printf("Msg(hessianMCAxHD): I'm sorry, memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Initialization
	for(int n = 0; n < maxbox; n++)
	{
		isready[n] = false;
		isdone[n] = false;
		Ut[n] = NULL;
		Ur[n] = NULL;
	}

	// 2.3 Some memory allocation, initialization and file opening (if necessary)
	//
	FILE *f_file = NULL;
	if(file != NULL) // the Hessian will be dumped into a HD disk file
	{
		if( !(f_file = fopen(file, "w")) ) // if file creation error...
		{
			printf("Msg(hessianMCAxHD): Sorry, Hessian-Matrix file creation failed!\nForcing exit!\n");
			exit(1);
		}
		double dsize = (double) size;
		fwrite(&dsize,sizeof(double),1,f_file); // first number is the matrix rank
	}
	if(!justfile) // the Hessian will be dumped into RAM.
	{
		// Hessian memory allocation
		if(*p_hess_matrix == NULL)
		{
			if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
			{
				printf("Msg(hessianMCAxHD): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
				exit(1);
			}
			*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
		}
		else
			hess_matrix = *p_hess_matrix;

		// Hessian initialization
		if(debug)
			printf("Msg(hessianMCAxHD): hess_matrix size= %ld bytes (size= %ld)\n",(long int) sizeof(floating) * size*(size+1)/2, size);
		for(i=0; i< size*(size+1)/2; i++) // "size" is "long int"
			hess_matrix[i]=0.0;
	}

	// 2.4 Some mutex related stuff...
	int jobs_done = 0; // The number of jobs done.
	pthread_mutex_t mtx_isready; // Mutex to update "isready" matrix.
	pthread_cond_t cond_isready; // Condition to check "isready"

	hessianCApar_data *threads_data; // Array to store threads data
	pthread_t *threads; // Array of threads handles
	int return_code; // Integer return code for thread creation

	// Initialize mutex and condition variable objects
	pthread_mutex_init(&mtx_isready, NULL);
	pthread_cond_init (&cond_isready, NULL);

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------

	// The multi-threaded routine should be called here...
	// 0. Get Input: "a0","b0","boxsize","asipas", "isready", "Ut", "Ur", "hess_matrix", "undh", "size".
	// 1. Allocate Box's Uij (it is smaller).
	// 2. Store properly "top" (Ut) and "right" (Ur) elements of neighboring Boxes.
	// 3. Compute Hessian
	// 4. Once finished, check "isready" elements belonging to neighboring Boxes.
	// 5.a. MON's way: Each child thread launches new box computation (left and down).
	// 5.b. PABLO's way: Check list of "ready" boxes, and run the first one ready that you find.

	// Allocating threads data
	if( !(threads_data = (hessianCApar_data *) malloc(sizeof(hessianCApar_data) * nthreads)) )
	{
		fprintf(stderr,"Msg(hessianMCAxHD_par): I'm sorry, thread memory allocation failed!\nForcing exit!\n");
		exit(1);
	}
	if( !(threads = (pthread_t *) malloc(sizeof(pthread_t) * nthreads)) )
	{
		fprintf(stderr,"Msg(hessianMCAxHD_par): I'm sorry, thread allocation failed!\nForcing exit!\n");
		exit(1);
	}

	// Setting top-right box as ready (this seeds the whole process)
	isready[(nbox-1)*(nbox)/2] = true;

	// CREATING THREADS...
	if(debug2)
		fprintf(stderr,"\nMsg(hessianMCAxHD_par): Launching threads...\n");
	for(int i = 0; i < nthreads; i++)
	{
		// Loading initial data for threads...
		threads_data[i].k = i; // Top-right row box.
		threads_data[i].l = nbox-1; // Top-right column box.
		threads_data[i].icbox = icbox; // Box size.
		threads_data[i].asipas = asipas;
		threads_data[i].nasipas = nasipas; // number of sipas within each box
		threads_data[i].nipa = x_nipa;
		threads_data[i].isready = isready;
		threads_data[i].isdone = isdone;
		threads_data[i].Ut = Ut;
		threads_data[i].Ur = Ur;
		threads_data[i].hess_matrix = hess_matrix;
		threads_data[i].size = size;
		threads_data[i].undh = undh;
		threads_data[i].erx = erx; // the "erx" array
		threads_data[i].jobs_done = &jobs_done;
		threads_data[i].coord = coordCA; // Current atomic model array of cartesian coordinates.
		threads_data[i].mtx_isready = &mtx_isready;
		threads_data[i].cond_isready = &cond_isready;

		// Creating threads...
		fprintf(stderr,"Msg(hessianMCAxHD_par): Creating thread %2d (%ld)...",i,pthread_self());
		return_code = pthread_create( &threads[i], NULL, hessianMCAxHD_thread, (void *) &threads_data[i] );
		fprintf(stderr,"done!\n");
		if(return_code )
		{
			fprintf(stderr,"Msg(hessianMCAxHD_par): ERROR; return code from pthread_create() is %d for thread %d\n",return_code,i);
			exit(-1);
		}
	}

	pthread_mutex_lock(&mtx_isready);
	while(jobs_done < maxbox)
	{
		pthread_cond_wait(&cond_isready, &mtx_isready); // Waits until some thread sends a signal
		if(debug2)
		{
			fprintf(stderr,"Control signal received, checking jobs: jobs_done= %d\n",jobs_done);
			fprintf(stderr,"Showing \"isdone\" matrix:\n");
			for(int x=0; x < nbox; x++)	 // x --> box (row)
			{
				for(int y=0; y < nbox; y++) // y --> box (column)
					if(y>=x)
						if(isdone[x + y*(y+1)/2])
							fprintf(stderr,"1 ");
						else
							fprintf(stderr,"0 ");
					else
						fprintf(stderr,"  ");
				fprintf(stderr,"\n");
			}
		}
	}
	if(debug2)
		fprintf(stderr,"All jobs completed! jobs_done= %d\n",jobs_done);
	pthread_mutex_unlock(&mtx_isready);

	// Destroy mutex related stuff...
	pthread_mutex_destroy(&mtx_isready);
	pthread_cond_destroy(&cond_isready);

	free(threads_data);
	free(threads);

	// Freeing common to all threads memory...
	for(int i=0;i<6;i++)
		free( erx[i] ); // erx[6][size]
	free( erx );
	free( undh );

	if(file != NULL) // the Hessian will be dumped into a HD disk file
		fclose(f_file); // closing file

}


// Multi-threaded routine to compute a "Hessian-Box"
void *hessianMCAxHD_thread(void *targ)
{
	// Variables definition
	bool debug = true; // Dump some debugging info...
	int return_code;
	bool readyfound;
	int cbox; // current box-index
	int tbox; // top box-index
	int rbox; // right box-index
	int lbox; // left box-index
	int bbox; // bottom box-index
	int nbox; // Maximum number of boxes per matrix dimension
	int icbox; // Internal coordinates boxsize.
	int k,l; // Current box indices
	int ut_low,ut_high; // Current top unit min and max indices
	int ur_low,ur_high; // Current right unit min and max indices
	int ut_size,ur_size; // Current top and right unit sizes
	int col, col_buf;
	long int a,b,ai,bi,af,bf; // Some IC indices
	long int size; // Total number of ICs
	int isipas = 0; // "Sorted ipas" array index
	twid ***asipas; // Array with pointers to an array of "sorted ipas" belonging only to a single Hessian box.
	int *nasipas; // number of sipas whithin each box
	int nipa; // Number of "ipas".
	twid **sipas; // array of "sorted ipas" belonging only to a single Hessian box.
	bool *isready; // This bool array tells if box(k,l) is ready to be computed.
	bool *isdone; // This bool array tells if box(k,l) was already computed.
	double ****Ut,****Ur; // Ut = U-top, Ur = U-right. Both matrices will store corresponding U interaction elements fore every Hessian-Box.
	double *hess_matrix; // The Hessian matrix.
	double **erx; // the "erx" array
	double temp;
	double R[6];
	double **T; // T matrix (accounts for the interaction between units)
	// double T[6][6]; // T matrix (accounts for the interaction between units)
	float *coord; // atomic model cartesian coordinates array
	int *undh; // This returns the first unit-index on the left side of the dihedral.
	int num_units; // Maximum number of units
	int *available_cores; // The number of free processing cores.
	int *jobs_done; // The number of jobs done.
	pthread_mutex_t *mtx_isready; // Mutex to update "isready" matrix.
	pthread_cond_t *cond_isready; // Condition to check "jobs_done"

	// Loading data from "thread argument"
	k = ((hessianCApar_data *) targ)->k;
	l = ((hessianCApar_data *) targ)->l;
	icbox = ((hessianCApar_data *) targ)->icbox;
	asipas = ((hessianCApar_data *) targ)->asipas;
	nasipas = ((hessianCApar_data *) targ)->nasipas;
	nipa = ((hessianCApar_data *) targ)->nipa;
	isready = ((hessianCApar_data *) targ)->isready;
	isdone = ((hessianCApar_data *) targ)->isdone;
	Ut = ((hessianCApar_data *) targ)->Ut;
	Ur = ((hessianCApar_data *) targ)->Ur;
	hess_matrix = ((hessianCApar_data *) targ)->hess_matrix;
	size = ((hessianCApar_data *) targ)->size;
	undh = ((hessianCApar_data *) targ)->undh;
	erx =  ((hessianCApar_data *) targ)->erx;
	jobs_done = ((hessianCApar_data *) targ)->jobs_done;
	coord = ((hessianCApar_data *) targ)->coord;
	mtx_isready = ((hessianCApar_data *) targ)->mtx_isready;
	cond_isready = ((hessianCApar_data *) targ)->cond_isready;

	nbox = (long int) ceil((double)size/icbox); // number of boxes per box-matrix dimension ("ceiled")
	num_units = undh[size-1] + 2; // Inter-unit variables (rot-trans) don't increase units number

	fprintf(stderr,"nbox= %d  icbox= %d  num_units= %d\n",nbox,icbox,num_units);

	//
	// T memory allocation (6x6 matrix)
	if( (T = (double **) malloc( sizeof(double *) * 6 ) ) )
	{
		for(int n=0; n<6; n++)
			if( !(T[n] = (double *) malloc( sizeof(double) * 6 )) )
			{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an T-matrix element!\n"); exit(1); }
	}
	else
	{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an T-matrix element!\n"); exit(1); }

	// TASKS:
	// 1. Allocate Box's Uij (now it is smaller).
	// 2. Compute Hessian while storing properly "top" (Ut) and "right" (Ur) elements.
	// 3. Once finished, check "isready" elements belonging to neighboring Boxes.
	// 4. Check list of "ready" boxes, and process the first one ready you find.

	do {
		// Screening current "jobs-list" searching for a "ready-job".
		readyfound = false;
		while( !readyfound )
		{
			pthread_mutex_lock(mtx_isready); // Locking "isready"
			for(int x=0; x < nbox && !readyfound; x++)	 // x --> box (row)
				for(int y=nbox-1; y >= x && !readyfound; y--) // y --> box (column)
				{
					cbox = x + y*(y+1)/2; // box index
					if(isready[cbox]) // If current box is ready, current thread can do the job.
					{
						isready[cbox] = false; // so un-check it and work...
						readyfound = true;
						if(isdone[cbox]) // If box was already computed... then this is just a mistake!
						{
							fprintf(stderr,"Sorry thread (%d,%d) parallelization logic error in (%d,%d)...\nForcing exit!\n",k,l,x,y);
							exit(-1);
						}
						k = x;
						l = y;
					}
				}

			if(!readyfound) // If there is not any box to be computed yet... then waiting...
			{
				fprintf(stderr,"Thread %ld waiting ...\n",pthread_self());
				pthread_cond_wait(cond_isready, mtx_isready); // Waits until some thread sends a signal
				fprintf(stderr,"Condition signal received by thread %ld, re-checking for ready boxes.\n",pthread_self());
			}
			pthread_mutex_unlock(mtx_isready); // Un-locking "isready"
		}
		// At this point a new valid box (k,l updated) is going to be processed...

		// Computing the Hessian (doing the job)...
		if(false)
		{
			// Some fake job for testing...
			float rantime;
			rantime = 1 + 3*(float)rand()/RAND_MAX;
			fprintf(stderr,"Hi, I'm thread %ld (%d,%d). Job: waiting for some time... %f\n",pthread_self(),k,l,rantime);
			usleep(rantime*1000);
			fprintf(stderr,"Thread %ld (%d,%d) finished!\n",pthread_self(),k,l);
		}
		else
		{
			// BEGIN MEAT ------------------------------------------------------------------------------------------------------------
			cbox = k + l*(l+1)/2; // current box index
			tbox = (k-1) + l*(l+1)/2; // top box index
			rbox = k + (l+1)*(l+2)/2; // right box index
			lbox = k + (l-1)*l/2; // left box index
			bbox = (k+1) + l*(l+1)/2; // bottom box index

			// Translating box limits from IC (a,b) to units (i,j). And computing sizes.
			a = k*icbox; // initial IC row (top)
			b = (size-1) - (nbox-(l+1))*icbox; // initial IC col (boxes are taken into account from most right column)
			if(b-icbox >= 0)
				ut_low = undh[b-icbox]; // Lowest top unit index for IC b. (most right side)
			else
				ut_low = undh[0]; // Lowest top unit index for IC b. (most right side column)
			if(b >= 0)
				ut_high = undh[b]; // Highest top unit index for IC b. (most left side column)
			else
				ut_high = undh[0]; // Highest top unit index for IC b. (most left side column)

			if(a < size-1)
				ur_low = undh[a]; // Lowest right unit index for IC a.
			else
				ur_low = undh[size-1]; // Lowest right unit index for IC a.
			if(a+icbox < size)
				ur_high = undh[a+icbox]; // Highest right unit index for IC a.
			else
				ur_high = undh[size-1]; // Highest right unit index for IC a.
			ut_size = ut_high - ut_low + 1; // Top unit size
			ur_size = ur_high - ur_low + 2; // Right unit size (+1 to account for (a-1,b+1) elements)
			sipas = asipas[cbox]; // Current box "sipas"

			fprintf(stderr,"Thread %ld (%d,%d) start! \n",pthread_self(),k,l);
			fprintf(stderr,"a=%ld b=%ld ut_low=%d ut_high=%d ur_low=%d ur_high=%d ut_size=%d ur_size=%d\n",a,b,ut_low,ut_high,ur_low,ur_high,ut_size,ur_size);

			// 1. Allocating current box Uij, U-top and U-right), and T matrix.
			//
			// Uabmn minimal-memory allocation (full-square + 2 cols. + Tij "on the fly" computation)
			// Now indexing is a bit different
			int uij1x,uij1y,uij2x,uij2y,uij3x,uij3y,uijx,uijy;
			int col,col_buf;
			// Now we will allocate just a (2 cols. x num_units rows.) matrix
			// with 6x6 elements on each position.
			double ****Uij;
			double ***Uitemp; // Column temporal buffer
			if( (Uij = (double ****) malloc( sizeof(double ***) * 2 )) ) // Just 2 columns!
				for(int i=0; i<2; i++)
				{
					if( (Uij[i] = (double ***) malloc( sizeof(double **) * ur_size)) ) // rows
						for(int j=0; j < ur_size; j++)
							if( (Uij[i][j] = (double **) malloc( sizeof(double *) * 6 )) )
								for(int m=0;m<6;m++)
								{
									if( (Uij[i][j][m] = (double *) malloc( sizeof(double) * 6 )) )
										for(int n=0; n<6; n++)
											Uij[i][j][m][n] = 0.0; // initialization
									else
									{ fprintf(stderr,"Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
								}
							else
							{ fprintf(stderr,"Msg(hessianM): Unable to allocate an Uij-matrix element!\n"); exit(1); }
					else
					{ fprintf(stderr,"Msg(hessianM): Unable to allocate an Uij-matrix column!\n"); exit(1); }
				}
			else
			{ fprintf(stderr,"Msg(hessianMCAxHD): Unable to allocate Uij-matrix (2 x num_units) !\n"); exit(1);	}
			//
			// Current thread is who allocates memory for top and right U's.
			//
			// Bottom-box's U-top...
			if(k < nbox-1) // Only there exists bottom-box for current box when current box is not in the last row (k=nbox-1).
			{
				// Allocating memory for left-box's U-top
				if( !(Ut[bbox] = (double ***) malloc(sizeof(double **) * ut_size ) ) )
				{
					fprintf(stderr,"Msg(hessianMCAxHD_thread): Sorry, memory allocation error for: Ut. Forcing exit!\n");
					exit(1);
				}
				for(int n=0; n < ut_size; n++)
					if( (Ut[bbox][n] = (double **) malloc( sizeof(double *) * 6 ) ) )
						for(int m=0;m<6;m++)
							if( (Ut[bbox][n][m] = (double *) malloc( sizeof(double) * 6 ) ) )
								for(int x=0; x<6; x++)
									Ut[bbox][n][m][x] = 0.0; // initialization
							else
							{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an Uab-matrix element!\n"); exit(1); }
					else
					{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an Uij-matrix element!\n"); exit(1); }
			}
			//
			// Left-box's U-right...
			if(l > 0) // Only there exists left-boxes for current box when current box is not in the most-left column (l=0).
			{
				// Allocating memory for bottom-box's U-right
				if( !(Ur[lbox] = (double ***) malloc(sizeof(double **) * ur_size ) ) )
				{
					fprintf(stderr,"Msg(hessianMCAxHD_thread): Sorry, memory allocation error for: Ur. Forcing exit!\n");
					exit(1);
				}
				for(int n=0; n < ur_size; n++)
					if( (Ur[lbox][n] = (double **) malloc( sizeof(double *) * 6 ) ) )
						for(int m=0;m<6;m++)
							if( (Ur[lbox][n][m] = (double *) malloc( sizeof(double) * 6 ) ) )
								for(int x=0; x<6; x++)
									Ur[lbox][n][m][x] = 0.0; // initialization
						else
						{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an Uab-matrix element!\n"); exit(1); }
				else
				{ fprintf(stderr,"Msg(hessianMCAxHD_thread): Unable to allocate an Uij-matrix element!\n"); exit(1); }
			}

			// 2. Compute Hessian while storing properly "top" (Ut) and "right" (Ur) elements.

			// Uab Computation Main-Loop (Computing Tij elements "on the fly") <--- (EASILY PARALLELIZABLE)
			// (Which "unipa" belongs to the "outer" units corresponding to "a" and "b" dihedrals)

			// col_buf = undh[size-1]+1; // chieck this! It is senseless... Initialize col=0.
			// col_buf = undh[(l*icbox)-1]+1; // should col=0 ???
			//col = 0;
			// From box first row (ai) to the bottom most part of it (af), i.e. from top to bottom.
			ai = k*icbox; // initial IC row (top)
			af = ((k+1)*icbox)-1; // final IC row (bottom)
			// From box right column (bi) to the left most part of the box (bf), i.e. from right to left.
			bi = (size-1) - ((nbox-(l+1))*icbox); // initial IC col (most right)
			bf = size - ((nbox-l)*icbox); // final IC col (most left)

			col = undh[bi]+1; // CHECK THIS!

			// Loading current-box U-right into Uij... (something better would have been done here?)
			//
			// Given in CA-model there are two consecutive dihedrals between units, sometimes happens that a U column has been already computed
			// when changing from one box into another (in IC). In this case, no Uij initialization or column switching is needed!
			// Taking into account this, if initial "b" (bi) and previous "b" (bi-1) belong to the same unit, then column switching is not needed
			// and Ur should be inserted in Uij[0]. On the contrary, if initial and previous "b" belong to different units, Ur should be inserted
			// in Uij[1] and Uij[0] should be initialized with 0's.
			//
			if(l < nbox-1) // if there exist U-right in current-box
			{
				if( bi < size-1 ) // if not right most IC column...
				{
					col = undh[bi+1]+1; // previous "b" unit (at right-side of bond "b")

					if( undh[bi] == undh[bi+1] ) // If initial and previous "b" belong to the same unit
					{
						// Loading current-box U-right into Uij[0]
						for(int n=0; n < ur_size; n++)
							for(int m=0;m<6;m++)
								for(int x=0; x<6; x++)
									Uij[0][n][m][x] = Ur[cbox][n][m][x]; // copy Ur into Uij[0]
					}
					else  // If initial and previous "b" belong to different unit
					{
						// Loading current-box U-right into Uij[1] (Uij[0] has been already initialized)
						for(int n=0; n < ur_size; n++)
							for(int m=0;m<6;m++)
								for(int x=0; x<6; x++)
									Uij[1][n][m][x] = Ur[cbox][n][m][x]; // copy Ur into Uij[1]
					}
				}
//				else
//					col = undh[size-1]+1; // right most unit index
			}
			fprintf(stderr,"ai=%ld bi=%ld af=%ld bf=%ld col=%d\n",ai,bi,af,bf,col);

			// SCREENING COLUMNS...
			// for(long int b=size-1; b >= 0; b--) // b --> dihedrals (column)
			for(long int b = bi; b >= bf; b--) // screen "b" ICs (columns) within current box range
			{
				// Loading current-box U-top into Uij...
				//
//				if(k > 0 && undh[b] != undh[b+1] ) // if there exist U-top in current-box, and if initial and previous "b" belong to different unit.
//				{
//					for(int m=0;m<6;m++)
//						for(int x=0; x<6; x++)
//							Uij[0][0][m][x] = Ut[cbox][ undh[bi]-undh[b] ][m][x]; // copy Ut[cbox][curr.col.] into Uij[0][0]
//				}

				col_buf = col; // Uij column buffer
				col = undh[b]+1; // Current Uij column

				// if( col != col_buf && b < size-1) // If we changed Uij column (i.e. unit)
				if( col != col_buf && b < size-1) // If we had changed Uij column (i.e. unit) and we were not in the right most column of current box.
				{
					fprintf(stderr,"Changed Uij column and not in current right most column\n");

					fprintf(stderr,"Deleting previous column...\n");
					// Displacing Uij one position to the right...
					// Last element initialization (or erasing)
					for(int j = 0; j < ur_size; j++)
						for(int m=0; m<6; m++)
							for(int n=0; n<6; n++)
								Uij[1][j][m][n] = 0.0; // initialization

					fprintf(stderr,"Swapping column...\n");
					// Second, displacing 1 position to the right... (Column swapping)
					Uitemp = Uij[1]; // buffer
					Uij[1] = Uij[0]; // current element displacement
					Uij[0] = Uitemp;

//					// Copy current-box Ut-element (undh[bi]-undh[b]) into Uij[0] only when col. changes
//					if(a > 0)
//						if(undh[a] != undh[a-1]) // if "a" belongs to a different unit than "a-1", then Ut-element should go to Uij[0][0]
//							for(int m=0;m<6;m++)
//								for(int x=0; x<6; x++)
//									Uij[0][0][m][x] = Ut[cbox][undh[bi]-undh[b]][m][x];
//						else // Otherwise, it should go to Uij[0][1]. (This way it will not be computed again)
//							for(int m=0;m<6;m++)
//								for(int x=0; x<6; x++)
//									Uij[0][1][m][x] = Ut[cbox][undh[bi]-undh[b]][m][x];

				}

				// Copy current-box Ut-element (undh[bi]-undh[b]) into Uij[0] only when "a" changes unit.
				if(k > 0)
					if(undh[ai] != undh[ai-1]) // if "a" belongs to a different unit than "a-1", then Ut-element should go to Uij[0][0]
						for(int m=0;m<6;m++)
							for(int x=0; x<6; x++)
								Uij[0][0][m][x] = Ut[cbox][undh[bi]-undh[b]][m][x];
					else // Otherwise, it should go to Uij[0][1]. (This way it will not be computed again)
						for(int m=0;m<6;m++)
							for(int x=0; x<6; x++)
								Uij[0][1][m][x] = Ut[cbox][undh[bi]-undh[b]][m][x];


				fprintf(stderr,"BEGIN column b=%ld\n",b);

				// SCREENING ROWS...
				// for(long int a=0; a <= b; a++)	 // a --> dihedrals (row)
				for(long int a = ai; a <= af && a <= b; a++)	 // a --> dihedrals (row)
				{ 								     // it fills upper triangular part (including diagonal)
					fprintf(stderr,"\nBegin row a=%ld (col. b=%ld) ",a,b);

					// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
					//
					//   a0   a1   a2
					// O----O----O----O
					// u0   u1   u2   u3
					//
					// U(a,b)
					// uijx = undh[a];
					uijx = undh[a] - undh[ai] + 1; // (current row). The "+1" offset will account further for "(a-1)" units.
					uijy = 0; // current col.
					fprintf(stderr,"uijx= %d  uijy= %d ",uijx,uijy);
					fprintf(stderr,"Uij[uijy][uijx][0][0]= %f, ",Uij[uijy][uijx][0][0]);

					// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
					// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
					// the Uij for the 4 ICs combination should be computed only once!
					if( Uij[uijy][uijx][0][0] == 0.0 ) // If "Uij" was not computed yet (compute Uij only once!)
					{
						// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
						// Computing valid indices

						fprintf(stderr,"uijx=%d uijy=%d ",uijx,uijy);

						// U(a,b+1)
						if(undh[b]+1 < num_units-1) // if valid "uij1"
						{
							// uij1x = undh[a]; // i= undh[a] // (row)
							uij1x = undh[a] - undh[ai] + 1; // (row) (the "+1" is to account for the "(a-1)" unit)
							uij1y = 1 ; // j= (undh[b]+1)+1 // (col.)
							fprintf(stderr,"uij1x=%d uij1y=%d ",uij1x,uij1y);
						}

						// U(a-1,b)
						if( undh[a] > 0 ) // if valid "uij2"
						{
							// uij2x = undh[a]-1;
							uij2x = undh[a] - undh[ai];
							uij2y = 0;
							fprintf(stderr,"uij2x=%d uij2y=%d ",uij2x,uij2y);
						}

						// U(a-1,b+1)
						if( undh[a] > 0 && undh[b]+2 < num_units ) // if valid "uij3"
						{
							// uij3x = undh[a]-1;
							uij3x = undh[a] - undh[ai];
							uij3y = 1;
							fprintf(stderr,"uij3x=%d uij3y=%d ",uij3x,uij3y);
						}
//						fprintf(stderr,"uijx=%d uijy=%d uij1x=%d uij1y=%d uij2x=%d uij2y=%d uij3x=%d uij3y=%d\n",uijx,uijy,uij1x,uij1y,uij2x,uij2y,uij3x,uij3y);
						fprintf(stderr,"\n");

//						if(nasipas[cbox] > 0 && isipas < nasipas[cbox])
//							fprintf(stderr,"Current sipas[%d]->i= %d  sipas[%d]->j= %d  undh[a]= %d  undh[b]+1= %d\n",isipas,sipas[isipas]->i,isipas,sipas[isipas]->j,undh[a],undh[b]+1 );

						// Computing Uij element
//						if( nasipas[cbox] > 0 && isipas < nasipas[cbox] && sipas[isipas]->i == undh[a] && sipas[isipas]->j == undh[b]+1 ) // if "a"'s and "b"'s units interact
						if( isipas < nasipas[cbox] && sipas[isipas]->i == undh[a] && sipas[isipas]->j == undh[b]+1 ) // if "a"'s and "b"'s units interact
						{	// If "i and "j" units have ipas... then, "T" must be computed
							fprintf(stderr,"Current sipas[%d]->i= %d  sipas[%d]->j= %d  undh[a]= %d  undh[b]+1= %d\n",isipas,sipas[isipas]->i,isipas,sipas[isipas]->j,undh[a],undh[b]+1 );
							fprintf(stderr,"Units corresponding to a=%ld and b=%ld Interact (nasipas=%d): ",a,b,nasipas[cbox]);

							// Computing Tij element "on the fly" from "sorted ipas" array.
							// (ec.21) Noguti & Go 1983 pp.3685-90
							calcTij_new( T, sipas, &isipas, coord, nasipas[cbox]); // "isipas" is properly updated here!

							//					fprintf(stderr,"\nT-matrix: %d vs. %d\n",sipas[isipas]->i,sipas[isipas]->j);
							//					for(int n=0;n<6;n++)
							//					{
							//						for(int m=0;m<6;m++)
							//							fprintf(stderr," %f",T[n][m]);
							//						fprintf(stderr,"\n");
							//					}

							// (ec.23) Noguti & Go 1983 pp.3685-90
							if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
							{
								fprintf(stderr,"middle (using: uij1 uij2 uij3)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
										                        + Uij[uij2y][uij2x][n][m]
										                        - Uij[uij3y][uij3x][n][m]
										                        + T[n][m];
							}
							else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row
							{
								fprintf(stderr,"top-row (using: uij1)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
										                        + T[n][m];
							}
							else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column
							{
								fprintf(stderr,"right-column (using: uij2)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m]
										                        + T[n][m];
							}
							else // a == 0 && b == size-1  --> (0,n)
							{
								fprintf(stderr,"top-left-box (using: none)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = T[n][m];
							}
						}
						else // if "a" and "b" 's units don't interact
						{	 // (no ipas...) --> then, "T" computation is not necessary!
							fprintf(stderr,"Units corresponding to a=%ld and b=%ld Don't interact!\n",a,b);

							if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
							{
								fprintf(stderr,"middle (using: uij1 uij2 uij3)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
										                        + Uij[uij2y][uij2x][n][m]
										                        - Uij[uij3y][uij3x][n][m];
							}
							else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row, but (0,n) (computed above)
							{
								fprintf(stderr,"top-row (using: uij1)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m];
							}
							else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column, but (0,n)
							{
								fprintf(stderr,"right-column (using: uij2)\n");
								for(int n=0;n<6;n++)
									for(int m=0;m<6;m++)
										Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m];
							}
							// a == 0 && b == size-1
							else // a == 0 && b == size-1  --> (0,n)
								fprintf(stderr,"top-left-box (using: none)\n");
							// if there is no interaction, then not Tij addition to Uab needed!
						}
					} // At this point, the Uij element is fully computed (if needed).
					// valgrind: This causes some "invalid reads"
					// fprintf(stderr,"a=%ld, b=%ld, uijy=%d, uijx=%d, Uij[uijy][uijx][0][0]= %f, T[0][0]= %f\n",a,b,uijy,uijx,Uij[uijy][uijx][0][0],T[0][0]);

					// *********************************************************
					// Rab Computation ( the same Single- or Multi- Chain)
					// Table I. Noguti & Go 1983 pp.3685-90
					// (In CA-model, only "backbone vs. backbone" --> Rab = Uab)
					// *********************************************************

//					if(a>=58 && a<=61 && b==943)
//					{
//						fprintf(stderr,"ZERO - This is the full current Uij:\n");
//						for(int ii=0; ii<6; ii++)
//						{
//							for(int jj=0; jj<6; jj++)
//								fprintf(stderr,"%f ",Uij[uijy][uijx][ii][jj]);
//							fprintf(stderr,"\n");
//						}
//					}

					// (ec.16) Noguti & Go 1983 pp.3685-90
					//			fprintf(stderr,"\nR-vector: uijy= %d  uijx= %d\n",uijy,uijx);
					for(int n=0; n<6; n++)
					{
						R[n] = 0.0;
						for(int m=0; m<6; m++)
						{
							R[n] += erx[m][a] * Uij[uijy][uijx][m][n]; // era x Rab
							//					fprintf(stderr," %f*%f",erx[m][a],Uij[uijy][uijx][m][n]);
						}
						//				fprintf(stderr," %f",R[n]);
					}
					//			fprintf(stderr,"\n");

					temp = 0.0;
					for(int n=0; n<6; n++)
						temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

					// if(!justfile) // the Hessian will be dumped into RAM.
					hess_matrix[a + b*(b+1)/2] = temp; // Rab' x erb = Hessian

					fprintf(stderr,"a=%ld b=%ld  uijx=%d uijy=%d  hess= %f\n",a,b,uijy,uijx,hess_matrix[a + b*(b+1)/2]);


//					if(a>=58 && a<=61 && b==943)
//					{
//						fprintf(stderr,"This is the full current Uij:\n");
//						for(int ii=0; ii<6; ii++)
//						{
//							for(int jj=0; jj<6; jj++)
//								fprintf(stderr,"%f ",Uij[uijy][uijx][ii][jj]);
//							fprintf(stderr,"\n");
//						}
//
//						fprintf(stderr,"This is the full current erx[a]:\n");
//						for(int ii=0; ii<6; ii++)
//							fprintf(stderr,"%f ",erx[ii][a]);
//						fprintf(stderr,"\n");
//
//						fprintf(stderr,"This is the full current erx[b]:\n");
//						for(int ii=0; ii<6; ii++)
//							fprintf(stderr,"%f ",erx[ii][b]);
//						fprintf(stderr,"\n");
//					}

					// if(file != NULL) // the Hessian will be dumped into a HD disk file
					//	fwrite(&temp,sizeof(double),1,f_file);
					//					printf("%5d  %21.18e\n",a+1,temp);
				}

				// Filling bottom-box's U-top before deleting current-box Uij
				if(k < nbox-1) // Only there exists U-top for the bottom-box of current box when current box is not in the last row (k=nbox-1).
//					if( col != col_buf && b < size-1) // If we had changed Uij column (i.e. unit) and we were not in the right most column of current box.
					{
						for(int m=0;m<6;m++)
							for(int x=0; x<6; x++)
								Ut[bbox][undh[bi]-undh[b]][m][x] = Uij[0][ur_size-1][m][x]; // last row element from Uij[0]
					}

			} // loop end

			// Filling the left-box's U-right before deleting current-box Uij.
			if(l > 0) // Only there exists U-right for the left-boxes of current box when current box is not in the most-left column (l=0).
			{
				for(int n=0; n < ur_size; n++)
					for(int m=0;m<6;m++)
						for(int x=0; x<6; x++)
							Ur[lbox][n][m][x] = Uij[0][n][m][x];
			}

			// Free "Uij"
			for(int i=0; i<2; i++)
			{
				for(int j=0; j < ur_size; j++)
				{
					for(int m=0;m<6;m++)
						free(Uij[i][j][m]);
					free(Uij[i][j]);
				}
				free(Uij[i]);
			}
			free(Uij);

			fprintf(stderr,"Thread %ld (%d,%d) finished! \n",pthread_self(),k,l);

//			if(k==0 && l==nbox-3)
			if(k==0 && l==nbox-3)
			{
				show_matrix(hess_matrix, size, "Hessian Matrix (F) hess_tr:");
				double *bench_matrix=NULL;
				int size1;
				// bench_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 );
				fprintf(stderr,"Reading hess_ref.bin ...\n");
				read_matrixB(&bench_matrix,&size1,"hess_ref.bin");
				fprintf(stderr,"Matrix hess_ref.bin read!\n");

				double *diff_matrix;
				diff_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 );
				// Compute matrix difference (much easier visualization)
				for(int i=0; i < size*(size+1)/2; i++)
					diff_matrix[i] = bench_matrix[i] - hess_matrix[i]; // computing matrix difference
				show_matrix(diff_matrix, size, "Hessian Difference Matrix (D) diff_matrix:");
				fflush(stdout);
				exit(0);
			}

		}  // END MEAT ---------------------------------------------------------------------------------------------------------------

		// At this point the box should have been already computed and Ur and Ut properly updated...
		pthread_mutex_lock(mtx_isready); // Locking "isready"
		isdone[k + l*(l+1)/2] = true; // setting current-box as done.
		(*jobs_done)++; // Updating the number of finished jobs

		// Checking whether LEFT-box should be set "ready", i.e. if top-left (k-1,l-1) is done.
		if( l > 0 && (k == 0 || (k > 0 && isdone[k-1 + (l-1)*l/2])) )
			if( !isdone[k + (l-1)*l/2] ) // if left-box has not been done yet by any other thread...
				isready[k + (l-1)*l/2] = true; // then set left-box as ready.

		// Checking whether BOTTOM-box should be set "ready", i.e. if bottom-right (k+1,l+1) is done.
		if( k < nbox-1 && (l == nbox-1 || (l < nbox-1 && isdone[k+1 + (l+1)*(l+2)/2])) )
			if( !isdone[k+1 + l*(l+1)/2] ) // if bottom-box has not been done yet by any other thread...
				isready[k+1 + l*(l+1)/2] = true; // then set bottom-box as ready.

		pthread_cond_broadcast(cond_isready); // sending broadcast signal to check new "jobs_done" value...
		pthread_mutex_unlock(mtx_isready); // Un-locking "isready"
	} while(l >= 0); // Do - while the row is not over...

	for(int i=0;i<6;i++)
		free( T[i] );
	free( T );

	fprintf(stderr,"FORCING EXIT, DEBUGGING....\n");
}


//// (05/01/2011)
//// Minimal memory requirements: It outputs Hessian matrix directly into FILE.
//// Fast Hessian Matrix computation // (Multi-Chain & CA-model)
//// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
////        -Single row (N,CA,C)-PDB coordinates,
//// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
//// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
//// Noguti & Go (1983) pp. 3685-90
////void hessianMCAxF(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,char *file,tri *props, int *unat, bool *fix)
//void hessianMCAx_test(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,floating **p_hess_matrix,tri *props, int *unat, bool *fix)
//{
//	bool debug = false;
//	double prod,prod1;
//	int k,l,ind2,ind3,ls,ks;
//	double *dummy;
//	int prin, fin,i, j, j2, m, k1, k2, buff;
//	double r_alpha[3],r_beta[3],r[3],e[3];
//	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
//	double c,d,temp;
//	double v[6];
//	int resn,num_res,num_atoms,num_seg;
//	Residue *res;
//	pdbIter *iter = new pdbIter( mol ); // Iterator to screen atoms
//	pdbIter *iter_frag; // Iterator to screen fragments
//	pdbIter *iter_seg; // Iterator to screen segments
//	num_atoms = mol->get_num_atoms();
//	num_res = mol->get_num_fragments();
//	int index_res = 0;
//	floating *hess_matrix;
//	Htimer ht_timer; // timer
//	TMOL fragtype;
//
//	// File related variables
//	FILE *f_file; // file handle
//	char *file="hessian.bin"; // file name
//
//	// ********************************************
//	// Storing==> (ea, ea x ra) (== (eb, eb x rb) )
//	// ********************************************
//	double **erx;
//	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
//	for(int i=0;i<6;i++)
//		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
//
//	int *undh; // returns the first unit-index on the left side of the dihedral
//	int un_index = 0;
//	undh = (int *) malloc( sizeof(int) * size );
//
//	bool any_rotrans = false;
//
//	// ****************************************************************************
//	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
//	// ****************************************************************************
//	j = 0; // current dihedral index
//	j2 = 0; // mobile dihedral index
//	k1 = 0;
//	int k0 = 0; // first-NH + CA + last-CO model index
//
//	int n_inter; // number of mobile intersegment variables
//	int num_units=0;
//	int n_fixed=0; // counter of fixed internal coordinates
//	int fix_offset=0; // unit-offset (joins units)
//	int natom=0;
//	bool break_unit = false;
//	bool already_broken = false;
//
//	// ---------------------------------------------------
//	// 1. BUILDING "erx" array and other auxiliars: "undh"
//	// ---------------------------------------------------
//
//	Segment * seg;
//	Atom * atom;
//	// Screening segments
//	iter_seg = new pdbIter( mol ); // Iterator to screen segments
//	num_seg = iter_seg->num_segment();
//	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
//	{
//		seg = ( Segment * ) iter_seg->get_segment();
//		iter_frag = new pdbIter( seg );
//		num_res = iter_frag->num_fragment(); // gets number of residues of current segment (to detect ending)
//
//		if(iter_seg->pos_segment != 0) // non-first segment
//		{
//			n_inter=0; // number of mobile intersegment variables
//			for(int n=0;n<6;n++)
//				if(fix == NULL || fix[j+n]) // if mobile
//					n_inter++;
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if( fix == NULL || (fix[j] || fix[j+1] || fix[j+2] || fix[j+3] || fix[j+4] || fix[j+5]) )
//			{
//				any_rotrans = true;
//				// Defining rigid units
//			}
//			else
//				any_rotrans = false;
//
//			// ********************
//			// Adding 3 TRASLATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (0)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					// ea x ra (Psi)  (-gamma_v)
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//					erx[3+axis][j2] = -1.0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 trans added
//
//			// ********************
//			// Adding 3 ROTATIONS
//			// ********************
//			for(int axis=0; axis<3; axis++) // 3 traslations
//			{
//				if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//				{
//					// ea (Phi)  (delta_v)
//					erx[0][j2] = 0;
//					erx[1][j2] = 0;
//					erx[2][j2] = 0;
//					erx[axis][j2] = 1.0;
//					// ea x ra (Psi)  (delta_v  x  u)  it should be 0 as well...
//					erx[3][j2] = 0;
//					erx[4][j2] = 0;
//					erx[5][j2] = 0;
//
//					undh[j2] = un_index;
//
//					j2++;
//				}
//				else
//				{
//					fix_offset++;
//					n_fixed++;
//				}
//
//				j++; // adding traslational degree of freedom
//			} // 3 rots added
//
//			// Check whether any of the 6D inter-segment variables are fixed
//			if(any_rotrans)
//			{
//				un_index++; // count units (due to segment ending)
//			}
//		}
//
//		// Screen ALL fragments
//		for ( iter_frag->pos_fragment = 0; !iter_frag->gend_fragment(); iter_frag->next_fragment() )
//		{
//			res = ( Residue * ) iter_frag->get_fragment();
//			resn = resnum_from_resname( res->getName() );
//			fragtype = res->getMolType();
//
//			// PROTEIN FRAGMENT
//			if( fragtype == tmol_protein )
//			{
//
//				if(iter_frag->pos_fragment == 0) // has NH
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//				// NOT FIRST, NOT LAST, NOT PRO -> PHI
//				// ("PHI-bond")
//				if ( !(iter_frag->pos_fragment == 0 || iter_frag->pos_fragment == num_res-1 || resn == PRO) )
//				{
//					if( fix == NULL || fix[j] ) // Check whether current variable is mobile
//					{
//						//			  printf ("%d Not first not PRO -> PHI \n",j);
//						// ("k2" is updated below (for the 1st not PRO -> PHI)
//
//						// y_lambda (NH pos 0)
//						e[0] = -coord[k1 * 3];
//						e[1] = -coord[k1 * 3 + 1];
//						e[2] = -coord[k1 * 3 + 2];
//						// CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// e_lambda
//						e[0] += r[0]; // NH --> CA
//						e[1] += r[1];
//						e[2] += r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] ); // normalization factor
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] index
//				}  // NOT FIRST NOT PRO
//
//				if(break_unit)
//				{
//					un_index++;
//					break_unit = false; // by default not breaking unit
//				}
//
//				unat[k0] = un_index;
//				k0++;
//
//				already_broken = false;
//
//				// NOT LAST, NOT FIRST RESIDUE--> PSI
//				// ("PSI-bond")
//				if( iter_frag->pos_fragment != 0 && iter_frag->pos_fragment != num_res - 1 )
//				{
//					//			 printf("%d Not last residue -> PSI\n",j);
//					if( fix == NULL || fix[j] ) // Check whether current variable it's mobile
//					{
//						// get CA pos 1
//						r[0] = coord[(k1+1) * 3];
//						r[1] = coord[(k1+1) * 3 + 1];
//						r[2] = coord[(k1+1) * 3 + 2];
//						// get C pos 2 (C=O in 3BB2R) or (C in Full-Atom)
//						e[0] = coord[(k1+2) * 3];
//						e[1] = coord[(k1+2) * 3 + 1];
//						e[2] = coord[(k1+2) * 3 + 2];
//						// e_lambda ==> CA --> C (unit vector)
//						e[0] -= r[0];
//						e[1] -= r[1];
//						e[2] -= r[2];
//
//						temp = sqrt( e[0] * e[0] + e[1] * e[1] + e[2] * e[2] );
//						for ( m = 0; m < 3; m++ )
//							e[m] /= temp; // Unit vector normalization
//
//						// ea
//						erx[0][j2] = e[0];
//						erx[1][j2] = e[1];
//						erx[2][j2] = e[2];
//						// ea x ra
//						erx[3][j2] = e[1] * r[2] - e[2] * r[1];
//						erx[4][j2] = e[2] * r[0] - e[0] * r[2];
//						erx[5][j2] = e[0] * r[1] - e[1] * r[0];
//
//						undh[j2] = un_index; // which unit has the "j" dihedral seen from the right (left side)
//
//						if(!already_broken)
//						{
//							break_unit = true;
//							already_broken = true;
//						}
//
//						j2++;
//					}
//					else
//						n_fixed++;
//
//					j++; // der[] dihedral index
//				}  // NOT LAST
//
//				if(iter_frag->pos_fragment == num_res-1) // has CO
//				{
//					unat[k0] = un_index; // bb unit index
//					k0++;
//				}
//
//			} // END PROTEIN FRAGMENT
//			else // NOT FOUND MOL-TYPE
//			{
//				printf("Msg(hessianMCAx): Sorry, Molecular-type (TMOL) not identified!\nForcing exit!\n");
//				exit(2);
//			}
//
//			k1+=3; // residue index // first residue atom index
//			index_res++; // counts residues (due to "props")
//		}
//		delete(iter_frag);
//	}
//	delete(iter_seg);
//
//	// ---------------------------------------------------------------------------
//	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
//	// ---------------------------------------------------------------------------
//	// Tij related pre-computations (to get direct Tij computation)
//
//	num_units = un_index + 1; // Inter-unit variables (rot-trans) dont increase units number
//
//	if(debug)
//		printf("Msg(hessianMCAx): Number of units: %d (%d)\n",num_units,un_index);
//
//	int unipa_size=0;
//	int unipa_index;
//	int *p_int;
//	int **unipa;
//
//	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
//	// (triangular packing storage)
//	unipa = (int **) malloc( sizeof(int *) * num_units*(num_units+1)/2 );
//	unipa_size += sizeof(int *) * num_units*(num_units+1)/2;
//
//	// "unipa" initialization
//	for(int i=0; i<num_units*(num_units+1)/2; i++)
//		unipa[i] = NULL;
//
//	if(debug)
//		printf("Msg(hessianMCAx): Triangular unipa_size= %d (nipa=%d)\n",unipa_size,nipa);
//
//	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
//	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
//	{
//		k = decint[index].k; // k
//		l = decint[index].l; // l
//		i = unat[k]; // "k" atom unit index
//		j = unat[l]; // "l" atom unit index
//		if(i != j) // the same unit is rigid, so inner contacts must not be taken into account!
//		{		   // (as well as the non contacting pairs!)
//			// The following must be allways true:  (k < l) && (i < j)
//			if( i > j )
//			{
//				buff = j;
//				j = i;
//				i = buff;
//			}
//			// translates from "squared" to "triangular" matrix elements
//			unipa_index = i + j*(j+1)/2;
//
//			// Allocating memory for the new element (memmory efficient)
//			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
//			{
//				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa[unipa_index][0] = 2;
//				unipa_size +=  2 * sizeof(int);
//			}
//			else // If already allocated, then increases it size one int.
//			{
//				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
//				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
//				{
//					printf("Msg(hessianMCAx): Memmory allocation failed in realloc()!\nForcing exit!\n");
//					exit(1);
//				}
//				unipa[unipa_index] = p_int;
//				unipa_size +=  sizeof(int);
//			}
//
//			// Storing ipa index for k,l interacting pair of atoms
//			unipa[unipa_index][ unipa[unipa_index][0]-1 ] = index;
//		}
//	}
//
//	if(debug)
//		printf("Msg(hessianMCAx): Final UNIPA size = %d bytes (%.3f Mb)\n",unipa_size,(float)unipa_size/1e6);
//
//
//	// ---------------------------------------------------------------------------
//	// 2'. Sorting "ipas" instead of building "unipa" (mandatory for huge systems!)
//	// ---------------------------------------------------------------------------
//
//	// Sorting "ipas" to facilitate Hessian computation (mandatory for huge systems!)
//	sort_ipas(decint,&nipa,unat,num_units);
//
//	// -----------------------------------------------------------------------------------
//	// 3. HESSIAN BUILDING from Uab elements:
//	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
//	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//	// -----------------------------------------------------------------------------------
//
//	// Hessian memory allocation
//	if(*p_hess_matrix == NULL)
//	{
//		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
//		{
//			printf("Msg(hessianMCAx): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
//			exit(1);
//		}
//		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
//	}
//	else
//		hess_matrix = *p_hess_matrix;
//
//	// Hessian initialization
//	for(i=0; i<size*(size+1)/2; i++)
//		hess_matrix[i]=0.0;
//
//	// Uabmn minimal-memory allocation (full-square + 2 cols. + Tij "on the fly" computation)
//	int Usize=0;
//
//	// Now indexing is a bit different
//	int uij_index,uij1x,uij1y,uij2x,uij2y,uij3x,uij3y,uijx,uijy;
//	int col, col_buf;
//
//	// Now we will allocate just a (6 cols. x num_units rows.) matrix
//	// with 6x6 elements on each position.
//	double ****Uij;
//	double ***Uitemp;
//	if( Uij = (double ****) malloc( sizeof(double ***) * 2 ) ) // Just 2 columns!
//		for(int i=0; i<2; i++)
//		{
//			if( Uij[i] = (double ***) malloc( sizeof(double **) * num_units) ) // rows
//				for(int j=0; j<num_units; j++)
//					if( Uij[i][j] = (double **) malloc( sizeof(double *) * 6 ) )
//					{
//						for(int m=0;m<6;m++)
//						{
//							if( Uij[i][j][m] = (double *) malloc( sizeof(double) * 6 ) )
//								for(int n=0; n<6; n++)
//									Uij[i][j][m][n] = 0.0; // initialization
//							else
//							{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
//						}
//						Usize += sizeof(double *)*6 + sizeof(double)*36;
//					}
//					else
//					{ printf("Msg(hessianM): Unable to allocate an Uij-matrix element!\n"); exit(1); }
//			else
//			{ printf("Msg(hessianM): Unable to allocate an Uij-matrix column!\n"); exit(1); }
//		}
//	else
//	{
//		printf("Msg(hessianMCAx): Unable to allocate Uij-matrix (6 x size) !\n");
//		exit(1);
//	}
//	Usize += ( sizeof(double **) * num_units * 2 ) + ( sizeof(double ***) * 2 );
//	if(debug)
//		printf("Msg(hessianMCAx): Final Uij size = %d bytes (%.3f Mb)\n",Usize,(float)Usize/1e6);
//
//	// T memory allocation (6x6 matrix)
//	double **T;
//	T = (double **) malloc( sizeof(double *) * 6 );
//	for(int i=0; i<6; i++)
//		T[i] = (double *) malloc( sizeof(double) * 6 );
//
//	col_buf = undh[size-1]+1;
//	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
//	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
//	for(int b=size-1; b >= 0; b--) // b --> dihedrals (column)
//	{
//		col_buf = col; // Uij column buffer
//		col = undh[b]+1; // Current Uij column
//
//		if( col != col_buf && b < size-1) // If we changed Uij column (i.e. unit)
//		{
//			// Displacing Uij one position to the right...
//			// First, last element initialization
//			for(int j=0; j<num_units; j++)
//				for(int m=0;m<6;m++)
//					for(int n=0; n<6; n++)
//						Uij[1][j][m][n] = 0.0; // initialization
//			// Second, displacing 1 position to the right...
//			Uitemp = Uij[1]; // buffer
//			for(int j=0; j>=0; j--)
//				Uij[j+1] = Uij[j]; // current element displament
//			Uij[0] = Uitemp;
//		}
//
//		for(int a=0; a <= b; a++)	 // a --> dihedrals (row)
//		{ 							 // it fills upper triangular part (including diagonal)
//			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
//			//
//			//   a0   a1   a2
//			// O----O----O----O
//			// u0   u1   u2   u3
//			//
//			// U(a,b)
//			uij_index = undh[a] + (undh[b]+1)*(undh[b]+2)/2; // needed by "unipa"
//			uijx = undh[a];
//			uijy = 0;
//
//			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
//			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
//			// the Uij for the 4 ICs combination should be computed only once!
//			if( Uij[uijy][uijx][0][0] == 0.0 ) // If "Uij" was not computed yet (compute Uij only once!)
//			{
//				// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
//				// Computing valid indices
//
//				// U(a,b+1)
//				if(undh[b]+1 < num_units-1 ) // if valid "uab1"
//				{
//					uij1x = undh[a]; // i= undh[a] // (row)
//					uij1y = 1 ; // j= (undh[b]+1)+1 // (col.)
//				}
//
//				// U(a-1,b)
//				if( undh[a] > 0 )
//				{
//					uij2x = undh[a]-1;
//					uij2y = 0;
//				}
//
//				// U(a-1,b+1)
//				if( undh[a] > 0 && undh[b]+2 < num_units )
//				{
//					uij3x = undh[a]-1;
//					uij3y = 1;
//				}
//
//				// Computing Uij element
//				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
//				{								 // (whether they have ipas...) --> then, "T" must be computed
//					// Computing Tij element "on the fly"
//					// (ec.21) Noguti & Go 1983 pp.3685-90
//					calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );
//
//					// (ec.23) Noguti & Go 1983 pp.3685-90
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
//													   + Uij[uij2y][uij2x][n][m]
//													   - Uij[uij3y][uij3x][n][m]
//													   + T[n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
//													   + T[n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m]
//													   + T[n][m];
//					}
//					else // a == 0 && b == size-1  (0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = T[n][m];
//					}
//				}
//				else // if "a" and "b" 's units don't interact
//				{	 // (no ipas...) --> then, "T" computation is not necessary!
//					if( undh[a] > 0 && undh[b]+1 < num_units-1 ) // In the middle... (most of them)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m]
//													   + Uij[uij2y][uij2x][n][m]
//													   - Uij[uij3y][uij3x][n][m];
//					}
//					else if( undh[a] == 0 && undh[b]+1 < num_units-1 ) // In the top row (execepting: 0,n  already computed above)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij1y][uij1x][n][m];
//					}
//					else if( undh[a] > 0 && undh[b]+1 == num_units-1 ) // In the right column (execepting: 0,n)
//					{
//						for(int n=0;n<6;n++)
//							for(int m=0;m<6;m++)
//								Uij[uijy][uijx][n][m] = Uij[uij2y][uij2x][n][m];
//					}
//					// if no interaction, then no Tij addition to Uab needed!
//				}
//				// at this point, the Uij element is fully computed!
//			}   // only if not computed yet!
//
//			// ****************************************************
//			// Rab - Matrix Computation
//			// Noguti & Go (1983) pp. 3685-90
//			// ****************************************************
//			// Rab Computation ( the same Single- or Multi- Chain)
//			// Table I. Noguti & Go 1983 pp.3685-90
//			// backbone vs. backbone --> Rab = Uab
//
//			// (ec.16) Noguti & Go 1983 pp.3685-90
//			for(int n=0; n<6; n++)
//			{
//				R[n] = 0.0;
//				for(int m=0; m<6; m++)
//					R[n] += erx[m][a] * Uij[uijy][uijx][m][n]; // era x Rab
//			}
//			temp = 0.0;
//			for(int n=0; n<6; n++)
//				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian
//			hess_matrix[a + b*(b+1)/2] = temp; // Rab' x erb = Hessian
//		}
//	} // loop ended
//
//	// Free "Uij"
//	for(int i=0; i<1; i++)
//	{
//		for(int j=0; j<num_units; j++)
//		{
//			for(int m=0;m<6;m++)
//				free(Uij[i][j][m]);
//			free(Uij[i][j]);
//			Usize -= sizeof(double *)*6 + sizeof(double)*36;
//		}
//		free(Uij[i]);
//		Usize -= sizeof(double **)*num_units;
//	}
//	free(Uij);
//	Usize -= sizeof(double ***)*2;
//	if(debug)
//		printf("Msg(hessianMCAx): Final freed Uij size = %d bytes (should be =0)\n",Usize);
//
//	// Free "unipa"
//	for(int i=0; i<num_units*(num_units+1)/2; i++)
//		free( unipa[i] );
//	free( unipa );
//
//	for(int i=0;i<6;i++)
//	{
//		free( erx[i] ); // erx[6][size]
//		free( T[i] );
//	}
//	free( erx );
//	free( T );
//	free( undh );
//}

// Hessian matrix computation by "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
void hess_naive(twid *decint,trd *der,int nipa,int size,float *coord,double **p_hess_matrix)
{
	double prod,prod1;
	float r[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	double *hess_matrix;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
	{
		prod1 = decint[index].C / pow(decint[index].d,2);
		k=decint[index].k; // k
		l=decint[index].l; // l
		ks=k*size;
		ls=l*size;
		r[0] = coord[l * 3]-coord[k * 3];
		r[1] = coord[l * 3 + 1]-coord[k * 3 + 1];
		r[2] = coord[l * 3 + 2]-coord[k * 3 + 2];

		for(int i=0;i<size;i++) // screening dihedrals
		{
			// Does this screen diagonally the K-matrix ??
			ind2=ls+i; // i-dihedral and l-pseudo-atom derivative index
			ind3=ks+i; // i-dihedral and k-pseudo-atom derivative index
			// eigval is used as dummy var
			dummy[i] = r[0]*(der[ind2].x-der[ind3].x) +
			r[1]*(der[ind2].y-der[ind3].y) +
			r[2]*(der[ind2].z-der[ind3].z);

		}

		// fill upper triangular part (including diagonal)
		for(int i=0 ;i<size; i++)
		{
			prod = prod1 * dummy[i]; // eigval is used as dummy var
			for(int j=i;j<size;j++)
				hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
		}
	}
	free(dummy);
}


// Hessian matrix computation by "naive" method O(n^4)
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
void hess_naiveVW_old(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix)
{
	double prod,prod1;
	float r[3],rk[3],rl[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	double *hess_matrix;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
	{
		prod1 = decint[index].C / pow(decint[index].d,2);
		k=decint[index].k; // k
		l=decint[index].l; // l
		rk[0] = coord[k*3];
		rk[1] = coord[k*3+1];
		rk[2] = coord[k*3+2];
		rl[0] = coord[l*3];
		rl[1] = coord[l*3+1];
		rl[2] = coord[l*3+2];
		r[0] = rl[0]-rk[0];
		r[1] = rl[1]-rk[1];
		r[2] = rl[2]-rk[2];

		for(int i=0;i<size;i++) // screening dihedrals
		{
			if( body1[i][k] && !body1[i][l] ) // if first 1st body and second 2nd body
				dummy[i] = derVW(V[i],W[i],rk,rl);
			else if( !body1[i][k] && body1[i][l] ) // first 2nd body and second 1st body (due to CHI angles)
				dummy[i] = derVW(V[i],W[i],rl,rk);
			else
				dummy[i] = 0.0; // both in the same body
		}

		// fill upper triangular part (including diagonal)
		for(int i=0 ;i<size; i++)
		{
			prod = prod1 * dummy[i]; // eigval is used as dummy var
			for(int j=i;j<size;j++)
				hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
		}
	}
	free(dummy);
}

// Hessian matrix computation by "naive" method O(n^4)
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
void hess_naiveVW(float *coord,twid *decint,double ***V,double ***W,int **body1,int nipa,int size,double **p_hess_matrix,int model)
{
	double prod,prod1;
	float r[3],rk[3],rl[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	double *hess_matrix;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	if( model == 1 || model == 2 ) // C5/HA models
		for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
		{
			prod1 = decint[index].C / pow(decint[index].d,2);
			k=decint[index].k; // k
			l=decint[index].l; // l
			rk[0] = coord[k*3];
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			rl[0] = coord[l*3];
			rl[1] = coord[l*3+1];
			rl[2] = coord[l*3+2];
			r[0] = rl[0]-rk[0];
			r[1] = rl[1]-rk[1];
			r[2] = rl[2]-rk[2];

			for(int i=0;i<size;i++) // screening dihedrals
			{
				// If "k" belongs to body1 and "l" to body2
				if(
						((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) && // body1
						((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) )    // body2
				)
					dummy[i] = derVW(V[i],W[i],rk,rl);
				else if(
						((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) ) && // body1
						((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) )    // body2
				)
					dummy[i] = derVW(V[i],W[i],rl,rk);
				else
					dummy[i] = 0.0; // both in the same body (this saves a lot CPU)

			}

			// fill upper triangular part (including diagonal)
			for(int i=0 ;i<size; i++)
			{
				prod = prod1 * dummy[i]; // eigval is used as dummy var
				for(int j=i;j<size;j++)
					hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
			}
		}
	else // CA-model (or NCAC-model)
		for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
		{
			prod1 = decint[index].C / pow(decint[index].d,2);
			k=decint[index].k; // k
			l=decint[index].l; // l
			rk[0] = coord[k*3];
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			rl[0] = coord[l*3];
			rl[1] = coord[l*3+1];
			rl[2] = coord[l*3+2];
			r[0] = rl[0]-rk[0];
			r[1] = rl[1]-rk[1];
			r[2] = rl[2]-rk[2];

			for(int i=0;i<size;i++) // screening dihedrals
			{
				// If "k" belongs to body1 and "l" to body2
				if( k < body1[i][0] && l >= body1[i][0] )
					dummy[i] = derVW(V[i],W[i],rk,rl);
				else if( l < body1[i][0] && k >= body1[i][0] ) // now "l" and "k" (reversed order)
					dummy[i] = derVW(V[i],W[i],rl,rk);
				else
					dummy[i] = 0.0; // both in the same body (this saves a lot CPU)
			}

			// fill upper triangular part (including diagonal)
			for(int i=0 ;i<size; i++)
			{
				prod = prod1 * dummy[i]; // eigval is used as dummy var
				for(int j=i;j<size;j++)
					hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
			}
		}

	free(dummy);
}

// Hessian matrix computation by "naive" method O(n^4) (Valid for all "type")
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Avoids zero derivatives computation (Mon, 19/02/2009) (zero-tolerance not needed with V/W!)
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
//void hess_naiveMFAVW(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix)
void hess_naiveVW0_old(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix)
{
	double prod,prod1;
	float r[3],rk[3],rl[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	double *hess_matrix;

	// (Mon, 19/02/2009)
	int lim_l,lim_r;
	bool liml_set, limr_set;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
	{
		prod1 = decint[index].C / pow(decint[index].d,2);
		k=decint[index].k; // k
		l=decint[index].l; // l
		rk[0] = coord[k*3];
		rk[1] = coord[k*3+1];
		rk[2] = coord[k*3+2];
		rl[0] = coord[l*3];
		rl[1] = coord[l*3+1];
		rl[2] = coord[l*3+2];
		r[0] = rl[0]-rk[0];
		r[1] = rl[1]-rk[1];
		r[2] = rl[2]-rk[2];

		liml_set = false;
		limr_set = false;

		lim_l = 0;
		lim_r = -1;

		for(int i=0;i<size;i++) // screening dihedrals
		{
			if( body1[i][k] && !body1[i][l] ) // k-body1 & l-body2
				dummy[i] = derVW(V[i],W[i],rk,rl);
			else if( body1[i][l] && !body1[i][k] ) // l-body1 & k-body2 (order reversal)
				dummy[i] = derVW(V[i],W[i],rl,rk);
			else
				dummy[i] = 0.0; // both in the same body

			if(lim_r == -1)
				// if null dummy
				if( dummy[i] == 0.0 )
					lim_l=i; //
				else
					lim_r=i;
			else
				// if not-null dummy
				if( dummy[i] != 0.0 ) // Mon added (1/12/2009)
					lim_r=i;
		}

		// fill upper triangular part (including diagonal)
		for(int i=lim_l ;i<=lim_r; i++)
		{
			prod = prod1 * dummy[i]; // eigval is used as dummy var
			for(int j=i;j<=lim_r;j++)
				hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
		}
	}
	free(dummy);
}

// Hessian matrix computation by "naive" method O(n^4) (Valid for all "type")
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Avoids zero derivatives computation, 10-20x faster (Mon, 19/02/2009) (zero-tolerance not needed with V/W!)
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
//void hess_naiveMFAVW(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix)
void hess_naiveVW0(float *coord,twid *decint,double ***V,double ***W,int **body1,int nipa,int size,double **p_hess_matrix,int model,int row)
{
	double prod,prod1;
	float r[3],rk[3],rl[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	double *hess_matrix;

	// (Mon, 19/02/2009)
	int lim_l,lim_r;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (double *) malloc( sizeof(double) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		// Hessian initialization
		for(int i=0; i<size*(size+1)/2; i++)
			hess_matrix[i]=0.0;

		// Outputs Hessian Matrix
		*p_hess_matrix = hess_matrix;
	}
	else
		hess_matrix = *p_hess_matrix;

	if(row < 0) // Full matrix computation (standard mode)
	{
		// Dummy memory allocation & initialization
		if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
		{
			printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
			exit(1);
		}
		for(int i=0; i<size; i++)
			dummy[i]=0.0; // initialization
	}

	if( model == 1 || model == 2 ) // C5/HA models
		for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
		{
			prod1 = decint[index].C / pow(decint[index].d,2);
			k=decint[index].k; // k
			l=decint[index].l; // l
			rk[0] = coord[k*3];
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			rl[0] = coord[l*3];
			rl[1] = coord[l*3+1];
			rl[2] = coord[l*3+2];
			r[0] = rl[0]-rk[0];
			r[1] = rl[1]-rk[1];
			r[2] = rl[2]-rk[2];

			if(row < 0) // Full matrix computation (standard mode)
			{
				lim_l = 0;
				lim_r = -1;
				for(int i=0;i<size;i++) // screening dihedrals
				{
					if(
							((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) && // body1
							((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) )    // body2
					) // If "k" belongs to body1 and "l" to body2
						dummy[i] = derVW(V[i],W[i],rk,rl);
					else if(
							((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) ) && // body1
							((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) )    // body2
					) // now "l" and "k" (reversed order)
						dummy[i] = derVW(V[i],W[i],rl,rk);
					else
						dummy[i] = 0.0; // both in the same body (this saves a lot CPU)

					if(lim_r == -1)
						// if null dummy
						if( dummy[i] == 0.0 )
							lim_l=i; //
						else
							lim_r=i;
					else
						// if not-null dummy
						if( dummy[i] != 0.0 ) // Mon added (1/12/2009)
							lim_r=i;
				}

				// fill upper triangular part (including diagonal)
				for(int i=lim_l ;i<=lim_r; i++)
				{
					prod = prod1 * dummy[i]; // eigval is used as dummy var
					for(int j=i;j<=lim_r;j++)
						hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
				}
			}
			else // Single row computation mode (for further paralellization)
			{
				if(		((k < body1[row][0] || (body1[row][1]>=0 && k > body1[row][1])) && (k != body1[row][2]) ) && // body1
						((l < body1[row][0] || (body1[row][1]>=0 && l > body1[row][1])) && (l != body1[row][2]) )    // body2
				) 		// If "k" belongs to body1 and "l" to body2
				{
					// "prod" computation ("prod" is needed first!)
					prod = prod1 * derVW(V[row],W[row],rk,rl);
					for(int i=row;i<size;i++) // screening dihedrals
						if(		((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) && // body1
								((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) )    // body2
						) 		// If "k" belongs to body1 and "l" to body2
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rk,rl); // fill just selected row
						else if(
								((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) ) && // body1
								((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) )    // body2
						) 		// now "l" and "k" (reversed order)
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rl,rk); // fill just selected row
						// both in the same body (this saves a lot CPU) --> Then do nothing!
				}
				else if(
						((l < body1[row][0] || (body1[row][1]>=0 && l > body1[row][1])) && (l != body1[row][2]) ) && // body1
						((k < body1[row][0] || (body1[row][1]>=0 && k > body1[row][1])) && (k != body1[row][2]) )    // body2
				)		// now "l" and "k" (reversed order)
				{
					// "prod" computation ("prod" is needed first!)
					prod = prod1 * derVW(V[row],W[row],rk,rl);
					for(int i=row;i<size;i++) // screening dihedrals
						if(		((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) ) && // body1
								((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) )    // body2
						) 		// If "k" belongs to body1 and "l" to body2
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rk,rl); // fill just selected row
						else if(
								((l < body1[i][0] || (body1[i][1]>=0 && l > body1[i][1])) && (l != body1[i][2]) ) && // body1
								((k < body1[i][0] || (body1[i][1]>=0 && k > body1[i][1])) && (k != body1[i][2]) )    // body2
						) 		// now "l" and "k" (reversed order)
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rl,rk); // fill just selected row
					// both in the same body (this saves a lot CPU) --> Then do nothing!
				}
				// both in the same body (this saves a lot CPU) --> Then do nothing!
			}
		}
	else // CA-model
		for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
		{
			prod1 = decint[index].C / pow(decint[index].d,2);
			k=decint[index].k; // k
			l=decint[index].l; // l
			rk[0] = coord[k*3];
			rk[1] = coord[k*3+1];
			rk[2] = coord[k*3+2];
			rl[0] = coord[l*3];
			rl[1] = coord[l*3+1];
			rl[2] = coord[l*3+2];
			r[0] = rl[0]-rk[0];
			r[1] = rl[1]-rk[1];
			r[2] = rl[2]-rk[2];

			if(row < 0) // Full matrix computation (standard mode)
			{
				lim_l = 0;
				lim_r = -1;
				for(int i=0;i<size;i++) // screening dihedrals
				{
					if( k < body1[i][0] && l >= body1[i][0] ) // If "k" belongs to body1 and "l" to body2
						dummy[i] = derVW(V[i],W[i],rk,rl);
					else if( l < body1[i][0] && k >= body1[i][0] ) // now "l" and "k" (reversed order)
						dummy[i] = derVW(V[i],W[i],rl,rk);
					else
						dummy[i] = 0.0; // both in the same body (this saves a lot CPU)

					if(lim_r == -1)
						// if null dummy
						if( dummy[i] == 0.0 )
							lim_l=i; //
						else
							lim_r=i;
					else
						// if not-null dummy
						if( dummy[i] != 0.0 ) // Mon added (1/12/2009)
							lim_r=i;
				}

				for(int i=lim_l ;i<=lim_r; i++)
				{
					prod = prod1 * dummy[i]; // eigval is used as dummy var
					for(int j=i;j<=lim_r;j++) // fill upper triangular part (including diagonal)
						hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
				}
			}
			else // Single row computation mode (for further paralellization)
			{
				if( k < body1[row][0] && l >= body1[row][0] ) // If "k" belongs to body1 and "l" to body2
				{
					// "prod" computation ("prod" is needed first!)
					prod = prod1 * derVW(V[row],W[row],rk,rl); // alpha
					for(int i=row;i<size;i++) // screening dihedrals (beta)
						if( k < body1[i][0] && l >= body1[i][0] ) // If "k" belongs to body1 and "l" to body2
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rk,rl); // fill just selected row
						else if( l < body1[i][0] && k >= body1[i][0] ) // now "l" and "k" (reversed order)
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rl,rk); // fill just selected row
						// both in the same body (this saves a lot CPU) --> Then do nothing!
				}
				else if( l < body1[row][0] && k >= body1[row][0] ) // now "l" and "k" (reversed order)
				{
					// "prod" computation ("prod" is needed first!)
					prod = prod1 * derVW(V[row],W[row],rl,rk);
					for(int i=row;i<size;i++) // screening dihedrals
						if( k < body1[i][0] && l >= body1[i][0] ) // If "k" belongs to body1 and "l" to body2
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rk,rl); // fill just selected row
						else if( l < body1[i][0] && k >= body1[i][0] ) // now "l" and "k" (reversed order)
							hess_matrix[ row + i*(i+1)/2 ] += prod * derVW(V[i],W[i],rl,rk); // fill just selected row
						// both in the same body (this saves a lot CPU) --> Then do nothing!
				}
				// both in the same body (this saves a lot CPU) --> Then do nothing!
			}
		}

	if(row < 0)
		free(dummy);
}


// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// Avoids zero derivatives computation (Mon, 1/12/2009) (warning zero-tolerance=0.0000000001)
void hess_naive0(twid *decint,trd *der,int nipa,int size,float *coord,floating **p_hess_matrix)
{
	double prod,prod1;
	float r[3];
	int k,l,ind2,ind3,ls,ks;
	double *dummy;
	floating *hess_matrix;

	// (Mon added, 1/12/2009)
	double tol=0.0000000001;
	int lim_l,lim_r;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (double *) malloc(size * sizeof(double)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
	{
		prod1 = decint[index].C / pow(decint[index].d,2);
		k=decint[index].k; // k
		l=decint[index].l; // l
		ks=k*size;
		ls=l*size;
		r[0] = coord[l * 3]-coord[k * 3];
		r[1] = coord[l * 3 + 1]-coord[k * 3 + 1];
		r[2] = coord[l * 3 + 2]-coord[k * 3 + 2];

		lim_l = 0;
		lim_r = -1;
		for(int i=0;i<size;i++) // screening dihedrals
		{
			// Does this screen diagonally the K-matrix ??
			ind2=ls+i; // i-dihedral and l-pseudo-atom derivative index
			ind3=ks+i; // i-dihedral and k-pseudo-atom derivative index
			// eigval is used as dummy var
			dummy[i] = r[0]*(der[ind2].x-der[ind3].x) +
			r[1]*(der[ind2].y-der[ind3].y) +
			r[2]*(der[ind2].z-der[ind3].z);

			if(lim_r == -1)
				// if null dummy
				if( dummy[i] < tol && dummy[i] > -tol ) // Mon added (1/12/2009)
					lim_l=i; //
				else
					lim_r=i;
			else
				// if not-null dummy
				if( !(dummy[i] < tol && dummy[i] > -tol) ) // Mon added (1/12/2009)
					lim_r=i;
		}

		// fill upper triangular part (including diagonal)
		for(int i=lim_l ;i<=lim_r; i++)
		{
			prod = prod1 * dummy[i]; // eigval is used as dummy var
			for(int j=i;j<=lim_r;j++)
				hess_matrix[ i + j*(j+1)/2 ] += (floating) prod * dummy[j];
		}
	}
	free(dummy);
}

// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// Avoids zero derivatives computation (Mon, 1/12/2009) (warning zero-tolerance=0.0000000001)
void hess_naive0_double(twid *decint,trd *der,int nipa,int size,float *coord,floating **p_hess_matrix)
{
	floating prod,prod1;
	float r[3];
	int k,l,ind2,ind3,ls,ks;
	floating *dummy;
	floating *hess_matrix;

	// (Mon added, 1/12/2009)
	floating tol=0.0000000001;
	int lim_l,lim_r;

	if(*p_hess_matrix == NULL)
	{
		if( !(hess_matrix = (floating *) malloc( sizeof(floating) * size*(size+1)/2 )) ) // Triangular
		{
			printf("Msg(hess_naive): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
			exit(1);
		}
		*p_hess_matrix = hess_matrix; // outputs Hessian Matrix
	}
	else
		hess_matrix = *p_hess_matrix;

	// Hessian initialization
	for(int i=0; i<size*(size+1)/2; i++)
		hess_matrix[i]=0.0;

	// Dummy memory allocation & initialization
	if( !( dummy  = (floating *) malloc(size * sizeof(floating)) ) )
	{
		printf("Sorry, unable to allocate Dummy memory!!!\nForcing exit!!!\n\n");
		exit(1);
	}
	for(int i=0; i<size; i++)
		dummy[i]=0.0; // initialization

	for(int index=0;index<nipa;index++) // screens contacts (k vs. l)
	{
		prod1 = decint[index].C / pow(decint[index].d,2);
		k=decint[index].k; // k
		l=decint[index].l; // l
		ks=k*size;
		ls=l*size;
		r[0] = coord[l * 3]-coord[k * 3];
		r[1] = coord[l * 3 + 1]-coord[k * 3 + 1];
		r[2] = coord[l * 3 + 2]-coord[k * 3 + 2];

		lim_l = 0;
		lim_r = -1;
		for(int i=0;i<size;i++) // screening dihedrals
		{
			// Does this screen diagonally the K-matrix ??
			ind2=ls+i; // i-dihedral and l-pseudo-atom derivative index
			ind3=ks+i; // i-dihedral and k-pseudo-atom derivative index
			// eigval is used as dummy var
			dummy[i] = r[0]*(der[ind2].x-der[ind3].x) +
			r[1]*(der[ind2].y-der[ind3].y) +
			r[2]*(der[ind2].z-der[ind3].z);

			if(lim_r == -1)
				// if null dummy
				if( dummy[i] < tol && dummy[i] > -tol ) // Mon added (1/12/2009)
					lim_l=i; //
				else
					lim_r=i;
			else
				// if not-null dummy
				if( !(dummy[i] < tol && dummy[i] > -tol) ) // Mon added (1/12/2009)
					lim_r=i;
		}

		// fill upper triangular part (including diagonal)
		for(int i=lim_l ;i<=lim_r; i++)
		{
			prod = prod1 * dummy[i]; // eigval is used as dummy var
			for(int j=i;j<=lim_r;j++)
				hess_matrix[ i + j*(j+1)/2 ] += prod * dummy[j];
		}
	}
	free(dummy);
}
