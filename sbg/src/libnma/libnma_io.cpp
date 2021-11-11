/************************************************************************
 *                     LIBRARY: libnma_io                                *
 *************************************************************************
 * Program is part of the ADP package URL: http://sbg.cib.csic.es        *
 * (c) Jose Ramon Lopez-Blanco (Mon), Jose Ignacio Garzon and            *
 * Pablo Chacon. CSIC-CIB's Structural Bioinformatics Group (2004-2010)  *
 *************************************************************************
 *                                                                       *
 *   Input-Output library.                                               *
 *                                                                       *
 *************************************************************************
 * This library is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation; either version 2 of the License, or     *
 * (at your option) any later version.                                   *
 ************************************************************************/

#include <libnma_io.h>
#include <sstream>
#include <iostream>

// Memory allocator
void *allocate(size_t size, char *message)
{
	void *pointer;
	if( !(pointer = (void *) malloc(size) ) ) // Allocating memory with malloc
	{
		printf("Memory allocation error...\n%s\nForcing exit!\n",message);
		exit(1);
	}
	return pointer;
}

// Memory re-allocator
void *reallocate(void *pointer, size_t size, char *message)
{
	if( !(pointer = (void *) realloc(pointer,size) ) ) // Allocating memory with malloc
	{
		printf("Memory allocation error...\n%s\nForcing exit!\n",message);
		exit(1);
	}
	return pointer;
}

// Saves & Normalizes eigenpairs according to "ptraj" format (normalizes by default)
void save_ptraj_modes (char * name, int size, int start_mode, int end_mode, floating *eigval, floating *eigvect,bool norm_switch)
{
	floating current_ampl, vectorrmsd,maximumamplitude,norm;
	FILE *out;
	int i,j,nmodes;
	nmodes=end_mode-start_mode;
	if( !(out = fopen(name, "w")) )
	{
		printf("Msg(save_ptraj_modes): Unable open to write %s\nForcing exit!\n\n",name);
		exit(1);
	}

	// Header line
	fprintf(out," Eigenvector file: COVAR\n %d %d Contains %d eigenvectors\n",size,size,nmodes);

	for(i=start_mode; i<end_mode; i++)
	{
		// Normalization
		if(norm_switch)
		{
			vectorrmsd = 0.0;
			for(j=i*size; j<(i+1)*size; j++)
			{ vectorrmsd += eigvect[j]* eigvect[j]; }
			norm = sqrt(vectorrmsd);
		}

		fprintf(out,"****\n%d     %18.10e\n",i-start_mode+1,eigval[i]); // (12/5/2009) Avoids precision waste
		for(j=0; j<size; j++)
		{
			// fprintf(out,"%11.5lf",rounder( eigvect[i*size+j]/norm, 5) ); // precision: 5 decimals
			if(norm_switch) // NORMALIZE
				// fprintf(out,"%11.5lf",eigvect[i*size+j]/norm); // precision: 5 decimals
				fprintf(out,"%18.10e",eigvect[i*size+j]/norm); // (12/5/2009) Avoids precision waste
			else // NOT-NORMALIZE
				// fprintf(out,"%11.5lf",eigvect[i*size+j]); // precision: 5 decimals
				fprintf(out,"%18.10e",eigvect[i*size+j]); // (12/5/2009) Avoids precision waste
			if ((j+1) % 7 ==0)  fprintf(out,"\n");
			// if ((j+1) % 1 ==0)  fprintf(out,"\n");
		}
		if (j % 7 !=0)  fprintf(out,"\n");
		// if (j % 1 !=0)  fprintf(out,"\n");

	}

	fclose(out);
}

// Saves & Normalizes eigenpairs with "ptraj" format (normalizes by default) SINGLE PRECISION
void save_ptraj_modes (char * name, int size, int start_mode, int end_mode, float *eigval, float *eigvect,bool norm_switch)
{
	double vectorrmsd,norm;
	FILE *out;
	int i,j,nmodes;
	nmodes=end_mode-start_mode;
	if( !(out = fopen(name, "w")) )
	{
		printf("Msg(save_ptraj_modes): Unable open to write %s\nForcing exit!\n\n",name);
		exit(1);
	}

	// Header line
	fprintf(out," Eigenvector file: COVAR\n %d %d Contains %d eigenvectors\n",size,size,nmodes);

	for(i=start_mode; i<end_mode; i++)
	{
		// Normalization
		if(norm_switch)
		{
			vectorrmsd = 0.0;
			for(j=i*size; j<(i+1)*size; j++)
			{ vectorrmsd += (double) eigvect[j]* eigvect[j]; }
			norm = sqrt(vectorrmsd);
		}

		fprintf(out,"****\n%d     %18.10e\n",i-start_mode+1,eigval[i]); // (12/5/2009) Avoids precision waste
		for(j=0; j<size; j++)
		{
			// fprintf(out,"%11.5lf",rounder( eigvect[i*size+j]/norm, 5) ); // precision: 5 decimals
			if(norm_switch) // NORMALIZE
				// fprintf(out,"%11.5lf",eigvect[i*size+j]/norm); // precision: 5 decimals
				fprintf(out,"%18.10e",eigvect[i*size+j]/norm); // (12/5/2009) Avoids precision waste
			else // NOT-NORMALIZE
				// fprintf(out,"%11.5lf",eigvect[i*size+j]); // precision: 5 decimals
				fprintf(out,"%18.10e",eigvect[i*size+j]); // (12/5/2009) Avoids precision waste
			if ((j+1) % 7 ==0)  fprintf(out,"\n");
			// if ((j+1) % 1 ==0)  fprintf(out,"\n");
		}
		if (j % 7 !=0)  fprintf(out,"\n");
		// if (j % 1 !=0)  fprintf(out,"\n");

	}

	fclose(out);
}

// Reverse eigenvalues/eigenvectors ordering
// eval --> eigenvalues array
// evec --> eigenvector matrix
// size --> number of coordinates
// nevs --> number of eigenpairs
void reverse_ptraj(double *eval,double *evec,int size,int nevs)
{
	int i,j,k;
	double val; // value buffer
	double *vec; // vector buffer

	if(evec!=NULL)
		if(!(vec=(double *)malloc(sizeof(double)*size)))
		{
			printf("Msg(reverse_ptraj): Sorry, memory allocation failed!\nForcing exit.\n");
			exit(1);
		}

	for(i=0,j=nevs-1; i<(int)nevs/2; i++,j--)
	{
//		fprintf(stderr,"Swaping %d-th by %d-th eigenpairs. nevs/2=%d\n",i,j,(int)nevs/2);
		// Swapping i-th and j-th values
		if(eval!=NULL)
		{
			val = eval[i]; // Buffering i-th value
			eval[i] = eval[j];
			eval[j] = val;
		}

		// Swapping i-th and j-th vectors
		if(evec!=NULL)
			for(k=0; k<size; k++)
			{
				vec[k] = evec[size*i + k]; // Buffering i-th vector
				evec[size*i + k] = evec[size*j + k];
				evec[size*j + k] = vec[k];
			}
	}

	if(evec!=NULL)
		free(vec);
}

// Mass Weight/Un-weight eigenvectors.
// 	evec --> eigenvector matrix
// 	mass --> masses array (length: number of atoms)
// 	size --> number of coordinates
// 	nevs --> number of eigenpairs
// 	weight --> true: weights (multiply by sqrt(mi)), false: un-weights (divide by sqrt(mi))
void weight_ptraj(double *evec,double *mass,int size,int nevs, bool weight)
{
	int i,k;
	double val; // value buffer

	if(weight)
		for(i=0; i<nevs; i++)
			for(k=0; k<size; k++)
				evec[size*i + k] *= mass[(int)k/3];
	else
		for(i=0; i<nevs; i++)
			for(k=0; k<size; k++)
				evec[size*i + k] /= mass[(int)k/3];
}

// Reads an entire ptraj file allocating memory
// Warning: if only 1 mode provided, then --> n_components=-1
int read_ptraj( char * file, double **evect, double **evals, int *numatoms, int *n_vectors, int *n_components  )
{
	FILE *f_file;
	char word[LINE_LENGTH];
	bool correct;
	int n_index=0;
	double *vector,*values;
	float nulo;
	*n_vectors = 0;
	*n_components = 0;

	// oppening nms files
	if( !( f_file = fopen(file, "r") ) )
	{
		printf("\nMsg (read_ptraj): Input file (%s) reading failed!\n\n",file);
		exit(1);
	}

	// Reading all NMs in "nmac_ptraj" format
	correct = false;
	while( fscanf(f_file,"%s",word) != EOF )
	{
		if( strncmp(word,"Eigenvector",11) == 0 )
		{
			correct = true; // If false then Not ptraj file ---> exit
			break;
		}
	}

	if( !correct )
	{ printf("\nMsg (ptraj_read): Not a COVAR/ptraj file detected!!!\n\tForcing exit...\n\n"); exit(1); }

	// Locates the reading pointer at the first eigenvalue component
	while( fscanf(f_file,"%s",word) != EOF )
		if( strncmp(word,"****",4) == 0 ) break;
	(*n_vectors)++;

	while( fscanf(f_file,"%s",word) != EOF )
	{ // Counting number of vectors
		if( strncmp(word,"****",4) == 0 )
		{ (*n_vectors)++; }

		if( *n_vectors == 1 ) // whether it's the first vector...
		{ (*n_components)++; }
	}
//	(*n_components) -= 3; // the first is "****" and then came the two first numbers (the index and the eigenvalue)
	(*n_components) -= 2; // the first is "****" and then came the two first numbers (the index and the eigenvalue)

	*numatoms = *n_components/3; // taking num of atoms from the 3D eigenvectors number (ONLY VALID IF CARTESIAN MODES SUPPLIED)
	//	printf("Msg (ptraj_read): The %s file has: %d atoms, %d vectors and %d components\n",file,*numatoms,*n_vectors,*n_components);
	rewind( f_file );

	// Allocating memory
	vector = (double *) malloc( (*n_vectors) * sizeof(double) * (*n_components) );
	values = (double *) malloc( (*n_vectors) * sizeof(double) );
	// Reading values
	while( fscanf(f_file,"%s",word) != EOF )
	{
		if( strncmp(word,"****",4) == 0 )//	if( !mod_switch ) {
		{ // Reads NM index
			fscanf(f_file,"%d",&n_index);
			n_index -= 1; // "n_index" is the vector index
			fscanf(f_file,"%lf",values+n_index);
			for(int i=0; i< (*n_components); i++)
			{ // Reads NM	vector coordinates
				fscanf(f_file,"%lf",vector + n_index * (*n_components) + i);
			}
		}
	}
	*evect = vector;
	*evals = values;
	return ( 0 );
}

// Reads an entire ptraj file (it allocates memory)
// (This one returns a "trs" structure array)
trs *read_ptraj( char * file, int numatoms, int *n_vectors  )
{
	FILE *f_file;
	char word[FILE_NAME];
	bool correct;
	int n_index=0;
	trs *vector;
	float nulo;
	*n_vectors = 0;
	// oppening nms files
	if( !( f_file = fopen(file, "r") ) )
	{
		printf("\nMsg (ptraj_read): Input file (%s) reading failed!\n\n",file);
		exit(1);
	}

	// Reading all NMs in "nmac_ptraj" format
	correct = false;
	while( fscanf(f_file,"%s",word) != EOF )
	{
		if( strncmp(word,"Eigenvector",11) == 0 )
		{
			correct = true; // If false then Not ptraj file ---> exit
			break;
		}
	}

	if( !correct )
	{ printf("\nMsg (ptraj_read): Not a COVAR/ptraj file detected!!!\n\tForcing exit...\n\n"); exit(1); }
	while( fscanf(f_file,"%s",word) != EOF )
	{ // Counting number of vectors
		if( strncmp(word,"****",4) == 0 )
			(*n_vectors)++;
	}
	rewind( f_file );

	// Allocating memory
	vector = (trs *) malloc( (*n_vectors) * sizeof(trs) * (numatoms) );

	// Reading First file values
	printf("Msg (ptraj_read): Reading file vector values\n");
	while( fscanf(f_file,"%s",word) != EOF )
	{
		if( strncmp(word,"****",4) == 0 )//	if( !mod_switch ) {
		{ // Reads NM index
			fscanf(f_file,"%d %f",&n_index,&nulo);
			n_index -= 1; // "n_index" is the vector index
			for(int i=0; i< numatoms; i++)
			{ // Reads NM	vector coordinates
				fscanf(f_file,"%f %f %f",&(vector + n_index * numatoms + i)->x
						,&(vector + n_index * numatoms + i)->y
						,&(vector + n_index * numatoms + i)->z);
			}
		}
	}
	printf("Msg (ptraj_read): Reading succeded !!!\n");
	return (vector);
}

// Secondary Structure Input (allocates table)
// In the input file, the residue index now begins with 0.
void read_ss(char *file_ss, char **table, int *p_num)
{
	bool debug = false;
	FILE *f_ss;
	int nline=0;
	int index;
	char ss;

	if(!(f_ss=fopen(file_ss,"r")))
	{
		printf("\nMsg(read_ss): Unable to read %s\nForcing Exit!\n\n",file_ss);
		exit(1);
	}

	if(debug)
		printf("Msg(read_ss): Reading Secondary Structure File from: %s\n",file_ss);

	while( fscanf(f_ss,"%d %c",&index,&ss) != EOF )
		nline++;
	*table=(char *)malloc(sizeof(char)*nline); // Allocate memory
	rewind(f_ss);
	nline=0;
	while( fscanf(f_ss,"%d %c",&index,&ss) != EOF )
	{
		if(index<1)
		{
			printf("\nMsg(read_ss): Bad residue index found in %s!\nForcing Exit!\n\n",file_ss);
			exit(1);
		}
		*((*table)+index)=ss;
		nline++;
	}
	fclose(f_ss);

	if(debug)
		printf("Msg(read_ss): Reading Success! (Readed %d residues SS)\n",nline);

	if(p_num != NULL)
		*p_num = nline;
}

// Reads input TS functions
void read_TSfunc(char *file, TSfunc **funcs, int *n_func)
{
	bool debug = false;
	FILE *f_file;
	char *line = NULL;
	int n_line=0;
	size_t len;
	ssize_t read;
	char caracter;
	TSfunc *p_funcs;

	if(debug)
		printf("read_TSfunc> Reading Input File (%s):\n",file);
	if( !(f_file=fopen(file,"r")) )
	{
		printf("\nread_TSfunc> You should include a file with the (initial) functional parameters\n\n");
		exit(1);
	}

	// Counting the number of lines (equal to the number of functions)

	while ((read = getline(&line, &len, f_file)) != -1)
	{
		sscanf(line,"%c",&caracter);
		// Whether the line doesn't begin with "#" or "H" (header info)
		if( !(caracter=='#' || caracter=='*') )
			n_line++; // counts functions
	}

//	std::string lineF;
//	while(std::getline(f_file,lineF))
//	{
//	    if(lineF.size() && (lineF[0] =='#' || lineF[0]=='*') )
//	    	n_line++; // counts functions
//	}


	rewind(f_file);
	if(debug)
		printf("read_TSfunc> Reading %d functions from %s .\n",n_line,file);

	// Allocating memory to store functional parameters
	if( !( p_funcs = (TSfunc *) malloc( n_line * sizeof(TSfunc) ) ) )
	{
		printf("\nread_TSfunc> Unable to allocate memory\n\n");
		exit(1);
	}

	n_line=0;
	while ((read = getline(&line, &len, f_file)) != -1)
	{
		sscanf(line,"%c",&caracter);
		// Whether the line doesn't begin with "#"
		if( !(caracter=='#' || caracter=='*') )
		{
			sscanf(line,"%c%c %d %f %f %f", &(p_funcs[n_line].i), &(p_funcs[n_line].j), &(p_funcs[n_line].t), &(p_funcs[n_line].a), &(p_funcs[n_line].b), &(p_funcs[n_line].c) );
			n_line++;
		}
	}
	if(debug)
		printf("read_TSfunc> %d functions readed.\n",n_line);

	// Showing readed functions
	if(debug)
		for(int i=0; i<n_line; i++)
			printf("read_TSfunc> %c%c %d %f %f %f\n", p_funcs[i].i, p_funcs[i].j, p_funcs[i].t, p_funcs[i].a, p_funcs[i].b, p_funcs[i].c );

	*n_func=n_line;
	*funcs = p_funcs; // outputing...
}

// Read Force constants file (Kfile) allocating memory
// Warning: if "coord" not provided, distances ".d" will not be updated!
void read_Kfile(twid **p_ipas,int *p_nipa,char *file, float *coord)
{
	bool debug=false;
	FILE *f_Ks;
	int nipa=0,i,j,cont;
	int block=cont=1000; // Memory will be allocated in "cont" elements blocks
	twid *decint=NULL;
	float cte;

	// Reading contacts form a force constants file (Klist)
	if(debug)
		printf("read_Kfile> Reading force constants file: %s ",file);
	if( !(f_Ks = (FILE *)fopen(file,"r") ) )
	{
		printf("read_Kfile> Sorry, unable to open contacts file: %s\nCheck input!!!\n",file);
		exit(1);
	}

	decint = ( twid * ) malloc( block * sizeof( twid ) ); // first contact-list structure allocation

//	fprintf(stderr,"read_Kfile(): Breakpoint 0\n");

	while( fscanf(f_Ks,"%d %d %e",&i,&j,&cte) != EOF )
	{
//		printf("i= %d  j= %d  cte= %f\n",i,j,cte);
		i--; // indices run from 1,...,N
		j--;
//		fprintf(stderr,"read_Kfile(): Breakpoint 1\n");

		if( j > i )
		{
			nipa++; // Counts Number of Interacting Pairs of Atoms
			if(cont==0)
			{
//				if(!(decint = ( twid * ) realloc( decint, block * nipa * sizeof( twid ) ))) // resizes contact-list structure
				if(!(decint = ( twid * ) realloc( decint, (block + nipa) * sizeof( twid ) ))) // resizes contact-list structure
				{
					fprintf(stderr,"Sorry, realloc failed!\n");
					exit(1);
				}
				cont=block;
			}
			cont--;
//			decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact-list structure
			decint[nipa - 1].k = i; // k-pseudo-atom index (i-atom index)
			decint[nipa - 1].l = j; // l-pseudo-atom index (j-atom index)
			if(coord!=NULL)
				decint[nipa - 1].d = sqrt( pow(coord[3*i]-coord[3*j],2) +
						pow(coord[3*i+1]-coord[3*j+1],2) +
						pow(coord[3*i+2]-coord[3*j+2],2) ); // set distance
			decint[nipa - 1].C = cte; // force constant
		}
//		fprintf(stderr,"read_Kfile(): Breakpoint 2\n");
	}
	if(debug)
		printf("read_Kfile> (%d nipas)\n",nipa);
	fclose( f_Ks );
	*p_nipa = nipa;
	*p_ipas = decint;
//	fprintf(stderr,"read_Kfile(): Breakpoint 3\n");
}

// Reads a file with the residue indices whose internal coordinates will be fixed.
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int fixFrag( char *file, Macromolecule *mol, tri *props,bool *fixed)
{
	// Mobile Internal Coordinates selection
	int num_fix=0;
	int num_ic = 0; // counts the mobile IC's
	int ic_index = 0;
	int res_index = 0;
	int num_res;

	FILE *ffix;
	if( !(ffix = fopen(file,"r")) )
	{
		printf("Msg(fselectIC): Reading fix-file: %s, Failed!\n",file);
		exit(2);
	}

	Segment *seg;
	pdbIter *iter_seg = new pdbIter( mol ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues
	num_res = mol->get_num_fragments();

	fscanf(ffix,"%d",&num_fix); // gets first fixed residue number
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		// Inter-segment Internal Coordinates are allways mobile
		if(iter_seg->pos_segment !=0)
		{
			for(int i=0; i<6; i++,ic_index++)
			{
				fixed[ic_index] = true; // mobile
				num_ic++; // counts the mobile IC's
			}
		}

		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			for(int k=0; k < props[res_index].nan; k++)
			{
				if( res_index == num_fix )
				{
					fixed[ic_index] = false; // fixed
				}
				else
				{
					fixed[ic_index] = true; // mobile
					num_ic++; // counts the mobile IC's
				}
				ic_index++;
			}
			if( res_index == num_fix && num_fix < num_res)
				fscanf(ffix,"%d",&num_fix);
			res_index++;
		}
		delete iter_frags;
	}
	delete(iter_seg);
	fclose(ffix);
	return num_ic;
}

// Reads a file with the fixed internal coordinates indices.
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int read_fixIC( char *file, Macromolecule *mol, tri *props,bool *fixed)
{
	bool debug=false;
	// Mobile Internal Coordinates selection
	int num_fix=0;
	int num_ic = 0; // counts the mobile IC's
	int ic_index = 0;
	int res_index = 0;
	int num_res;

	FILE *ffix;
	if( !(ffix = fopen(file,"r")) )
	{
		printf("Msg(fselectIC): Reading fix-file: %s, Failed!\n",file);
		exit(2);
	}

	Segment *seg;
	pdbIter *iter_seg = new pdbIter( mol ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues
	num_res = mol->get_num_fragments();

	// First, setting all mobile...
	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		// Inter-segment Internal Coordinates are allways mobile
		if(iter_seg->pos_segment !=0)
		{
			for(int i=0; i<6; i++,ic_index++)
				fixed[ic_index] = true; // mobile
		}

		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			for(int k=0; k < props[res_index].nan; k++,ic_index++)
				fixed[ic_index] = true; // mobile
			res_index++;
		}
		delete iter_frags;
	}
	delete(iter_seg);

	// Now fixing
	while( fscanf(ffix,"%d",&num_fix) != EOF )
	{
		fixed[num_fix] = false; // fixed
		if(debug)
			printf("Fixed %4d!\n",num_fix);
		num_ic++;
	}

	fclose(ffix);
	return ic_index-num_ic; // Returns the number of mobile variables
}

// Reads a file with the residue indices whose internal coordinates will be fixed.
// FORMAT (protein):      "[residue-index] [phi] [chi] [psi]"
// FORMAT (nucleic acid): "[fragment-index] [alpha] [beta] [gamma] [chi] [epsilon] [zeta]"
//   [angle]= (0=fix / 1=mobile) All allways required! (according to fragment ID)
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int read_fixDH( char *file, Macromolecule *mol, tri *props,bool *fixed, int type, int model, bool *addrot)
{
	bool debug=false;
	// Mobile Internal Coordinates selection
	int num_fix=0;
	int num_ic = 0; // counts the mobile IC's
	int ic_index = 0;
	int index_res = 0;
	int num_res;
	int phi,chi,psi,rotrans,resn;
	int iang_max;
	int angle[6]; // allocates 0/1 (fixed/mobile) array
	TMOL fragtype;

	FILE *ffix;
	if( !(ffix = fopen(file,"r")) )
	{
		printf("Msg(read_fixDH): Reading fix-file: %s, Failed!\n",file);
		exit(2);
	}

	Segment *seg;
	Residue *res;
	pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues

	int seg_atoms_old,interchain=0;
	int seg_atoms = 0; // should be initialized
	bool may_have_rot = false;

	if(debug)
		printf("Msg(fselectDH): Reading fix-file: %s:\n",file);

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		num_res = iter_frags->num_fragment();

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms;
		seg_atoms = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		// Inter-segment Internal Coordinates are always mobile
		if(iter_seg->pos_segment !=0)
		{
			if(seg_atoms != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
				interchain = 6;
			else
				interchain = 3;

			for(int i=0; i<interchain; i++,ic_index++)
			{
				fixed[ic_index] = true; // mobile
				num_ic++; // counts the mobile IC's
				fscanf(ffix,"%d %d",&num_fix,&rotrans);
				if(debug)
					printf("seg= %2d  -->  %d %d\n",iter_seg->pos_segment,num_fix,rotrans);

				if(num_fix==index_res)
				{
					if(rotrans == 0)
					{
						fixed[ic_index] = false; // fix
						num_ic--;
					}
				}
				else
				{
					printf("Msg(read_fixDH): I'm sorry index mismatch in file: %s\n",file);
					exit(2);
				}
			}
		}

		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			res = ( Residue * ) iter_frags->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				fscanf(ffix,"%d %d %d %d",&num_fix,&phi,&chi,&psi);
				if(debug)
					printf("%d %d %d %d  type= %d\n",num_fix,phi,chi,psi,fragtype);

				// PHI
				if( !(iter_frags->pos_fragment == 0 || resn == PRO
						|| (iter_frags->pos_fragment == num_res-1 && (model == 0 || model == 3) ) ) )
				{
					if(phi!=0)
					{
						if(debug)
							printf("Mobile PHI ic= %d\n",ic_index);
						fixed[ic_index] = true; // mobile
						num_ic++; // counts the mobile IC's
					}
					else
					{
						if(debug)
							printf("Fixing PHI ic= %d\n",ic_index);
						fixed[ic_index] = false; // fixed
					}
					ic_index++;
				}

				// CHI
				if( type == 2 && model != 0 && model != 3)
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frags->pos_fragment==0 ||
											(iter_frags->pos_fragment==num_res-1 && model==1) ) ) )
					{
						if(chi!=0)
						{
							if(debug)
								printf("Mobile CHI ic= %d\n",ic_index);
							fixed[ic_index] = true; // mobile
							num_ic++; // counts the mobile IC's
						}
						else
						{
							if(debug)
								printf("Fixing CHI ic= %d\n",ic_index);
							fixed[ic_index] = false; // fixed
						}

						ic_index++;
					}

				// PSI
				if( !(iter_frags->pos_fragment == num_res-1
						|| (iter_frags->pos_fragment == 0 && (model==0 || model==3) ) )
						|| model==2 ) // Full-Atom allways has PSI
				{
					if(psi!=0)
					{
						if(debug)
							printf("Mobile PSI ic= %d\n",ic_index);
						fixed[ic_index] = true; // mobile
						num_ic++; // counts the mobile IC's
					}
					else
					{
						if(debug)
							printf("Fixing PSI ic= %d\n",ic_index);
						fixed[ic_index] = false; // fixed
					}
					ic_index++;
				}
			}
			// if RNA/DNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA ||
					 resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				fscanf(ffix,"%d %d %d %d %d %d %d",&num_fix,&angle[0],&angle[1],&angle[2],&angle[3],&angle[4],&angle[5]);
				if(debug)
					printf("%d %d %d %d %d %d %d\n",num_fix,angle[0],angle[1],angle[2],angle[3],angle[4],angle[5]);

				if(iter_frags->pos_fragment != num_res-1)
					iang_max=6; // not last fragment
				else
					iang_max=4; // last fragment

				if(debug)
					printf("fixed:");
				for(int iang=0; iang<iang_max; iang++)
				{
					if( !(iang == 3 && type != 2) ) // if not (is CHI and type!=2)
					{
						if(angle[iang] != 0)
						{
							fixed[ic_index] = true; // mobile
							num_ic++; // counts the mobile IC's
						}
						else
							fixed[ic_index] = false; // fixed
						if(debug)
							printf(" %d",fixed[ic_index]);
						ic_index++;
					}
				}
				if(debug)
					printf("\n");
			}
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype != tmol_smol )
			{
				printf("Msg(read_fixDH): Sorry, unknown fragment type! Forcing exit!\n");
				exit(1);
			}

			index_res++;
		}
		delete iter_frags;
	}
	iter_seg->clean_virtual();
	delete iter_seg;
	fclose(ffix);
	return num_ic;
}

// Reads the "fix" file and outputs a "number of residues sized" array of "ints" with the indices
// of those "rigid bodies" as a function of the residues. It searchs for chunks of continuously fixed
// sequences of residues and assigns correlative indices to the corresponding "rigid bodies".
//
// File FORMAT (protein):      "[residue-index] [phi] [chi] [psi]"
// File FORMAT (nucleic acid): "[fragment-index] [alpha] [beta] [gamma] [chi] [epsilon] [zeta]"
//   [angle]= (0=fix / 1=mobile) All always required! (according to fragment ID)
// Memory for the "fixed" array must be already allocated (for the total number of residues)
// Returns the number of rigid bodies
int read_fixRes( char *file, Macromolecule *mol, int *fixed)
{
	bool debug=false;
	int rigid_body = 0; // Rigid body index
	bool isfixed = false; // Does it have any inter-segment rotational/translational DoF?
	int num_fix=0;
	int index_res = 0;
	int num_res;
	int phi,chi,psi,rotrans,resn;
	int iang_max;
	int angle[6]; // allocates 0/1 (fixed/mobile) array
	TMOL fragtype;

	FILE *ffix;
	if( !(ffix = fopen(file,"r")) )
	{
		printf("Msg(read_fixRes): Reading fix-file: %s, Failed!\n",file);
		exit(2);
	}

	Segment *seg;
	Residue *res;
	pdbIter *iter_seg = new pdbIter( mol, true, true, true, true ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues

	int seg_atoms_old,interchain=0;
	int seg_atoms = 0; // should be initialized
	bool may_have_rot = false;
	bool last_frag = false; // true if last fragment was flexible

	if(debug)
		printf("Msg(read_fixRes): Reading fix-file: %s:\n",file);

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		num_res = iter_frags->num_fragment();

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms;
		seg_atoms = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		// Inter-segment Internal Coordinates
		if(iter_seg->pos_segment !=0)
		{
			if(seg_atoms != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
				interchain = 6;
			else
				interchain = 3;
			isfixed = false; // by default, false
			for(int i=0; i<interchain && !isfixed; i++)
			{
				fscanf(ffix,"%d %d",&num_fix,&rotrans);
				if(debug)
					printf("seg= %2d  -->  %d %d\n",iter_seg->pos_segment,num_fix,rotrans);

				if(num_fix==index_res)
				{
					if(rotrans == 0)
						isfixed = true;
				}
				else
				{
					printf("Msg(read_fixDH): I'm sorry index mismatch in file: %s\n",file);
					exit(2);
				}
			}
			if(isfixed && !last_frag) // Only if last fragment was not considered mobile so far.
				rigid_body++; // Split current rigid-body

			last_frag = false; // should be reseted after use
		}

		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			res = ( Residue * ) iter_frags->get_fragment();

			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			fixed[index_res] = rigid_body;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				fscanf(ffix,"%d %d %d %d",&num_fix,&phi,&chi,&psi);
				if(debug)
					printf("%d %d %d %d  type= %d\n",num_fix,phi,chi,psi,fragtype);

//				if( !(phi == 0 && chi == 0 && psi == 0) )
//				if( !(phi == 0 && psi == 0) ) // If PHI or PSI mobile
//				{
//					rigid_body++; // Split current rigid-body
//					if(iter_frags->pos_fragment == num_res-1)
//						last_frag = true; // if last fragment is flexible
//				}
				if(iter_frags->pos_fragment == num_res-1) // if last fragment
				{
					if( !(phi == 0 && psi == 0) ) // If PHI or PSI mobile (if flexible)
						last_frag = true; // if last fragment is flexible
				}
				else // if not-last fragment
					if( !(phi == 0 && psi == 0) ) // If PHI or PSI mobile (if flexible)
						rigid_body++; // Split current rigid-body

			}
			// if RNA/DNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA ||
					 resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				fscanf(ffix,"%d %d %d %d %d %d %d",&num_fix,&angle[0],&angle[1],&angle[2],&angle[3],&angle[4],&angle[5]);
				if(debug)
					printf("%d %d %d %d %d %d %d\n",num_fix,angle[0],angle[1],angle[2],angle[3],angle[4],angle[5]);

				if(iter_frags->pos_fragment != num_res-1)
					iang_max=6; // not last fragment
				else
					iang_max=4; // last fragment

				if(debug)
					printf("fixed:");
				isfixed = false; // by default, false
				for(int iang=0; iang<iang_max && !isfixed; iang++)
				{
					if(angle[iang] != 0)
						isfixed = true;
				}
				if(debug)
					printf("\n");

				if(iter_frags->pos_fragment == num_res-1) // if last fragment
				{
					if(isfixed) // if flexible
						last_frag = true; // if last fragment is flexible
				}
				else // if not-last fragment
					if(isfixed) // if flexible
						rigid_body++; // Split current rigid-body
			}
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype != tmol_smol )
			{
				printf("Msg(read_fixRes): Sorry, unknown fragment type! Forcing exit!\n");
				exit(1);
			}

			index_res++;
		}
		delete iter_frags;
	}
	iter_seg->clean_virtual();
	delete iter_seg;
	fclose(ffix);
	return (rigid_body+1);
}

// Writes a fixation-file with the residue indices and which internal coordinates are fixed.
// FORMAT (protein):      "[residue-index] [phi] [chi] [psi]"
// FORMAT (nucleic acid): "[fragment-index] [alpha] [beta] [gamma] [chi] [epsilon] [zeta]"
//   [angle]= (0=fix / 1=mobile) All always required! (according to fragment Mol-type: protein/nucleic-acid)
// "fixed" array must be already made! If fixed=NULL, then a fully mobile file will be written.
void write_fixDH( char *file, Macromolecule *mol, tri *props,bool *fixed, int type, int model)
{
	bool debug=false;
	// Mobile Internal Coordinates selection
	int ic_index = 0;
	int index_res = 0;
	int num_res;
	int resn;
	int iang_max;
	int angle[6]; // allocates 0/1 (fixed/mobile) array
	TMOL fragtype;

	FILE *ffix;
	if( !(ffix = fopen(file,"w")) )
	{
		printf("Msg(write_fixDH): Writing fix-file: %s, Failed!\n",file);
		exit(2);
	}

	Segment *seg;
	Residue *res;
	pdbIter *iter_seg = new pdbIter( mol, true, true ,true, true ); // Iterator to screen segments
	pdbIter *iter_frags; // iter residues

	int seg_atoms_old,interchain=0;
	int seg_atoms = 0; // should be initialized
	bool may_have_rot = false;

	if(debug)
		printf("Msg(write_fixDH): Writing fix-file: %s:\n",file);

	for( iter_seg->pos_segment = 0; !iter_seg->gend_segment(); iter_seg->next_segment() )
	{
		seg = iter_seg->get_segment();
		iter_frags = new pdbIter( seg );
		num_res = iter_frags->num_fragment();

		// This fix the bug... (for a nucleotide-residue "getMolType()" does not work ¿?
		fragtype = seg->getMolType();

		// Checking single atom segments... SMOL's
		seg_atoms_old = seg_atoms;
		seg_atoms = seg->num_atoms();
		if(seg_atoms_old > 1)
			may_have_rot = true;

		if(debug)
			printf("seg= %d  num_res= %d  seg_atoms= %d (old=%d) may_have_rot= %d\n",iter_seg->pos_segment,num_res,seg_atoms,seg_atoms_old,may_have_rot);

		// Inter-segment Internal Coordinates are always mobile
		if(iter_seg->pos_segment !=0)
		{
// MON: Revise this,... shomething may be removed... "may_have_rot"
			if(seg_atoms != 1 && (seg_atoms_old != 1 || may_have_rot) ) // NON single atom segments
				interchain = 6;
			else
				interchain = 3;

			for(int i=0; i<interchain; i++,ic_index++)
			{
				if(fixed!=NULL)
				{
					fprintf(ffix,"%d %d",index_res,fixed[ic_index]);
					if(debug)
						printf("%d %d  inter= %d\n",index_res,fixed[ic_index],interchain);
				}
				else
				{
					fprintf(ffix,"%d %d",index_res,1);
					if(debug)
						printf("%d %d  inter= %d\n",index_res,1,interchain);
				}
				fprintf(ffix,"\n"); // line jump
			}
		}

		for( iter_frags->pos_fragment = 0; !iter_frags->gend_fragment(); iter_frags->next_fragment() )
		{
			res = ( Residue * ) iter_frags->get_fragment();
//			fragtype = res->getMolType();
			if(fragtype != tmol_smol)
				resn = resnum_from_resname( res->getName() );
			else
				resn = 666;

			// if Protein fragment
			if( fragtype == tmol_protein )
			{
				fprintf(ffix,"%d",index_res); // fragment index
				if(debug)
					printf("%d",index_res);

				// PHI
				if( !(iter_frags->pos_fragment == 0 || resn == PRO
						|| (iter_frags->pos_fragment == num_res-1 && (model == 0 || model == 3) ) ) )
				{
					if(fixed!=NULL)
					{
						fprintf(ffix," %d",fixed[ic_index]);
						if(debug)
							printf(" %d (phi)\n",fixed[ic_index]);
					}
					else
					{
						fprintf(ffix," 1");
						if(debug)
							printf(" 1 (phi)\n");
					}
					ic_index++;
				}
				else
				{   // hasn't got PHI
					fprintf(ffix," 0");
					if(debug)
						printf(" 0 (phi)");
				}

				// CHI
				if( type == 2 && model != 0 && model != 3)
				{
					if( props[index_res].nan==3 ||
							(props[index_res].nan==2 &&
									(iter_frags->pos_fragment==0 ||
											(iter_frags->pos_fragment==num_res-1 && model==1) ) ) )
					{   // has CHI
						if(fixed!=NULL)
						{
							fprintf(ffix," %d",fixed[ic_index]);
							if(debug)
								printf(" %d (chi)\n",fixed[ic_index]);
						}
						else
						{
							fprintf(ffix," 1");
							if(debug)
								printf(" 1 (chi)\n");
						}
						ic_index++;
					}
					else
					{   // hasn't CHI
						fprintf(ffix," 0");
						if(debug)
							printf(" 0 (chi)");
					}
				}
				else
				{   // hasn't CHI
					fprintf(ffix," 0");
					if(debug)
						printf(" 0 (chi)");
				}

				// PSI
				if( !(iter_frags->pos_fragment == num_res-1
						|| (iter_frags->pos_fragment == 0 && (model==0 || model==3) ) )
						|| model==2 ) // Full-Atom allways has PSI
				{
					if(fixed!=NULL)
					{
						fprintf(ffix," %d",fixed[ic_index]);
						if(debug)
							printf(" %d (psi)\n",fixed[ic_index]);
					}
					else
					{
						fprintf(ffix," 1");
						if(debug)
							printf(" 1 (psi)\n");
					}
					ic_index++;
				}
				else
				{   // hasn't PSI
					fprintf(ffix," 0");
					if(debug)
						printf(" 0 (psi)");
				}
				fprintf(ffix,"\n"); // line jump
				if(debug)
					printf("\n");
			}
			// if RNA fragment
			else if( resn == GUA || resn == ADE || resn == CYT || resn == URA ||
				 resn == DGUA || resn == DADE || resn == DCYT || resn == DTHY )
			{
				fprintf(ffix,"%d",index_res); // fragment index
				if(debug)
					printf("%d",index_res);

				if(iter_frags->pos_fragment != num_res-1)
					iang_max=6; // not last fragment
				else
					iang_max=4; // last fragment

				for(int iang=0; iang<iang_max; iang++)
				{
					if(iang != 3) // if not CHI
					{
						if(fixed!=NULL)
							fprintf(ffix," %d",fixed[ic_index]);
						else
							fprintf(ffix," %d",1);
						if(debug)
							printf(" %d",fixed[ic_index]);
						ic_index++;
					}
					else // if CHI
						if(type == 2) // has CHI
						{
							if(fixed!=NULL)
								fprintf(ffix," %d",fixed[ic_index]);
							else
								fprintf(ffix," %d",1);
							if(debug)
								printf(" %d (chi)",fixed[ic_index]);
							ic_index++;
						}
						else // has not CHI
						{
							fprintf(ffix," %d",0);
							if(debug)
								printf(" %d (chi)",0);
						}
				}

				if(iter_frags->pos_fragment == num_res-1)
				{
					fprintf(ffix," %d %d",0,0);
					if(debug)
						printf(" %d %d",0,0);
				}

				fprintf(ffix,"\n"); // line jump
				if(debug)
					printf("\n");
			}
			// For SMOL - Small MOLecules (ligands) --> Nothing "flexible" should be done here.
			else if( fragtype == tmol_smol )
			{
				if(debug)
					printf("%d (SMOL)\n",index_res);
			}
			else
			{
				printf("Msg(write_fixDH): Sorry, unknown fragment type! Forcing exit!\n");
				exit(1);
			}

			index_res++;
		}
		delete iter_frags;
	}
	iter_seg->clean_virtual();
	delete(iter_seg);
	fclose(ffix);
}

// It writes an IC coordinates fixation mask (fixIC-format)
void write_fixIC(char *text,bool *fixed,int size)
{
	FILE *f_fix;
	if( !(f_fix=(FILE *)fopen(text,"w") ) )
	{
		printf("Msg(write_fixIC): Sorry, unable to write FIX FILE: %s\n",text);
		exit(1);
	}
	for(int i=0; i<size; i++)
	{
		if(!fixed[i]) // if DoF is fixed
			fprintf(f_fix,"%d\n",i);
	}
	fclose(f_fix);
}

// Write Force constants file (Kfile)
void write_Kfile(twid *ipas,int nipa,char *file)
{
	FILE *f_Kout;
	if( !(f_Kout = (FILE *)fopen(file,"w") ) )
	{
		printf("Msg(write_Kfile): Sorry, unable to open Output contacts file: %s\nCheck input!!!\n",file);
		exit(1);
	}

	for(int i=0; i<nipa; i++)
		fprintf(f_Kout,"%d %d %E\n",ipas[i].k+1,ipas[i].l+1,ipas[i].C);

	fclose(f_Kout);
}

// Write Secondary Structure file
void write_SSfile(char *ss,int nss,char *file)
{
	FILE *f_SS;
	if( !(f_SS = (FILE *)fopen(file,"w") ) )
	{
		printf("write_SSfile> Sorry, unable to open output SS file: %s\nCheck input!!!\n",file);
		exit(1);
	}

	for(int i=0; i<nss; i++)
		fprintf(f_SS,"%5d %c\n",i+1,ss[i]);

	fclose(f_SS);
}

// Shows a triangular packed matrix (standard output)
void show_matrix(floating *matrix, int size, char *name)
{
	int i=0;
	int j=0;
	printf("%s\n",name);
	for ( i = 0; i < size; i++ )
	{
		for ( j = 0; j < size; j++ )
		{
			if(j>=i)
				printf("%12.5f ",matrix[i + j*(j+1)/2]);
			else
				printf("%12.5f ",0.0);
		}
		printf("\n");
	}
}

// Shows a triangular packed matrix (standard output) with single precision
void show_matrix(float *matrix, int size, char *name)
{
	int i=0;
	int j=0;
	printf("%s\n",name);
	for ( i = 0; i < size; i++ )
	{
		for ( j = 0; j < size; j++ )
		{
			if(j>=i)
				printf("%10.5f ",matrix[i + j*(j+1)/2]);
			else
				printf("%10.5f ",0.0);
		}
		printf("\n");
	}
}

// Saving a triangular packed matrix (into file) (Huge)
// type --> 0= default, 1= save indices: [i] [j] [covariance], 2= Gnuplot matrix format (row and col wise)
// save_sym --> true, save full matrix (apply symmetry to obtain an squared symmetric matrix)
void save_matrix(floating *matrix, int size, char *name, int type, bool save_sym)
{
	long int i=0;
	long int j=0; // Needed for huge systems!
	FILE *f_file;

	f_file = fopen(name, "w");
	if(type != 2)
		fprintf(f_file,"%d\n",size);
	switch(type)
	{
	case 1:
		if(save_sym)
		{
			for ( i = 0; i < size; i++ )
				for ( j = i+1; j < size; j++ ) // diagonal not included
					fprintf(f_file,"%ld %ld %22.15e\n%ld %ld %22.15e\n",i,j,matrix[i + j*(j+1)/2],j,i,matrix[i + j*(j+1)/2]);
		}
		else
		{
			for ( i = 0; i < size; i++ )
				for ( j = i+1; j < size; j++ ) // diagonal not included
					fprintf(f_file,"%ld %ld %22.15e\n",i,j,matrix[i + j*(j+1)/2]);
		}
		for ( i = 0; i < size; i++ ) // diagonal
			fprintf(f_file,"%ld %ld %22.15e\n",i,i,matrix[i + i*(i+1)/2]);
		break;
	case 2:
		for ( i = 0; i < size; i++ )
		{
			for ( j = 0; j < size; j++ )
				if(j >= i)
					fprintf(f_file,"%f ",matrix[i + j*(j+1)/2]);
				else
					fprintf(f_file,"%f ",matrix[j + i*(i+1)/2]);
			fprintf(f_file,"\n");
		}
		break;
	case 0:
		for ( i = 0; i < size; i++ )
			for ( j = i; j < size; j++ ) // diagonal included
				fprintf(f_file,"%22.15e\n",matrix[i + j*(j+1)/2]);
		break;
	}

	fclose(f_file);
}

// Saving a triangular packed matrix (into BINARY file) (Huge)
void save_matrixB(floating *matrix, int size, char *name)
{
	int i=0;
	long int j=0; // Needed for huge systems!
	FILE *f_file;
	floating dsize;
	dsize = (floating) size; // casting integer into floating

	f_file = fopen(name, "w");
	fwrite(&dsize,sizeof(floating),1,f_file); // first number is the matrix rank
	for ( i = 0; i < size; i++ )
		for ( j = i; j < size; j++ ) // diagonal included
			fwrite(&matrix[i + j*(j+1)/2],sizeof(floating),1,f_file);
	fclose(f_file);
}

// Save array into a plain-text file (DOUBLE)
void save_array(double *array,int size,char *file)
{
	FILE *f_file;
	if( !( f_file = fopen(file,"w") ) )
	{
		printf("Msg(save_array): I'm sorry, unable to open: %s\n",file);
		exit(1);
	}

	for(int i=0; i<size; i++)
		fprintf(f_file,"%22.15e\n",array[i]);

	fclose(f_file);
}

// Save array into a plain-text file (SINGLE)
void save_array(float *array,int size,char *file)
{
	FILE *f_file;
	if( !( f_file = fopen(file,"w") ) )
	{
		printf("Msg(save_array): I'm sorry, unable to open: %s\n",file);
		exit(1);
	}

	for(int i=0; i<size; i++)
		fprintf(f_file,"%22.15e\n",array[i]);

	fclose(f_file);
}

// Saving a rectangular matrix (into file)
void save_matrix_rec(floating *matrix, int vec, int size, char *name)
{
	int i=0;
	int j=0;
	FILE *f_file;

	f_file = fopen(name, "w");
	fprintf(f_file,"%d %d\n",vec,size);
	for ( i = 0; i < vec; i++ )
		for ( j = 0; j < size; j++ ) // diagonal included
			fprintf(f_file,"%22.15e\n",matrix[j + i*size]);
	fclose(f_file);
}

// Saving a rectangular matrix (into BINARY file)
void save_matrix_recB(floating *matrix, int vec, int size, char *name)
{
	int i=0;
	int j=0;
	FILE *f_file;
	floating dvec,dsize;
	dvec = (floating) vec; // casting integer into floating
	dsize = (floating) size; // casting integer into floating

	f_file = fopen(name, "w");
	fwrite(&dvec,sizeof(floating),1,f_file); // first number is the number of vectors
	fwrite(&dsize,sizeof(floating),1,f_file); // second number is the matrix rank
	for ( i = 0; i < vec; i++ )
		for ( j = 0; j < size; j++ ) // diagonal included
			fwrite(&matrix[j + i*size],sizeof(floating),1,f_file);
	fclose(f_file);
}

// Saving a rectangular matrix (into file) with single precision
void save_matrix_rec(float *matrix, int vec, int size, char *name)
{
	int i=0;
	int j=0;
	FILE *f_file;

	f_file = fopen(name, "w");
	fprintf(f_file,"%d %d\n",vec,size);
	for ( i = 0; i < vec; i++ )
		for ( j = 0; j < size; j++ ) // diagonal included
			fprintf(f_file,"%18.10e\n",matrix[j + i*size]);
	fclose(f_file);
}

// Reading a triangular packed matrix (into memory) (Huge)
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrix(floating **p_matrix, int *p_size, char *name)
{
	int i=0;
	long int j=0;
	long int size;
	FILE *f_file;
	floating *matrix;

	if( !(f_file = fopen(name, "r") ) )
	{
		printf("Msg(read_matrix): I'm sorry, unable to open %s file\nForcing exit\n",name);
		exit(1);
	}
	fscanf(f_file,"%d",p_size); // getting size
	size = *p_size;
//	printf("Msg(read_matrix): %d rank matrix detected\n",size);

	// allocating memory
	if(*p_matrix == NULL)
	{
		if( !(matrix = (floating *) malloc( sizeof( floating ) * size*(size+1)/2) ) )
		{
			printf("Msg(read_matrix): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( floating ) * size*(size+1)/2);
			exit(1);
		}
		else
			*p_matrix = matrix;
	}

	for ( i = 0; i < size; i++ )
		for ( j = i; j < size; j++ ) // diagonal included
			fscanf(f_file,"%lf",&matrix[i + j*(j+1)/2]);
	fclose(f_file);
}

// Reading a triangular packed matrix (into memory) from BINARY file. (Huge)
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrixB(floating **p_matrix, int *p_size, char *name)
{
	int i=0;
	long int j=0;
	long int size;
	FILE *f_file;
	floating *matrix;
	floating dsize;
	fprintf(stderr,"Hi I'm inside read_matrixB ...\n");

	if( !(f_file = fopen(name, "r") ) )
	{
		printf("Msg(read_matrixB): I'm sorry, unable to open %s file\nForcing exit\n",name);
		exit(1);
	}
	fread(&dsize,sizeof(floating),1,f_file); // getting size
	*p_size = (int) dsize; // casting double into integer
	size = *p_size;
//	printf("Msg(read_matrix): %d rank matrix detected\n",size);
//	fflush(stdout);

	// allocating memory
	if(*p_matrix == NULL)
	{
		if( !(matrix = (floating *) malloc( sizeof( floating ) * size*(size+1)/2) ) )
		{
			printf("Msg(read_matrixB): I'm sorry, unable to allocate %ld bytes\nForcing exit\n",(long int) sizeof( floating ) * size*(size+1)/2);
			exit(1);
		}
		else
			*p_matrix = matrix;
	}

	fprintf(stderr,"Hi I've allocated memory (inside read_matrixB) size= %ld ...\n",size);

	for ( i = 0; i < size; i++ )
		for ( j = i; j < size; j++ ) // diagonal included
			fread(&matrix[i + j*(j+1)/2],sizeof(floating),1,f_file);
	fprintf(stderr,"Hi I've finished (inside read_matrixB) ...\n");
	fclose(f_file);
	fprintf(stderr,"Hi I've closed file (inside read_matrixB) ...\n");
}

// Reading a triangular packed matrix (into memory) from BINARY file (Reverse-row-wise)
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrixB2(floating **p_matrix, int *p_size, char *name)
{
	int i=0;
	long int j=0;
	long int size;
	FILE *f_file;
	floating *matrix;
	floating dsize;

	if( !(f_file = fopen(name, "r") ) )
	{
		printf("Msg(read_matrixB2): I'm sorry, unable to open %s file\nForcing exit\n",name);
		exit(1);
	}
	fread(&dsize,sizeof(floating),1,f_file); // getting size
	*p_size = (int) dsize; // casting double into integer
	size = *p_size;
//	printf("Msg(read_matrix): %d rank matrix detected\n",size);
//	fflush(stdout);

	// allocating memory
	if(*p_matrix == NULL)
	{
		if( !(matrix = (floating *) malloc( sizeof( floating ) * size*(size+1)/2) ) )
		{
			printf("Msg(read_matrixB2): I'm sorry, unable to allocate %d bytes\nForcing exit\n",sizeof( floating ) * size*(size+1)/2);
			exit(1);
		}
		else
			*p_matrix = matrix;
	}

//	for ( i = 0; i < size; i++ )
//		for ( j = i; j < size; j++ ) // diagonal included
//			fread(&matrix[i + j*(j+1)/2],sizeof(floating),1,f_file);

	for(j=size-1; j >= 0; j--) // b --> dihedrals (column)
		for(i=0; i <= j; i++)	 // i --> dihedrals (row)
			fread(&matrix[i + j*(j+1)/2],sizeof(floating),1,f_file);

	fclose(f_file);
}


// *******************************************************************************
// PABLO's OLD IO-FUNCTIONS (DEFPROT-back compatibility)
// *******************************************************************************
void save_ascii_modes (int size, int start_mode, int end_mode, float *eigval, float *eigvect)
{
	float current_ampl, vectorrmsd,maximumamplitude;
	char file_name[8];
	FILE *out;
	int i,j;

	printf( "nmac>\nnmac>         Saving  modes in ascii files    \nnmac>\n" );
	printf( "nmac>                  Eigenvalue    File\nnmac>\n" );


	/* first six must be zero or close to, not saved */
	printf("nmac>\nnmac>          MODE    EIGENVALUE\nnmac>\n");
	for(i=0; i<start_mode; i++) {
		printf("nmac>   Mode %4d %15.5E  not saved\n",i+1,eigval[i]);
	}

	for(i=start_mode; i<end_mode; i++) {
		sprintf(file_name,"mod%d",i+1);
		out = fopen(file_name, "w");
		vectorrmsd = 0.0;maximumamplitude=-10000;
		for(j=0; j<size; j+=3) {
			fprintf(out,"%d %12.5E %12.5E %12.5E\n",j/3+1,
					eigvect[i*size+j+0],
					eigvect[i*size+j+1],
					eigvect[i*size+j+2]);
			current_ampl =  eigvect[i*size+j+0]* eigvect[i*size+j+0] +
			eigvect[i*size+j+1]* eigvect[i*size+j+1] +
			eigvect[i*size+j+2]* eigvect[i*size+j+2];
			vectorrmsd+=current_ampl;
			if (maximumamplitude < current_ampl )
				maximumamplitude = current_ampl;
		}
		fclose(out);
		maximumamplitude = sqrt(maximumamplitude);
		vectorrmsd = sqrt(vectorrmsd/((double) (size/3.0)));
		// if (i+start_mode<20) printf("nmac>   Mode %4d %15.5E  %s\n",i+1,eigval[i],file_name);
	}
	fprintf(stderr,"nmac>\nnmac>     Max amplitude %f Vector rmsd %f\nnmac>\n",
			vectorrmsd,maximumamplitude);
}

void save_ptraj_modes (int size, int start_mode, int end_mode, float *eigval, float *eigvect)
{
	float current_ampl, vectorrmsd,maximumamplitude,norm;
	char filename[20];
	FILE *out;
	int i,j,nmodes;
	nmodes=end_mode-start_mode;
	sprintf(filename,"nmac_ptraj.evec");
	out = fopen(filename, "w");
	fprintf(out," Eigenvector file: COVAR\n %d %d Contains %d eigenvectors\n",size,size,nmodes);
	/* first six must be zero or close to, not saved */
	printf("nmac>\nnmac>          MODE    EIGENVALUE\nnmac>\n");
	for(i=0; i<start_mode; i++) {
		printf("nmac>   Mode %4d %15.5E  not saved\n",i+1,eigval[i]);
	}
	for(i=start_mode; i<end_mode; i++) {
		// Max amplitude
		vectorrmsd = 0.0;maximumamplitude=-10000;
		for(j=0; j<size; j+=3) {
			current_ampl =  eigvect[i*size+j+0]* eigvect[i*size+j+0] +
			eigvect[i*size+j+1]* eigvect[i*size+j+1] +
			eigvect[i*size+j+2]* eigvect[i*size+j+2];
			vectorrmsd+=current_ampl;
			if (maximumamplitude < current_ampl )
				maximumamplitude = current_ampl;
		}
		// printf(" %11.5f\n",vectorrmsd);
		maximumamplitude = sqrt(maximumamplitude);
		norm=vectorrmsd;
		vectorrmsd = sqrt(vectorrmsd/((double) (size/3.0)));
		//getchar();

		// mode
		fprintf(out,"****\n%d     %11.5f\n",i-start_mode+1,eigval[i]);
		for(j=0; j<size; j++) {
			fprintf(out,"%11.5f",eigvect[i*size+j]/sqrt(norm));
			if (fmod((float)j+1,(float)7)==0)  fprintf(out,"\n");
		}
		if (fmod((float)j,(float)7)!=0)  fprintf(out,"\n");

		// if (i+start_mode<20) printf("nmac>   Mode %4d %15.5E \n",i+1,eigval[i]);

	}
	//printf("nmac>\nnmac>     Max amplitude %f Vector rmsd %f norm %f\nnmac>\n", vectorrmsd,maximumamplitude,norm);
	fclose(out);
}
