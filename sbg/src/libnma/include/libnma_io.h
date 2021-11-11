#ifndef LIBNMA_IO_H_
#define LIBNMA_IO_H_

#include "nma.h"

// Memory allocator
void *allocate(size_t size, char *message);
// Memory re-allocator
void *reallocate(void *pointer, size_t size, char *message);

// Saves & Normalizes eigenpairs according to "ptraj" format (normalizes by default)
void save_ptraj_modes (char *, int, int, int, floating *, floating *, bool norm_switch=true);
// Saves & Normalizes eigenpairs with "ptraj" format (normalizes by default) SINGLE PRECISION
void save_ptraj_modes (char * name, int size, int start_mode, int end_mode, float *eigval, float *eigvect,bool norm_switch);

// Reverse eigenvalues/eigenvectors ordering
// eval --> eigenvalues array
// evec --> eigenvector matrix
// size --> number of coordinates
// nevs --> number of eigenpairs
void reverse_ptraj(double *eval,double *evec,int size,int nevs);

// Mass Weight/Un-weight eigenvectors.
// 	evec --> eigenvector matrix
// 	mass --> masses array (length: number of atoms)
// 	size --> number of coordinates
// 	nevs --> number of eigenpairs
// 	weight --> true: weights (divide by sqrt(mi)), false: un-weights (multiply by sqrt(mi))
void weight_ptraj(double *evec,double *mass,int size,int nevs, bool weight=true);

// Reads an entire ptraj file allocating memory
int read_ptraj( char * file, double **evect, double **evals, int *numatoms, int *n_vectors, int *n_components  );

// USED in "nmaview" (18/3/2010)
// Reads an entire ptraj file allocating memory
// (This one returns a "trs" structure array)
trs *read_ptraj( char * file, int numatoms, int *n_vectors  );

// USED in "nmafit" (18/3/2010)
// Secondary Structure Input (allocates table)
// In the input file, the residue index now begins with 0.
void read_ss(char *file_ss, char **table, int *p_num=NULL);

// USED in "nmafit" (18/3/2010)
// Reads input TS functions
void read_TSfunc(char *file, TSfunc **funcs, int *n_func);

// Read Force constants file (Kfile) allocating memory
// Warning: if "coord" not provided, distances ".d" will not be updated!
void read_Kfile(twid **p_ipas,int *p_nipa,char *file, float *coord=NULL);

// Reads a file with the residue indices whose internal coordinates will be fixed.
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int fixFrag( char *file, Macromolecule *mol, tri *props,bool *fixed);

// Reads a file with the fixed internal coordinates indices.
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int read_fixIC( char *file, Macromolecule *mol, tri *props,bool *fixed);

// Reads a file with the residue indices whose internal coordinates will be fixed.
// FORMAT: "[residue-index] [0/1](phi) [0/1](chi) [0/1](psi)"
//          (0=fix / 1=mobile) (phi, chi and psi allways required!)
// "fixed" array must be already allocated for full- "size" elements!
// Returns the number of mobile variables
int read_fixDH( char *file, Macromolecule *mol, tri *props,bool *fixed, int type, int model, bool *addrot=NULL);

// Reads the "fix" file and outputs a "number of residues sized" array of "ints" with the indices
// of those "rigid bodies" as a function of the residues. It searchs for chunks of continuously fixed
// sequences of residues and assigns correlative indices to the corresponding "rigid bodies".
//
// File FORMAT (protein):      "[residue-index] [phi] [chi] [psi]"
// File FORMAT (nucleic acid): "[fragment-index] [alpha] [beta] [gamma] [chi] [epsilon] [zeta]"
//   [angle]= (0=fix / 1=mobile) All always required! (according to fragment ID)
// Memory for the "fixed" array must be already allocated (for the total number of residues)
// Returns the number of rigid bodies
int read_fixRes( char *file, Macromolecule *mol, int *fixed);

//// Parses a Clustal's ".aln" file to determine residue correspondence between two sequences
//// type = 1, 2 or 3 (match *, *: or *:. ,respectively)
//void read_aln(char *file, bool **p_mask1, int *nres1, bool **p_mask2, int *nres2, int *nres, int type=3);

// Writes a fixation-file with the residue indices and which internal coordinates are fixed.
// FORMAT (protein):      "[residue-index] [phi] [chi] [psi]"
// FORMAT (nucleic acid): "[fragment-index] [alpha] [beta] [gamma] [chi] [epsilon] [zeta]"
//   [angle]= (0=fix / 1=mobile) All always required! (according to fragment Mol-type: protein/nucleic-acid)
// "fixed" array must be already made! If fixed=NULL, then a fully mobile file will be written.
void write_fixDH( char *file, Macromolecule *mol, tri *props,bool *fixed, int type, int model);

// It writes an IC coordinates fixation mask (fixIC-format)
void write_fixIC(char *text,bool *fixed,int size);

// Write Force constants file (Kfile)
void write_Kfile(twid *ipas,int nipa,char *file);

// Write Secondary Structure file
void write_SSfile(char *ss,int nss,char *file);

// Shows a triangular packed matrix (standard output)
void show_matrix(floating *matrix, int size, char *name="Matrix");
// Shows a triangular packed matrix (standard output) with single precision
void show_matrix(float *matrix, int size, char *name);

// Saving a triangular packed matrix (into file)
// type --> 0= default, 1= save indices: [i] [j] [covariance], 2= Gnuplot matrix format (row and col wise)
// save_sym --> true, save full matrix (apply symmetry to obtain an squared symmetric matrix)
void save_matrix(floating *matrix, int size, char *name, int type=0, bool save_sym=false);
// Saving a triangular packed matrix (into BINARY file)
void save_matrixB(floating *matrix, int size, char *name);
// Save array into a plain-text file
void save_array(double *array,int size,char *file);
// Save array into a plain-text file
void save_array(float *array,int size,char *file);

// Saving a rectangular matrix (into file)
void save_matrix_rec(floating *matrix, int vec, int size, char *name);
// Saving a rectangular matrix (into file) with single precision
void save_matrix_rec(float *matrix, int vec, int size, char *name);
// Saving a rectangular matrix (into BINARY file)
void save_matrix_recB(floating *matrix, int vec, int size, char *name);

// Reading a triangular packed matrix (into memory)
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrix(floating **p_matrix, int *p_size, char *name);
// Reading a triangular packed matrix (into memory) from BINARY file
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrixB(floating **p_matrix, int *p_size, char *name);
// Reading a triangular packed matrix (into memory) from BINARY file (Reverse-row-wise)
// Allocates memory (if *p_matrix == NULL) and returns the matrix size (in "*p_size")
void read_matrixB2(floating **p_matrix, int *p_size, char *name);

// *******************************************************************************
// PABLO's OLD IO-FUNCTIONS (DEFPROT-back compatibility)
// *******************************************************************************
// save modes in separate ascii files
void save_ascii_modes (int, int , int , float *, float *);
// save modes in a ptraj ascii file
void save_ptraj_modes (int, int , int , float *, float *);

#endif /*LIBNMA_IO_H_*/


// USED in "nmamon" (18/3/2010)
// Reads input functions
//void read_funcSS(char *file, float **functions, int *n_func, int *n_parm);
// Saves contacts file
//void save_contacts(char *name, int *offset, double *parm, float *cont_matrix, double *dist_matrix, int n);
// Reads a list from a text-file (First string per line considered only!)
//void read_list(char *file, char ***list, int *n_elements);
// Saves parameters file (for functional parameters fitting)
//void save_parms(int n_parm, char *name, int *offset, double *parm, double *eps, int n);
// Reads a Hessian Matrix from a "Force_Constants_File"
// (Allocates memory if *p_hessian == NULL)
//void read_hessian(char *file, myfloat **p_hessian, int *n);
// Reads the Contacts Matrix from a ".k" file.
// Contact list = (i,j,k,dist)
//void read_contacts(char *file, double **p_contacts, int *n_contacts);
// Reads the hydrogen bonds according to a HBplus file.
//void read_hbplus( char *file_hb2, int **list, int *n_bonds );
// Reads Simplices from file (.sx) (Simplices with ordered atoms and t-numbers) (from: simp2k.pl)
//void read_simplices(char *file, int **simplices, int *n_simplices);
// save & show the eigenvalues
//void save_eigenvalues(char *f_name, int start_mode, int end_mode, myfloat *eigval);
