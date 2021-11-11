#ifndef LIBNMA_HESSIAN_H_
#define LIBNMA_HESSIAN_H_

#include "nma.h"

// Computing Tij element "on the fly"
// (T 6x6 matrix must be already allocated!)
void calcTij( double **T, int *index, int n, twid *decint, float *coord );

// Computing Tij element "on the fly" from "sorted ipas" array.
// (ec.21) Noguti & Go 1983 pp.3685-90
void calcTij_new(double **T, twid **sipas, int *index, float *coord, int nipa);

// Cartesian Coordinate Space (CCS) Hessian Building
// (Mass-weighted CCS if mass!=NULL; i.e. including kinetic energy matrix into the hessian)
// Triangular matrix storage and Hessian memory allocation (if *p_hess_matrix==NULL)
void hessianC(double **p_hess_matrix,double *mass,float *coord,double *dist_matrix,double *cont_matrix, int size);

// Dihedral Angle Space (DAS) Hessian Building (Multi-chain)
// Adds torsional springs to a previously built hessian matrix (triangular packing)
// (See second term of ec.5 from: Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471)
void hessianMDHx(Macromolecule *mol,tri *props,floating *hess,double omega,int size,bool *fix,bool *addrot=NULL);

// Dihedral Angle Space (DAS) Hessian Building (Multi-chain/Protein/DNA/RNA/SMOL)
// Adds torsional springs to a previously built hessian matrix from file (triangular packing)
// (See second term of ec.5 from: Lu, Poon and Ma. J.Chem. Theory Comput. 2006, 2, 464-471)
void hessianMDHxHD(Macromolecule *mol,tri *props,char *file,double omega,bool *fix,bool *addrot);

// Fast Hessian Matrix computation (Multi-Chain & Full-Atom / 3BB2R)
// (Very Efficient memory allocation)
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
// Noguti & Go (1983) pp. 3685-90
void hessianMFAx(Macromolecule *mol, twid *decint,int nipa,int size,float *coord,floating **p_hess_matrix,tri *props, int *unat, int type, int model, bool *fix, bool *addrot=NULL);

// Fast Hessian Matrix computation // (Multi-Chain & CA-model)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// (Very Efficient memory allocation)
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// Noguti & Go (1983) pp. 3685-90
// (3/2/2010)
void hessianMCAx(Macromolecule *mol,twid *decint,int nipa,int size,float *coord,float *coordCA,floating **p_hess_matrix,tri *props, int *unat, bool *fix);

// Fast Hessian Matrix computation: Multi-Chain & fix & CA-model. (4/4/2012)
// Minimal memory requirements. (Efficient "Uab" computation, i.e. not using "unipa")
// Input: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// [Ref.] Noguti & Go (1983) pp. 3685-90
void hessianMCAxHD(Macromolecule *mol,twid *decint,int nipa,long int size,float *coord,float *coordCA,
		floating **p_hess_matrix,tri *props, int *unat, bool *fix, char *file=NULL, bool justfile=false);

// Fast Hessian Matrix computation: Multi-Chain & fix & CA-model. (7/11/2012)
// Minimal memory requirements. (Efficient "Uab" computation, i.e. not using "unipa")
// Input: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Hessian-Matrix if p_hess_matrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "nthreads" --> Number of threads to be used.
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// [Ref.] Noguti & Go (1983) pp. 3685-90
void hessianMCAxHD_par(Macromolecule *mol,twid *decint,int nipa,long int size,float *coord,float *coordCA,
		floating **p_hess_matrix,tri *props, int *unat, bool *fix, int nthreads, char *file=NULL, bool justfile=false);

// Multi-threaded routine to compute a "Hessian-Box"
void *hessianMCAxHD_thread(void *threadarg);

// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
void hess_naive(twid *decint,trd *der,int nipa,int size,float *coord,double **p_hess_matrix);
// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// ********* void hess_naiveVW(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix);
void hess_naiveVW(float *coord,twid *decint,double ***V,double ***W,int **body1,int nipa,int size,double **p_hess_matrix,int model);
// Hessian matrix computation by Kovack's "naive" method O(n^4) (Valid for all "type")
// The same as hess_naive(), but computing derivatives "on the fly" via V/W arrays
// Avoids zero derivatives computation (Mon, 19/02/2009) (zero-tolerance not needed with V/W!)
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// ******** void hess_naiveVW0(float *coord,twid *decint,double ***V,double ***W,bool **body1,int nipa,int size,double **p_hess_matrix);
void hess_naiveVW0(float *coord,twid *decint,double ***V,double ***W,int **body1,int nipa,int size,double **p_hess_matrix,int model,int row=-1);
// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// Avoids zero derivatives computation (Mon, 1/12/2009) (warning zero-tolerance=0.0000000001)
void hess_naive0(twid *decint,trd *der,int nipa,int size,float *coord,floating **p_hess_matrix);
// Hessian matrix computation by Kovack's "naive" method O(n^4)
// The same as hess_nipa_old(), but with Triangular matrix packing
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_hess_matrix == NULL --> automatic memory allocation!)
// Avoids zero derivatives computation (Mon, 1/12/2009) (warning zero-tolerance=0.0000000001)
void hess_naive0_double(twid *decint,trd *der,int nipa,int size,float *coord,double **p_hess_matrix);

#endif /*LIBNMA_HESSIAN_H_*/
