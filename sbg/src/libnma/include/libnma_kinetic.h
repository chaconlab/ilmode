#ifndef LIBNMA_KINETIC_H_
#define LIBNMA_KINETIC_H_

#include "nma.h"

// Calculates: elMY, emMY, el, em vectors. Needed for Fast Kinetic Energy computation.
// Table I. and (ec. 44) from Braun et al. (1984)
void calcKine(double mb1, double *Y1, double mb3, double *Y3, double **I1, double **I3, double *Phi, double *Psi, double *elMY, double *emMY, double *el, double *em);

// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain)
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naive( Macromolecule *mol, floating **p_mass_matrix, trd *der, int size);
// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain)
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naive_double( Macromolecule *mol, double **p_mass_matrix, trd *der, int size);
// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain) (for all models and types)
// The same as kinetic_naive(), but computing derivatives "on the fly" via V/W arrays
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naiveVW( Macromolecule *mol, float *coord, floating **p_mass_matrix,double ***V,double ***W,bool **body1, int size);
// Computes the Kinetic Energy Matrix (H) (Single- and Multi-Chain) (for all models and types)
// The same as kinetic_naive(), but computing derivatives "on the fly" via V/W arrays
// Noguti & Go (1983) pp. 3283-8 (ec. 6)
// Valid for any CG-model, just introduce the appropriate "der"
// Upper Triangular matrix packing:  i + j*(j+1)/2  -->  for 0<=i<=j<size (begining with 0)
// (if *p_mass_matrix == NULL --> automatic memory allocation!)
void kinetic_naiveVW_double( Macromolecule *mol, float *coord, double **p_mass_matrix,double ***V,double ***W,bool **body1, int size);

// Computes the Kinetic Energy Matrix (H) (Multi-Chain & FULL-ATOM & Protein/RNA/DNA/SMOL)
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMFAx( Macromolecule *mol, float *coord, tri *props, int size, floating **p_Hmatrix, int type, int model, bool *fix, bool *addrot=NULL);

// Computes the Kinetic Energy Matrix (H) (Multi-Chain & CA-model)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMCAx( Macromolecule *mol, float *coord, tri *props, int size, floating **p_Hmatrix, int model, bool *fix);

// Computes the Kinetic Energy Matrix (H) (Multi-Chain & CA-model) (Huge-systems)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
// Allocates Kinetic-Energy Matrix if p_Hmatrix=NULL (triangular packing)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "file" --> Hessian matrix file name
// "justfile" --> true, it only outputs the disk file (no memory allocation)
// Noguti & Go (1983) pp. 3283-8 (ec. 26)
// Warning: Coordinates should be already centered! (Center of Mass)
void kineticMCAxHD( Macromolecule *mol, float *coord, tri *props, int size, floating **p_Hmatrix, int model, bool *fix, char *file=NULL, bool justfile=false);

#endif /*LIBNMA_KINETIC_H_*/
