#ifndef LIBNMA_DYDQ_H_
#define LIBNMA_DYDQ_H_

#include "nma.h"

// Computes derivatives for an interaction pair of atoms given "v"/"w" elements
double derVW(double **v, double **w,float *ri,float *rj);

// K-matrix coeficients computation
// Ec's 24 & 25, form Noguti & Go, 1983
//     computes the 4 vectors (from Go's paper) in terms of which the derivatives
//     of the cartesian coordinates w/r/t the dihedral coordinates are expressed as:
//     Dr_b/Dq_j = v[0] + w[0] x r_b  (for body 1)
//     Dr_c/Dq_j = v[1] + w[1] x r_c  (for body 2)
// Computes "v" and "w" (Multi-Chain) given Phi & Psi, etc...
void calcoefM( double mb1, double mb2, double phi[3], double psi[3], double Y1[3], double I1[3] [3], double I2[3] [3],
	      double J[3] [3], double v[2] [3], double w[2] [3] );

// Compute derivatives of cartesian coordinates w/r/t dihedrals (Multi-Chain & FULL-ATOM)
// NEEDS: -Macromolecule (pseudo-atom model), -Single row PDB coordinates,
// -Properties structure array (props[]), -Derivatives matrix **reference
// INTERNAL COORDS. MODEL: type = 0 --> phi,psi,... type = 1 & 2 --> phi,chi,psi...
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
// "ictype" & "icmodel" --> Derivatives IC-model and type
// "addrot" --> bool array. "true", if 3 additional rotations should be added due to fixing.
// WARNING: "size" should be set according to the IC-model (not the atomic model)
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
void dydqMFAx(Macromolecule *molr,float *coord, tri *props, trd **p_der, int type, int model, int size, bool *fix=NULL, bool *addrot=NULL);

// Compute derivatives of cartesian coordinates w/r/t dihedrals (Multi-Chain)
// NEEDS: -Macromolecule (N,CA,CA model) (3-atoms model),
//        -Single row (N,CA,C model) coordinates,
//        -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void dydqMCA3x(Macromolecule *mol,float *coord, trd **p_der, int size, bool *fix);

// Compute derivatives of cartesian coordinates w/r/t dihedrals (Multi-Chain + Fixation)
// NEEDS: -Macromolecule (first-NH + CA-atom + last-CO, per segment, model),
//        -Single row (N,CA,C)-PDB coordinates,
//        -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
void dydqMCAx(Macromolecule *mol,float *coord, trd **p_der, int size, bool *fix);
// Compute V/W-arrays - Memory Efficient (see ec. 40 Braun et al.) (Multi-Chain + Fixation)
// NEEDS: -Macromolecule (N,CA,CA model) (3-atoms model),
//        -Single row (N,CA,C)-PDB coordinates,
//        -Arrays V/W by reference.
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: CA-only (phi,psi)
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void vwMCA3x(Macromolecule *mol,float *coord, double ****p_V, double ****p_W, int ***p_body1, int size, bool *fix, int model=3);
// Compute V/W-arrays (see ec. 40 Braun et al.) (Multi-Chain + Fixation)
// (Multi-Chain & FULL-ATOM & Protein/RNA/DNA)
// NEEDS: Macromolecule, Single row PDB coordinates,
// -Properties structure array (props[]), -Derivatives matrix **reference
// WARNING: "coord" must be centered, ie. CoM = (0,0,0)
// INTERNAL COORDS. MODEL: type = 0 --> phi,psi,... type = 1 & 2 --> phi,chi,psi...
// "fix" --> bool array (number of Internal Coords sized) with the mask of mobile IC's.
void vwMFAx(Macromolecule *molr,float *coord, tri *props, double ****p_V, double ****p_W, int ***p_body1, int type, int model, int size, bool *fix, bool *addrot=NULL);

// Dihedral Eigenvector ---> Cartesian Eigenvector (*mass = NULL)
// Dihedral Eigenvector ---> Mass-Weighted Cartesian Eigenvector (*mass != mass-array)
// (evec==NULL --> allocates memory)
// hess_matrix = dihedral modes
void di2cart(floating *hess_matrix,trd *der,floating *evec,int size,int num_atoms,int modes_saved, double *mass=NULL);

// Translate from IC-modes to Mass-weighted/unweighted Cartesian modes
// If *evec != NULL --> Unweighted Cartesian modes modes output in "*evec"
void di2cartVW(floating *hess_matrix,float *coord,double ***V,double ***W,int **body1,floating *evec,int size,int num_atoms,int modes_saved, int model);

// Multi-thread routine to convert CC modes into IC modes in parallel
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void *di2cartVW_thread(void *threadarg);

// Parallel routine to convert CC modes into IC modes
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void di2cartVW_par(floating *hess_matrix,float *coord,double ***V,double ***W,int **body1,floating *evec,int size,int num_atoms,int modes_saved, int model, int nthreads);

// Translates from Cartesian normal modes to Mass-Weighted ones
// (mass <-> mass square root)
void cart2wcart(floating *evec,int num_atoms,int modes_saved,double *mass);

#endif /*LIBNMA_DYDQ_H_*/
