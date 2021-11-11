#ifndef LIBNMA_MOVE_H_
#define LIBNMA_MOVE_H_

#include "nma.h"

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix);

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain): N-CA-C-model (for CA you need 3 atoms!)
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
// "iters" array of residue-level (pdbIter *)'s
void move_dihedralMCAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix, pdbIter **iters);

// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain): 3BB2R & Full-Atom models
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
// If fix=NULL, no fixed dihedrals (equal to previous-version)
void move_dihedralMFAx(Macromolecule *mol, double *uu, tri *props, double step, int type, int model, bool *fix=NULL, bool *addrot=NULL);

// It applies a Cartesian Coordinates Normal Mode with "step" amplitude to an
// input structure "moli" and outputs a different "molf" structure.
// "uu"-array has 3N elements, where N= number of atoms. (Valid for all atomic models)
void move_cart(Macromolecule *mol, double *uu, double step);

// It applies many consecutive CCS Normal Modes with their amplitudes stored
// in "amp"-array to an input structure "mol".
// "uu"-array has 3N elements, where N= number of atoms. (Valid for all atomic models)
void move_cart(Macromolecule *mol, double *uu, double *amp, int nmodes);

//// Moves atoms given an Internal Coordinates Normal Mode and a Step (CA atomic model)
//// Using V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
//void move_vwMCAx(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1 );
//
//// Multi-thread routine to move atoms given an Internal Coordinates Normal Mode and a Step (CA atomic model)
//// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
//void *move_vwMCAx_thread(void *threadarg);
//
//// Moves atoms given an Internal Coordinates Normal Mode and a Step (CA atomic model)
//// Using V/W-arrays and multi-threaded parallel schedule (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
//void move_vwMCAx_par(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1, int nthreads);

// Moves atoms given an Internal Coordinates Normal Mode and a Step (C5 & HA atomic models)
// Using V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void move_VW(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1, int model );

// Multi-thread routine to move atoms given an Internal Coordinates Normal Mode and a Step
// It uses V/W-arrays (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void *move_VWthread(void *threadarg);

// Moves atoms given an Internal Coordinates Normal Mode and a Step in Parallel
// Using V/W-arrays and multi-threaded parallel schedule (Memory efficient version of Jacobian) (Multi-Chain/Prot/DNA/RNA/SMOL)
void move_VWpar(Macromolecule *mol, double step, double *evec, int index, int size, double ***V, double ***W, int **body1, int nthreads, int model);

#endif /*LIBNMA_MOVE_H_*/
