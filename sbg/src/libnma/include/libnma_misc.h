#ifndef LIBNMA_IC_H_
#define LIBNMA_IC_H_

#include "nma.h"
extern CRandomMersenne *rg; // Mersenne Twister global object

// Creates contacts list (IPA) directly from a Macromolecule ("mol")
void make_ipas_new(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff);

// Creates contacts list (IPA) directly from a Macromolecule ("mol")
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
void make_ipas_new0(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff);

// Multiplies by "factor" every ipa force constant ("C") if both atoms belong to different molecules.
void modify_intermolec_ipas(Macromolecule *mol, twid *decint, int nipa, float factor);

// Creates contacts list (IPA) directly from a Macromolecule ("mol")
// Contacts according to a predefined set of functions (depending on Topology and Secondary Structure)
void make_ipasTS(Macromolecule *mol, twid **p_decint, int *p_nipa, float cutoff, TSfunc *funcs, int nfuncs, char *ss);

// "MIXED" contacting method. (WARNING: only fully valid for single chain and CA models)
// Creates contacts list (IPA) and assigns force constants given a Macromolecule ("mol")
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
void make_ipas_mix0(Macromolecule *mol, twid **p_decint, int *p_nipa);

// "MIXED" contacting method. (WARNING: only fully valid for single chain and CA models)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_mix(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa);

// Yang,...,Jernigan's "pfENM". PNAS (2009)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_pfENM(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa);

// Kovaks's "inverse exponential" (Rueda et al.)
// Fills a contact_matrix with force constants given a Macromolecule ("mol")
void cont_matrix_Kovacs(Macromolecule *mol, double *cont_matrix, int num_atoms, int *p_nipa);

// Computes NIPAs distances
void dist_nipa(twid *decint,float *coord, int nipa);

// Computes some residue properties (Full-Atom & Multi-Chain & RNA)
// Creates an array to relate atoms and the rigid units which atoms belongs to (unat).
void properMFA(Macromolecule *molr,tri **p_props,int **p_unat,int type,int model);

// Computes some residue properties (CA-model & Multi-Chain)
// "p_props" stores the properities.
// "p_unat" stores the unit which an atom belongs to.
void properCA(Macromolecule *molr,tri **p_props,int **p_unat);

// Creates two auxiliar arrays with segment properties (needed due to fixing):
//   addrot[#seg] --> true, if 3 additional rotations should be added due to fixing.
//   effseg[#seg] --> <int>, with the number of "effective segment" for #seg.
// (Allocates memory itself, if it's needed)
int seg_props(Macromolecule *mol, tri *props, bool *fix, int model, int type, bool **p_addrot, int **p_effseg);

// Creates a random mask of fixed internal coordinates according to a probability "prob"
// (Given a random number [0:1), an IC will be selected if "random>prob")
// "fix" array must be already allocated for full- "size" elements!
int fixRand(bool *fix,int size,float prob=0.5);

// Creates a random mask of fixed DiHedral coordinates according to a probability "prob"
// (Inter-chain Rotational/Translational coordinates will be always kept mobile)
// (Given a random number [0:1), an IC will be selected if "random>prob")
// "fix" array must be already allocated for full- "size" elements!
int fixRandDH(Macromolecule *mol,tri *props,bool *fix,int type,int model,float prob=0.5);

// Creates a mask of fixed DiHedral coordinates according to Secondary Structure
// (Inter-chain Rotational/Translational coordinates will be always kept mobile)
// "fix" array must be already allocated for full- "size" elements!
// "ss_table" is a num_res sized array with the molecule single char SS identifiers.
// "fix_ss" is an string with the residue char SS identifiers which will be fixed.
int fixSS(Macromolecule *mol,tri *props,bool *fix,char *ss_table,char *fix_ss,int type,int model);

// Extracts Dihedral-Angle components directly from the Hessian Matrix raw modes
// Allocates the mode array, only for the first time.
//void dihedral_comps(Macromolecule *mol,double *evec,int nev,int size,tri *props,double **pp_mode);
void dihedral_comps(double *evec,int nev,int size,double **pp_mode);

// Change normal modes (internal coordinates) from "i" model into "f" model
//void changemodelIC(Macromolecule *moli, Macromolecule *molf, floating *eveci, floating *evecf, int sizei, int sizef, int model, int type, int modelf, int typef);
void changemodelIC(Macromolecule *moli, tri *propsi, tri *propsf, floating *eveci, floating *evecf, bool *fixi, bool *fixf, int nmodes, int sizei, int sizef,
		int model, int type, int modelf, int typef, bool inverse=false);
// Change fix-arrays (internal coordinates) from a "i" model into a "f" model.
// WARNING: "i" is ALLWAYS the finer-grained model, "f" is ALLWAYS the coarser-grained model.
//          (no memory allocation is performed)
// inverse = false (default): Matching "f" DoFs will be set to "i"'s ones (all "f" DoFs will be filled)
// inverse = true: Matching "i" DoFs will be set to "f"'s ones, (not-matching "i" DoFs will be = zero)
void changefixIC(Macromolecule *moli, tri *propsi, tri *propsf, bool *fixi, bool *fixf,
		int model, int type, int modelf, int typef, bool inverse=false);

// Pablo's "inverse exponential function"
double inv_exp(double k, double x, double x0, double power);

// Normalices eigenvectors
int norm_evec(double *evec,int nevec,int size);

// Selects the "n_chain"th Chain from a Macromolecule
// *****   Watch out with NUCLEIC ACIDS and SMOL's ******
Macromolecule *select_chain(Macromolecule *in, int n_chain);

// Selects the "n_seg"th Segment from a Macromolecule
// *****   Watch out with NUCLEIC ACIDS and SMOL's ******
Macromolecule *select_segment(Macromolecule *in, int n_seg);

// Number of Degrees of Freedom in a Macromolecule "mol" (given "props")
int num_dofs(Macromolecule *mol,tri *props, int *size_dh=NULL, int *n_chain=NULL, int *n_seg=NULL);

// Fills a vector with random values from "low" to "high" [low:high)
void rand_vector(double *v, int size, double low, double high);

twid **sort_ipas(twid *ipas, int *pnipas, int *unat, int num_units); // Quick-sort for ipas
void quicksort0(twid **arr, int low, int high, int size); // Quick-sort recursive routine
void quicksort(twid **arr, int low, int high, int size); // Quick-sort recursive routine (ready for big systems)
void quicksort2(twid **arr, long int low, long int high, int size); // Quick-sort recursive routine

//
void swapSPmatrix(char *file, long int mem);

// Computing the Radius of Gyration (Rg) of a Macromolecule.
// See: Seeliger and de Groot. PLOS (2010).
float radius_gyration(Macromolecule *mol);

// Length of one vector
double length(double *v1);

// Dot product between two vectors
double dotp(double *v1, double *v2);

#endif /*LIBNMA_IC_H_*/
