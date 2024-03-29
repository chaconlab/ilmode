//============================================================================
// Name        : ilmode.cpp
// Author      : Mon
// Version     : 1.4
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <random>
#include <chrono>
#include <limits>

#include "cmdl/CmdLine.h" // Command line parser
#include "libnma/include/libnma_io.h" // Mon's NMA library (Dihedral included)
#include "libnma/include/libnma_cg.h" // Mon's NMA Coarse-Graining library
#include "libnma/include/libnma_matrix.h" // Some vector and matrix algebra routines
#include "libpdb/include/Rotamer.h"
#include "libnma/include/libnma_misc.h" // Mon's NMA Internal Coordinates related library

char version[]="1.5"; // version code
char prog[]="ilmode"; // program name


// Input data simple structure
typedef struct
{
	char fname[FILE_NAME]; // File name
	int start;			   // Start residue number of the loop
	int end;			   // End residue number of the loop
	char chain;			   // Chain ID (A,B,C,D...  )
} inputdata;

#define ptr_check(p) { \
		if ((p) == NULL) { \
			fprintf(stderr, "%s> Memory allocation failure. Exiting.\n", prog); \
			exit(1); \
		} \
}

#define SWAPPING(_a,_b,_type) \
		{\
	_type _tmp;\
	\
	_tmp = (_a);\
	(_a) = (_b);\
	(_b) = _tmp;\
		}


/*==============================================================================================*/


// Input variables
char file_pdb[FILE_NAME]; // Initial PDB
char file_pdb2[FILE_NAME]; // Target PDB (final)
char file_loops[FILE_NAME]; // Loops Multi-PDB
char file_align[FILE_NAME];	// Movie Multi-PDB structure file name

char file_movie[FILE_NAME];	// Movie Multi-PDB structure file name
char file_final[FILE_NAME];	// Morphed PDB structure file name
char file_log[FILE_NAME]; // Log file
char name[FILE_NAME];
char base[FILE_NAME];
char text[FILE_NAME*10];
char dummy[FILE_NAME*10];
char saved_files[2000]; // string to buffer screen output until the program end.

// INPUT PARAMETERS (PRE-DEFINED)
//==============================================================================================
unsigned int seed; // Random number generator seed for Mersenne Twister
int verb = 0; // Verbose level (0= none, 1= low, 2= high, ...)
int saved_files_len; // store the "saved_files" array current length
int nevec = -1; // number of eigenvectors to be computed (set by parser)
int imod = 0; // Index of mode to be shown (0,1,...,N-1)
int strategy = 0; // Sampling strategy (0=raw sampling for selected mode "i", 1= hyper-square sampling, etc...
int nsamples = 100;   // Number of samples for selected strategy
int max_loops_save = 0;   // Number of samples for selected strategy

int nsteps   = 1000;   // Number of samples for selected strategy
float rmsd_conv = 1e-4; // RMSD convergence threshold
float nevec_fact = -1; // % of eigenvectors to be computed (set by parser)
bool debug = false; // Debug mode
bool parse_verb = false;
bool saveformat_switch = false; // = true --> save formated input PDB
bool savemodel_switch = false; // = true --> save current CG-model PDB
bool ss_switch = false; // the NMA force constants will be set according to SS rules
bool func_switch = false; // input file with function coefficients (Topology and SS)
bool linear_switch = false; // =true --> linear motion
bool morph_switch = false; // =true --> Enable morphing stuff
bool mr_switch = false; // =true --> Enable Monte-Carlo stuff (random mode amplitudes)
bool loops_switch = false; // =true --> Enable multiple loops processing
bool anchors_present = true; // =true --> Anchors are present in loops Multi-PDB, otherwise just mobile residues present
int traji=0; // traj count;

// INPUT PARAMETERS (CUSTOMIZABLE by parser)
//==============================================================================================
int model = 0; // Model Coarse-Graining
int potential = 0; // Potential method
int contact = 0; // Contact method
int scoring = 1; // Scoring method
int type = 0; // Internal Coordinates Coarse-Graining type
int start = -1; // N-terminal loop residue index, first mobile residue (PDB residue index)
int loop_end = -1; // C-terminal loop residue index, last mobile residue (PDB residue index)
int nflanks = 0; // Number of loop flanking residues (e.g. "2" means two flanking residues at Nt and other two at Ct, "0" means disabled)
bool norm_modes = false; // true = Normalizes modes (norm=1)
bool delHydrogens_switch = true; // Delete hydrogens
bool delHeteros_switch = false; // Delete heteroatoms
bool delWaters_switch = true; // Delete waters
bool already_seed=false; // = true if seed has been defined by user
// bool loop_startend = true; // = true if Start and End residue indices were parsed
bool skip_missingatoms = false; // Set true to skip missing atoms
bool ali_flanks = false; // Set true to enable flanks alignment
float cte_k0 = 1.0; // Force constant factor for Inverse Exponential
float cte_k02 = 1.0; // Force constant factor for Inverse Exponential
float cte_k1 = 1.0; // Force constant factor for Tirion's simple cutoff method
float x0 = 3.8; // Inflexion point of the Inverse Exponential
float power = 6.0; // Power term of the Iverse Exponential
float cutoff_k0 = 10; // Distance cutoff for the Inverse Exponential
float cutoff_k1 = 10; // Distance cutoff for the Tirion's simple cutoff method
float target_rmsd = 2.0; // Target RMSD for mode-following strategies
float delta_rmsd = 999999; // Delta RMSD in mode-following strategies (high value means disabled)
float flank_rmsd = 0; //
double maxang = 0.0; // Maximum angular increment to normalize modes [deg]
double dc; // Characteristic distance factor for Deformability computations (see: hardy and compute_def functions in libnma_def.cpp)
char chain = '*'; // Chain ID for initial PDB
char chain2 = '*'; // Chain ID for target PDB
int id_chain;
int id_chain2;
int check_model;

FILE *f;





extern "C" {

//SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
//$                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
//*  -- LAPACK driver routine (version 3.2) --
//CHARACTER          JOBVL, JOBVR
//INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
//DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
//$                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
//$                   VR( LDVR, * ), WORK( * )
void dggev_(char *jobvl, char *jobvr, int *sizex, double *hess_matrix,
		int *lda, double *mass_matrix, int *ldb, double *alphar,
		double *alphai, double *beta, double *vl, int *ldvl, double *vr,
		int *ldvr, double *work, int *lwork, int *info);

//void dsygvx_(int *ITYPE, char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, double *B, int *LDB, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z, int *LDZ, double *WORK, int *LWORK, int *IWORK, int *IFAIL, int *INFO);
}

// Compute the derivatives (der) of Cartesian coordinates (r) wrt dihedral angles (q) for loops:  der = dr/dq
// OUTPUT:
//  Returns the derivatives matrix (automatic memory allocation)
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  props --> Residue properties array of the whole macromolecule
//  ifr   --> Index of First loop Residue (internal numeration)
//  ilr   --> Index of Last loop Residue (internal numeration)
//  nla   --> Number of Loop Atoms (only mobile atoms)
//  size  --> Number of dihedral angles (the mobile variables)
//  model --> Coarse-Grained model
trd *drdqC5x(float *coord, tri *props, int ifr, int ilr, int nla, int size, int model);
inline void drdqC5x(float *coord, tri *props, int ifr, int ilr, int nla, int size, int model, trd *del);

// Get array of atomic masses from a Macromolecule iterator (num_atoms is for cross-checking purposes)
float *get_masses(pdbIter *iter, int num_atoms);

// Compute the Kinetic Energy matrix (masses matrix, M) for loops NMA
// OUTPUT:
//  Returns the Kinetic Energy matrix (automatic memory allocation)
// INPUT:
//  der    --> Derivatives
//  masses --> Array of masses, one mass per (pseudo)atom
//  props  --> Residue properties array of the whole macromolecule
//  ifr    --> Index of First loop Residue (internal numeration)
//  nla    --> Number of Loop Atoms (only mobile atoms)
//  size   --> Number of dihedral angles (the mobile variables)
//  nco    --> Number of Constraints
double *kineticC5x(trd *der, float *masses, tri *props, int ifr, int nla, int size, int nco);
inline void kineticC5x(trd *der, float *masses, tri *props, int ifr, int nla, int size, int nco, double *mass_matrix);

// Pablo's "inverse exponential function"
inline double Inv_exp(double k, double x, double x0, double power)
{
	return ( k / ( 1.0 + pow( x/x0, power ) ) );
}



// Normalices eigenvectors
inline int Norm_evec(double *evec,int nevec,int size)
{
	double norm = 0.0;
	int offset;
	// Normalizing eigenvectors
	for(int i=0; i<nevec; i++)
	{
		norm = 0.0;
		offset = i * size;
		for(int j=0; j<size; j++) norm += pow(evec[offset+j],2); // computes the norm
		norm = sqrt(norm);
		for(int j=0; j<size; j++) evec[offset+j] /= norm ; // this normalizes
	}
	return( 0 ); // normal exit status
}

// Compute the Hessian matrix (2nd derivatives of the potential energy, H) for loops NMA
// OUTPUT:
//  Returns the Heesian matrix (automatic memory allocation)
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  der    --> Derivatives
//  props  --> Residue properties array of the whole macromolecule
//  decint --> List of Interacting Pairs of Atoms (contacts list)
//  nipa   --> Number of Interacting Pairs of Atoms (number of contacts)
//  ifr    --> Index of First loop Residue (internal numeration)
//  nla    --> Number of Loop Atoms (only mobile atoms)
//  size   --> Number of dihedral angles (the mobile variables)
//  nco    --> Number of Constraints
double *hessianC5x(float *coord, trd *der, tri *props, twid *decint, int nipa, int ifr, int nla, int size, int nco);
inline void hessianC5x(float *coord, trd *der, tri *props, twid *decint, int nipa, int ifr, int nla, int size, int nco,double *rdr, double *hess_matrix );
double *hessianFast(float *coord, float *coordCA, trd *der, tri *props, twid *decint, int nipa, int ifr, int ilr, int size, int nco, int num_atoms);
inline void  hessianFast(float *coord, float *coordCA, trd *der, tri *props, twid *decint, int nipa,  int ifr, int ilr, int size, int nco, int num_atoms, double *hess_matrix);

// Computing Tij element "on the fly"
// (ec.21) Noguti & Go 1983 pp.3685-90
// (T 6x6 matrix must be already allocated!)
void calcTij( double **T, int *index, int num, twid *decint, float *coord );

int diag_dggev(double *eigval, double *eigvect, double *mass_matrix, double *hess_matrix, int size, int nco, int *neig);
inline int diag_dggev(double *eigval, double *eigvect, double *mass_matrix, double *hess_matrix, int size, int nco, int *neig, double *alphar, double *alphai, double *beta, double *vr,  double *vl,  double *work);



//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx(double *hess_matrix,double *mass_matrix,double *eigval, int size, int il, int iu);

// Scale eigenvectors so that the maximum value of the components is "maxang"
void scale_vectors(double *eigvect, int size, int neig, double maxang);

// Compute N-dimensional vector modulus (i.e. vector length)
double vector_modulus(double *v, int size);

// Compute N-dimensional vector modulus (i.e. vector length)
float vector_modulus(float *v, int size);

// Element-wise difference between vectors "v" and "w" (d = v - w)
void vector_diff(float *v, float *w, double *d, int size);

// Element-wise difference between dihedral angle [deg] arrays "v" and "w" (d = v - w) taking into account rotation
void dihedrals_diff(float *v, float *w, double *d, int size);

// Compute N-dimensional dot-product between two vectors (i.e. the modulus of the projection between both vectors)
double dotprod(double *v, double *w, int size);

// Compute the normalized (0,1) N-dimensional dot-product between two vectors (i.e. the cosine of the angle between both vectors)
double dotprodnorm(double *v, double *w, int size, double mv = 0.0, double mw = 0.0);

// Compute the normalized (0,1) N-dimensional dot-product between two vectors (i.e. the cosine of the angle between both vectors)
float dotprodnorm(float *v, float *w, int size, float mv = 0.0, float mw = 0.0);

// Convert IC (dihedral) modes into Cartesian modes
//	masses_sqrt --> Array of the square root of the atomic masses (one per loop atom)
inline double *ic2cart(double *eigvect, int nevec, trd *der, int size, int nla, float *masses_sqrt = NULL);
inline void *ic2cart(double *eigvect, int nevec, trd *der, int size, int nla, double *Aop, float *masses_sqrt = NULL);

void move_loop_linear_steps(char *file, Macromolecule *mol, double *cevec, int imod, int ifpa, int nla, int steps=10, float factor=1.0);

void move_loop_linear(pdbIter *iter, double *cevec, int imod, int ifpa, int nla);
void move_loop_linear_factor(pdbIter *iter, double *cevec, int imod, int ifpa, int nla, float factor);


inline void move_loop_dihedral(pdbIter *iter, int ifr, int ilr, tri *props, double *uu, int size, int model, float step);

float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla);

float rmsd_loop_residue(pdbIter *iter, pdbIter *iter2, int ifpr, int ifr2, int nla);



float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, int ifpa2);

float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta);



// Get dihedral angles (an array) from a continuous Macromolecule loop
// WARNING! Dihedral angles order: Psi_NtAnchor, Phi_1st, Psi_1st, Phi_2nd, Psi_2nd, ... , Phi_Nth, Psi_Nth, Phi_CtAnchor
//	iter --> Macromolecule iterator
//	ifr  --> Index (internal) of First Residue
//	nlr  --> Number of (mobile) loop residues
//	p_dihedrals --> Pointer to the array of dihedrals (=NULL forces automatic memory allocation)
void loop_dihedrals(pdbIter *iter, int ifr, int nlr, float **p_dihedrals);

// Compute just the loop dihedral angles considered in Loop-NMA
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  props --> Residue properties array of the whole macromolecule
//  ifr   --> Index of First loop Residue (internal numeration)
//  ilr   --> Index of Last loop Residue (internal numeration)
//  size  --> Number of dihedral angles (the mobile variables)
//  model --> Coarse-Grained (CG) model
//  type  --> Chi type (0= No-chi, 2= Phi,Chi,Psi)
// OUTPUT:
//	**p_dhs --> Pointer to dihedral angles array according to current CG model (automated memory allocation if *NULL)
void loop_dihedrals(float *coord, tri *props, int ifr, int ilr, int size, int model, int type, float **p_dhs);

// Compute the RMSD between two arrays
float vector_rmsd(float *array, float *array2, int size);

// Compute the RMSD for a given array of increments
float vector_rmsd(double *array, int size);

// Compute the Dihedral angles RMSD between two protein loops
float rmsd_dihedral_loop(pdbIter *iter, pdbIter *iter2, int ifr, int ifr2, int nlr, float *dhs=NULL, float *dhs2=NULL);

void delta_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta);

void delta_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta, int ifpa2);

float rmsd_flank(pdbIter *iter, pdbIter *iter2, int iffa, int ifa, int ila, int ilfa, int iffa2);

// Does some loop atom clash with its environment?
// 	iter,iter2 --> 2 differerent iterators to the same macromolecule
//	ifpa       --> Index of first (pseudo)atom
//	nla        --> Number of loop atoms
//	cut2       --> cutoff distance squared (2.0 A by default, i.e. 4.0)
bool clashed_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, float cut2);

// Does some loop atom clash with its environment?
// 	iter,iter2 --> 2 differerent iterators to the same macromolecule
//	ifpa       --> Index of first (pseudo)atom
//	nla        --> Number of loop atoms
//	cut2       --> cutoff distance squared (2.0 A by default, i.e. 4.0)
//  ifpa2      --> Index of first (pseudo)atom in 2nd PDB
bool clashed_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, float cut2, int ifpa2);

// Creates contacts list (IPA) for some loop (intra-loop + loop vs. environment) from two iterators pointing to the same Macromolecule.
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
inline void make_ipas_loop(pdbIter *iterA, pdbIter *iterB, int ifpa, int nla, float cutoff, twid **p_decint, int *p_nipa);
inline void make_ipas_loop(pdbIter *iterA, pdbIter *iterB, int ifpa, int nla, float cutoff, twid *decint, int *p_nipa);


// Get the array of atomic masses for some loop from one macromolecule iterator that points to a Macromolecule.
//	sqrt --> Set "true" to compute the square root of the masses (to mass-weight Cartesian modes), otherwise it just gets the masses
float *get_masses_loop(pdbIter *iter, int ifpa, int nla, bool sqrt = false);

// Measure loop-closure quality at the C-t end
//  iter   --> Iterator of first Macromolecule (e.g. reference)
//  iter2  --> Iterator of moving Macromolecule (e.g. target)
//	ilr    --> Index of last residue (internal numeration)
//  p_dist --> OUTPUT distance between C atom of last loop residue and the N atom of Ct-anchor
//  p_ang  --> OUTPUT bond angle between C atom of last loop residue and the N and C atoms of Ct-anchor
void anchor_drift(pdbIter *iter, pdbIter *iter2, tri *props, int ilr, float *p_dist, float *p_ang);

// Mode following routine
//	refmode0 --> Reference mode to drive motion (mode following)
//	mol --> Input macromolecule
//	model --> Atomic model
//	type --> ICs type (with/without chi)
//	props --> Properties structure
//	masses --> Atomic masses of the whole protein
//	masses_loop --> Atomic masses of the loop
//	ifa --> Index of first atom of the loop
//	ifr --> Index of first residue of the loop
//	ilr --> Index of the last residue of the loop
//	na --> Number of atoms of the loop
//	size --> Number of ICs of the loop
//	nco --> Number of constraints (typically 6)
//	cutoff --> Distance cutoff for elastic network definition
//	nsamples --> Number of samples (frames)
//	eigval --> IC eigenvalues, automatic memory allocation if NULL (OUTPUT)
//	eigvect --> IC eigenvectors, automatic memory allocation if NULL (OUTPUT)
//	rmsd_conv --> (OPTIONAL) Motion amplitude [RMSD] (motion will stop upon "nsamples" or if "rmsd_conv" RMSD is reached)
//	iterini --> (OPTIONAL) Reference pdb iterator to compute RMSD
//	file_movie --> (OPTIONAL) Output Multi-PDB file
//	delta_rmsd --> (OPTIONAL) Delta RMSD in mode-following strategies (high value means disabled)
//	p_fi --> (OPTIONAL) Frame index (model index in Multi-PDB file). Automatically updated inside.
//	chain --> (OPTIONAL) Chain-ID for output Multi-PDB
//	update --> (OPTIONAL) if "true", the "refmode" will be updated by current most overlapping mode
//	RETURN --> Number of non-null eigenpairs
inline int follow_mode(double *refmode0, Macromolecule *mol, int model, int type, tri *props, float *masses, float *masses_loop, int ifa, int ifr, int ilr,
		int na, int size, int nco, double maxang, float cutoff, int nsamples, int max_loops_save, double *eigval, double *eigvect, float rmsd_conv = 99999, pdbIter *iterini = NULL,
		char *file_movie = NULL, float delta_rmsd = 99999, int *p_fi = NULL, char chain = 'A', bool update=false);


inline int MC_mode(double *refmode0, Macromolecule *mol, int model, int type, tri *props, float *masses, float *masses_loop, int ifa, int ifr, int ilr,
		int na, int size, int nco, double maxang, float cutoff, int nsamples, double *eigval, double *eigvect, float rmsd_conv = 99999, pdbIter *iterini = NULL,
		char *file_movie = NULL, float delta_rmsd = 99999, int *p_fi = NULL, char chain = 'A', bool update=false);



// Loop NMA routine to just compute the eigenvectors/values given some macromolecular loop
//	mol --> Input macromolecule
//	model --> Atomic model
//	type --> ICs type (with/without chi)
//	props --> Properties structure
//	masses --> Atomic masses array
//	ifa --> Index of first atom of the loop
//	ifr --> Index of first residue of the loop
//	ilr --> Index of the last residue of the loop
//	na --> Number of atoms of the loop
//	size --> Number of ICs of the loop
//	nco --> Number of constraints (typically 6)
//	cutoff --> Distance cutoff for elastic network definition
//	eigval --> IC eigenvalues, automatic memory allocation if NULL (OUTPUT)
//	eigvect --> IC eigenvectors, automatic memory allocation if NULL (OUTPUT)
//	der --> Derivatives (dr/dq), always allocated here (free elsewhere) (OUTPUT)
//	RETURN --> Number of non-null eigenpairs
int nma_loop(Macromolecule *mol, int model, int type, tri *props, float *masses, int ifa, int ifr, int ilr,
		int na, int size, int nco, float cutoff, double *eigval, double *eigvect, trd **p_der);

// Get loop coordinates into a pre-allocated "coord" array (intended for copy & paste or coordinates backup)
void get_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord);

// Set loop coordinates from a pre-allocated "coord" array (intended for copy & paste or coordinates backup)
void set_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord);

inline void update_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord);

// Show square matrix
void show_matrix(double *matrix, int size, char *name = "Showing matrix:", char *fmt = " %6.2f");

// Show cartesian normal mode in VMD
void show_cartmode(float *coord, double *Aop, tri *props, int ifr, int nla, char *text, int imod, float factor=1.0);

// Show a matrix vector-wise...
void show_vectors(FILE *f, double *v, int size, int n, char *name = "Showing matrix vector-wise:", char *fmt = " %5.2e");

// Show a single vector (double)
void show_vector(FILE *f, double *v, int size, char *name = "Showing some vector:", char *fmt = " %5.2e", bool newline=true, bool newline2=true);

// Show a single vector (float)
void show_vector(FILE *f, float *v, int size, char *name = "Showing some vector:", char *fmt = " %5.2e", bool newline=true, bool newline2=true);

// Reads all lines from a text file and returns the number of lines read (automatic memory allocation)
//  file       --> File name
//  p_lines    --> Pointer to the lines array (it will be allocated automatically)
//  linelength --> Length of each line (number of characters)
int readTextLines(char *file, char ***p_lines, int linelength = 50); // Reading all rows from textfile

// Compute IC "mode" from "eigvect" IC eigenvectors and "alpha" and "delta" arrays. Requires:
//	"nevec"  --> Number of eigenvectors
//	"size"   --> Number of dihedral coordinates
void make_mode(double *alpha, double *delta, double *eigvect, double *mode, int nevec, int size);

// Compute element-wise "p" power of "in" vector into "out" vector, both of same length "len"
void pow_vector(double *in, double *out, int len, int p);

// Adds the "len" elements in vector "in"
double sum_vector(double *in, int len);

// Get the index of the maximum value in the array
int get_max_index(double *v, int n);


//==============================================================================================
int main( int argc, char * argv[] )
{
	fprintf( stdout, "%s>\n%s> Welcome to the Internal coordinates Loops Modal analysis tool v%s\n%s>\n",prog,prog,version,prog);
	std::string temp;

	auto t1 = std::chrono::steady_clock::now();


	// COMMAND-LINE PARSER:
	using namespace TCLAP;
	CmdLine cmd(argv[0],"Some text here... Shown with -h option", version );
	try {

		//		fprintf(stdout,"EXAMPLES: \n"
		//				"  Example-1: %s DHV15_3agoA2.pdb 29 40 -t DHV15_3agnA1_ali.pdb -m 2 -i 1 -a 1 -C 1 -o helixchi -x\n"
		//				"  Example-2: %s DHV15_3agoA2.pdb 29 40 -t DHV15_3ahwA3_ali.pdb -m 2 -i 1 -a 1 -C 1 -o morphchi -x\n"
		//				"  Example-3: %s DHV15_3agoA2.pdb 29 40 -m 2 -n 3 -i 1 -a 1 -C 1 -o kk -x\n"
		//				"  Example-4: %s DHV15_3agoA2.pdb 29 40 -m 2 -i 1 -a 1 -s 2 --ns 100 --rmsd 2 -o kk\n",prog,prog,prog,prog);

		// Define required arguments no labeled (mandatory)
		// ---------------------------------------------------------------------------------------------------------------

		UnlabeledValueArg<std::string> Input("PDB","Single-run: Input PDB file (Receptor PDB for docking)\n"
				"Multi-run: Plain text file (.txt) with one basename of the PDB and loops per row, e.g. 1tca_1 corresponds to PDB file 1tca_1.pdb and loops file 1tca_1<loops_suffix>","default","pdb");
		cmd.add( Input );

		UnlabeledValueArg<int> Start("start","N-terminal loop residue index, first mobile residue (PDB residue index).",1,"start");
		cmd.add( Start );

		UnlabeledValueArg<int> End("end","C-terminal loop residue index, last mobile residue (PDB residue index).",2,"end");
		cmd.add( End );

		//		ValueArg<int> Start("","start","N-terminal loop residue index, first mobile residue (PDB residue index).",false,-1,"int");
		//		cmd.add( Start );
		//
		//        ValueArg<int> End("","end","C-terminal loop residue index, last mobile residue (PDB residue index).",false,-1,"int");
		//		cmd.add( End );


		// Optional arguments
		// ---------------------------------------------------------------------------------------------------------------

		ValueArg<int> Verb("","verb", "Verbose level (0=low, 1=medium, 2=high) (default=0).",false,0,"int");
		cmd.add( Verb );

		ValueArg<unsigned int> Seed("","seed", "Pre-define the random number generator SEED (Mersenne Twister) (default=random-seed from /dev/urandom)",false,386,"unsigned int");
		cmd.add( Seed );

		ValueArg<double> DC("","dc", "Characteristic Distance of the Hardy's Quadric Interpolation used in Deformability computations (default=15). It should be > 0.",false,15,"double");
		cmd.add( DC );

		SwitchArg SkipMissingAtoms("","skip_missingatoms", "Disable missing atoms check. This way you can use N,CA,C,O (or others) Coarse-Grained models as long as dihedrals can be defined. (default=disabled)", false);
		cmd.add( SkipMissingAtoms );

		SwitchArg DelHeteros("","delete_heteros", "Delete Hetero-atoms, including waters (default=disabled).", true);
		cmd.add( DelHeteros );

		SwitchArg KeepWaters("","keep_waters", "Disables Water molecules deletion (default=disabled).", true);
		cmd.add( KeepWaters );

		SwitchArg KeepHydrogens("","keep_hydrogens", "Disables Hydrogen atoms deletion (default=disabled).", true);
		cmd.add( KeepHydrogens );

		SwitchArg Norm("", "norm","Enables (norm=1) eigenvector normalization. "
				"Note this does not change vector direction (default=disabled).", true);
		cmd.add( Norm );

		SwitchArg NoTors("","notors", "Disables extra torsional potential (default=disabled).", true);
		cmd.add( NoTors );

		ValueArg<double> MaxAng("a", "maxang","Maximum angular increment to normalize modes [deg] (default=1.0).",false,1.0,"float");
		cmd.add( MaxAng );

		ValueArg<float> k1_Cte("", "k1_k","Tirion's method stiffness constant (default=1.0).",false,1.0,"float");
		cmd.add( k1_Cte );

		ValueArg<float> k1_Cutoff("","k1_c","Tirion's method distance cutoff (default=10A).", false, 10,"float");
		cmd.add( k1_Cutoff );

		ValueArg<float> k0_Power("","k0_p", "Sigmoid function power term (default=6).",false, 6.0,"float");
		cmd.add( k0_Power);

		ValueArg<float> k0_X0("", "k0_x0","Sigmoid function inflexion point (default=3.8A).",false, 3.8,"float");
		cmd.add( k0_X0 );

		ValueArg<float> k0_Cte("", "k0_k","Sigmoid function stiffness constant (default=1000.0).",false, 1000.0,"float");
		cmd.add( k0_Cte );

		ValueArg<float> k0_Cutoff("","k0_c","Sigmoid function distance cutoff (default=10A).", false, 10,"float");
		cmd.add( k0_Cutoff );

		ValueArg<int> IndexMode("i","imod", "Index to select some mode (1,2,...,nevs) for computations (default=1).",false,1,"int");
		cmd.add( IndexMode );

		SwitchArg Linear("","linear", "Enable linear motion using 1st order dy/dq derivatives, otherwise exact dihedral angles rotation (default=disabled).", true);
		cmd.add( Linear );

		SwitchArg Chi("x","chi", "Considers CHI dihedral angle (default=disabled).", true);
		cmd.add( Chi );

		ValueArg<float> Nevs("n","nevs", "Used modes range, either number [1,N] <integer>, or ratio [0,1) <float>. 0=All modes, 1= 1 mode. (default=All).",false,0.0,"int/float");
		cmd.add( Nevs );

		ValueArg<int> Scoring("","score", "Scoring model (only utilized in morphing): (default=1)\n"
				"  1= Atomic RMSD between current and target structures.\n"
				"  2= Dihedral angles RMSD between current and target structures.\n",false,1,"int");
		cmd.add( Scoring );

		ValueArg<int> Contact("C","contact", "Contact rule: (default=1)\n"
				"  0= Consider all intra-loop contacts, i.e. all contacts involving only loop (pseudo)atoms (loop vs. loop).\n"
				"  1= Consider all contacts involving any loop (pseudo)atom (both loop vs. loop and loop vs. neighborhood).\n",false,1,"int");
		cmd.add( Contact );

		ValueArg<int> Potential("P","potential", "Pairwise interaction potential: (default=0)\n"
				"  0= Sigmoid function (= k/(1+(x/x0)^p), if x < c, else k=0).\n"
				"  1= Tirion's cutoff (= k, if x < c, else k=0).\n"
				"  4= edNMA formalism (CA-model only).\n"
				"  By default an extra torsional potential will be added.",false,0,"int");
		cmd.add( Potential );

		SwitchArg NoAnchors("","noanchors", "Anchors not present in loops Multi-PDB. By default anchors should be found in loops Multi-PDB. (default=disabled)", false);
		cmd.add( NoAnchors );

		ValueArg<std::string> Loops("l","loops", "Loops Multi-PDB suffix. Either use \"_rasp.pdb\" suffix to work with RCD generated loops in Multi-run mode, or loops Multi-PDB file name for Single-run mode. Each loop coordinates will be pasted into the main PDB file to perform calculations. (default=none).",false,"none","string");
		cmd.add( Loops );

		SwitchArg AliFlanks("","aliflanks", "Enable flanks alignment. (default=disabled)", false);
		cmd.add( AliFlanks );

		ValueArg<int> NFlanks("","flanks", "Number of loop flanking residues. E.g. \"2\" means two flanking residues at Nt and other two at Ct, \"0\" means disabled (default=0).",false,0,"int");
		cmd.add( NFlanks );

		ValueArg<float> dRMSD("","drmsd", "RMSD increment (Delta-RMSD) to save a sample in mode-following strategies (default disabled).",false,999999,"float");
		cmd.add( dRMSD );

		ValueArg<float> RMSD("r","rmsd", "Target RMSD for mode-following strategies (default=2).",false,2.0,"float");
		cmd.add( RMSD );

		ValueArg<int> NSamples("","nr", "Number of runs/frames for -s 2  (default=1000)",false,1000,"int");
		cmd.add( NSamples );

		ValueArg<int> NSteps("","ns", "Max number of sampling steps (default=1000).",false,1000,"int");
		cmd.add( NSteps );

		ValueArg<int> Nloops("","nloops", "Exact number of loop saved",false,0,"int");
	    cmd.add( Nloops );



		ValueArg<int> Strategy("s","strategy","Sampling strategy direction (default=2)\n"
				"    0= Single mode direction defined by -i\n"
				"    1= Random contribution of all the modes\n"
				"    2= Mutiple sampling S=1 --ns times\n"
				"    2= Mode following if morphing\n"
				,false,0,"int");
		cmd.add( Strategy );

		ValueArg<std::string> Target("t","target", "Set target PDB file name to enable morphing (default=none).",false,"none","string");
		cmd.add( Target );

		SwitchArg MR("","mc", " Multiple random modal directions.", true);
		cmd.add( MR );

		ValueArg<std::string> Name("o","name", "Output files basename (default=ilmode).",false,prog,"string");
		cmd.add( Name );

		ValueArg<int> Model("m","model", "Coarse-Grained model: 0=CA, 1=C5, 2=Heavy-Atom, 3=N,CA,C (default=2)",false,2,"int");
		cmd.add( Model );

		ValueArg<char> Chain2("","chain2", "Chain ID for the target PDB (default=first chain)",false,1,"int");
		cmd.add( Chain2 );

		ValueArg<char> Chain("","chain", "Chain ID for the initial PDB (default=first chain of loops Multi-PDB)",false,1,"int");
		cmd.add( Chain );




		// Parse the command line.
		// ---------------------------------------------------------------------------------------------------------------
		cmd.parse(argc,argv);

		//  Start measuring time
		//auto start = std::chrono::high_resolution_clock::now();






		// Load variables from parser
		// ---------------------------------------------------------------------------------------------------------------

		strcpy(file_pdb,((temp=Input.getValue()).c_str()));

		start = Start.getValue();
		loop_end = End.getValue();

		if (loop_end - start  < 2  )
		{
			fprintf(stderr, "ilmode>  Error loop index (%d) < start index (%d) + 2 \n",loop_end, start);
			exit(1);
		}


		//		if(Start.isSet() && End.isSet())
		//			loop_startend = true; // = true if Start and End residue indices were parsed

		if(Chain.isSet())
			chain = Chain.getValue();
		fprintf(stdout,"%s> Initial loop selected from residue %d (Nt) to %d (Ct) (both included) of chain %c\n", prog, start, loop_end, chain);

		strcpy(base,((temp=Name.getValue()).c_str())); // Gets Basename

		if(Target.isSet())
		{
			strcpy(file_pdb2,((temp=Target.getValue()).c_str())); // Target PDB file name (INPUT)
			morph_switch = true; // =true --> Enable morphing stuff

			if(Chain2.isSet())
				chain2 = Chain2.getValue();
			fprintf(stdout,"%s> Target loop selected from residue %d (Nt) to %d (Ct) (both included) of chain %c\n", prog, start, loop_end, chain2);
		}
		if(Loops.isSet())
		{
			strcpy(file_loops,((temp=Loops.getValue()).c_str())); // Loops Multi-PDB file name
			loops_switch = true; // =true --> Enable multiple loops processing
		}

		anchors_present = !NoAnchors.isSet();

		// Contacting method parameters

		power = k0_Power.getValue();
		power = power/2;

		cte_k0 = k0_Cte.getValue();


		x0 = k0_X0.getValue()*k0_X0.getValue();

		cutoff_k0 = k0_Cutoff.getValue();

		cte_k1 = k1_Cte.getValue();

		cutoff_k1 = k1_Cutoff.getValue();

		maxang = MaxAng.getValue(); // Maximum angular increment to normalize modes [deg]

		// Number of eigenvectors to be computed
		nevec_fact = Nevs.getValue();
		if(nevec_fact < 0) // checking
		{
			fprintf(stdout,"%s> Error, invalid number of eigenvectors requested (%f)!\nForcing exit!\n",prog,nevec_fact);
			exit(1);
		}



		imod = IndexMode.getValue() - 1; // indices must begin in "0"

		if (imod < 0 ) // checking
		{
			fprintf(stdout,"%s>  Invalid mode number (%d) > 1\nForcing exit!\n",prog,imod);
			exit(1);
		}




		verb = Verb.getValue();
		fprintf(stdout,"%s> Verbose level: %d\n", prog, verb);

		if(Norm.isSet())
		{
			norm_modes = true;
			fprintf(stdout,"%s> Normal modes will be normalized into unit vectors.\n", prog);
		}

		// Setting model and chi
		model = Model.getValue();
		if(Chi.isSet())
			type = 2; // phi,chi,psi
		else
			type = 0; // phi,psi

		if (model==0 or model==3) {
			if (type == 2) {
				type = 0; // phi,psi
				fprintf(stdout,"%s> Chi angle incompatible with this model\n", prog);

			}
		}

		linear_switch = Linear.isSet(); // =true --> linear motion

		contact = Contact.getValue(); // Contact method

		scoring = Scoring.getValue(); // Scoring method

		potential = Potential.getValue(); // Potential method

		if( potential == 4 && model != 0 ) // ED-NMA only valid for CA-model
		{
			fprintf(stdout,"%s> At this moment, the edNMA potential is only valid for CA-model!\nForcing exit!\n", prog);
			exit(1);
		}

		mr_switch = MR.isSet(); // Enable mutliple random directions

		strategy = Strategy.getValue(); // Sampling strategy(0=raw sampling for selected mode "i", 1= hyper-square sampling, etc...

		if (strategy) {
			mr_switch = true;
		}

		if (morph_switch) {
			mr_switch = false;
		}



		nsamples = NSamples.getValue(); // Number of samples for selected strategy
		nsteps = NSteps.getValue(); // Number of samples for selected strategy
		max_loops_save = Nloops.getValue(); // Max Number of loops saved
		target_rmsd = RMSD.getValue(); // Target RMSD for mode-following strategies
		delta_rmsd = dRMSD.getValue(); // Delta RMSD in mode-following strategies

		ali_flanks = AliFlanks.isSet(); // Enables flanks alignment

		if( NFlanks.isSet() )
			nflanks = NFlanks.getValue(); // Number of loop flanking residues (e.g. "2" means two flanking residues at Nt, and other two at Ct)

		if(SkipMissingAtoms.isSet())
			skip_missingatoms = true; // Set true to skip missing atoms

		if(Seed.isSet()) // Fixed seed
		{
			seed = (unsigned int) Seed.getValue(); // Gets seed
			already_seed = true; // = true if seed has been defined by user
		}
		else // Random seed (time initialization)
		{     // Needed to avoid seed repetition between different runs.
			FILE *fp;
			unsigned char b[4];
			int l=0;
			if ((fp = fopen("/dev/urandom", "r")) == NULL)
			{
				fprintf(stderr, "%s> Error! Could not open /dev/urandom for read\n%s> Exit\n", prog, prog);
				exit(2);
			}
			fread(b,1,4,fp);
			l |= b[0] & 0xFF;
			l <<= 8;
			l |= b[1] & 0xFF;
			l <<= 8;
			l |= b[2] & 0xFF;
			l <<= 8;
			l |= b[3] & 0xFF;
			seed = (unsigned int) l;
			fclose(fp);
		}
		fprintf(stdout,"%s> Mersenne Twister's SEED: --seed = %u\n",prog, seed);

		if(DelHeteros.isSet())
			delHeteros_switch = true; // Delete heteroatoms

		if(KeepHydrogens.isSet())
			delHydrogens_switch = false; // Keep hydrogens

		if(KeepWaters.isSet())
			delWaters_switch = false; // Keep waters

		dc = DC.getValue();
		if(dc <= 0.0)
			dc = 0.001;
	}
	catch ( ArgException& e )
	{
		std::cout << "  Error->" << e.error() << " " << e.argId() << std::endl;
	}
	using namespace std;

	// Mersenne Twister Seed Initialization
	// Output random float number in the interval 0 <= x < 1, with Random()
	rg = new CRandomMersenne( seed );

	// Saving Input Log File
	FILE *f_log;

	// CA  Conditions
	Condition *calpha = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	calpha->add( " CA " );
	Conditions *calpha2 = new Conditions();
	calpha2->add( calpha );

	// NCAC  Conditions ( N-, CA-, C- selection)
	Condition *ncac = new Condition( -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 );
	ncac->add( " N  " );
	ncac->add( " CA " );
	ncac->add( " C  " );
	Conditions *ncac2 = new Conditions();
	ncac2->add( ncac );

	//	// Saving Input Command-Log File
	//	FILE *f_com;
	//	sprintf(text,"%s%s.log",prog,base);
	//	if( !(f_com=(FILE *)fopen(text,"w") ) )
	//	{
	//		fprintf(stdout,"Sorry, unable to open LOG FILE: %s\n",text);
	//		exit(1);
	//	}
	//	for(int i=0; i<argc; i++)
	//		fprintf(f_com,"%s ",argv[i]);
	//	if(!already_seed)
	//		fprintf(f_com,"--seed %u\n",seed); // This allows user to carry out again the same run.
	//	fclose(f_com);
	//
	//	sprintf(saved_files,"%s> Com file:                            %35s\n",prog,text);

	Residue *res;
	Atom *at;
	int nipa;

	// Initialize aminoacids and nucleotids
	init_aminoacids();

	inputdata *input = NULL;
	//	if( !(input = (inputdata*) malloc( sizeof(inputdata) * 1 ))) // Allocate just one element (later realloc)
	//	{
	//		fprintf(stdout,"Sorry, unable to allocate memory!!!\n");
	//		exit(1);
	//	}

	// Selecting Single- or Multiple- run modes
	int npdbs = 0; // Total number of cases to be processed (1 by default)
	//	char **pdbnames = NULL; // List of PDBs for Multiple Run mode
	//	char *currentloops; // Current Loops filename
	if(strncmp(file_pdb + strlen(file_pdb)-4, ".txt", 4) == 0)
	{
		fprintf(stdout,"MULTIPLE RUN MODE (%s has .txt extension)\n",file_pdb);

		//-------------- READING INPUT TEXT FILE FOR MULTI-LOOP ---------------//
		FILE *p_file = fopen(file_pdb, "r");		/*Read input file*/
		if(p_file == NULL)					/*If problem with opening file...*/
		{
			fprintf(stdout,"rcd> Input file %s could not be opened...\nrcd>\n", file_pdb);
			exit(1);
		}									/*File opening check passed...*/

		char mystring[2];
		char myline[1024];
		int nscan = 0;
		while( fgets(myline, 1024, p_file) )
		{
			//fprintf(stdout,"%s\n",myline);

			if(myline[0] != '#')
			{
				// Memory allocation
				if( !(input = (inputdata*) realloc(input, sizeof(inputdata) * (npdbs+1) ))) // Reallocate elements
				{
					fprintf(stdout,"Sorry, unable to re-allocate memory for %d elements!!!\n", (npdbs+1));
					exit(1);
				}

				nscan = sscanf(myline,"%s %d %d %s\n", &input[npdbs].fname, &input[npdbs].start, &input[npdbs].end, mystring);

				if( nscan != 4  )
				{
					fprintf(stderr,"rcd> Please, check input text file... nscan=%d \nForcing exit!\n",nscan);
					exit(1);
				}

				input[npdbs].chain = mystring[0]; // with %c this does not work...

				fprintf(stdout,"rcd> %2d %8s %8d %8d %c\n",npdbs+1, input[npdbs].fname, input[npdbs].start, input[npdbs].end, input[npdbs].chain);

				npdbs++;
			}
		}

		fclose(p_file); // Close file upon reading input text file...
		for(int p = 0; p < npdbs; p++)
		{
			fprintf( stdout, "%s> %2d %8s %8d %8d %c\n", prog, p+1, input[p].fname, input[p].start, input[p].end, input[p].chain);
		}

	}
	else if(strncmp(file_pdb+strlen(file_pdb)-4,".pdb",4) == 0)
	{
		fprintf(stdout,"ilmode> Single run mode (%s has .pdb extension)\n",file_pdb);
		npdbs = 1; // Single-run

		// Memory allocation
		input = (inputdata*) realloc(input, sizeof(inputdata) ); // Allocate just one element

		strcpy(input[0].fname, file_pdb);
		input[0].start = start;
		input[0].end = loop_end;
		input[0].chain = chain;
	}
	else
	{
		fprintf(stdout,"%s has invalid extension: %s\n Forcing exit!!!\n", file_pdb, file_pdb+strlen(file_pdb)-4);
		exit(2);
	}


	// exit(0);


	// Process all input PDBs in single (npdbs=1) or multiple run (npdbs>1) modes
	for(int p = 0; p < npdbs; p++)
	{
		float delta0; // Initial Delta
		float rmsd0; // Initial RMSD
		float rmsd0_ncac; // Initial RMSD
		float rmsd; // Current RMSD

		// Generate basename for current "p"
		//		strcpy(name,base); // Overwrite previous "name" with the general "base"
		//		strncat(name, input[p].fname, strlen(input[p].fname) - 4); // Concatenate without extension and dot
		int len = strlen(input[p].fname) - 4;
		strncpy(name, input[p].fname, len); // Concatenate without extension and dot
		name[len] = '\0'; // mandatory...
		strcat(name,base); // Overwrite previous "name" with the general "base"
		fprintf( stdout, "%s> Current basename: %s\n", prog, name );

		// Copy required variables
		start = input[p].start;
		loop_end = input[p].end;
		chain = input[p].chain;

		// READING INPUT INITIAL PDB

		fprintf( stdout, "%s> Reading PDB file: %s\n", prog, input[p].fname );
		Macromolecule *molr = new Macromolecule( input[p].fname );
		molr->readPDB(input[p].fname);
		if(delHydrogens_switch)
		{
			if (verb > 1)
				fprintf( stdout, "%s> Deleting Hydrogen atoms (if any)...\n",prog );
			molr->deleteHYDS();
		}
		if(delHeteros_switch || model == 0) // CA model can't deal with HETATM's... (TO DO)
		{
			if (verb > 1)
				fprintf( stdout, "%s> Deleting Hetero-atoms (if any)...\n",prog );
			molr->delete_heteros();
		}
		if(delWaters_switch)
		{
			if (verb > 1)
				fprintf( stdout, "%s> Deleting Water molecules (if any)...\n", prog );
			molr->delete_waters();
		}
		if (verb > 1)
			molr->info(stdout);

		// Formating Initial PDB first
		if(verb > 1)
			fprintf( stdout, "%s> Formatting residues order and checking for missing atoms\n", prog );

		if (model==0) check_model=0;
		if (model==3) check_model=0;
		if (model==1) check_model=1;
		if (model==2) check_model=2;

		// fprintf( stdout, "%s> Model  %d %d\n", prog, model, check_model );



		if(molr->format_residues(false,check_model) > 0)
		{
			if(skip_missingatoms) // skip missing atoms
				fprintf( stdout, "%s> Warning, missing atom(s) found! Be aware of wrong results!\n", prog );
			else
			{
				fprintf( stdout, "%s> Error, missing atom(s) found! Forcing exit!\n", prog );
				exit(1);
			}
		}

		if(debug)
		{
			sprintf(dummy,"%s_format.pdb",name);
			molr->writePDB( dummy );
			sprintf(text,"%s> Formated model PDB:                  %35s\n", prog, dummy);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files + saved_files_len, text);
		}

		Macromolecule *mol, *molNCAC, *molini, *molini2, *molr2, *mol2, *molNCAC2;

		if(morph_switch) // Morphing stuff
		{
			// READING TARGET PDB
			fprintf( stdout, "%s> Reading Target PDB file: %s\n", prog, file_pdb2 );
			molr2 = new Macromolecule( file_pdb2 );
			molr2->readPDB(file_pdb2);
			if(delHydrogens_switch)
			{
				if (verb > 1)
					fprintf( stdout, "%s> Deleting Hydrogen atoms (if any)...\n",prog );
				molr2->deleteHYDS();
			}
			if(delHeteros_switch || model == 0) // CA model can't deal with HETATM's... (TO DO)
			{
				if (verb > 1)
					fprintf( stdout, "%s> Deleting Hetero-atoms (if any)...\n",prog );
				molr2->delete_heteros();
			}
			if(delWaters_switch)
			{
				if (verb > 1)
					fprintf( stdout, "%s> Deleting Water molecules (if any)...\n", prog );
				molr2->delete_waters();
			}
			if (verb > 1)
				molr2->info(stdout);

			// Formating TARGET PDB first
			if(verb > 1)
				fprintf( stdout, "%s> Formatting residues order and checking for missing atoms\n", prog );


			//	fprintf( stdout, "%s> Model  %d %d\n", prog, model, check_model );


			if(molr2->format_residues(false,check_model) > 0)
			{
				if(skip_missingatoms) // skip missing atoms
					fprintf( stdout, "%s> Warning, missing atom(s) found! Be aware of wrong results!\n", prog );
				else
				{
					fprintf( stdout, "%s> Error, missing atom(s) found! Forcing exit!\n", prog );
					exit(1);
				}
			}
		}

		// Setting Coarse-Graining model ("mol" will hold current CG'ed Macromol.)
		switch(model)
		{
		case 0: // CA-IC model: CA + (NH and CO)-terminal
		{
			fprintf( stdout, "%s> Coarse-Graining model: CA-model\n", prog);
			fprintf( stdout, "%s> Not fully implemented Exit\n", prog);
			exit(1);

			// N,CA,C selection
			molNCAC = molr->select( ncac2 );

			// Creates a CA-model with first NH and last CO pseudo-atoms of each segment.
			// Warning, atoms are not copied, they're just pointers to the original atoms.
			// setmass = true --> adds masses to Occupancy and Bfactor, otherwise left unchanged.
			// nomass = true --> sets 1.0 masses to all CAs excepting NH and CO at segment endings
			//                   (unit mass is divided between NH or CO and their CAs)
			// nomass = false --> whole residue mass applied to all CAs excepting NH and CO at segment endings
			//                   ( 15, 28 and Residue_mass-(NH_mass or CO_mass) for NH, CO and its CAs)

			// Makes model (if it's needed)
			mol = cg_CA( molNCAC ); // masses will be computed
			//mol = cg_CA(molNCAC,  false,  false, false, false );

			if(morph_switch)
			{
				// N,CA,C selection
				molNCAC2 = molr2->select( ncac2 );
				// Makes model (if it's needed)
				mol2 = cg_CA( molNCAC2 ); // masses will be computed
			}
			break;
		}
		case 3: // N,CA,C-model
		{
			fprintf( stdout, "%s> Coarse-Graining model: N,CA,C-model (ideal for KORP integration)\n", prog);

			mol = molr->select( ncac2 ); // N,CA,C selection
			// mass_NCAC( mol ); // Makes model bug pablo 2022
			mass_NCAC( mol, false, true, false, false );

			if(morph_switch)
			{
				mol2 = molr2->select( ncac2 ); // N,CA,C selection
				// mass_NCAC( mol2 ); // Makes model bug pablo 2022
				mass_NCAC( mol2, false, true, false, false );

			}

			break;
		}
		case 1: // 3BB2R model
		{
			fprintf( stdout, "%s> Coarse-Graining model: 3BB2R\n", prog);
			mol = molr;

			// Makes 3BB2R reduced model
			//     Each residue is represented by 3 backbone atoms and 2 side-chain atoms.
			//     Thus, 3 dihedral angles are needed for each residue (phi, psi, chi).
			//     There are a few exceptions: for Ala, Gly and Pro,
			//     and for the 1st and last residues.
			fprintf( stdout, "%s> Creating 3BB2R reduced model:\n", prog);
			cg_3BBR2( mol );

			if(morph_switch)
			{
				mol2 = molr2;
				cg_3BBR2( mol2 );
			}

			break;
		}
		case 2: // Full-Atom
		{
			fprintf( stdout, "%s> Coarse-Graining model: All heavy atoms (no coarse-graining)\n", prog);
			mol = molr;

			// Add appropriate masses to the All-Heavy-Atoms model
			mass_FA( mol ); // sets masses

			if(morph_switch)
			{
				mol2 = molr2;
				mass_FA( mol2 ); // sets masses
			}

			break;
		}
		}

		// Create iterators for RMSD computations
		pdbIter *itermol = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule
		pdbIter *itermol2 = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule to speed up Elastic Network refresh

		// Copy initial molecule (needed further to convert normal modes into the appropriate atomic model and for RMSD computations)
		molini = new Macromolecule(mol); // Initial coarse-grained macromolecule copy
		pdbIter *iterini = new pdbIter( molini, true, true, true, true ); // Iterator to the initial (reference) macromolecule

		molini2 = new Macromolecule(mol); // Initial coarse-grained macromolecule copy
		pdbIter *iterini2 = new pdbIter( molini2, true, true, true, true ); // Iterator to the initial (reference) macromolecule


		pdbIter *itertar; // Iterator to target macromolecule
		if(morph_switch)
			itertar = new pdbIter( mol2, true, true, true, true ); // Iterator to the initial (reference) macromolecule



		// Initializing some internal coords. related residue properties
		// (# atoms, # units, # internal coordinates, etc...)
		tri *props,*props2;
		int *unat,*unat2;
		switch(model)
		{
		case 0:
		case 3:
			//						properCA(mol,&props,&unat);
			//						if(morph_switch)
			//							properCA(mol2,&props2,&unat2);
			//						break;
		case 1:
		case 2:
			properMFA(mol,&props,&unat,type,model);
			if(morph_switch)
				properMFA(mol2,&props2,&unat2,type,model);
			break;
		}

		// Saving current model
		if(savemodel_switch)
		{
			sprintf(dummy,"%s_model.pdb",name);
			mol->writePDB( dummy ); // NO-renumber the PDB
			sprintf(text,"%s> Model PDB:                           %35s\n", prog, dummy);
			fprintf(stdout,"%s> Model PDB: %s\n", prog, dummy);
			saved_files_len = strlen(saved_files);
			strcpy(saved_files + saved_files_len, text);
		}



		// Enable multiple loops processing
		Macromolecule *loopsr;
		Macromolecule *loops;
		pdbIter *iter_loops; // Loops iterator
		int ri,rf;
		int nloops = 1; // Number of loops per Loops Multi-PDB file
		if(loops_switch)
		{
			// Generate muli-loops file name in "dummy"
			if(npdbs > 1) // Multi-run mode
			{
				int len = strlen(input[p].fname) - 4;
				strncpy(dummy, input[p].fname, len); // Overwrite with current pdb name without extension
				dummy[len] = '\0'; // required for "strncpy"
				strcat(dummy, file_loops); // Concatenate without extension and dot
				fprintf( stdout, "%s> Current Loops Multi-PDB file name: %s (file_loops= %s)\n", prog, dummy, file_loops );
			}
			else // Single-run
			{
				strcpy(dummy, file_loops); // Copy loops Multi-PDB file name
			}

			if(debug)
				fprintf( stdout, "%s> Reading Loops Multi-PDB: %s\n", prog, dummy );
			loopsr = new Macromolecule(dummy); // Reading loops Multi-PDB into a Macromolecule (each loop will be a molecule)
			loopsr->readPDB(dummy);
			if(debug)
				fprintf( stdout, "%s> Deleting Hydrogen atoms (if any)...\n", prog );
			loopsr->deleteHYDS();
			if(debug)
				fprintf( stdout, "%s> Deleting Hetero-atoms (if any)...\n", prog );
			loopsr->delete_heteros();
			if(debug)
				fprintf( stdout, "%s> Deleting Water molecules (if any)...\n", prog );
			loopsr->delete_waters();
			//		fprintf( stdout, "%s> Deleting duplicate atoms within residue (if any)...\n", prog );
			//		loopsr->delete_duplicates(); // Remove duplicate atoms

			// Formating residues of all Loops at once
			if(debug)
				fprintf( stdout, "%s> Formatting Loops Multi-PDB residues order and checking for missing atoms, model= %d\n", prog, model );



			if(loopsr->format_residues(false,check_model) > 0)
			{
				if(skip_missingatoms) // skip missing atoms
					fprintf( stdout, "%s> Warning, missing atom(s) found in loops Multi-PDB! Be aware of wrong results!\n", prog );
				else
				{
					fprintf( stdout, "%s> Error, missing atom(s) found in Loops Multi-PDB! (according to %d CG-model). Forcing exit!\n", prog, model );
					exit(1);
				}
			}
			else
				if(debug)
					fprintf(stderr, "%s> No missing atom(s) in Loops Multi-PDB!\n", prog );

			if(debug)
			{
				fprintf( stdout, "%s> Written formated Loops Multi-PDB\n", prog );
				loopsr->writePDB( "loopsformatted.pdb" );
			}

			iter_loops = new pdbIter(loopsr); // Loops iterator
			nloops = iter_loops->num_molecule(); // Get number of loops in Loops Multi-PDB
			// fprintf(stderr,"nloops= %d\n", nloops);
			// exit(0);

			// If not provided, it gets chain ID from first loop
			Chain *ch;
			if(chain=='*' || npdbs > 1) // If no-chain was specified or in multi-pdb run mode, then get chain-ID from Loops
			{
				iter_loops->pos_chain = 0;
				ch = iter_loops->get_chain();
				chain = ch->getName()[0];
			}
			//			fprintf(stderr,"chain= %c\n",chain);
			//			exit(0);

			// Get anchor residue indices from first loop in the Multi-PDB file
			pdbIter *iter_seg;
			Segment *seg;
			iter_loops->pos_segment=0;
			seg = iter_loops->get_segment();
			iter_seg = new pdbIter(seg);
			iter_seg->pos_fragment=0;
			ri = (iter_seg->get_fragment())->getIdNumber(); // get Nt anchor residue number (PDB) <-- First residue in PDB
			iter_seg->pos_fragment=iter_seg->num_fragment()-1;
			rf = (iter_seg->get_fragment())->getIdNumber(); // get Ct anchor residue number (PDB), i.e. the last that moves... <-- Last residue in PDB
			delete iter_seg;

			if(!anchors_present) // If anchors are not present in loops, indices must be modified
			{
				// converting into anchor indices
				ri--;
				rf++;
			}
			fprintf( stdout, "%s> Anchor residues obtained from Multi-PDB loop %s --> Nt %d, Ct %d (PDB numeration) Chain_id= \"%c\"\n", prog, file_loops, ri, rf, chain );

			// Mobile loop conditions
			if(anchors_present) // Anchors must be removed if present
			{
				Condition *mobile;
				Conditions *mobile2 = new Conditions();
				mobile = new Condition(-1,-1,-1,-1,-1,-1,ri+1,rf-1,-1,-1); // get residues from Nt+1 to Ct-1, i.e. only those mobile...
				mobile2->add(mobile);

				loops = loopsr->select_cpy(mobile2); // select only mobile residues
				delete loopsr;
				delete iter_loops;
				iter_loops = new pdbIter(loops); // Loops iterator
				loopsr = loops; // Now "loopsr" has only mobile residues (all atoms read)
			}

			if(debug)
				loops->writePDB("loopsmobile.pdb");

			// Cross-checking input
			if( start >= 0 && loop_end >= 0 && (ri+1 != start || rf-1 != loop_end) )
			{
				fprintf( stdout, "%s> Error, loop residue indices mismatch between parser or Loops-Multi-pdb (ri=%d rf=%d) and loops file (start=%d end=%d). Forcing exit!\n", prog, ri, rf, start, loop_end);
				exit(1);
			}

			// Setting indices of first and last mobile residues of loop
			start = ri+1; // Residue PDB index of first mobile residue
			loop_end = rf-1; // Residue PDB index of last mobile residue
			//			if(anchors_present) // If anchors are not present in loops, indices must be modified
			//			{
			//			}
			//			else
			//			{
			//				start = ri; // Residue PDB index of first mobile residue
			//				loop_end = rf; // Residue PDB index of last mobile residue
			//			}
		}
		// END of loops multi-pdb stuff...


		// Compute internal residue indices of initial loop boundaries (NOTE: CONSIDER CHAIN-ID SOME DAY...)
		int ifr = -1; // internal index of the first mobile residue of the initial loop
		int ilr = -1; // internal index of the last mobile residue of the initial loop
		Chain *ch;
		pdbIter *iter_ch = new pdbIter( mol, true, true, true, true ); // iter to screen fragments (residues)
		int size = 0; // NOT considering Nt-anchor PSI angle (its O atom should move if it were considered)
		int num_atoms_loop=0; // Number of loop pseudo-atoms
		int res_index = 0;
		int resid;
		for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() )
		{
			ch = iter_ch->get_chain();
			pdbIter *iter = new pdbIter( ch, true, true, true, true ); // iter to screen fragments (residues)

			if(chain == ch->getName()[0])
			{
				fprintf( stdout, "%s> Sequence: ",prog);
				id_chain=iter_ch->pos_chain ;
				for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
				{
					res = (Residue *) iter->get_fragment();
					resid = res->getIdNumber(); // Get residue index in PDB numeration

					if(resid >= start && resid <= loop_end)
					{
						num_atoms_loop += props[res_index].nat;
						size += props[res_index].nan;
						fprintf(stdout,"%c",AA[resnum_from_resname( res->getName() )].aa_name1);
					}

					// fprintf(stdout,"%d  nid= %d\n",iter->pos_fragment,res->getIdNumber());
					if(ifr < 0)
					{
						if(res->getIdNumber() == start) // ifr found!
							ifr = res_index;
					}
					else // ifr already found!
					{
						if(res->getIdNumber() == loop_end) // ilr found!
						{
							ilr = res_index;
							break; // exit for loop
						}
					}
					// fprintf(stderr, "ch= %c  res= %4d  ifr= %4d  ilr= %4d  props= %d\n", ch->getName()[0], res_index, ifr, ilr, props[res_index].k1);

					res_index++; // update absolute residue index
				}

			}
			else
				res_index += iter->num_fragment();

			delete iter;
		}
		delete iter_ch;

		fprintf( stdout, "\n");

		fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) mobile residues of the initial loop (chain %c)\n", prog, ifr, ilr, chain);

		// Some checking
		if(ifr < 0 || ilr < 0)
		{
			fprintf(stderr,"%s> Sorry, internal indices of first (%d) or last (%d) mobile residues of the initial loop (chain %c) not found! "
					"Forcing exit!\n", prog, ifr, ilr, chain);
			exit(1);
		}
		// Get the internal indices of first or last mobile atoms of the initial loop (required to get "just loop contacts")
		int ifa; // internal index of the first mobile atom of the loop
		int ila; // internal index of the last mobile atom of the loop
		ifa = props[ifr].k1;
		ila = props[ilr+1].k1 - 1;
		fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) mobile atoms of the initial loop\n", prog, ifa, ila);


		// Compute internal residue indices of target loop boundaries (NOTE: CONSIDER CHAIN-ID SOME DAY...)
		int ifr2 = -1; // internal index of the first mobile residue of the target loop
		int ilr2 = -1; // internal index of the last mobile residue of the target loop
		if(morph_switch)
		{
			iter_ch = new pdbIter( mol2, true, true, true, true ); // iter to screen fragments (residues)
			int size2 = 0; // NOT considering Nt-anchor PSI angle (its O atom should move if it were considered)
			int num_atoms_loop2=0; // Number of loop pseudo-atoms
			int res_index2 = 0;
			int resid2;
			bool chain_found = false; // safety check
			for( iter_ch->pos_chain = 0; !iter_ch->gend_chain(); iter_ch->next_chain() )
			{
				ch = iter_ch->get_chain();
				pdbIter *iter = new pdbIter( ch, true, true, true, true ); // iter to screen fragments (residues)
				fprintf( stdout, "%s> Sequence: ",prog);
				if(chain2 == ch->getName()[0] || chain2 == '*' || npdbs > 1) // If no-chain was specified or in multi-pdb run mode, then get chain-ID from Loops
				{
					chain_found = true; // safety check
					id_chain2=iter_ch->pos_chain ;

					for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() ) // screen residues
					{
						res = (Residue *) iter->get_fragment();
						resid2 = res->getIdNumber(); // Get residue index in PDB numeration

						if(resid2 >= start && resid <= loop_end)
						{
							num_atoms_loop2 += props[res_index].nat;
							size2 += props[res_index].nan;
							fprintf(stdout,"%c",AA[resnum_from_resname( res->getName() )].aa_name1);

						}

						// fprintf(stdout,"%d  nid= %d\n",iter->pos_fragment,res->getIdNumber());
						if(ifr2 < 0)
						{
							if(res->getIdNumber() == start) // ifr found!
								ifr2 = res_index2;
						}
						else // ifr already found!
						{
							if(res->getIdNumber() == loop_end) // ilr found!
							{
								ilr2 = res_index2;
								break; // exit for loop
							}
						}
						// fprintf(stderr, "ch= %c  res= %4d  ifr= %4d  ilr= %4d  props= %d\n", ch->getName()[0], res_index, ifr, ilr, props[res_index].k1);

						res_index2++; // update absolute residue index
					}
				}
				else
					res_index2 += iter->num_fragment();

				delete iter;
			}
			fprintf( stdout, "\n");

			delete iter_ch;

			if(!chain_found) // safety check
			{
				fprintf(stderr, "%s> Sorry, chain %c not found in target PDB! Forcing exit!\n\n", prog, chain2);
				exit(1);
			}

			fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) mobile residues of the target loop (chain %c)\n", prog, ifr2, ilr2, chain);

			// Some checking
			if(ifr2 < 0 || ilr2 < 0)
			{
				fprintf(stderr,"%s> Sorry, internal indices of first (%d) or last (%d) mobile residues of the target loop (chain %c) not found! "
						"Forcing exit!\n", prog, ifr2, ilr2, chain2);
				exit(1);
			}
		}


		int iffa; // internal index of the first flanking atom of the loop
		int ilfa; // internal index of the last flanking atom of the loop
		if(nflanks > 0)
		{
			iffa = props[ifr-nflanks].k1; // first
			ilfa = props[ilr+1+nflanks].k1 - 1; // last
			fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) flanking atoms of the initial loop\n", prog, iffa, ilfa);
		}

		// Saving Input Log File
		sprintf(file_log,"%s_traj.log", name);
		if( !(f_log=(FILE *)fopen(file_log,"w") ) )
		{
			fprintf(stdout,"Sorry, unable to open LOG FILE: %s\n",file_log);
			exit(1);
		}
		fprintf(f_log,"# %s> Welcome to %s v%s\n%s> COMMAND: ",prog,prog,version,prog);
		bool already_seed=false;
		for(int i=0; i<argc; i++)
		{
			fprintf(f_log,"%s ",argv[i]);
			if(strcmp(argv[i],"--seed") == 0 )
				already_seed = true;
		}
		if(!already_seed)
			fprintf(f_log,"--seed %u\n",seed); // This allows user to carry out again the same run.
		else
			fprintf(f_log,"\n");

		sprintf(saved_files,"%s> Log file:                            %35s\n",prog,file_log);



		int ifa2; // internal index of the first mobile atom of the target loop
		int ila2; // internal index of the last mobile atom of the target loop
		int iffa2; // internal index of the first flanking atom of the target loop
		int ilfa2; // internal index of the last flanking atom of the target loop
		if(morph_switch)
		{
			// Get the internal indices of first and last mobile atoms of the target loop (required to get "just loop contacts")
			ifa2 = props2[ifr2].k1;
			ila2 = props2[ilr2+1].k1 - 1;
			fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) mobile atoms of the target loop\n", prog, ifa2, ila2);

			// Get the internal indices of first and last flanking atoms of the target loop
			//			if(nflanks > 0)
			//			{
			//				iffa2 = props2[ifr2-nflanks].k1;
			//				ilfa2 = props2[ilr2+1+nflanks].k1 - 1;
			//				fprintf( stdout, "%s> Internal indices of first (%d) or last (%d) flanking atoms of the target loop\n", prog, iffa2, ilfa2);
			//
			//				// Computing RMSD of flanks
			//				flank_rmsd = rmsd_flank(itermol, itertar, iffa, ifa, ila, ilfa, iffa2);
			//				// fprintf(stderr,"Flanks RMSD: %f\n", flank_rmsd);
			//				// sprintf(dummy, "%s> Flanks RMSD with %d residues= %8f\n", prog, flank_rmsd, nflanks);
			//				sprintf(dummy, "%s> Flanks RMSD with %d residues= %8.2f ", prog, nflanks, flank_rmsd);
			//				fprintf(f_log, "%s", dummy); // Dump log info
			//				fprintf(stdout,"%s",dummy);
			//			}

			float flank_rmsd_min = 0.0; // Flanks RMSD upon alignment, if any

			// Flanks rigid body alignment using Kabsch & Sander stuff...

			if(nflanks > 0)
			{


				Conditions *conds = new Conditions();
				Condition *cond = new Condition(-1,-1,-1,-1,id_chain,id_chain,start-nflanks,start-1,-1,-1);
				Condition *condB = new Condition(-1,-1,-1,-1,id_chain,id_chain,loop_end+1, loop_end+nflanks,-1,-1);
				cond->add(" N  ");
				cond->add(" CA ");
				cond->add(" C  ");
				condB->add(" N  ");
				condB->add(" CA ");
				condB->add(" C  ");
				//cond->add(" O  ");
				// cond->add(" CB ");
				conds->add(condB);
				conds->add(cond);
				Macromolecule *loopInit = mol->select(conds);
				//loopInit->writePDB("loopI.pdb");



				Conditions *conds2 = new Conditions();
				Condition *cond2 = new Condition(-1,-1,-1,-1,id_chain2,id_chain2,start-nflanks,start-1,-1,-1);
				Condition *condB2 = new Condition(-1,-1,-1,-1,id_chain2,id_chain2,loop_end+1, loop_end+nflanks,-1,-1);
				cond2->add(" N  ");
				cond2->add(" CA ");
				cond2->add(" C  ");
				condB2->add(" N  ");
				condB2->add(" CA ");
				condB2->add(" C  ");
				conds2->add(condB2);
				conds2->add(cond2);
				Macromolecule *loopTarget = mol2->select(conds2);
				//loopTarget->writePDB("loopT.pdb");



				sprintf(dummy, "%s> Flanks RMSD with %d residues= %8.2f ", prog, nflanks, loopInit->rmsd(loopTarget));
				fprintf(f_log, "%s", dummy); // Dump log info
				fprintf(stdout,"%s",dummy);



				if(ali_flanks) {

					sprintf(file_align,"%s_align.pdb", name);
					f= fopen(file_align,"w");
					fclose(f);

					float matrix4[4][4];

					flank_rmsd_min = loopInit->minRmsd(loopTarget, matrix4);

					M4Rot *matrix4_op = new M4Rot(matrix4);
					mol2->applyAtoms(matrix4_op); // superpose
					//loopTarget->applyAtoms(matrix4_op); // superpose
					//loopTarget->writePDB("loopTA.pdb");

					delete matrix4_op;
					mol2->writeMloop(file_align, 1, ifr2-nflanks, ilr2+nflanks, chain2);



					sprintf(dummy, " Flanks Min_RMSD (%d residues)= %8.2f\n", nflanks, flank_rmsd_min);
					fprintf(f_log, "%s", dummy); // Dump log info
					fprintf(stdout,"%s",dummy);


				}

				/*
				bool *mask, *mask2; // Atomic masks for alignment

				// Allocate masks memory
				mask = (bool *) malloc( sizeof(bool) * mol->get_num_atoms() );
				mask2 = (bool *) malloc( sizeof(bool) * mol2->get_num_atoms() );

				// Reset masks
				for(int i = 0; i < mol->get_num_atoms(); i++)
					mask[i] = false;
				for(int i = 0; i < mol2->get_num_atoms(); i++)
					mask2[i] = false;

				// Activate flanks in input PDB mask
				for(int i = iffa; i < ifa; i++)
					mask[i] = true;

				for(int i = ila+1; i <= ilfa; i++)
					mask[i] = true;

				// Activate flanks in target PDB mask
				for(int i = iffa2; i < ifa2; i++)
					mask2[i] = true;
				for(int i = ila2+1; i <= ilfa2; i++)
					mask2[i] = true;

				fprintf(stdout, "%s> ---> Flanks %d %d %d %d \n", prog, iffa, ila, iffa2, ila2);


				// Compute transformation of target PDB and RMSD evaluation
				float matrix4[4][4];
				flank_rmsd_min = mol->minRmsd(mol2, matrix4, mask, mask2);
				mol2->writePDB("prealiflanks.pdb");

				// Alignment of target PDB (minRMSD) (computed with flanks selection, but applied to full model)
				M4Rot *matrix4_op = new M4Rot(matrix4);
				mol2->applyAtoms(matrix4_op); // superpose
				delete matrix4_op;

				mol2->writePDB("aliflanks.pdb");

				flank_rmsd = rmsd_flank(itermol, itertar, iffa, ifa, ila, ilfa, iffa2);
				fprintf(stderr,"Flanks RMSD after align %f\n", flank_rmsd);
				 */

			}


		}


		// Computing Total number of degrees of freedom (hessian matrix rank)
		int num_res; // Number of residues (mol)
		int num_atoms; // Number of (pseudo)atoms (mol)
		pdbIter *iter = new pdbIter( mol, true, true, true, true ); // iter to screen fragments (residues)
		num_res = iter->num_fragment();
		num_atoms = iter->num_atom();


		int reglen = loop_end - start + 1; // Number of loop residues (mobile)
		if(props[ifr + reglen].nan != 1) // If not Proline
			size++; // considering Ct-anchor Phi angle? (If Ct-anchor is Proline it does not have Phi)

		// sprintf(dummy, "%s> Residues: %d  Dihedrals(DoFs): %d \n", prog, ilr-ifr+1, size);

		int nco = 6; // Number of (scalar) constraints
		sprintf( dummy , "%s> Selected model residues: %d\n", prog, num_res );
		sprintf( dummy + strlen(dummy), "%s> Selected model (pseudo)atoms: %d\n", prog, num_atoms );
		sprintf( dummy + strlen(dummy), "%s> Number of residues in loop: %d\n", prog, reglen);
		sprintf( dummy + strlen(dummy), "%s> Number of pseudo-atoms in loop: %d\n", prog, num_atoms_loop);
		sprintf( dummy + strlen(dummy), "%s> Number of constraints: %d\n", prog, nco);
		sprintf( dummy + strlen(dummy), "%s> Number of free variables in loop: %d\n", prog, size);
		fprintf(f_log, "%s", dummy);
		fprintf(stdout,"%s",dummy);

		if (imod +1 >  size-nco ) {
			fprintf(stderr, "%s> Error selected mode %d must <=  %d\n%s> Exit\n%s>\n",prog, imod+1, size-nco, prog, prog );
			exit(1);
		}



		// Get masses array for eigenvector normalization
		float *masses_loop = NULL;
		masses_loop = get_masses_loop(itermol, ifa, num_atoms_loop + 3, true); // "true" to compute the square root of the masses for mass-weighting

		// Allocate "decint" to store the contacts list (ipas-list)
		twid *decint = NULL; // Contacts data structure
		twid *decint2 = NULL; // Contacts data structure

		// Number of eigenvectors to be computed (we need to know "size" first)

		if(nevec_fact == 0.0) // all modes requested
		{
			nevec = size - nco;
			fprintf( stdout, "%s> All non-trivial modes will be computed: %d\n", prog, nevec);
		}
		else if(nevec_fact >= 1.0) // number of modes
		{
			nevec = (int) nevec_fact;
			fprintf( stdout, "%s> Range of computed modes: 1-%d\n", prog, nevec);
		}
		else
		{
			nevec = (int) (nevec_fact * size);
			fprintf( stdout, "%s> Range of computed modes: 1-%d (%.0f%)\n", prog, nevec, nevec_fact*100);
		}
		// Checking
		if(nevec > size)
		{
			fprintf( stdout, "%s> Sorry, more eigenvectors requested (%d) than available (%d), forcing maximum.\n", prog, nevec, size);
			nevec = size;
		}
		else if(nevec <= 0) // checking
		{
			fprintf( stdout, "%s> Error, invalid number of eigenvectors requested %d (%f)!\nForcing exit!\n", prog, nevec, nevec_fact);
			exit(1);
		}

		double *eigval; // Eigenvalues
		eigval  = (double *) malloc( (size + nco) * sizeof(double) );
		ptr_check(eigval);

		double *eigvect; // Eigenvectors
		eigvect = (double *) malloc( (size + nco) * size * sizeof(double) );
		ptr_check(eigvect);

		double *mass_matrix; // Kinetic energy matrix
		double *hess_matrix; // Hessian matrix

		// Get array of atomic masses from a Macromolecule iterator
		float *masses; // masses array
		masses = get_masses(iter,num_atoms);

		// Set name for trajectory file
		//		if(linear_switch)
		//			sprintf(file_movie,"%s_linear%02d.pdb", name, imod+1);
		//		else
		//			sprintf(file_movie,"%s_dihedral%02d.pdb", name, imod+1);
		sprintf(file_movie,"%s_traj.pdb", name);
		sprintf(file_final,"%s_final.pdb", name); // Morphed PDB structure file name

		// Delete previous trajectory file (if any)


		f= fopen(file_movie,"w");
		fclose(f);
		f= fopen(file_final,"w");
		fclose(f);


		int ncomps = 3*(num_atoms_loop+3); // Number of components (x,y,z for each atom of mobile loop + 3 Ct-anchor atoms)

		double *refmode; // Reference mode to prevent random sign reversal upon sucesive NMAs and to drive motion
		refmode = (double *) malloc(ncomps * sizeof(double)); // Reference mode to drive motion

		// Alpha overlaps array
		double *alpha = (double *) malloc( (size+nco) * sizeof(double)); // allocate memory for the maximum possible

		// Total delta overlaps array
		double *tdelta = (double *) malloc( (size+nco) * sizeof(double)); // allocate memory for the maximum possible

		for(int i=0; i<size+nco; i++)
			tdelta[i] = 0.0; // Initialize

		// Delta elements array
		double *vdelta = (double *) malloc( (size+nco) * sizeof(double)); // allocate memory for the maximum possible
		double *vdeltaS = (double *) malloc( (size+nco) * sizeof(double)); // allocate memory for the maximum possible

		// Merged mode
		double *mode = (double *) malloc( ncomps * sizeof(double)); // current merged mode for motion

		// Current loop coordinates
		float *coordx = (float *) malloc( sizeof(float) * num_atoms_loop * 3 ); // for future copy & paste loop coordinates

		float *coord2; // (pseudo)atomic coordinates single row vector
		//		if(morph_switch && scoring == 2)
		if(morph_switch)
		{
			if(verb > 1)
				fprintf(stdout, "%s> Getting Target PDB coordinates in single row format (pseudo-atom model)\n", prog);
			mol2->coordMatrix( &coord2 );
		}


		for(int l=0; l<nloops; l++) // Iter loops in Loops Multi-PDB
		{
			// fprintf(stdout,"Loop stuff %d of %d, num_atoms_loop= %d, ifr= %d  ilr= %d  chain= %c\n", l, nloops, num_atoms_loop, ifr, ilr, chain);

			//			mol->writePDB("mol0.pdb");
			//			molini->writePDB("molini0.pdb");
			//			loops->writePDB("loops0.pdb");

			if(loops_switch)
			{
				// Get current loop coordinates
				get_loop_coords(iter_loops, l * num_atoms_loop, num_atoms_loop, coordx);
				//				get_loop_coords(iterini, ifa, num_atoms_loop, coordx);

				// Set initial loop coordinates to prevent unwanted distortions
				set_loop_coords(itermol, ifa, num_atoms_loop, coordx);
				set_loop_coords(iterini, ifa, num_atoms_loop, coordx); // "iterini comes from molini, a different Macromolecule...
			}

			// Sampling strategy (0=raw sampling for selected mode "i", 1= hyper-square sampling, etc...)
			int neig; // Number of valid eigenvectors (typically --> size-nco)
			int info; // Info output parameter for diagonalization routine
			float ddrift = 0.0; // Distance drift at Ct end
			float adrift = 0.0; // Angle (bond angle) drift at Ct end
			float adist = 0.0; // Anchor distance
			float aang = 0.0; // Anchor angle
			float adist0 = 0.0; // Initial anchor distance
			float aang0 = 0.0; // Initial anchor angle

			// Arrays required to compute Dihedrals RMSD
			// WARNING! Dihedral angles order: Psi_NtAnchor, Phi_1st, Psi_1st, Phi_2nd, Psi_2nd, ... , Phi_Nth, Psi_Nth, Phi_CtAnchor
			int nlrs = ilr-ifr+1; // Number of mobile residues in loop
			float *dhs = NULL;
			float *dhs2 = NULL;

			// Measure anchor drift
			anchor_drift(iterini, itermol, props, ilr, &adist0, &aang0);
			fprintf( stdout, "%s> Initial anchor distance and angle: %f A and %f deg\n", prog, adist0, aang0);

			//			mol->writePDB("mol.pdb");
			//			molini->writePDB("molini.pdb");
			// exit(0);
			float rmsd_old; // Old RMSD
			float last_rmsd = 0; // RMSD of last saved structure


			if(morph_switch)
			{
				rmsd_old = 999999; // some high value required
			}
			else
				rmsd_old = 0.0; // zero value required




			switch (strategy)
			{
			//
			// CASE 0  "i"-th mode-following sampling
			//
			case 0: //   single random direction at f=0
			case 1: // Multiple random direction at f=0
			{
				double modref; // Reference vector modulus
				double modcurr; // Current vector modulus
				traji=0;
				if(verb > 1)
				fprintf(stdout,"> Initial structure dumped into Muli-PDB\n");
				mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);

				float *coord; // (pseudo)atomic coordinates single row vector

				for(int f = 0; f < nsteps; f++) // Generate N-samples (frames)
				{
					// Initialize Eigenvalues
					for(int i=0; i<size+nco; i++)
						eigval[i] = 0.0;

					// Initialize Eigenvectors
					for(int i=0; i<(size+nco)*size; i++)
						eigvect[i] = 0.0;

					if(verb > 1)
						fprintf(stdout,"%s> Getting coordinates single row (pseudo-atom model)\n", prog);

					if (f==0)
						mol->coordMatrix( &coord );
					else {   // just update coord instead initializes
						//fprintf(stderr, "---> %d %d %d  %d\n", props[ifr].k1, num_atoms_loop +3, props[ilr].k1,props[ilr+1].k1-1 );

						update_loop_coords(itermol, props[ifr].k1, num_atoms_loop +3, coord);
					}


					trd *der; // Derivatives
					if(verb > 1)
						fprintf(stdout,"%s> Computing derivatives...\n", prog);
					der = drdqC5x(coord, props, ifr, ilr, num_atoms_loop, size, model);

					if(verb > 1)
						fprintf(stdout, "%s> Computing Kinetic Energy matrix (masses matrix)...\n", prog);
					mass_matrix = kineticC5x(der, masses, props, ifr, num_atoms_loop, size, nco);

					// Updating Elastic network
					if(decint != NULL)
						free(decint);

					if( !(decint = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )  // Required for "realloc"
					{
						fprintf(stdout,"Sorry, \"decint\" memory allocation failed!\n");
						exit(1);
					}

					switch(potential)
					{
					case 0: // INVERSE EXPONENTIAL (power of distance for contact matrix)
					{
						// Making Interacting Pair of (non-virtual) Atoms (ipas)
						make_ipas_loop(itermol, itermol2, ifa, num_atoms_loop, cutoff_k0, &decint, &nipa);
						// fprintf( stdout, "%s> Inverse Exponential (%d nipas) cutoff= %.1f, k= %f, x0= %.1f ", prog, nipa, cutoff_k0, cte_k0, x0);
						for(int i=0; i<nipa; i++)
							decint[i].C = Inv_exp( cte_k0, decint[i].d, x0, power); // setting Force Constants
						break;
					}

					case 1: // DISTANCE CUTOFF METHOD
					{
						// Making Interacting Pair of (non-virtual) Atoms (ipas)
						make_ipas_loop(itermol, itermol2, ifa, num_atoms_loop, cutoff_k1, &decint, &nipa);
						// fprintf( stdout, "%s> Cutoff Distance (%d nipas) cutoff= %.1f, k= %f ", prog, nipa, cutoff_k1, cte_k1);
						for(int i=0; i<nipa; i++)
							decint[i].C = cte_k1; // setting Force Constants
						break;
					}
					}

					// Prune contacts, just considering loop contacts
					int nipa2 = 0;
					switch(contact)
					{
					case 0: // loop vs. loop
					{
						if( !(decint2 = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )
						{
							fprintf(stdout,"Sorry, \"decint2\" memory allocation failed!\n");
							exit(1);
						}

						// Only store contacts that involve loop atoms
						for(int i=0; i<nipa; i++)
							if( (decint[i].k >= ifa && decint[i].k <= ila) && (decint[i].l >= ifa && decint[i].l <= ila) )
							{
								nipa2++; // Counts number of Interacting Pairs of Atoms
								decint2 = ( twid * ) realloc( decint2, nipa2 * sizeof( twid ) ); // resizes contact list-structure
								decint2[nipa2 - 1].k = decint[i].k; // k-pseudo-atom index (i-atom index)
								decint2[nipa2 - 1].l = decint[i].l; // l-pseudo-atom index (j-atom index)
								decint2[nipa2 - 1].d = decint[i].d; // set distance
								decint2[nipa2 - 1].C = decint[i].C; // force constant will be set in the future
							}
						// fprintf( stdout, "%s> %d initial contacts pruned to just %d\n", prog, nipa, nipa2);

						// Just store the requested contacts
						free(decint);

						decint = decint2;
						nipa = nipa2;
						break;
					}
					case 1:
						break;
					default:
						fprintf( stdout, "%s> Please, introduce a valid Contact method to continue!!!\n\nForcing exit!\n\n", prog);
						exit(1);
						break;
					}

					// IPAs checking
					if(verb > 2 ) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
						for(int i=0; i<nipa; i++)
							fprintf(stdout,"ipa %4d: k= %d  l= %d  d= %f  C= %f\n",i,decint[i].k,decint[i].l,decint[i].d,decint[i].C);


					// HESSIAN
					if(verb > 1)
						fprintf(stdout,"%s> Computing Hessian matrix (potential energy matrix)...\n", prog);
					hess_matrix = hessianC5x(coord, der, props, decint, nipa, ifr, num_atoms_loop, size, nco);
					//					if(scoring != 2) // "coord" required for Dihedrals RMSD
					//						free(coord);
					//hess_matrix = hessianFast(coord, coord, der, props, decint, nipa, ifr, ilr, size, nco, num_atoms);

					if(verb > 1)
					{
						// Show Hessian matrix
						show_matrix(hess_matrix, size + nco, "Hessian:", " %7.2f");
						// Show Kinetic Energy matrix
						show_matrix(mass_matrix, size + nco, "Kinetic:", " %7.0f");
					}
					// COMPUTING THE EIGENVECTORS AND EIGENVALUES
					info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig);
					free(mass_matrix);
					free(hess_matrix);



					// MON: check this, seems unnecessary...
					if(neig < nevec)
					{
						fprintf(stderr,"Warning more eigenvectors requested (%d) than available (%d)\n",nevec,neig);
						// exit(1);
					}

					// Some checking...
					if( info ) // if info != 0
					{
						fprintf(stderr,"\n%s> An error occurred in the matrix diagonalization: %d\n", prog, info);
						exit(1);
					}

					if(verb > 1)
					{
						fprintf( stdout, "%s> Eigensolver successfully finished!!! (neig=%d)\n", prog, neig);
						show_vector(stdout,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
						show_vectors(stdout,eigvect,size,neig,"Dumping Raw Eigenvectors:"," %5.2e");
					}

					// Scale eigenvectors so that the maximum value of the components is "maxang"
					if(maxang != 0.0)
					{
						scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
						if(verb > 1)
							show_vectors(stdout,eigvect,size,neig,"Dumping Scaled Eigenvectors:"," %5.2e");
					}

					// Compute the Cartesian eigenvectors from the Internal Coordinates eigenvectors
					double *cevec; // Cartesian eigenvectors

					if(scoring != 2) // Not angular scoring method
					{
						// MON: neig should be nevec, shoudn't it?
						cevec = ic2cart(eigvect, neig, der, size, num_atoms_loop + 3, masses_loop);
						free(der); // Free obsolete derivatives

						//			cevec = ic2cart(eigvect, neig, der, size, num_atoms_loop + 3);
					}

					if ((f==0)&&(verb > 1))  {
						mol->writePDB("ptraj.pdb"); // Write the morphed loop together with the complete PDB structure
						save_ptraj_modes("ptraj.evec", (num_atoms_loop + 3)*3, 0, 10, eigval, cevec, false);
						//show_vectors(stdout,cevec,ncomps,neig,"Dumping Raw CC Eigenvectors:", " %5.2e");

						//exit(1);
					}
					// show_vectors(stdout,cevec,ncomps,neig,"Dumping Raw CC Eigenvectors:", " %5.2e");

					// fprintf(stdout,"size= %d  nlrs= %d %d\n",size,nlrs, f);
					//					exit(0);

					// Store reference mode "refmode" to maintain the requested direction
					if(morph_switch) // Morphing
					{
						if(f==0) {
							rmsd0 = rmsd_loop(itermol, itertar, ifa, num_atoms_loop, ifa2); // store initial RMSD
							rmsd0_ncac = rmsd_loop_residue(itermol, itertar, ifr, ifr2, nlrs );
						}

						switch(scoring) // Scoring method
						{
						case 1: // Atomic RMSD scoring
						{
							//							// MON: testing
							//							loop_dihedrals(coord, props, ifr, ilr, size, model, type, &dhs); // Get dihedral angles arrays
							//							show_vector(stderr, dhs, size, "dhs: ", " %6.1f");
							//							loop_dihedrals(coord2, props2, ifr2, ilr2, size, model, type, &dhs2);
							//							show_vector(stderr, dhs2, size, "dhs2:", " %6.1f");
							//							// Element-wise difference between dihedral angle [deg] arrays "v" and "w" (d = v - w) taking into account rotation
							//							dihedrals_diff(dhs2, dhs, refmode, size); // Element-wise difference between vectors "v" and "w" (d = v - w)
							//							modref = vector_modulus( refmode, size ); // Reference vector modulus

							// Use the atomic coordinates "delta" vector as reference
							delta_loop(itermol, itertar, ifa, num_atoms_loop +3 , refmode, ifa2); // Compute "delta" vector between both conformations
							//show_vector(stderr, refmode, ncomps, "ref: ", " %6.1f");

							//fprintf(stderr,"%d %d\n",num_atoms_loop,ncomps);
							rmsd = vector_rmsd(refmode, size);
							modref = vector_modulus( refmode, ncomps ); // Reference vector modulus
							//fprintf(stderr,"dihedral_diff_module= %f  dihedral_rmsd= %f\n",modref,rmsd);

							break;
						}

						case 2: // Dihedral angles RMSD scoring
						{
							//							loop_dihedrals(itermol, ifr, nlrs, &dhs); // Get dihedral angles arrays
							//							loop_dihedrals(itertar, ifr2, nlrs, &dhs2);
							//							vector_diff(dhs, dhs2, refmode, nlrs*2); // Element-wise difference between vectors "v" and "w" (d = v - w)
							//							modref = vector_modulus( refmode, nlrs*2 ); // Reference vector modulus

							loop_dihedrals(coord, props, ifr, ilr, size, model, type, &dhs); // Get dihedral angles arrays
							show_vector(stderr, dhs, size, "dhs: ", " %6.1f");
							loop_dihedrals(coord2, props2, ifr2, ilr2, size, model, type, &dhs2);
							show_vector(stderr, dhs2, size, "dhs2:", " %6.1f");

							// vector_diff(dhs, dhs2, refmode, size); // Element-wise difference between vectors "v" and "w" (d = v - w)

							// Element-wise difference between dihedral angle [deg] arrays "v" and "w" (d = v - w) taking into account rotation
							dihedrals_diff(dhs2, dhs, refmode, size); // Element-wise difference between vectors "v" and "w" (d = v - w)
							show_vector(stderr, refmode, size, "ref: ", " %6.1f");

							modref = vector_modulus( refmode, size ); // Reference vector modulus
							break;
						}
						}
					}  // end morph switch
					else if(mr_switch) //  multiple mode following.
					{

						//				// Generate some random vector in CC
						//				for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
						//					refmode[ k ] = 2*rg->Random() - 1.0; // rg->Random() --> Output random float number in the interval 0 <= x < 1

						//				// Randomly choose some mode as reference
						//				int ranm = rg->IRandom(0, nevec-1); // Output random integer in the interval min <= x <= max
						//				for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
						//					refmode[ k ] = cevec[ranm * ncomps + k ];
						//				modref = vector_modulus( refmode, ncomps ); // Reference vector modulus

						// Use selected random combinations of mode as reference vector (only on the first iteration!)
						if(f==0)
						{
							// Reset "refmode"
							for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
								refmode[ k ] = 0.0;

							// Update "refmode"
							for(int n=0; n<nevec; n++)
							{
								// One random mode amplitude for each eigenvector
								double rand_factor = 2*rg->Random() - 1.0; // rg->Random() --> Output random float number in the interval 0 <= x < 1

								for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
									refmode[ k ] += rand_factor * cevec[n*ncomps + k ];
							}

							modref = vector_modulus( refmode, ncomps ); // Reference vector modulus
						}
					}
					else // Single mode following
					{

						// Use selected mode as reference vector (only on the first iteration!)
						if(f==0)
						{
							for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
								refmode[ k ] = cevec[imod*ncomps + k ];
							//fprintf(stdout,"hola2\n");

							modref = vector_modulus( refmode, ncomps ); // Reference vector modulus
							//fprintf(stdout,"hola3\n");

						}

					}


					if (verb > 0) {
						sprintf(text,"%s> %4d ", prog, f);
						sprintf(dummy,"");
					}

					// Compute modes amplitudes into "alpha" array
					double delta = 0.0;
					double ddot;
					switch(scoring) // Scoring method
					{
					case 1: // Atomic RMSD scoring
					{
						for(int n=0; n<nevec; n++)
						{
							ddot = dotprodnorm( refmode, cevec + n*ncomps, ncomps, modref );
							alpha[n] = ddot;
							vdelta[n] = ddot*ddot; // Compute vectors with delta components
							tdelta[n] += vdelta[n]; // Compute vectors with delta components
							delta += vdelta[n];
							if (verb > 0) {
								sprintf(dummy, " %7.5f", vdelta[n]);
								strcat(text, dummy);
							}
						}

						if ((f==0)&&(verb > 1)) {
							float vdump;
							float vdump2=0;

							double *cevec2  = (double *) malloc( ncomps * sizeof(double));


							for(int i=0; i<ncomps; i++)
								cevec2[i]=0;

							//show_cartmode(coord, refmode, props, ifr, num_atoms_loop + 3, "refmode.txt", 0);
							float sample0=100.0;

							double ga=0, siga=0;
							for(int j=0; j<sample0; j++) {

								// generate random eigenvectors
								for(int n=0; n<nevec; n++)
									for(int i=0; i<ncomps-6; i++)
										cevec2[i]=rg->Random()-0.5;

								// calculate random gammas
								double ga0=0;
								for(int n=0; n<nevec; n++) {
									ddot = dotprodnorm( refmode, cevec2, ncomps, modref );
									ga0 += ddot*ddot;
								}
								// stats
								ga += ga0;
								siga += ga0*ga0;
								// fprintf( stdout, "%d %7.3f\n", j, ga0);

							}

							ga/=sample0;
							siga = sqrt(siga/sample0 - ga*ga );
							fprintf(stdout, " avg %7.3f sig %7.3f %.0f\n", ga,siga,sample0);

							//show_cartmode(coord, cevec2, props, ifr, num_atoms_loop + 3, "evec.txt", 0);

							//getchar();
							// for(int n=0; n<nevec; n++) {
							// ddot = dotprodnorm( refmode, cevec + n*ncomps, ncomps, modref );
							// / Compute vectors with delta components


							for(int n=0; n<nevec; n++) {
								vdump2+=1/eigval[n];
							}
							fprintf( f_log, "%s> Variance  ", prog );
							for(int n=0; n<nevec; n++) {
								fprintf( f_log, " %7.3f", 1/eigval[n]/vdump2);
							}
							fprintf( f_log, " \n");

							fprintf( stdout, "%s> Overlap %7.3f Modes ", prog, delta );
							fprintf( f_log, "%s> Overlap %7.3f Modes ", prog, delta );

							for(int n=0; n<nevec; n++) {
								// fprintf( stdout, "% 7.3f",vdelta[n]);
								vdeltaS[n]=vdelta[n];
							}
							// fprintf( stdout, "\n");
							//							for(int i=0; i<nevec; i++)
							//							 for(int j=i+1; j<nevec; j++)
							//								 if (vdeltaS[i]<vdeltaS[j]) {
							//									 vdump=vdeltaS[i];
							//									 vdeltaS[i]=vdeltaS[j];
							//									 vdeltaS[j]=vdump;
							//								 }
							vdump=0;
							vdump2=0;
							float variance = 0;
							for(int n=0; n<nevec; n++) {

								//fprintf( stdout, " variance %7.3f %d",variance, n);

								fprintf( stdout, " %7.3f",vdeltaS[n]);
								fprintf( f_log, " %7.3f",vdeltaS[n]);
								//fprintf( f_log, " %7.3f",vdeltaS[n]-vdelta2[n])/(sq_sum[n]);

								vdump+=vdeltaS[n];
								//							   if ((variance+=1/eigval[n]/vdump2) > 0.98) {
								//								   fprintf( stdout, " var %7.3f %d ",variance-1/eigval[n]/vdump2, n);
								//								   fprintf( f_log, " var %7.3f %d ",variance-1/eigval[n]/vdump2, n);
								//								   break;
								//
								//							   }
							}
							fprintf( stdout, " sum %7.3f Z %7.3f\n", vdump, (vdump-ga)/(siga));
							fprintf( f_log, " sum %7.3f Z %7.3f\n", vdump,  (vdump-ga)/(siga));


						}

						//show_vector(stderr, vdelta, nevec, "", " %7.5f", false, false); //
						break;
					}

					case 2: // Dihedral angles RMSD scoring
					{
						for(int n=0; n<nevec; n++)
						{
							// show_vector(stderr, eigvect + n*size, size, "evn:", " %6.3f");
							ddot = dotprodnorm( refmode, eigvect + n*size, size, modref );
							alpha[n] = ddot;
							delta += ddot*ddot;
							if (verb > 0) {
								sprintf(dummy, " %7.5f", ddot*ddot);
								strcat(text, dummy);
							}
						}
						if (verb > 2)
							show_vector(stderr, alpha, nevec, "alph", " %6.3f");
						break;
					}
					}

					//					fprintf(stderr,"nevec/2 = %d  nevec= %d\n",nevec/2, nevec);
					//					show_vector(stderr, vdelta, nevec, "", " %7.5f", false, false); //

					// WARNING: these "delta" elements are not "Square-rooted" as "delta"
					delta = sqrt(delta);
					if(f==0)
						delta0 = delta; // store initial delta

					if (verb > 0) {
						// Total Delta
						sprintf(dummy, "  %4.2f %4.2f", sum_vector(vdelta, nevec/2), (sum_vector(vdelta + (nevec/2), nevec - (nevec/2) )));
						strcat(text, dummy);
						sprintf(dummy, "  %7.5f", delta); // Total Delta
						strcat(text, dummy);
					}



					// Show Cartesian normal mode in VMD
					// sprintf(text, "%s_mode%02d.vmd", name, imod+1);
					// show_cartmode(coord, cevec, props, ifr, num_atoms_loop + 3, text, imod);

					// Initialize current mode
					for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
						mode[ k ] = 0.0; // zero initialization

					for(int n=0; n<nevec; n++)
					{
						double ampn = pow(alpha[n],2); // Amplitude of n-th mode
						// fprintf(stderr, "n= %d  ampn= %f\n", n, ampn);

						if(alpha[n] < 0.0) // Reverse n-th mode
							ampn *= -1.0; // Reverse current CC mode to maintain initial direction

						// Multiple modes merging
						if(linear_switch) // linear motion
						{
							// Generate "merged" Cartesian mode for motion
							for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
								mode[ k ] += cevec[ n * ncomps + k ] * ampn;
						}
						else // dihedral motion
						{
							// Generate "merged" Dihedral-coordinates mode for motion
							for(int k = 0; k < size; k++) // Mon added the +3, watch out!
								mode[ k ] += eigvect[ n * size + k ] * ampn;
						}
					}

					// Generate trajectory
					if(linear_switch) // linear motion
					{
						if(maxang != 0.0)
							scale_vectors(mode, 3 * (num_atoms_loop + 3), 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode

						if (true) {
							move_loop_linear(itermol, mode, 0, ifa, num_atoms_loop + 3); // including Ct-anchor (+3) for checking purposes...
						}
						else  {

							// apply linear in steps
							int num_steps=1000;

							move_loop_linear_factor(itermol, mode, 0, ifa, num_atoms_loop + 3,1.0/num_steps); // including Ct-anchor (+3) for checking purposes...

						    update_loop_coords(itermol, props[ifr].k1, num_atoms_loop+3, coord);

							if (f==0)
								mol->writeMloop("linear_traj.pdb", 1, ifr-1, ilr+1, chain);

							for(int s = 2; s < num_steps; s++)
							{
								der = drdqC5x(coord, props, ifr, ilr, num_atoms_loop, size, model);
								free(cevec);
								cevec = ic2cart(eigvect, neig, der, size, num_atoms_loop + 3, masses_loop);
								free(der); // Free obsolete derivatives
								// normalization
								for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
									mode[ k ] = 0.0; // zero initialization
								for(int n=0; n<nevec; n++)
								{
									double ampn = pow(alpha[n],2);
									if(alpha[n] < 0.0) // Reverse n-th mode
										ampn *= -1.0; // Reverse current CC mode to maintain initial direction
									for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
										mode[ k ] += cevec[ n * ncomps + k ] * ampn;
								}
								scale_vectors(mode, 3 * (num_atoms_loop + 3), 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode
								move_loop_linear_factor(itermol, mode, 0, ifa, num_atoms_loop +3 ,1.0/num_steps); // including Ct-anchor (+3) for checking purposes...
								//update_loop_coords(itermol, props[ifr].k1, num_atoms_loop +3 , coord);
								if (f==0)
							 	mol->writeMloop("linear_traj.pdb", s, ifr-1, ilr+1, chain);

							}
						}
						// move_loop_linear_steps(file_movie, mol, mode, 0, ifa, num_atoms_loop + 3, 1, 1.0); // including Ct-anchor (+3) for checking purposes...
					}
					else // dihedral motion
					{
						if(maxang != 0.0)
							scale_vectors(mode, size, 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode

						move_loop_dihedral(itermol, ifr, ilr, props, mode, size, model, 1.0);
						// mol->writeMPDB(file_movie, f+1); // Dump current conformation (frame) to a Multi-PDB file

						// fprintf(stdout,"follow_mode> Structure dumped into Muli-PDB dRMSD = %8f > %8f\n", rmsd - last_rmsd, delta_rmsd);
						//						mol->writeMloop(file_movie, f+1, ifr-1, ilr+1, chain);
					}

					// Does some loop atom clash with its environment?


					// Compute anchor drift
					anchor_drift(iterini, itermol, props, ilr, &adist, &aang);
					ddrift = adist - adist0; // Distance increment wrt. initial distance
					adrift = aang - aang0; // Angle increment wrt. initial angle
					if (verb > 0) {
						sprintf(dummy, " %6.3f %5.2f", ddrift, adrift);
						strcat(text, dummy);
					}

					if(morph_switch) // Morphing protocol
					{
						switch(scoring)
						{
						case 1: // Cartesian coordinates RMSD
							rmsd = rmsd_loop(itermol, itertar, ifa, num_atoms_loop, ifa2);
							break;

						case 2: // Dihedral angles RMSD
							// rmsd = rmsd_dihedral_loop(itermol, itertar, ifr, ifr2, ilr-ifr+1, dhs, dhs2);
							// rmsd = vector_rmsd(dhs, dhs2, size);
							rmsd = vector_rmsd(refmode, size);
							// fprintf(stderr,"Dihedral_RMSD= %6.2f\n", rmsd);
							break;
						}


						if(verb > 0) {
							//sprintf(dummy, " %d", clashed_loop( itermol, itermol2, ifa, num_atoms_loop, 1.0, ifa2) );
							//strcat(text, dummy);
							sprintf(dummy, " %7.4f", rmsd);
							strcat(text, dummy);
							fprintf(stdout, "%s\n", text); // Dump all output for current frame
							fprintf(f_log, "%s\n", text); // Dump log info

						}


						// Morphing Convergence Test
						if(rmsd_old - rmsd < rmsd_conv)
						{

							sprintf(dummy, "%s> Morphing convergence reached! dRMSD = %8f < %8f\n", prog, rmsd_old-rmsd, rmsd_conv);
							fprintf(f_log, "%s", dummy); // Dump log info
							fprintf(stdout,"%s",dummy);

							break;
						}
					}
					else // Not-morphing protocol
					{
						rmsd = rmsd_loop(iterini, itermol, ifa, num_atoms_loop);


						if(verb > 0) {
							sprintf(dummy, " %7.4f", rmsd);
							strcat(text, dummy);
							fprintf(stdout,"%s\n", text); // Dump all output for current frame
							fprintf(f_log, "%s\n", text); // Dump log info

						}


						// Convergence test
						//				if(!mr_switch && (rmsd - rmsd_old < rmsd_conv || clashed_loop( itermol, itermol2, ifa, num_atoms_loop, 1.0)))

						if(!mr_switch && rmsd - rmsd_old < rmsd_conv)
						{
							sprintf(dummy, "%s> Motion convergence reached! dRMSD = %8f < %8f\n", prog, rmsd-rmsd_old, rmsd_conv);
							fprintf(f_log, "%s", dummy); // Dump log info
							if (verb > 0) fprintf(stdout,"%s",dummy);
							break;
						}

						if(rmsd > target_rmsd)
						{
							sprintf(dummy, "%s> Motion convergence reached! RMSD = %8f > %8f target_rmsd\n", prog, rmsd, target_rmsd);
							fprintf(f_log, "%s", dummy); // Dump log info
							if (verb > 0) fprintf(stdout,"%s",dummy);
							break;
						}
					}

					if(f==0) // only the first time
						last_rmsd = rmsd;

					// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
					// if(delta_rmsd < rmsd - last_rmsd)
					if(delta_rmsd < fabsf(rmsd - last_rmsd) )
					{
						if(verb > 0)
							fprintf( stdout, "%s> Structure dumped into Muli-PDB dRMSD = %8f > %8f\n", prog, rmsd - last_rmsd, delta_rmsd);
						mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);
						last_rmsd = rmsd; // Keep last saved RMSD
					}

					rmsd_old = rmsd;

					if(scoring == 1) // Cartesian RMSD scoring
						free(cevec); // free obsolete Cartesian modes
				}

				if(delta_rmsd != 999999)
					fprintf( stdout, "%s> Saving %s and %s dRMSD = %8f > %8f\n", prog, file_movie, file_final, fabsf(rmsd - last_rmsd), delta_rmsd);
				else
					fprintf( stdout, "%s> Final structure dumped into Muli-PDBs dRMSD = %8f\n", prog, rmsd );

				mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);

				if(morph_switch) // Morphing protocol
				{
					mol->writePDB(file_final); // Write the morphed loop together with the complete PDB structure
					float rmsd_ncac = rmsd_loop_residue(itermol, itertar, ifr, ifr2, nlrs );

					if (verb>1) {
						scale_vectors(tdelta, nevec, 1, 1.0);
						// Dump "tdelta"
						fprintf(stdout, "%s> %5s", prog, "Total"); // Dump log info
						show_vector(stdout, tdelta, nevec, "", " %7.5f", false, true);
						fprintf(f_log, "%s> %5s", prog, "Total"); // Dump log info
						show_vector(f_log, tdelta, nevec, "", " %7.5f", false, true);
					}
					sprintf(dummy, "%s> Morphing:  Initial_RMSD= %-8.2f NCAC %-8.2f Final_RMSD= %-8.2f NCAC %-8.2f Delta_RMSD= %-8.2f  Initial_Delta= %-8.2f  Motion= %-8.2f\n",
							prog, rmsd0, rmsd0_ncac, rmsd, rmsd_ncac, rmsd0 - rmsd, delta0, (rmsd0 - rmsd)/ rmsd0);
					fprintf(f_log, "%s", dummy); // Dump log info
					fprintf(stdout,"%s",dummy);
					free(coord2);
				}
				free(coord);

			}
			break;


			//
			// CASE 2   Mode-following walks from initial random merged-modes.
			//

			case 2:
			{
				trd *der;

				// Compute the eigenvectors/values for some macromolecular loop (Required to define the initial "refmode" each iteration)
				neig = nma_loop(mol, model, type, props, masses, ifa, ifr, ilr, num_atoms_loop, size, nco, cutoff_k0, eigval, eigvect, &der);
				//fprintf(stderr,"model= %d  type= %d  ifa= %d  ifr= %d  ilr= %d  na= %d  size= %d  nco= %d  cutoff= %f  neig= %d\n", model, type, ifa, ifr, ilr, num_atoms_loop, size, nco, cutoff_k0, neig);
				// exit(0);

				if(verb > 1)
					fprintf( stdout, "%s> Eigensolver successfully finished!!! (neig=%d)\n", prog, neig);
				// show_vector(stdout,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
				// show_vectors(stdout,eigvect,size,neig,"Dumping Raw IC Eigenvectors:", " %5.2e");

				// Scale eigenvectors so that the maximum value of the components is "maxang"
				if(maxang != 0.0)
				{
					scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
					if(verb > 1)
						show_vectors(stdout,eigvect,size,neig,"Dumping Scaled Eigenvectors:", " %5.2e");
				}

				// Compute the Cartesian eigenvectors from the Internal Coordinates eigenvectors
				double *cevec; // Cartesian eigenvectors
				//				masses_loop = NULL;
				cevec = ic2cart(eigvect, neig, der, size, num_atoms_loop + 3, masses_loop);
				free(der);

				//show_vectors(stdout,cevec,ncomps,neig,"Dumping Raw CC Eigenvectors:", " %5.2e");

				double *xnm; // Coordinate in Normal Modal space
				xnm = (double *) malloc( sizeof(double) * neig); // Allocate NM coordinate
				for(int i=0; i<neig; i++)
					xnm[i] = 0.0; // Initialize NM coordinate

				// Get initial loop coordinates
				float *coordini = (float *) malloc( sizeof(float) * num_atoms_loop * 3 ); // for future copy & paste loop coordinates
				get_loop_coords(iterini, ifa, num_atoms_loop, coordini);

				// Write just the indicated loop (from Nt anchor (ifr-1) to Ct anchor (ilr+1)) into a Multi-PDB
				int fi = 1; // Frame index
				mol->writeMloop(file_movie, fi++, ifr-1, ilr+1, chain); // MON: use parser option...

				// srand(1867);

				int indx[neig];
				for(int i=0; i<neig; i++) indx[i]=i;

				double sdrift=0, sadrift=0;

				for(int f = 0; f < nsamples; f++) // Generate N-samples (frames)
				{
					// Reset reference mode

					if (max_loops_save!=0) {
				    if (fi>=max_loops_save) break;
					if ((f==nsamples-1)&& (fi<max_loops_save)) nsamples++;
					}

					for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
						refmode[ k ] = 0.0;


					random_shuffle(&indx[0],&indx[neig]);

					int ranm = rg->IRandom(1, neig);

					//					fprintf( stdout, "%s> %d %d\n", prog, neig, ranm);
					//					for(int i=0; i<ranm; i++)
					//					fprintf( stdout, " %d", indx[i]);
					//					getchar();

					// ranm=neig;

					// Generate a random direction (in NM coordinates) based on initial conformation Normal Modes
					double sumar=0.0;
					while (sumar<=0.01) {
						sumar=0;
						for(int i=0; i<ranm; i++) {
							double randa=2*rg->Random() - 1.0;
							//double randa=rg->Random();
							xnm[indx[i]] = randa;
							sumar+=fabs(randa);
						}
					}

					//					for(int i=0; i<neig; i++)
					//					{
					//                        if (xnm[i]==0) {
					//						if (rg->Random()<0.5) xnm[i]=0.05;
					//						else  xnm[i]=-0.05;
					//                        } else if (xnm[i]<0) xnm[i]-=0.05;
					//                        	else  xnm[i]+=0.05;
					//					}


					//					double sense;
					//					if (rg->Random()>=0.5) sense=1.0;
					//					else sense=-1.0;

					for(int i=0; i<ranm; i++)
					{
						// Generate a Reference Mode (refmode) from the NM coordinate and the Cartesian eigenvectors
						for(int k = 6; k < ncomps-9; k++) // Mon added the +3, watch out!
							refmode[ k ] += xnm[indx[i]] * cevec[indx[i]*ncomps + k ];

					}

					//					show_vectors(stdout,cevec,ncomps,neig,"Dumping Raw CC Eigenvectors:", " %5.2e");
					//					getchar();

					//					for(int i=0; i<neig; i++)
					//					{
					//						// Generate a random direction (in NM coordinates) based on initial conformation Normal Modes
					//						xnm[i] = 2*rg->Random() - 1.0; // Random double in [-1,1) interval --> rg->Random() outputs a random float number in the interval 0 <= x < 1
					//
					//						//if (xnm[i]<0) xnm[i]-=0.05;
					//						//else  xnm[i]+=0.05;
					//
					//						// Generate a Reference Mode (refmode) from the NM coordinate and the Cartesian eigenvectors
					//						for(int k = 6; k < ncomps-9; k++) // Mon added the +3, watch out!
					//							refmode[ k ] += xnm[i] * cevec[i*ncomps + k ];
					//
					//					}





					if (verb > 0)
						show_vector(stdout,xnm,neig,"xnm:"," %6.3f");
					//show_vector(stdout, refmode+12, 12, "refmode0:", " %5.2f", false, true);

					// // Show refmode
					// for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
					//	fprintf(stderr," %5.2f",refmode[k]);
					//	fprintf(stderr,"\n");

					// MON: check the "1000" below... use some PARSER input???

					// Mode following routine with target RMSD

					double target_rmsd_ef;
					target_rmsd_ef= target_rmsd;
					// target_rmsd_ef=rg->Random()*target_rmsd;
					// if (target_rmsd_ef<delta_rmsd) target_rmsd_ef=delta_rmsd;

					follow_mode(refmode, mol, model, type, props, masses, masses_loop, ifa, ifr, ilr, num_atoms_loop, size, nco, maxang, cutoff_k0, nsteps, max_loops_save,
							eigval, eigvect, target_rmsd_ef, iterini, file_movie, delta_rmsd, &fi, chain, false);

					// mol->writeMPDB(file_movie, f+2); // Dump current conformation (frame) to a Multi-PDB file

					// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
					// mol->writeMloop(file_movie, f+2, ifr-1, ilr+1, chain);

					anchor_drift(iterini, itermol, props, ilr, &adist, &aang);
					sdrift += fabs(adist - adist0); // Distance increment wrt. initial distance
					sadrift += fabs(aang - aang0); // Angle increment wrt. initial angle


					// Set initial loop coordinates to prevent unwanted distortions
					set_loop_coords(itermol, ifa, num_atoms_loop, coordini);
					//get_loop_coords(itermol, ifa, num_atoms_loop, coordini);
					//set_loop_coords(iterini, ifa, num_atoms_loop, coordini);


				}

				fprintf( stdout, "%s> Drifts  Distance %5.3f Angle %5.3f\n", prog, sdrift/nsamples, sadrift/nsamples);

				if (delta_rmsd==999999)
					fprintf( stdout, "%s> Saved %d conformation in %s\n", prog, nsamples, file_movie);
				else
					fprintf( stdout, "%s> Saved %d conformation in %d runs in %s every %f Δrmsd\n", prog, fi, nsamples, file_movie, delta_rmsd);


				free(coordini);
				free(cevec);
				free(xnm);

			}
			break;
			//
			// CASE 3   MC
			//

			case 3:
			{
				trd *der;

				traji=0;
						if(verb > 1)
							fprintf(stdout,"> Initial structure dumped into Muli-PDB\n");
						mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);

				anchor_drift(iterini, itermol, props, ilr, &adist0, &aang0);
				fprintf( stdout, "%s> Initial anchor distance and angle: %f A and %f deg\n", prog, adist0, aang0);


				// Compute the eigenvectors/values for some macromolecular loop (Required to define the initial "refmode" each iteration)
				neig = nma_loop(mol, model, type, props, masses, ifa, ifr, ilr, num_atoms_loop, size, nco, cutoff_k0, eigval, eigvect, &der);
				//fprintf(stderr,"model= %d  type= %d  ifa= %d  ifr= %d  ilr= %d  na= %d  size= %d  nco= %d  cutoff= %f  neig= %d\n", model, type, ifa, ifr, ilr, num_atoms_loop, size, nco, cutoff_k0, neig);
				// exit(0);

				if(verb > 1)
					fprintf( stdout, "%s> Eigensolver successfully finished!!! (neig=%d)\n", prog, neig);
				// show_vector(stdout,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
				// show_vectors(stdout,eigvect,size,neig,"Dumping Raw IC Eigenvectors:", " %5.2e");

				// Scale eigenvectors so that the maximum value of the components is "maxang"
				if(maxang != 0.0)
				{
					scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
					if(verb > 1)
						show_vectors(stdout,eigvect,size,neig,"Dumping Scaled Eigenvectors:", " %5.2e");
				}

				// Compute the Cartesian eigenvectors from the Internal Coordinates eigenvectors
				double *cevec; // Cartesian eigenvectors
				//				masses_loop = NULL;
				cevec = ic2cart(eigvect, neig, der, size, num_atoms_loop + 3, masses_loop);
				free(der);

				//show_vectors(stdout,cevec,ncomps,neig,"Dumping Raw CC Eigenvectors:", " %5.2e");

				double *xnm; // Coordinate in Normal Modal space
				xnm = (double *) malloc( sizeof(double) * neig); // Allocate NM coordinate
				for(int i=0; i<neig; i++)
					xnm[i] = 0.0; // Initialize NM coordinate

				// Get initial loop coordinates
				float *coordini0 = (float *) malloc( sizeof(float) * num_atoms_loop * 3 ); // for future copy & paste loop coordinates
				float *coordini = (float *) malloc( sizeof(float) * num_atoms_loop * 3 ); // for future copy & paste loop coordinates

				get_loop_coords(iterini, ifa, num_atoms_loop, coordini0);
				get_loop_coords(iterini, ifa, num_atoms_loop, coordini);

				// Write just the indicated loop (from Nt anchor (ifr-1) to Ct anchor (ilr+1)) into a Multi-PDB
				int fi = 0; // Frame index
				//mol->writeMloop(file_movie, fi++, ifr-1, ilr+1, chain); // MON: use parser option...

				// srand(1867);

				int indx[neig];
				for(int i=0; i<neig; i++) indx[i]=i;



				// Allocating scaling factor array (normal-mode space)
				double *scaling = (double *) malloc( sizeof(double) * size );
				// Allocating absolute conformation (normal-mode space)
				double *xconf = (double *) malloc( sizeof(double) * size );
				// Allocating Energy (normal-mode space)
				double *ener = (double *) malloc( sizeof(double) * size );

				// Initialization

				for(int i=0; i< size; i++)
				{
					xconf[i] = 0.0;
					ener[i] = 0.0;
				}



				double  Ef=100;
				double rfactor = 3.889087297*2; // 7.778174594
				double Kb = 0.00198717; // Kboltz (in Angstroms)
				double Ta = 300;
				double KbT = Kb * Ta; // setting thermal-energy



				for(int i=0; i< size; i++)

				{
					eigval[i] *= Ef; // Applying first the user-introduced Energy-factor (another scaling term)
					scaling[i] = rfactor * sqrt( KbT/eigval[i] ); // Computing scale factor
					// scaling[i] = rfactor * sqrt( KbT/1.0 ); // Computing scale factor

				}



				// Scale eigenvectors so that the maximum value of the components is "maxang"

				if(maxang != 0.0)

				{
					scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
				//	show_vectors(stdout,eigvect,size,neig,"Dumping Scaled Eigenvectors:", " %5.2e");
				}








				// sampling

				int accepted=0;

				double enertot=0;
				if(maxang != 0.0)
					scale_vectors(mode, size, 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode


				get_loop_coords(itermol, ifa, num_atoms_loop, coordini);


				for(int f = 0; f < nsamples; f++) // Generate N-samples (frames)

				{

					int sel_mode = rg->IRandom(0, neig-2); // [0:1) Playing dice with Mersenne!

					// Choosing random mode
					// Choosing random displacement
					double dx = ( (double) rg->Random() ) - 0.5; // random displacement (-0.5:0.5) (in normal-mode space)
					// Scaling displacement increment (anisotropic step)

					double ds = scaling[sel_mode] * dx; // "ds" is the scaled amplitude of the motion (in normal-mode space)
					// Energy increment due to "ds" in "mode"
					double enertry = 0.5 * eigval[sel_mode] * pow(xconf[sel_mode] + ds,2) - ener[sel_mode];


					//fprintf(stdout, "sel_mode %d\n", sel_mode);
					// Metropolis et al. acceptance test
					if( enertry < 0 || exp(-enertry/KbT) > (double) rg->Random() )

					{


						accepted++;
						enertot += enertry; // updating total energy
						ener[sel_mode] += enertry; // updating mode energy

				        // Initialize current "merged" mode

						xconf[sel_mode] += ds; // updating conformation (in normal-mode space)

						for(int k = 0; k < size; k++) // Mon added the +3, watch out!
							mode[ k ] =  eigvect[sel_mode*size + k ]*xconf[k];; // zero initialization


						 // move_loop_dihedral(itermol, ifr, ilr, props, mode, neig, model, 1.0);

						 set_loop_coords(iterini2, ifa, num_atoms_loop, coordini);

						 follow_mode(mode, molini2, model, type, props, masses, masses_loop, ifa, ifr, ilr, num_atoms_loop, size, nco, maxang, cutoff_k0, nsteps, max_loops_save,
						 												eigval, eigvect, delta_rmsd, itermol, file_movie, 10000, &fi, chain, false);




						rmsd = rmsd_loop(iterini, iterini2, ifa, num_atoms_loop);
						anchor_drift(iterini, iterini2, props, ilr, &adist, &aang);
						ddrift = adist - adist0; // Distance increment wrt. initial distance
						adrift = aang - aang0; // Angle increment wrt. initial angle


						if (rmsd<target_rmsd  && fabs(ddrift)<0.02 && fabs(adrift)<5) {
							fprintf(stdout, "val  %d rmsd %7.4f drift %6.3f %5.2f\n", f, rmsd, ddrift, adrift);
							if(f % 500 == 0)
							mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);

							get_loop_coords(iterini2, ifa, num_atoms_loop, coordini);
						    set_loop_coords(itermol, ifa, num_atoms_loop, coordini);




						} else {
							//fprintf(stdout, "fail %d rmsd %7.4f drift %6.3f %5.2f\n", f, rmsd, ddrift, adrift);
							set_loop_coords(itermol, ifa, num_atoms_loop, coordini);
							//xconf[sel_mode] -= ds;
						}


					}
				}

















				get_loop_coords(itermol, ifa, num_atoms_loop, coordini);


				for(int f = 0; f < nsamples*0.0; f++) // Generate N-samples (frames)
				{
					// Reset reference mode
					for(int k = 0; k < ncomps; k++) // Mon added the +3, watch out!
						refmode[ k ] = 0.0;


					random_shuffle(&indx[0],&indx[neig]);

					int ranm = rg->IRandom(1, neig-2);


					// Generate a random direction (in NM coordinates) based on initial conformation Normal Modes
					double sumar=0.0;
					while (sumar<=0.01) {
						sumar=0;
						for(int i=0; i<ranm; i++) {
							double randa=2*rg->Random() - 1.0;
							//double randa=rg->Random();
							xnm[indx[i]] = randa;
							sumar+=fabs(randa);
						}
					}

					//					for(int i=0; i<neig; i++)
					//					{
					//                        if (xnm[i]==0) {
					//						if (rg->Random()<0.5) xnm[i]=0.05;
					//						else  xnm[i]=-0.05;
					//                        } else if (xnm[i]<0) xnm[i]-=0.05;
					//                        	else  xnm[i]+=0.05;
					//					}


					//					double sense;
					//					if (rg->Random()>=0.5) sense=1.0;
					//					else sense=-1.0;
					for(int i=0; i<1; i++)
					{
						// Generate a Reference Mode (refmode) from the NM coordinate and the Cartesian eigenvectors
						for(int k = 6; k < ncomps-9; k++) // Mon added the +3, watch out!
							refmode[ k ] += xnm[indx[i]] * cevec[indx[i]*ncomps + k ];

					}




					double target_rmsd_ef;
					target_rmsd_ef= 0.1;
					// target_rmsd_ef=rg->Random()*target_rmsd;
					// if (target_rmsd_ef<delta_rmsd) target_rmsd_ef=delta_rmsd;

					set_loop_coords(iterini2, ifa, num_atoms_loop, coordini);

					follow_mode(refmode, molini2, model, type, props, masses, masses_loop, ifa, ifr, ilr, num_atoms_loop, size, nco, maxang, cutoff_k0, nsteps, max_loops_save,
												eigval, eigvect, delta_rmsd, itermol, file_movie, 10000, &fi, chain, false);

					// mol->writeMPDB(file_movie, f+2); // Dump current conformation (frame) to a Multi-PDB file

					// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
					// mol->writeMloop(file_movie, f+2, ifr-1, ilr+1, chain);

					// Set initial loop coordinates to prevent unwanted distortions
					//set_loop_coords(itermol, ifa, num_atoms_loop, coordini);

					anchor_drift(iterini, iterini2, props, ilr, &adist, &aang);
					ddrift = adist - adist0; // Distance increment wrt. initial distance
					adrift = aang - aang0; // Angle increment wrt. initial angle

					rmsd = rmsd_loop(iterini, iterini2, ifa, num_atoms_loop);


					if (fabs(ddrift)<0.01 && fabs(adrift)<5) {
					get_loop_coords(iterini2, ifa, num_atoms_loop, coordini);
					set_loop_coords(itermol, ifa, num_atoms_loop, coordini);
					mol->writeMloop(file_movie, (traji++), ifr-1, ilr+1, chain);

					fprintf(stdout, "val  %d rmsd %7.4f drift %6.3f %5.2f\n", f, rmsd, ddrift, adrift);

					} else {
						set_loop_coords(itermol, ifa, num_atoms_loop, coordini);

					}





				}

					fprintf( stdout, "%s> Saved %d conformation in %s\n", prog, traji, file_movie);
					fprintf( stdout, "%s> Saved %d conformation in %d runs in %s every %f Δrmsd\n", prog, fi, traji, file_movie, delta_rmsd);


				free(coordini);
				free(cevec);
				free(xnm);

			}
			break;

			default:
				fprintf(stderr,"ERROR. Please, introduce a valid sampling \"strategy\" (strategy= %d).Forcing exit!\n", strategy);
				exit(1);
			}
		}

		if(loops_switch)
			delete iter_loops;


		free(coordx);
		free(eigvect);
		free(eigval);

		delete iterini;
		delete itermol;
		delete itermol2;

		if(morph_switch)
			delete itertar;

		delete molini;


		delete molr;

		//		if(morph_switch)
		//		{
		//			for(int k=0; k<nevec; k++)
		//
		//			sprintf(dummy, "%s> Morphing: Initial_RMSD= %8f  Final_RMSD= %8f  Delta_RMSD= %8f\n", prog, rmsd0, rmsd, rmsd0 - rmsd);
		//			fprintf(f_log, "%s", dummy); // Dump log info
		//			fprintf(stdout,"%s",dummy);
		//			free(coord2);
		//		}

		fclose(f_log); // close log file
	}





	auto t2 = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = t2-t1;


	fprintf(stdout,"%s> Processed in %.3f seconds.\n%s>\n", prog, elapsed_seconds.count(), prog);


	return 0;
}

// Compute just the loop dihedral angles considered in Loop-NMA
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  props --> Residue properties array of the whole macromolecule
//  ifr   --> Index of First loop Residue (internal numeration)
//  ilr   --> Index of Last loop Residue (internal numeration)
//  size  --> Number of dihedral angles (the mobile variables)
//  model --> Coarse-Grained (CG) model
//  type  --> Chi type (0= No-chi, 2= Phi,Chi,Psi)
// OUTPUT:
//	**p_dhs --> Pointer to dihedral angles array according to current CG model (automated memory allocation if *NULL)
void loop_dihedrals(float *coord, tri *props, int ifr, int ilr, int size, int model, int type, float **p_dhs)
{
	bool verb = false;
	float *dhs; // Dihedrals array
	Tcoor atN, atCA, atC, atCB, atNn, atCp, atXG; // atNn (Next Nt atom), atCp (Previous Ct atom), atXg (atom at Gamma position)

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)

	int CBindex = 4;
	if(model == 1) // cg_3BB2R exchanges O (position 3) by CB (position 4) for convenience, i.e. CB is now at 3 and O at 4 positions.
		CBindex = 3;

	if(verb)
		fprintf(stderr,"loop_dihedrals> ");

	// Allocate dihedrals array
	if(*p_dhs == NULL)
		*p_dhs = (float *) malloc( size * sizeof(float) );
	dhs = *p_dhs;

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue

	for(i = ifr-1; i <= ifr+reglen; i++) // Screen loop residues, i.e. from Nt-anchor residue to Ct-anchor, inclusive.
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1;

		// fprintf(stderr,"i= %4d  k1= %4d\n",i,k1);

		if(i == ifr-1) // Previous-to-Nt residue (Nt anchor)
		{
			for(c=0; c<3; c++)
			{
				atN[c] = coord[ 3 * k1 + c]; // k1+0 --> N atom
				atCA[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom
				atCp[c] = atC[c] = coord[ 3 * (k1 + 2) + c]; // k1+2 --> C atom (and C from previous residue at next iteration)
			}
		}
		else
		{
			for(c=0; c<3; c++)
				atNn[c] = coord[ 3 * k1 + c]; // k1+0 --> N atom

			// Compute current Psi angle
			if(i != ifr + reglen) // if not Ct-anchor, i.e. if it has mobile PSI
			{
				dhs[ j ] = dihedral_bk( atN, atCA, atC, atNn );
				if(verb)
					fprintf(stderr," %5.1f:Psi",dhs[j]);
				j++;
			}

			for(c=0; c<3; c++)
			{
				atN[c] = atNn[c]; // k1+0 --> N atom
				atCA[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom
				atC[c] = coord[ 3 * (k1 + 2) + c]; // k1+2 --> C atom (and C from previous residue at next iteration)
			}

			// Compute current Phi angle
			if(props[i].nan != 1) // NOT Proline, i.e. if it has mobile PHI
			{
				dhs[ j ] = dihedral_bk( atCp, atN, atCA, atC );
				if(verb)
					fprintf(stderr," %5.1f:Phi",dhs[j]);
				j++;
			}

			for(c=0; c<3; c++)
				atCp[c] = atC[c]; // Copy C into Cp

			// Compute current Chi-1 angle
			if(type == 2 && i != ifr+reglen && props[i].nan == 3) // If it has CHI and it is not the Ct anchor
			{
				for(c=0; c<3; c++)
				{
					atCB[c] = coord[ 3 * (k1 + CBindex) + c]; // k1+CBindex --> CB atom
					atXG[c] = coord[ 3 * (k1 + CBindex + 1) + c]; // k1+CBindex+1 --> Gamma atom
				}

				dhs[ j ] = dihedral_bk( atN, atCA, atCB, atXG );
				if(verb)
					fprintf(stderr," %5.1f:Chi",dhs[j]);
				j++;
			}
		}
	}

	if(verb)
		fprintf(stderr,"\nloop_dihedrals> Number of dihedral variables (Phi,Psi,Chi) %d\n", j);
}


// Compute the derivatives (der) of Cartesian coordinates (r) wrt dihedral angles (q) for loops:  der = dr/dq
// OUTPUT:
//  Returns the derivatives matrix (automatic memory allocation)
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  props --> Residue properties array of the whole macromolecule
//  ifr   --> Index of First loop Residue (internal numeration)
//  ilr   --> Index of Last loop Residue (internal numeration)
//  nla   --> Number of Loop Atoms (only mobile atoms)
//  size  --> Number of dihedral angles (the mobile variables)
//  model --> Coarse-Grained model
trd *drdqC5x(float *coord, tri *props, int ifr, int ilr, int nla, int size, int model)
{
	bool verb = false;
	double e[3]; // rotation axis (e_lambda)
	double y[3]; // center of rotation
	double r[3]; // position of current atom
	double temp; // dummy variable

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)
	int sized = nla + 3; // Size of Derivatives (number of atoms with derivatives), "+3" includes the 3 Ct-anchor atoms for constraints.
	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	int ilpa = props[ilr+1].k1 - 1; // Index of last (mobile) pseudo-atom of loop

	int CBindex = 4;
	if(model == 1)
		CBindex = 3;

	if(verb)
		fprintf(stdout,"drdqC5x> ifpa: %d  ilpa: %d  CBindex: %d\n",ifpa,ilpa,CBindex);

	// Allocate derivatives
	trd *der = (trd *) malloc( sized * size * sizeof(trd) );
	ptr_check(der);

	// Initialize derivatives (MON: consider removal)
	for(int i=0; i < sized*size; i++)
	{
		der[i].x=0.0;
		der[i].y=0.0;
		der[i].z=0.0;
	}

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue
	int k2 = 0; // Index of last atom of current residue

	for(i = ifr; i <= ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor, inclusive.
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1;
		k2 = k1 + props[i].nat-1;

		if(verb)
			fprintf(stdout,"drdqC5x> i=%d  nan=%d  nat=%d  k1= %d  k2= %d\n",i,props[i].nan,props[i].nat,k1,k2);

		// Derivatives of PHI dihedral angles
		// ----------------------------------
		if(props[i].nan != 1) // NOT Proline, i.e. if it has just one mobile angle (PHI).
		{   // q_j = phi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> PHI\n",i,j);

			// Compute:   e = r_CA - r_N    and    y = r_N
			for(c=0; c<3; c++)
			{
				// Get:  y = r_N
				y[c] = coord[ 3 * k1 + c]; // k1+0 --> N atom

				// Compute:  e = r_CA - r_N   (vector N-->CA)
				e[c] = coord[ 3 * (k1 + 1) + c] - y[c]; // k1+1 --> CA atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
			for(k = ifpa, kp = 0; k < k1 + 2; k++, kp++) // All before current residue plus the N or CA pseudo-atoms of current residue
			{
				der[kp*size+j].x = 0.0;
				der[kp*size+j].y = 0.0;
				der[kp*size+j].z = 0.0;

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
			for(k = k1 + 2; k <= ilpa + 3; k++, kp++) // All atoms after the CA pseudo-atom of current residue, up to the first 3 atoms of Ct-anchor
			{
				for(c=0; c<3; c++)
					r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

				// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
				der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
				der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
				der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			j++; // Update current dihedral variable index
		}

		// Derivatives of CHI dihedral angles
		// ----------------------------------
		if(i != ifr+reglen && props[i].nan == 3) // If it has CHI and it is not the Ct anchor
		{   // q_j = chi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> CHI\n",i,j);

			// Compute:   e = r_CB - r_CA    and    y = r_CA
			for(c=0; c<3; c++)
			{
				// Get:  y = r_CA
				y[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom

				// Compute:  e = r_CB - r_CA   (vector CB-->C)
				e[c] = coord[ 3 * (k1 + CBindex) + c] - y[c]; // k1+3 --> CB atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			for(k = ifpa, kp = 0; k <= ilpa; k++, kp++) // All current and previous residue atoms do not move with current PSI
			{
				if( k >= k1 + CBindex && k <= k2 ) // Only move the Side-chain atoms
				{
					// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
					for(c=0; c<3; c++)
						r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

					// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
					der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
					der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
					der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
					// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

					if(verb)
						fprintf(stdout,"i= %3d  j= %3d  CHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f  Side-chain!\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
				}
				else
				{
					// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
					der[kp*size+j].x = 0.0;
					der[kp*size+j].y = 0.0;
					der[kp*size+j].z = 0.0;

					if(verb)
						fprintf(stdout,"i= %3d  j= %3d  CHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
				}
			}

			j++; // Update current dihedral variable index
		}

		// Derivatives of PSI dihedral angles
		// ----------------------------------
		if(i != ifr + reglen) // if not Ct-anchor, i.e. if it has mobile PSI angle
		{   // q_j = psi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> PSI\n",i,j);

			// Compute:   e = r_C - r_CA    and    y = r_CA
			for(c=0; c<3; c++)
			{
				// Get:  y = r_CA
				y[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom

				// Compute:  e = r_C - r_CA   (vector CA-->C)
				e[c] = coord[ 3 * (k1 + 2) + c] - y[c]; // k1+2 --> C atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
			for(k = ifpa, kp = 0; k <= k2; k++, kp++) // All current and previous residue atoms do not move with current PSI
			{
				if(k == k1 + 3) // if current residue Oxygen atom
				{
					// Compute derivatives for current residue Oxygen, it belongs to body 2 (mobile) --> dr/dq = e x (r-y)
					for(c=0; c<3; c++)
						r[c] = coord[ 3 * k + c ]; // current pseudo-atom position
					// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
					der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
					der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
					der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				}
				else
				{   // Body 1 does not move...
					der[kp*size+j].x = 0.0;
					der[kp*size+j].y = 0.0;
					der[kp*size+j].z = 0.0;
				}

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PSI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
			for(k = k2+1; k <= ilpa + 3; k++, kp++) // All atoms after current current residue, including the first 3 atoms of Ct-anchor
			{
				for(c=0; c<3; c++)
					r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

				// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
				der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
				der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
				der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PSI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			j++; // Update current dihedral variable index
		}

	}

	if(verb)
		fprintf(stdout,"Number of dihedral variables (Phi,Psi,Chi) %d\n", j);

	return der; // Return derivatives
}

inline void drdqC5x(float *coord, tri *props, int ifr, int ilr, int nla, int size, int model, trd *der)
{
	bool verb = false;
	double e[3]; // rotation axis (e_lambda)
	double y[3]; // center of rotation
	double r[3]; // position of current atom
	double temp; // dummy variable

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)
	int sized = nla + 3; // Size of Derivatives (number of atoms with derivatives), "+3" includes the 3 Ct-anchor atoms for constraints.
	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	int ilpa = props[ilr+1].k1 - 1; // Index of last (mobile) pseudo-atom of loop

	int CBindex = 4;
	if(model == 1)
		CBindex = 3;

	if(verb)
		fprintf(stdout,"drdqC5x> ifpa: %d  ilpa: %d  CBindex: %d\n",ifpa,ilpa,CBindex);

	// Allocate derivatives
	// trd *der = (trd *) malloc( sized * size * sizeof(trd) );
	// ptr_check(der);

	// Initialize derivatives (MON: consider removal)
	for(int i=0; i < sized*size; i++)
	{
		der[i].x=0.0;
		der[i].y=0.0;
		der[i].z=0.0;
	}

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue
	int k2 = 0; // Index of last atom of current residue

	for(i = ifr; i <= ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor, inclusive.
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1;
		k2 = k1 + props[i].nat-1;

		if(verb)
			fprintf(stdout,"drdqC5x> i=%d  nan=%d  nat=%d  k1= %d  k2= %d\n",i,props[i].nan,props[i].nat,k1,k2);

		// Derivatives of PHI dihedral angles
		// ----------------------------------
		if(props[i].nan != 1) // NOT Proline, i.e. if it has just one mobile angle (PHI).
		{   // q_j = phi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> PHI\n",i,j);

			// Compute:   e = r_CA - r_N    and    y = r_N
			for(c=0; c<3; c++)
			{
				// Get:  y = r_N
				y[c] = coord[ 3 * k1 + c]; // k1+0 --> N atom

				// Compute:  e = r_CA - r_N   (vector N-->CA)
				e[c] = coord[ 3 * (k1 + 1) + c] - y[c]; // k1+1 --> CA atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
			for(k = ifpa, kp = 0; k < k1 + 2; k++, kp++) // All before current residue plus the N or CA pseudo-atoms of current residue
			{
				der[kp*size+j].x = 0.0;
				der[kp*size+j].y = 0.0;
				der[kp*size+j].z = 0.0;

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
			for(k = k1 + 2; k <= ilpa + 3; k++, kp++) // All atoms after the CA pseudo-atom of current residue, up to the first 3 atoms of Ct-anchor
			{
				for(c=0; c<3; c++)
					r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

				// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
				der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
				der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
				der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			j++; // Update current dihedral variable index
		}

		// Derivatives of CHI dihedral angles
		// ----------------------------------
		if(i != ifr+reglen && props[i].nan == 3) // If it has CHI and it is not the Ct anchor
		{   // q_j = chi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> CHI\n",i,j);

			// Compute:   e = r_CB - r_CA    and    y = r_CA
			for(c=0; c<3; c++)
			{
				// Get:  y = r_CA
				y[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom

				// Compute:  e = r_CB - r_CA   (vector CB-->C)
				e[c] = coord[ 3 * (k1 + CBindex) + c] - y[c]; // k1+3 --> CB atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			for(k = ifpa, kp = 0; k <= ilpa; k++, kp++) // All current and previous residue atoms do not move with current PSI
			{
				if( k >= k1 + CBindex && k <= k2 ) // Only move the Side-chain atoms
				{
					// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
					for(c=0; c<3; c++)
						r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

					// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
					der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
					der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
					der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
					// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

					if(verb)
						fprintf(stdout,"i= %3d  j= %3d  CHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f  Side-chain!\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
				}
				else
				{
					// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
					der[kp*size+j].x = 0.0;
					der[kp*size+j].y = 0.0;
					der[kp*size+j].z = 0.0;

					if(verb)
						fprintf(stdout,"i= %3d  j= %3d  CHI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
				}
			}

			j++; // Update current dihedral variable index
		}

		// Derivatives of PSI dihedral angles
		// ----------------------------------
		if(i != ifr + reglen) // if not Ct-anchor, i.e. if it has mobile PSI angle
		{   // q_j = psi_i
			if(verb)
				fprintf(stdout,"drdqC5x> residue_i= %3d  dihedral_j= %3d  --> PSI\n",i,j);

			// Compute:   e = r_C - r_CA    and    y = r_CA
			for(c=0; c<3; c++)
			{
				// Get:  y = r_CA
				y[c] = coord[ 3 * (k1 + 1) + c]; // k1+1 --> CA atom

				// Compute:  e = r_C - r_CA   (vector CA-->C)
				e[c] = coord[ 3 * (k1 + 2) + c] - y[c]; // k1+2 --> C atom
			}

			// "e" normalization --> |e| = 1.0
			temp = sqrt( e[0]*e[0] + e[1]*e[1] + e[2]*e[2] );
			for(c=0; c<3; c++)
				e[c] /= temp;

			// Compute derivatives for body 1 (fixed) --> dr/dq = 0.0
			for(k = ifpa, kp = 0; k <= k2; k++, kp++) // All current and previous residue atoms do not move with current PSI
			{
				if(k == k1 + 3) // if current residue Oxygen atom
				{
					// Compute derivatives for current residue Oxygen, it belongs to body 2 (mobile) --> dr/dq = e x (r-y)
					for(c=0; c<3; c++)
						r[c] = coord[ 3 * k + c ]; // current pseudo-atom position
					// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
					der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
					der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
					der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				}
				else
				{   // Body 1 does not move...
					der[kp*size+j].x = 0.0;
					der[kp*size+j].y = 0.0;
					der[kp*size+j].z = 0.0;
				}

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PSI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			// Compute derivatives for body 2 (mobile) --> dr/dq = e x (r-y)
			for(k = k2+1; k <= ilpa + 3; k++, kp++) // All atoms after current current residue, including the first 3 atoms of Ct-anchor
			{
				for(c=0; c<3; c++)
					r[c] = coord[ 3 * k + c ]; // current pseudo-atom position

				// Derivative of the rotation of current atom ("r") around the axis "e" and center "y" of rotation: dr/dq = e x (r-y)
				der[kp*size + j].x = e[1] * (r[2]-y[2]) - e[2] * (r[1]-y[1]);
				der[kp*size + j].y = e[2] * (r[0]-y[0]) - e[0] * (r[2]-y[2]);
				der[kp*size + j].z = e[0] * (r[1]-y[1]) - e[1] * (r[0]-y[0]);
				// These derivatives do not consider Eckart conditions (they are not needed since the loop anchors are fixed in space)

				if(verb)
					fprintf(stdout,"i= %3d  j= %3d  PSI  k= %3d  kp= %3d  der= %6.3f %6.3f %6.3f\n",i,j,k,kp,der[kp*size+j].x,der[kp*size+j].y,der[kp*size+j].z);
			}

			j++; // Update current dihedral variable index
		}

	}

	if(verb)
		fprintf(stdout,"Number of dihedral variables (Phi,Psi,Chi) %d\n", j);

}

// Compute the Kinetic Energy matrix (masses matrix, M) for loops NMA
// OUTPUT:
//  Returns the Kinetic Energy matrix (automatic memory allocation)
// INPUT:
//  der    --> Derivatives
//  masses --> Array of masses, one mass per (pseudo)atom
//  props  --> Residue properties array of the whole macromolecule
//  ifr    --> Index of First loop Residue (internal numeration)
//  nla    --> Number of Loop Atoms (only mobile atoms)
//  size   --> Number of dihedral angles (the mobile variables)
//  nco    --> Number of Constraints
double *kineticC5x(trd *der, float *masses, tri *props, int ifr, int nla, int size, int nco)
{
	int i = 0; // i-th dihedral angle
	int isi = 0; // i-th dihedral angle (for matrix increment)
	int j = 0; // j-th dihedral angle
	int jsi = 0; // j-th dihedral angle (for matrix increment)
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int ks = 0; // Atom index in loop context (for matrix increment)
	int c = 0; // Coordinate index
	int sizex = size + nco;
	int dks, dls; // some indices
	double mak; // mass of atom k
	double mdx, mdy, mdz; // some buffers

	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop

	double *mass_matrix; // kinetic energy matrix
	mass_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(mass_matrix);
	//for(i=0; i < sizex*sizex; i++)
	//	mass_matrix[i] = 0.0; // Mon: initialize the Bordered Mass matrix

	memset(mass_matrix, 0, sizex*sizex*sizeof(double));

	// Build mass matrix
	for(kp=0, ks=0; kp<nla; kp++, ks+=size) // Screen just Loop atoms
	{
		k = kp + ifpa;
		mak = masses[k]; // get mass for current atom

		for(i=0,isi=0; i<size; i++,isi+=sizex) // Screen i-th dihedrals
		{
			dks=ks+i;
			mdx = mak * der[dks].x;
			mdy = mak * der[dks].y;
			mdz = mak * der[dks].z;

			for(j=i;j<size;j++) // Screen j-th dihedrals
			{   // upper triangular part (including diagonal)
				dls = ks+j;
				mass_matrix[isi+j] += mdx * der[dls].x + mdy * der[dls].y + mdz * der[dls].z;
			}
		}
	}

	// Fill in lower triangular part
	for(i=0,isi=0; i<size; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			mass_matrix[isi+j] = mass_matrix[jsi+i];

	return mass_matrix;
}
inline void kineticC5x(trd *der, float *masses, tri *props, int ifr, int nla, int size, int nco, double *mass_matrix)
{
	int i = 0; // i-th dihedral angle
	int isi = 0; // i-th dihedral angle (for matrix increment)
	int j = 0; // j-th dihedral angle
	int jsi = 0; // j-th dihedral angle (for matrix increment)
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int ks = 0; // Atom index in loop context (for matrix increment)
	int c = 0; // Coordinate index
	int sizex = size + nco;
	int dks, dls; // some indices
	double mak; // mass of atom k
	double mdx, mdy, mdz; // some buffers

	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop


	//for(i=0; i < sizex*sizex; i++)
	//	mass_matrix[i] = 0.0; // Mon: initialize the Bordered Mass matrix

	memset(mass_matrix, 0, sizex*sizex*sizeof(double));

	// Build mass matrix
	for(kp=0, ks=0; kp<nla; kp++, ks+=size) // Screen just Loop atoms
	{
		k = kp + ifpa;
		mak = masses[k]; // get mass for current atom

		for(i=0,isi=0; i<size; i++,isi+=sizex) // Screen i-th dihedrals
		{
			dks=ks+i;
			mdx = mak * der[dks].x;
			mdy = mak * der[dks].y;
			mdz = mak * der[dks].z;

			for(j=i;j<size;j++) // Screen j-th dihedrals
			{   // upper triangular part (including diagonal)
				dls = ks+j;
				mass_matrix[isi+j] += mdx * der[dls].x + mdy * der[dls].y + mdz * der[dls].z;
			}
		}
	}

	// Fill in lower triangular part
	for(i=0,isi=0; i<size; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			mass_matrix[isi+j] = mass_matrix[jsi+i];

}



// Compute the Hessian matrix (2nd derivatives of the potential energy, H) for loops NMA
// OUTPUT:
//  Returns the Heesian matrix (automatic memory allocation)
// INPUT:
//  coord --> Coordinates array of the whole macromolecule (simple array of floats). WARNING: All (pseudo)atoms must be sorted ("formated")
//  der    --> Derivatives
//  props  --> Residue properties array of the whole macromolecule
//  decint --> List of Interacting Pairs of Atoms (contacts list)
//  nipa   --> Number of Interacting Pairs of Atoms (number of contacts)
//  ifr    --> Index of First loop Residue (internal numeration)
//  nla    --> Number of Loop Atoms (only mobile atoms)
//  size   --> Number of dihedral angles (the mobile variables)
//  nco    --> Number of Constraints
double *hessianC5x(float *coord, trd *der, tri *props, twid *decint, int nipa, int ifr, int nla, int size, int nco)
{
	int i = 0; // i-th dihedral angle
	int isi = 0; // i-th dihedral angle (for matrix increment)
	int j = 0; // j-th dihedral angle
	int jsi = 0; // j-th dihedral angle (for matrix increment)
	int k = 0; // k-th atom index
	int l = 0; // l-th atom index

	double r[3]; // inter-atomic vector
	double v[3]; // some vector
	double w[3]; // some vector
	double wcv[3]; // some vector
	double vw[3]; // some vector
	double prod, prod1, sum; // some buffers

	int kp = 0; // k-th atom index in loop context
	int lp = 0; // l-th atom index in loop context

	int ks = 0; // k-th atom index in loop context (for matrix increment)
	int ls = 0; // l-th atom index in loop context (for matrix increment)

	int c = 0; // Coordinate index
	int sizex = size + nco;
	int dks, dls; // some indices

	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	int ilpa = ifpa + nla - 1; // Index of last (mobile) pseudo-atom of loop

	double *rdr = (double *) malloc(size * sizeof(double)); // Auxiliar array
	ptr_check(rdr);

	double *hess_matrix; // Bordered Hessian matrix
	hess_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(hess_matrix);

	//for(i=0; i < sizex*sizex; i++)
	//	hess_matrix[i] = 0.0; // Initializing the Bordered Hessian matrix

	memset(hess_matrix, 0, sizex*sizex*sizeof(double));

	// compute Hessian matrix
	for(int index=0; index<nipa; index++)
	{
		// prod1 = decint[index].C / pow(decint[index].d,2);
		// prod1 = decint[index].C / pow(decint[index].d,2);
		prod1 = decint[index].C / decint[index].d;

		k=decint[index].k;
		l=decint[index].l;
		kp=k-ifpa;
		lp=l-ifpa;
		ks=kp*size;
		ls=lp*size;

		// l-->k vector
		for(c=0; c<3; c++)
			r[c] = coord[3*k + c] - coord[3*l + c];

		if(k < ifpa || k > ilpa) // if "k" (pseudo)atom is outside the loop
		{
			// Only "l" (pseudo)atom is mobile
			for(i=0; i<size; i++)
			{
				dls = ls+i;
				rdr[i]= r[0]*(-der[dls].x) + r[1]*(-der[dls].y) + r[2]*(-der[dls].z);
			}
		}
		else if(l < ifpa || l > ilpa) // if "l" (pseudo)atom is outside the loop
		{
			// Only "k" (pseudo)atom is mobile
			for(i=0; i<size; i++)
			{
				dks = ks+i;
				rdr[i]= r[0]*(der[dks].x) + r[1]*(der[dks].y) + r[2]*(der[dks].z);
			}
		}
		else
		{
			// Both (pseudo)atoms "k" and "l" are mobile (both belong to the mobile loop)
			for(i=0; i<size; i++)
			{
				dks = ks+i;
				dls = ls+i;
				rdr[i]= r[0]*(der[dks].x-der[dls].x) + r[1]*(der[dks].y-der[dls].y) + r[2]*(der[dks].z-der[dls].z);
			}
		}

		for(i=0,isi=0; i<size; i++, isi+=sizex) // screen rows (i)
		{
			prod  = prod1 * rdr[i];

			// fill upper triangular part (including diagonal)
			for(j=i;j<size;j++) // screen cols (j)
				hess_matrix[isi+j] += prod * rdr[j];
		}
	}
	free(rdr);
	// ************* Add entries corresponding to constraints ******************/

	// CA (pseudo)atom index of the Ct-anchor
	k = nla + ifpa + 1;

	// v = r_N - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		v[c] = coord[3*(k-1) + c] - coord[3*k + c];
		prod += pow(v[c], 2);
	}
	prod = sqrt(prod); // norm of r_CA - r_N

	for(c=0; c<3; c++)
		v[c] /= prod;

	// w = r_C - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		w[c] = coord[3*(k+1) + c] - coord[3*k + c];
		prod += pow(w[c], 2);
	}
	prod = sqrt(prod); // norm of r_C - r_CA

	for(c=0; c<3; c++)
		w[c] /= prod;

	// cosg = cos(v^w)
	double cosg = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];     // cos(ang(v,w))

	// wcv = w - cosg*v  (difference between w and the projection of w into v; thus wcv is orthogonal to both v and vw, see next...)
	for(c=0; c<3; c++)          // w - cosg*v
		wcv[c] = w[c] - cosg*v[c];

	// vw = v x w
	vw[0] = v[1]*w[2] - v[2]*w[1];
	vw[1] = v[2]*w[0] - v[0]*w[2];
	vw[2] = v[0]*w[1] - v[1]*w[0];    /*  v x w  */

	for(i=0,isi=0; i<size; i++,isi+=sizex) // Mon: Screen "size" dihedrals
	{
		// Mon: CA motion of Ct-anchor is constrained
		hess_matrix[isi+size  ] = der[(nla+1)*size+i].x; // CA pseudo-atom of Ct-anchor
		hess_matrix[isi+size+1] = der[(nla+1)*size+i].y;
		hess_matrix[isi+size+2] = der[(nla+1)*size+i].z;
		// Mon: The effect of all derivatives on changing the wcv vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * wcv[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * wcv[1];
		sum += der[nla*size+i].z * wcv[2];
		hess_matrix[isi+size+3] = sum; // Mon: der * wcv
		// Mon: Constraint the change in the displacement of the N atom of Ct-anchor vw The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * vw[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * vw[1];
		sum += der[nla*size+i].z * vw[2];
		hess_matrix[isi+size+4] = sum; // Mon: der * vw
		// Mon: The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[(nla+2)*size+i].x * vw[0]; // CO pseudo-atom of Ct-anchor
		sum += der[(nla+2)*size+i].y * vw[1];
		sum += der[(nla+2)*size+i].z * vw[2];
		hess_matrix[isi+size+5] = sum;
	}
	/***************************************************************/

	/* fill in lower triangular part */
	for(i=0,isi=0; i<sizex; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			hess_matrix[isi+j] = hess_matrix[jsi+i];

	//	// Mon: Tip-effect....
	//	for(i=0,isi=0; i<size; i++,isi+=sizex)
	//		hess_matrix[isi+i] += 1e-3;
	//
	//	hess_matrix[sizex*0+0] += 1e3;
	//	hess_matrix[sizex*1+1] += 1e3;
	//	hess_matrix[sizex*(size-2)+size-2] += 1e3;
	//	hess_matrix[sizex*(size-1)+size-1] += 1e3;

	return hess_matrix;
}
// This is without matrix
inline void hessianC5x(float *coord, trd *der, tri *props, twid *decint, int nipa, int ifr, int nla, int size, int nco, double *rdr, double *hess_matrix )
{
	int i = 0; // i-th dihedral angle
	int isi = 0; // i-th dihedral angle (for matrix increment)
	int j = 0; // j-th dihedral angle
	int jsi = 0; // j-th dihedral angle (for matrix increment)
	int k = 0; // k-th atom index
	int l = 0; // l-th atom index

	double r[3]; // inter-atomic vector
	double v[3]; // some vector
	double w[3]; // some vector
	double wcv[3]; // some vector
	double vw[3]; // some vector
	double prod, prod1, sum; // some buffers

	int kp = 0; // k-th atom index in loop context
	int lp = 0; // l-th atom index in loop context

	int ks = 0; // k-th atom index in loop context (for matrix increment)
	int ls = 0; // l-th atom index in loop context (for matrix increment)

	int c = 0; // Coordinate index
	int sizex = size + nco;
	int dks, dls; // some indices

	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	int ilpa = ifpa + nla - 1; // Index of last (mobile) pseudo-atom of loop

	//for(i=0; i < sizex*sizex; i++)
	//	hess_matrix[i] = 0.0; // Initializing the Bordered Hessian matrix

	memset(hess_matrix, 0, sizex*sizex*sizeof(double));

	// compute Hessian matrix
	for(int index=0; index<nipa; index++)
	{
		// prod1 = decint[index].C / pow(decint[index].d,2);
		// ojo Pablo remove the sqrt of d
		prod1 = decint[index].C / decint[index].d;

		k=decint[index].k;
		l=decint[index].l;
		kp=k-ifpa;
		lp=l-ifpa;
		ks=kp*size;
		ls=lp*size;

		// l-->k vector
		for(int c=0; c<3; c++)
			r[c] = coord[3*k + c] - coord[3*l + c];

		if(k < ifpa || k > ilpa) // if "k" (pseudo)atom is outside the loop
		{
			// Only "l" (pseudo)atom is mobile
			for(i=0; i<size; i++)
			{
				dls = ls+i;
				rdr[i]= r[0]*(-der[dls].x) + r[1]*(-der[dls].y) + r[2]*(-der[dls].z);
			}
		}
		else if(l < ifpa || l > ilpa) // if "l" (pseudo)atom is outside the loop
		{
			// Only "k" (pseudo)atom is mobile
			for(i=0; i<size; i++)
			{
				dks = ks+i;
				rdr[i]= r[0]*(der[dks].x) + r[1]*(der[dks].y) + r[2]*(der[dks].z);
			}
		}
		else
		{
			// Both (pseudo)atoms "k" and "l" are mobile (both belong to the mobile loop)
			for(i=0; i<size; i++)
			{
				dks = ks+i;
				dls = ls+i;
				rdr[i]= r[0]*(der[dks].x-der[dls].x) + r[1]*(der[dks].y-der[dls].y) + r[2]*(der[dks].z-der[dls].z);
			}
		}

		for(i=0,isi=0; i<size; i++, isi+=sizex) // screen rows (i)
		{
			//if (rdr[i]!=0.0)
			{
				prod  = prod1 * rdr[i];

				// fill upper triangular part (including diagonal)
				for(j=i;j<size;j++) // screen cols (j)
					hess_matrix[isi+j] += prod * rdr[j];
			}
		}
	}
	// free(rdr);
	// ************* Add entries corresponding to constraints ******************/

	// CA (pseudo)atom index of the Ct-anchor
	k = nla + ifpa + 1;

	// v = r_N - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		v[c] = coord[3*(k-1) + c] - coord[3*k + c];
		prod += pow(v[c], 2);
	}
	prod = sqrt(prod); // norm of r_CA - r_N

	for(c=0; c<3; c++)
		v[c] /= prod;

	// w = r_C - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		w[c] = coord[3*(k+1) + c] - coord[3*k + c];
		prod += pow(w[c], 2);
	}
	prod = sqrt(prod); // norm of r_C - r_CA

	for(c=0; c<3; c++)
		w[c] /= prod;

	// cosg = cos(v^w)
	double cosg = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];     // cos(ang(v,w))

	// wcv = w - cosg*v  (difference between w and the projection of w into v; thus wcv is orthogonal to both v and vw, see next...)
	for(c=0; c<3; c++)          // w - cosg*v
		wcv[c] = w[c] - cosg*v[c];

	// vw = v x w
	vw[0] = v[1]*w[2] - v[2]*w[1];
	vw[1] = v[2]*w[0] - v[0]*w[2];
	vw[2] = v[0]*w[1] - v[1]*w[0];    /*  v x w  */

	for(i=0,isi=0; i<size; i++,isi+=sizex) // Mon: Screen "size" dihedrals
	{
		// Mon: CA motion of Ct-anchor is constrained
		hess_matrix[isi+size  ] = der[(nla+1)*size+i].x; // CA pseudo-atom of Ct-anchor
		hess_matrix[isi+size+1] = der[(nla+1)*size+i].y;
		hess_matrix[isi+size+2] = der[(nla+1)*size+i].z;
		// Mon: The effect of all derivatives on changing the wcv vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * wcv[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * wcv[1];
		sum += der[nla*size+i].z * wcv[2];
		hess_matrix[isi+size+3] = sum; // Mon: der * wcv
		// Mon: Constraint the change in the displacement of the N atom of Ct-anchor vw The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * vw[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * vw[1];
		sum += der[nla*size+i].z * vw[2];
		hess_matrix[isi+size+4] = sum; // Mon: der * vw
		// Mon: The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[(nla+2)*size+i].x * vw[0]; // CO pseudo-atom of Ct-anchor
		sum += der[(nla+2)*size+i].y * vw[1];
		sum += der[(nla+2)*size+i].z * vw[2];
		hess_matrix[isi+size+5] = sum;
	}
	/***************************************************************/

	/* fill in lower triangular part */
	for(i=0,isi=0; i<sizex; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			hess_matrix[isi+j] = hess_matrix[jsi+i];

	//	// Mon: Tip-effect....
	for(i=0,isi=0; i<size; i++,isi+=sizex)
		hess_matrix[isi+i] += 1e-3;

	hess_matrix[sizex*0+0] += 1e3;
	hess_matrix[sizex*1+1] += 1e3;
	hess_matrix[sizex*(size-2)+size-2] += 1e3;
	hess_matrix[sizex*(size-1)+size-1] += 1e3;


}


void hessianFast(float *coord, float *coordCA, trd *der, tri *props, twid *decint, int nipa, int ifr, int ilr, int size, int nco, int num_atoms, double *hess_matrix)
{
	bool debug = false;
	double prod,prod1;
	int l,ind2,ind3,ls,ks;
	double *dummy;
	int prin, fin, j2, m, buff;
	double r_alpha[3],r_beta[3],r[3],e[3];
	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
	double d,temp;
	double v[6];
	int resn,num_res,num_seg;
	int index_res = 0;

	if(debug)
		printf("debug> WARNING! Using: hessianFast() nipa %d ifr %d ilr %d size %d nco %d num_atoms %d\n", nipa, ifr, ilr, size, nco, num_atoms);

	// Allocate Bordered Hessian
	int sizex = size + nco;
	//	double *hess_matrix; // Bordered Hessian matrix
	//	if( !(hess_matrix = (double *) malloc( sizex * sizex * sizeof(double)) ) ) // Square matrix
	//	{
	//		printf("Msg(hessianFast): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
	//		exit(1);
	//	}

	memset(hess_matrix, 0, sizex*sizex*sizeof(double));

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)
	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	//	int ilpa = ifpa + nla - 1; // Index of last (mobile) pseudo-atom of loop
	int ilpa = props[ilr+1].k1 - 1; // Index of last (mobile) pseudo-atom of loop
	int nla = ilpa - ifpa + 1; // Number of loop atoms (mobile) Check!

	if(debug)
		printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	//	size--;
	//	printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	// ********************************************
	// Storing ==> (ea, ea x ra) (== (eb, eb x rb) )
	// ********************************************
	//		double **erx;
	//		erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	//		for(int i=0;i<6;i++)
	//		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
	//
	//		int *undh; // returns the closest unit-index on the left side of the dihedral
	//		undh = (int *) malloc( sizeof(int) * size );
	double erx[6][size];
	int undh[size];

	int unat[num_atoms];

	for(int x=0; x<num_atoms; x++)
		unat[x] = 0; // all atoms belong to 0 unit by default (i.e. environment)

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	int k0 = 0; // first-NH + CA + last-CO model index

	for(int i=0; i<ifr; i++)
		k0 += props[i].nat;

	if(debug)
		printf("Msg(hessianFast): k0= %d\n", k0);

	int num_units = 0;
	int un_index = 0; // un_index=0 --> environment rigid unit (i.e. the initialized value), N-CA belong to environment too (they do not move at all)
	int natom = 0;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary: "undh"
	// ---------------------------------------------------

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue
	int k2 = 0; // Index of last atom of current residue

	//	for(i = ifr; i <= ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor, inclusive.
	for(i = ifr; i < ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor (without Ct)
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1; // index of 1st atom in residue "i"
		k2 = k1 + props[i].nat-1;

		if(debug)
			fprintf(stdout,"hessianFast> residue i=%d  nan=%d  nat=%d  k1= %d  k2= %d\n",i,props[i].nan,props[i].nat,k1,k2);

		unat[k0] = un_index;  // the unit N atom belongs to
		k0++; // update atom counter

		// Vectors for PHI dihedral angles
		// ----------------------------------
		if(props[i].nan != 1) // NOT Proline, i.e. if it has a mobile angle PHI
		{   // q_j = phi_i
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PHI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];


			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			un_index++; // update if it has PHI

			j++; // Update current dihedral variable index
		}

		//		if(i > ifr)  // Not-first-residue (1st residue N-CA is a single rigid unit)
		//			un_index++;

		unat[k0] = un_index;  // the unit CA atom belongs to
		k0++; // update atom counter

		// Vectors for PSI dihedral angles
		// ----------------------------------

		// if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
		{
			// q_j = psi_i  // all aminoacids have mobile PSI angle
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PSI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];

			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			j++; // Update current dihedral variable index
		}

		if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
			un_index++;
		//		un_index++;

		unat[k0] = un_index; // the unit C atom belongs to
		k0++; // update atom counter

	}

	// MON: PHI at Ct anchor dihedral ????
	undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)


	if(debug)
		fprintf(stderr,"hessianFast> undh, unat, (ea, ea x ra) and (eb, eb x rb) initialized!\n");


	if(debug)
	{
		fprintf(stderr,"hessianFast> (ea, ea x ra)\n");
		for(int x=0; x<j; x++)
		{
			for(int y=0; y<6; y++)
				fprintf(stderr," %f", erx[y][x]);
			fprintf(stderr,"\n");
		}


		fprintf(stderr,"hessianFast> undh=");
		for(int x=0; x<size; x++)
			fprintf(stderr," [%d]=%d",x,undh[x]);

		fprintf(stderr,"\nhessianFast> unat=");
		for(int x=0; x<num_atoms; x++)
		{
			if(x % 32 == 0)
				fprintf(stderr,"\n");
			fprintf(stderr," %2d",unat[x]);
		}
		fprintf(stderr,"\n");
	}

	// ---------------------------------------------------------------------------
	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
	// ---------------------------------------------------------------------------
	// Tij related pre-computations (to get direct Tij computation)

	//	num_units = un_index + 1 + 1; // Inter-unit variables (rot-trans) dont increase units number
	num_units = un_index + 1; // total number of units --> M = size + 1

	if(debug)
		printf("Msg(hessianFast): Number of units: %ld (%d)\n",num_units,un_index);

	int unipa_index;
	int *p_int;
	// int **unipa;

	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
	// (triangular packing storage)
	// unipa = (int **) malloc( sizeof(int *) * num_units * ( num_units + 1) / 2 );

	// "unipa" initialization
	//	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
	//		unipa[i] = (int *) malloc( sizeof(int) * 500); ;
	//
	//	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
	//		unipa[i][0]=0;

	int unipa[num_units * ( num_units + 1) / 2 ][500];

	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
		unipa[i][0]=0;

	if(debug)
		printf("Msg(hessianFast): nipa=%d\n",nipa);

	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
	{

		k = decint[index].k; // k atom index
		l = decint[index].l; // l atom index
		i = unat[k]; // unit index of "k" atom
		j = unat[l]; // unit index of "l" atom

		if(debug)
			printf("Msg(hessianFast): unipa[%d] --> k= %d  l= %d  i= %d  j= %d\n",index,k,l,i,j);

		if(i != j) 	// the same unit is rigid, inner contacts must not be taken into account!
		{		// (as well as the non contacting pairs!)
			// The following must be always true:  (k < l) && (i < j)
			if( i > j )
			{
				buff = j;
				j = i;
				i = buff;
			}

			// translates from "squared" to "triangular" matrix elements
			unipa_index = i + j*(j+1)/2;



			//			// Allocating memory for the new element (memory efficient)
			//			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
			//			{
			//				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
			//				{
			//					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
			//					exit(1);
			//				}
			//				unipa[unipa_index] = p_int;
			//				unipa[unipa_index][0] = 2;
			//			}
			//			else // If already allocated, then increases its size one int.
			//			{
			//				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
			//				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
			//				{
			//					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
			//					exit(1);
			//				}
			//				unipa[unipa_index] = p_int;
			//			}

			if (unipa[unipa_index][0] == 0)
				unipa[unipa_index][0] = 1;

			unipa[unipa_index][0]++;

			// Storing ipa index for k,l interacting pair of atoms
			unipa[unipa_index][ unipa[unipa_index][0] - 1 ] = index;
		}
	}

	if(debug)
	{
		fprintf(stderr, "Msg(hessianFast): UNIPA computed!\n");
	}

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------

	// Hessian initialization
	for(i=0; i < sizex * sizex; i++)
		hess_matrix[i] = 0.0;

	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
	int uij_index = 0;
	int uij1, uij2, uij3;

	// MON: check "+1"
	// +1 due to void element (when there are just Traslational/Rotational DoFs)
	double ***Uij;
	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
	{
		for(i=0; i<(num_units*(num_units+1)/2); i++)
			Uij[i] = NULL; // initialization
	}
	else
	{
		printf("Msg(hessianFast): Unable to allocate Uij-matrix!\n");
		exit(1);
	}

	if(debug)
		fprintf(stderr,"Uij size = %d\n", num_units*(num_units+1)/2);

	// T memory allocation (6x6 matrix)
	//double **T;
	//	T = (double **) malloc( sizeof(double *) * 6 );
	//	for(int i=0; i<6; i++)
	//	T[i] = (double *) malloc( sizeof(double) * 6 );
	double T[6][6];

	if(debug)
		fprintf(stderr,"Uab Computation Main-Loop (Computing Tij elements on-the-fly)\n");

	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
	for(int b= size-3, c= size+5; b >= 0; b--,c--) 	// b --> dihedrals (column)
	{							// c --> deleteable col. index
		for(int a=0; a <= b; a++)	 		// a --> dihedrals (row)
		{ 						// it fills upper triangular part (including diagonal)
			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
			uij_index = undh[a] + ( undh[b] + 1 ) * ( undh[b] + 2 ) / 2;
			// uij_index = undh[a] + undh[b] * ( undh[b] + 1 ) / 2;

			if(debug)
				fprintf(stderr,"a= %d  b= %d  undh[a]= %d  undh[b]+1= %d  uij_index= %d\n", a, b, undh[a], undh[b]+1, uij_index);

			//			fprintf(stderr,"%ld ",uij_index);

			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
			// the Uij for the 4 ICs combination should be computed only once!
			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
			{
				// "in situ" Uij memory allocation
				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
				{
					for(int m=0; m<6; m++)
					{
						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
						{
							for(int n=0; n<6; n++)
								Uij[uij_index][m][n] = 0.0;
						}
						else
						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
					}
				}
				else
				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }

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

				if(debug)
					fprintf(stderr,"\tUij: a=%d b=%d uij1(row=%d col=%d) uij2(row=%d col=%d) uij3(row=%d col=%d) \n",a,b,undh[a],undh[b]+2,undh[a]-1,undh[b]+1,undh[a]-1,undh[b]+2);

				// Computing Uij element
				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
				{				// (whether they have ipas...) --> then, "T" must be computed
					// Computing Tij element "on the fly"
					// (ec.21) Noguti & Go 1983 pp.3685-90
					//	calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );


					double r_alpha[3];
					double r_beta[3];
					double v[6];

					for(int i=0; i<6; i++)
						for(int j=0; j<6; j++)
							T[i][j] = 0.0; // T initialization

					for(int i=1; i < unipa[uij_index][0]; i++) // screens inter-unit contacts (k vs. l)
					{
						int in=unipa[uij_index][i];
						int ki = 3*decint[in].k;
						int li = 3*decint[in].l;

						// k-atom position retrieval (alpha)
						r_alpha[0] = coordCA[ki];
						r_alpha[1] = coordCA[ki + 1];
						r_alpha[2] = coordCA[ki + 2];

						// l-atom position retrieval (beta)
						r_beta[0] = coordCA[li];
						r_beta[1] = coordCA[li + 1];
						r_beta[2] = coordCA[li + 2];

						// D-Matrix
						v[0] = r_alpha[1] * r_beta[2] - r_alpha[2] * r_beta[1];
						v[1] = r_alpha[2] * r_beta[0] - r_alpha[0] * r_beta[2];
						v[2] = r_alpha[0] * r_beta[1] - r_alpha[1] * r_beta[0];
						v[3] = r_alpha[0] - r_beta[0];
						v[4] = r_alpha[1] - r_beta[1];
						v[5] = r_alpha[2] - r_beta[2];
						//double d1 = decint[in].C / pow(decint[in].d,2);
						double d1 = decint[in].C / decint[in].d;

						//fprintf(stderr,"%f %f %f %d %d\n", d1, decint[in].C , decint[in].d, ki, li);
						//getchar();

						for(int n=0;n<6;n++) // Diadic product
							for(int m=0;m<6;m++)
								T[n][m] += d1 * v[n] * v[m]; // Storing S-matrix
					}




					if(debug)
						fprintf(stderr,"\tTij computed\n");

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

				// at this point, the Uij element is fully computed!

				// Show Uij[]
				if(debug)
				{
					fprintf(stderr,"Uij[%d] a=%d b=%d\n", uij_index, a, b);
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							fprintf(stderr," %f", Uij[uij_index][n][m]);
					fprintf(stderr,"\n");
				}

			}   	// only if not computed yet!

			// ****************************************************
			// Rab - Matrix Computation
			// Noguti & Go (1983) pp. 3685-90
			// ****************************************************
			// Rab Computation ( the same Single- or Multi- Chain)
			// Table I. Noguti & Go 1983 pp.3685-90
			// backbone vs. backbone --> Rab = Uab

			// (ec.16) Noguti & Go 1983 pp.3685-90
			// fprintf(stderr,"\tRij= ");
			for(int n=0; n<6; n++)
			{
				R[n] = 0.0;
				for(int m=0; m<6; m++)
					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
				// fprintf(stderr,"%f ",R[n]);
			}
			// fprintf(stderr,"\n");

			temp = 0.0;
			for(int n=0; n<6; n++)
				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

			if(debug)
				fprintf(stderr,"\tHab= %f\n",temp);

			// hess_matrix[a+b*(b+1)/2] = temp; // Rab' x erb = Hessian
			hess_matrix[ a * sizex + b ] = temp; // Rab' x erb = Hessian (square matrix upper triangular part including diagonal?)

		}

		/*
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
				}
			}
		}
		 */
	}


	/*
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
				}
			}
		}
	free( Uij );
	 */

	//	for(i=0; i<num_units*(num_units+1)/2; i++)
	//		free( unipa[i] );
	//	free( unipa );

	//	for(i=0;i<6;i++)
	//	{
	//		free( erx[i] ); // erx[6][size]
	//		free( T[i] );
	//	}
	//	free( erx );
	//	free( T );
	//	free( undh );



	// ************* Add entries corresponding to constraints ******************
	// double v[3]; // some vector
	double w[3]; // some vector
	double wcv[3]; // some vector
	double vw[3]; // some vector
	double sum; // some buffers

	// CA (pseudo)atom index of the Ct-anchor
	k = nla + ifpa + 1;

	// v = r_N - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		v[c] = coord[3*(k-1) + c] - coord[3*k + c];
		prod += pow(v[c], 2);
	}
	prod = sqrt(prod); // norm of r_CA - r_N

	for(c=0; c<3; c++)
		v[c] /= prod;

	// w = r_C - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		w[c] = coord[3*(k+1) + c] - coord[3*k + c];
		prod += pow(w[c], 2);
	}
	prod = sqrt(prod); // norm of r_C - r_CA

	for(c=0; c<3; c++)
		w[c] /= prod;

	// cosg = cos(v^w)
	double cosg = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];     // cos(ang(v,w))

	// wcv = w - cosg*v  (difference between w and the projection of w into v; thus wcv is orthogonal to both v and vw, see next...)
	for(c=0; c<3; c++)          // w - cosg*v
		wcv[c] = w[c] - cosg*v[c];

	// vw = v x w
	vw[0] = v[1]*w[2] - v[2]*w[1];
	vw[1] = v[2]*w[0] - v[0]*w[2];
	vw[2] = v[0]*w[1] - v[1]*w[0];    //  v x w


	int isi, jsi;
	for(i=0, isi=0; i<size; i++, isi += sizex) // Mon: Screen "size" dihedrals
	{
		// Mon: CA motion of Ct-anchor is constrained
		hess_matrix[isi+size  ] = der[(nla+1)*size+i].x; // CA pseudo-atom of Ct-anchor
		hess_matrix[isi+size+1] = der[(nla+1)*size+i].y;
		hess_matrix[isi+size+2] = der[(nla+1)*size+i].z;
		// Mon: The effect of all derivatives on changing the wcv vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * wcv[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * wcv[1];
		sum += der[nla*size+i].z * wcv[2];
		hess_matrix[isi+size+3] = sum; // Mon: der * wcv
		// Mon: Constraint the change in the displacement of the N atom of Ct-anchor vw The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * vw[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * vw[1];
		sum += der[nla*size+i].z * vw[2];
		hess_matrix[isi+size+4] = sum; // Mon: der * vw
		// Mon: The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[(nla+2)*size+i].x * vw[0]; // CO pseudo-atom of Ct-anchor
		sum += der[(nla+2)*size+i].y * vw[1];
		sum += der[(nla+2)*size+i].z * vw[2];
		hess_matrix[isi+size+5] = sum;
	}
	//***************************************************************

	// fill in lower triangular part
	for(i=0,isi=0; i<sizex; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			hess_matrix[isi+j] = hess_matrix[jsi+i];

}

double *hessianFast_mon(float *coord, float *coordCA, trd *der, tri *props, twid *decint, int nipa, int ifr, int ilr, int size, int nco, int num_atoms)
{
	bool debug = false;
	double prod,prod1;
	int l,ind2,ind3,ls,ks;
	double *dummy;
	int prin, fin, j2, m, buff;
	double r_alpha[3],r_beta[3],r[3],e[3];
	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
	double d,temp;
	double v[6];
	int resn,num_res,num_seg;
	int index_res = 0;

	if(debug)
		printf("debug> WARNING! Using: hessianFast()\n");

	// Allocate Bordered Hessian
	int sizex = size + nco;
	double *hess_matrix; // Bordered Hessian matrix
	if( !(hess_matrix = (double *) malloc( sizex * sizex * sizeof(double)) ) ) // Square matrix
	{
		printf("Msg(hessianFast): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)
	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	//	int ilpa = ifpa + nla - 1; // Index of last (mobile) pseudo-atom of loop
	int ilpa = props[ilr+1].k1 - 1; // Index of last (mobile) pseudo-atom of loop
	int nla = ilpa - ifpa + 1; // Number of loop atoms (mobile) Check!

	if(debug)
		printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	//	size--;
	//	printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	// ********************************************
	// Storing ==> (ea, ea x ra) (== (eb, eb x rb) )
	// ********************************************
	double **erx;
	erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	for(int i=0;i<6;i++)
		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]

	int *undh; // returns the closest unit-index on the left side of the dihedral
	undh = (int *) malloc( sizeof(int) * size );
	int unat[num_atoms];

	for(int x=0; x<num_atoms; x++)
		unat[x] = 0; // all atoms belong to 0 unit by default (i.e. environment)

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	int k0 = 0; // first-NH + CA + last-CO model index

	for(int i=0; i<ifr; i++)
		k0 += props[i].nat;

	if(debug)
		printf("Msg(hessianFast): k0= %d\n", k0);

	int num_units = 0;
	int un_index = 0; // un_index=0 --> environment rigid unit (i.e. the initialized value), N-CA belong to environment too (they do not move at all)
	int natom = 0;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary: "undh"
	// ---------------------------------------------------

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue
	int k2 = 0; // Index of last atom of current residue

	//	for(i = ifr; i <= ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor, inclusive.
	for(i = ifr; i < ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor (without Ct)
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1; // index of 1st atom in residue "i"
		k2 = k1 + props[i].nat-1;

		if(debug)
			fprintf(stdout,"hessianFast> residue i=%d  nan=%d  nat=%d  k1= %d  k2= %d\n",i,props[i].nan,props[i].nat,k1,k2);

		unat[k0] = un_index;  // the unit N atom belongs to
		k0++; // update atom counter

		// Vectors for PHI dihedral angles
		// ----------------------------------
		if(props[i].nan != 1) // NOT Proline, i.e. if it has a mobile angle PHI
		{   // q_j = phi_i
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PHI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];


			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			un_index++; // update if it has PHI

			j++; // Update current dihedral variable index
		}

		//		if(i > ifr)  // Not-first-residue (1st residue N-CA is a single rigid unit)
		//			un_index++;

		unat[k0] = un_index;  // the unit CA atom belongs to
		k0++; // update atom counter

		// Vectors for PSI dihedral angles
		// ----------------------------------

		// if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
		{
			// q_j = psi_i  // all aminoacids have mobile PSI angle
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PSI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];

			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			j++; // Update current dihedral variable index
		}

		if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
			un_index++;
		//		un_index++;

		unat[k0] = un_index; // the unit C atom belongs to
		k0++; // update atom counter

	}

	// MON: PHI at Ct anchor dihedral ????
	undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)


	if(debug)
		fprintf(stderr,"hessianFast> undh, unat, (ea, ea x ra) and (eb, eb x rb) initialized!\n");


	if(debug)
	{
		fprintf(stderr,"hessianFast> (ea, ea x ra)\n");
		for(int x=0; x<j; x++)
		{
			for(int y=0; y<6; y++)
				fprintf(stderr," %f", erx[y][x]);
			fprintf(stderr,"\n");
		}


		fprintf(stderr,"hessianFast> undh=");
		for(int x=0; x<size; x++)
			fprintf(stderr," [%d]=%d",x,undh[x]);

		fprintf(stderr,"\nhessianFast> unat=");
		for(int x=0; x<num_atoms; x++)
		{
			if(x % 32 == 0)
				fprintf(stderr,"\n");
			fprintf(stderr," %2d",unat[x]);
		}
		fprintf(stderr,"\n");
	}

	// ---------------------------------------------------------------------------
	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
	// ---------------------------------------------------------------------------
	// Tij related pre-computations (to get direct Tij computation)

	//	num_units = un_index + 1 + 1; // Inter-unit variables (rot-trans) dont increase units number
	num_units = un_index + 1; // total number of units --> M = size + 1

	if(debug)
		printf("Msg(hessianFast): Number of units: %ld (%d)\n",num_units,un_index);

	int unipa_index;
	int *p_int;
	int **unipa;

	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
	// (triangular packing storage)
	unipa = (int **) malloc( sizeof(int *) * num_units * ( num_units + 1) / 2 );

	// "unipa" initialization
	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
		unipa[i] = NULL;

	if(debug)
		printf("Msg(hessianFast): nipa=%d\n",nipa);

	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
	{

		k = decint[index].k; // k atom index
		l = decint[index].l; // l atom index
		i = unat[k]; // unit index of "k" atom
		j = unat[l]; // unit index of "l" atom

		if(debug)
			printf("Msg(hessianFast): unipa[%d] --> k= %d  l= %d  i= %d  j= %d\n",index,k,l,i,j);

		if(i != j) 	// the same unit is rigid, inner contacts must not be taken into account!
		{		// (as well as the non contacting pairs!)
			// The following must be always true:  (k < l) && (i < j)
			if( i > j )
			{
				buff = j;
				j = i;
				i = buff;
			}

			// translates from "squared" to "triangular" matrix elements
			unipa_index = i + j*(j+1)/2;

			// Allocating memory for the new element (memory efficient)
			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
			{
				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
				{
					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
				unipa[unipa_index][0] = 2;
			}
			else // If already allocated, then increases its size one int.
			{
				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
				{
					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
					exit(1);
				}
				unipa[unipa_index] = p_int;
			}

			// Storing ipa index for k,l interacting pair of atoms
			unipa[unipa_index][ unipa[unipa_index][0] - 1 ] = index;
		}
	}

	if(debug)
	{
		fprintf(stderr, "Msg(hessianFast): UNIPA computed!\n");
	}

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------

	// Hessian initialization
	for(i=0; i < sizex * sizex; i++)
		hess_matrix[i] = 0.0;

	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
	int uij_index = 0;
	int uij1, uij2, uij3;

	// MON: check "+1"
	// +1 due to void element (when there are just Traslational/Rotational DoFs)
	double ***Uij;
	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
	{
		for(i=0; i<(num_units*(num_units+1)/2); i++)
			Uij[i] = NULL; // initialization
	}
	else
	{
		printf("Msg(hessianFast): Unable to allocate Uij-matrix!\n");
		exit(1);
	}

	if(debug)
		fprintf(stderr,"Uij size = %d\n", num_units*(num_units+1)/2);

	// T memory allocation (6x6 matrix)
	double **T;
	T = (double **) malloc( sizeof(double *) * 6 );
	for(int i=0; i<6; i++)
		T[i] = (double *) malloc( sizeof(double) * 6 );

	if(debug)
		fprintf(stderr,"Uab Computation Main-Loop (Computing Tij elements on-the-fly)\n");

	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
	for(int b= size-3, c= size+5; b >= 0; b--,c--) 	// b --> dihedrals (column)
	{							// c --> deleteable col. index
		for(int a=0; a <= b; a++)	 		// a --> dihedrals (row)
		{ 						// it fills upper triangular part (including diagonal)
			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
			uij_index = undh[a] + ( undh[b] + 1 ) * ( undh[b] + 2 ) / 2;
			// uij_index = undh[a] + undh[b] * ( undh[b] + 1 ) / 2;

			if(debug)
				fprintf(stderr,"a= %d  b= %d  undh[a]= %d  undh[b]+1= %d  uij_index= %d\n", a, b, undh[a], undh[b]+1, uij_index);

			//			fprintf(stderr,"%ld ",uij_index);

			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
			// the Uij for the 4 ICs combination should be computed only once!
			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
			{
				// "in situ" Uij memory allocation
				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
				{
					for(int m=0; m<6; m++)
					{
						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
						{
							for(int n=0; n<6; n++)
								Uij[uij_index][m][n] = 0.0;
						}
						else
						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
					}
				}
				else
				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }

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

				if(debug)
					fprintf(stderr,"\tUij: a=%d b=%d uij1(row=%d col=%d) uij2(row=%d col=%d) uij3(row=%d col=%d) \n",a,b,undh[a],undh[b]+2,undh[a]-1,undh[b]+1,undh[a]-1,undh[b]+2);

				// Computing Uij element
				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
				{				// (whether they have ipas...) --> then, "T" must be computed
					// Computing Tij element "on the fly"
					// (ec.21) Noguti & Go 1983 pp.3685-90
					calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );

					if(debug)
						fprintf(stderr,"\tTij computed\n");

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

				// at this point, the Uij element is fully computed!

				// Show Uij[]
				if(debug)
				{
					fprintf(stderr,"Uij[%d] a=%d b=%d\n", uij_index, a, b);
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							fprintf(stderr," %f", Uij[uij_index][n][m]);
					fprintf(stderr,"\n");
				}

			}   	// only if not computed yet!

			// ****************************************************
			// Rab - Matrix Computation
			// Noguti & Go (1983) pp. 3685-90
			// ****************************************************
			// Rab Computation ( the same Single- or Multi- Chain)
			// Table I. Noguti & Go 1983 pp.3685-90
			// backbone vs. backbone --> Rab = Uab

			// (ec.16) Noguti & Go 1983 pp.3685-90
			// fprintf(stderr,"\tRij= ");
			for(int n=0; n<6; n++)
			{
				R[n] = 0.0;
				for(int m=0; m<6; m++)
					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
				// fprintf(stderr,"%f ",R[n]);
			}
			// fprintf(stderr,"\n");

			temp = 0.0;
			for(int n=0; n<6; n++)
				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

			if(debug)
				fprintf(stderr,"\tHab= %f\n",temp);

			// hess_matrix[a+b*(b+1)/2] = temp; // Rab' x erb = Hessian
			hess_matrix[ a * sizex + b ] = temp; // Rab' x erb = Hessian (square matrix upper triangular part including diagonal?)

		}

		/*
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
				}
			}
		}
		 */
	}


	/*
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
				}
			}
		}
	free( Uij );
	 */

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



	// ************* Add entries corresponding to constraints ******************
	// double v[3]; // some vector
	double w[3]; // some vector
	double wcv[3]; // some vector
	double vw[3]; // some vector
	double sum; // some buffers

	// CA (pseudo)atom index of the Ct-anchor
	k = nla + ifpa + 1;

	// v = r_N - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		v[c] = coord[3*(k-1) + c] - coord[3*k + c];
		prod += pow(v[c], 2);
	}
	prod = sqrt(prod); // norm of r_CA - r_N

	for(c=0; c<3; c++)
		v[c] /= prod;

	// w = r_C - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		w[c] = coord[3*(k+1) + c] - coord[3*k + c];
		prod += pow(w[c], 2);
	}
	prod = sqrt(prod); // norm of r_C - r_CA

	for(c=0; c<3; c++)
		w[c] /= prod;

	// cosg = cos(v^w)
	double cosg = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];     // cos(ang(v,w))

	// wcv = w - cosg*v  (difference between w and the projection of w into v; thus wcv is orthogonal to both v and vw, see next...)
	for(c=0; c<3; c++)          // w - cosg*v
		wcv[c] = w[c] - cosg*v[c];

	// vw = v x w
	vw[0] = v[1]*w[2] - v[2]*w[1];
	vw[1] = v[2]*w[0] - v[0]*w[2];
	vw[2] = v[0]*w[1] - v[1]*w[0];    //  v x w


	int isi, jsi;
	for(i=0, isi=0; i<size; i++, isi += sizex) // Mon: Screen "size" dihedrals
	{
		// Mon: CA motion of Ct-anchor is constrained
		hess_matrix[isi+size  ] = der[(nla+1)*size+i].x; // CA pseudo-atom of Ct-anchor
		hess_matrix[isi+size+1] = der[(nla+1)*size+i].y;
		hess_matrix[isi+size+2] = der[(nla+1)*size+i].z;
		// Mon: The effect of all derivatives on changing the wcv vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * wcv[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * wcv[1];
		sum += der[nla*size+i].z * wcv[2];
		hess_matrix[isi+size+3] = sum; // Mon: der * wcv
		// Mon: Constraint the change in the displacement of the N atom of Ct-anchor vw The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * vw[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * vw[1];
		sum += der[nla*size+i].z * vw[2];
		hess_matrix[isi+size+4] = sum; // Mon: der * vw
		// Mon: The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[(nla+2)*size+i].x * vw[0]; // CO pseudo-atom of Ct-anchor
		sum += der[(nla+2)*size+i].y * vw[1];
		sum += der[(nla+2)*size+i].z * vw[2];
		hess_matrix[isi+size+5] = sum;
	}
	//***************************************************************

	// fill in lower triangular part
	for(i=0,isi=0; i<sizex; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			hess_matrix[isi+j] = hess_matrix[jsi+i];


	return( hess_matrix ); // outputs Hessian Matrix
}

double *hessianFast(float *coord, float *coordCA, trd *der, tri *props, twid *decint, int nipa,  int ifr, int ilr, int size, int nco, int num_atoms)
{
	bool debug = false;
	double prod,prod1;
	int l,ind2,ind3,ls,ks;
	double *dummy;
	int prin, fin, j2, m, buff;
	double r_alpha[3],r_beta[3],r[3],e[3];
	double S[6][6],R[6],C[6][6],D[6][6],U[6][6];
	double d,temp;
	double v[6];
	int resn,num_res,num_seg;
	int index_res = 0;

	if(debug)
		printf("debug> WARNING! Using: hessianFast() nipa %d ifr %d ilr %d size %d nco %d num_atoms %d\n", nipa, ifr, ilr, size, nco, num_atoms);

	// Allocate Bordered Hessian
	int sizex = size + nco;
	double *hess_matrix; // Bordered Hessian matrix
	if( !(hess_matrix = (double *) malloc( sizex * sizex * sizeof(double)) ) ) // Square matrix
	{
		printf("Msg(hessianFast): I'm sorry, Hessian-Matrix memory allocation failed!\nForcing exit!\n");
		exit(1);
	}

	int reglen = ilr - ifr + 1; // Loop length (all mobile residues)
	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop
	//	int ilpa = ifpa + nla - 1; // Index of last (mobile) pseudo-atom of loop
	int ilpa = props[ilr+1].k1 - 1; // Index of last (mobile) pseudo-atom of loop
	int nla = ilpa - ifpa + 1; // Number of loop atoms (mobile) Check!

	if(debug)
		printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	//	size--;
	//	printf("Msg(hessianFast): nipa= %d  size= %d\n", nipa, size);

	// ********************************************
	// Storing ==> (ea, ea x ra) (== (eb, eb x rb) )
	// ********************************************
	//		double **erx;
	//		erx = (double **) malloc( sizeof(double *) * 6); // erx[6]
	//		for(int i=0;i<6;i++)
	//		erx[i] = (double *) malloc( sizeof(double) * size); // erx[6][size]
	//
	//		int *undh; // returns the closest unit-index on the left side of the dihedral
	//		undh = (int *) malloc( sizeof(int) * size );
	double erx[6][size];
	int undh[size];

	int unat[num_atoms];

	for(int x=0; x<num_atoms; x++)
		unat[x] = 0; // all atoms belong to 0 unit by default (i.e. environment)

	// ****************************************************************************
	// * First, screening dihedrals to pre-compute (ea, ea x ra) and (eb, eb x rb) <-- diadic expresions
	// ****************************************************************************
	int k0 = 0; // first-NH + CA + last-CO model index

	for(int i=0; i<ifr; i++)
		k0 += props[i].nat;

	if(debug)
		printf("Msg(hessianFast): k0= %d\n", k0);

	int num_units = 0;
	int un_index = 0; // un_index=0 --> environment rigid unit (i.e. the initialized value), N-CA belong to environment too (they do not move at all)
	int natom = 0;

	// ---------------------------------------------------
	// 1. BUILDING "erx" array and other auxiliary: "undh"
	// ---------------------------------------------------

	int i = 0; // Residue index for current residue
	int j = 0; // Dihedral angle index
	int k = 0; // Atom index for current residue
	int kp = 0; // Atom index in loop context
	int c = 0; // Coordinate index
	int k1 = 0; // Index of first atom of current residue
	int k2 = 0; // Index of last atom of current residue

	//	for(i = ifr; i <= ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor, inclusive.
	for(i = ifr; i < ifr+reglen; i++) // Screen loop residues, i.e. from 1st loop residue to Ct-anchor (without Ct)
	{
		// 1st and last indices of atoms of residue
		k1 = props[i].k1; // index of 1st atom in residue "i"
		k2 = k1 + props[i].nat-1;

		if(debug)
			fprintf(stdout,"hessianFast> residue i=%d  nan=%d  nat=%d  k1= %d  k2= %d\n",i,props[i].nan,props[i].nat,k1,k2);

		unat[k0] = un_index;  // the unit N atom belongs to
		k0++; // update atom counter

		// Vectors for PHI dihedral angles
		// ----------------------------------
		if(props[i].nan != 1) // NOT Proline, i.e. if it has a mobile angle PHI
		{   // q_j = phi_i
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PHI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];


			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			un_index++; // update if it has PHI

			j++; // Update current dihedral variable index
		}

		//		if(i > ifr)  // Not-first-residue (1st residue N-CA is a single rigid unit)
		//			un_index++;

		unat[k0] = un_index;  // the unit CA atom belongs to
		k0++; // update atom counter

		// Vectors for PSI dihedral angles
		// ----------------------------------

		// if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
		{
			// q_j = psi_i  // all aminoacids have mobile PSI angle
			if(debug)
				fprintf(stdout,"hessianFast> residue_i= %3d  dihedral_j= %3d  --> PSI\n",i,j);

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
			erx[0][j] = e[0];
			erx[1][j] = e[1];
			erx[2][j] = e[2];
			// ea x ra
			erx[3][j] = e[1] * r[2] - e[2] * r[1];
			erx[4][j] = e[2] * r[0] - e[0] * r[2];
			erx[5][j] = e[0] * r[1] - e[1] * r[0];

			undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)

			j++; // Update current dihedral variable index
		}

		if(i < ifr+reglen-1)  // Not-last-residue (last residue CA-C can be considered a rigid unit)
			un_index++;
		//		un_index++;

		unat[k0] = un_index; // the unit C atom belongs to
		k0++; // update atom counter

	}

	// MON: PHI at Ct anchor dihedral ????
	undh[j] = un_index; // which unit has the "j" dihedral seen from the right (left side)


	if(debug)
		fprintf(stderr,"hessianFast> undh, unat, (ea, ea x ra) and (eb, eb x rb) initialized!\n");


	if(debug)
	{
		fprintf(stderr,"hessianFast> (ea, ea x ra)\n");
		for(int x=0; x<j; x++)
		{
			for(int y=0; y<6; y++)
				fprintf(stderr," %f", erx[y][x]);
			fprintf(stderr,"\n");
		}


		fprintf(stderr,"hessianFast> undh=");
		for(int x=0; x<size; x++)
			fprintf(stderr," [%d]=%d",x,undh[x]);

		fprintf(stderr,"\nhessianFast> unat=");
		for(int x=0; x<num_atoms; x++)
		{
			if(x % 32 == 0)
				fprintf(stderr,"\n");
			fprintf(stderr," %2d",unat[x]);
		}
		fprintf(stderr,"\n");
	}

	// ---------------------------------------------------------------------------
	// 2. BUILDING "unipa" (array of "ipa"s indices for interacting pair of units)
	// ---------------------------------------------------------------------------
	// Tij related pre-computations (to get direct Tij computation)

	//	num_units = un_index + 1 + 1; // Inter-unit variables (rot-trans) dont increase units number
	num_units = un_index + 1; // total number of units --> M = size + 1

	if(debug)
		printf("Msg(hessianFast): Number of units: %ld (%d)\n",num_units,un_index);

	int unipa_index;
	int *p_int;
	// int **unipa;

	// "unipa" stores "ipa"s indices of interacting atoms belonging to the same pair of units
	// (triangular packing storage)
	// unipa = (int **) malloc( sizeof(int *) * num_units * ( num_units + 1) / 2 );

	// "unipa" initialization
	//	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
	//		unipa[i] = (int *) malloc( sizeof(int) * 500); ;
	//
	//	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
	//		unipa[i][0]=0;

	int unipa[num_units * ( num_units + 1) / 2 ][500];

	for(i=0; i < num_units * ( num_units + 1 ) / 2; i++)
		unipa[i][0]=0;

	if(debug)
		printf("Msg(hessianFast): nipa=%d\n",nipa);

	// Screening "ipa"s to build "unipa" and "nunipa" (this will lead further to direct Tij computation)
	for(int index=0; index<nipa; index++) // screens contacts (k vs. l)
	{

		k = decint[index].k; // k atom index
		l = decint[index].l; // l atom index
		i = unat[k]; // unit index of "k" atom
		j = unat[l]; // unit index of "l" atom

		if(debug)
			printf("Msg(hessianFast): unipa[%d] --> k= %d  l= %d  i= %d  j= %d\n",index,k,l,i,j);

		if(i != j) 	// the same unit is rigid, inner contacts must not be taken into account!
		{		// (as well as the non contacting pairs!)
			// The following must be always true:  (k < l) && (i < j)
			if( i > j )
			{
				buff = j;
				j = i;
				i = buff;
			}

			// translates from "squared" to "triangular" matrix elements
			unipa_index = i + j*(j+1)/2;



			//			// Allocating memory for the new element (memory efficient)
			//			if(unipa[unipa_index] == NULL) // If not allocated yet... then it does it!
			//			{
			//				if( !(p_int = (int *) realloc( unipa[unipa_index], 2 * sizeof(int) ) ) )
			//				{
			//					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
			//					exit(1);
			//				}
			//				unipa[unipa_index] = p_int;
			//				unipa[unipa_index][0] = 2;
			//			}
			//			else // If already allocated, then increases its size one int.
			//			{
			//				unipa[unipa_index][0]++; // Counts number of Interacting Pairs of Atoms between (i,j) units
			//				if( !(p_int = (int *) realloc( unipa[unipa_index], unipa[unipa_index][0] * sizeof(int) ) ) )
			//				{
			//					printf("Msg(hessianFast): Memory allocation failed in realloc()!\nForcing exit!\n");
			//					exit(1);
			//				}
			//				unipa[unipa_index] = p_int;
			//			}

			if (unipa[unipa_index][0] == 0)
				unipa[unipa_index][0] = 1;

			unipa[unipa_index][0]++;

			// Storing ipa index for k,l interacting pair of atoms
			unipa[unipa_index][ unipa[unipa_index][0] - 1 ] = index;
		}
	}

	if(debug)
	{
		fprintf(stderr, "Msg(hessianFast): UNIPA computed!\n");
	}

	// -----------------------------------------------------------------------------------
	// 3. HESSIAN BUILDING from Uab elements:
	// Fastest Hessian Matrix building up... (n^2...) + Some additional optimizations...
	// Uab - Matrix Computation - (ec.23) Noguti & Go 1983 pp.3685-90
	// -----------------------------------------------------------------------------------

	// Hessian initialization
	for(i=0; i < sizex * sizex; i++)
		hess_matrix[i] = 0.0;

	// Uabmn minimal-memory allocation (full-square + 6 cols. + Tij "on the fly" computation)
	int uij_index = 0;
	int uij1, uij2, uij3;

	// MON: check "+1"
	// +1 due to void element (when there are just Traslational/Rotational DoFs)
	double ***Uij;
	if( Uij = (double ***) malloc( sizeof(double **) * (num_units*(num_units+1)/2) ) ) // Uab[size]
	{
		for(i=0; i<(num_units*(num_units+1)/2); i++)
			Uij[i] = NULL; // initialization
	}
	else
	{
		printf("Msg(hessianFast): Unable to allocate Uij-matrix!\n");
		exit(1);
	}

	if(debug)
		fprintf(stderr,"Uij size = %d\n", num_units*(num_units+1)/2);

	// T memory allocation (6x6 matrix)
	//double **T;
	//	T = (double **) malloc( sizeof(double *) * 6 );
	//	for(int i=0; i<6; i++)
	//	T[i] = (double *) malloc( sizeof(double) * 6 );
	double T[6][6];

	if(debug)
		fprintf(stderr,"Uab Computation Main-Loop (Computing Tij elements on-the-fly)\n");

	// Uab Computation Main-Loop (Computing Tij elements "on the fly")
	// (Which "unipa" belongs to the "outher" units corresponding to "a" and "b" dihedrals)
	for(int b= size-3, c= size+5; b >= 0; b--,c--) 	// b --> dihedrals (column)
	{							// c --> deleteable col. index
		for(int a=0; a <= b; a++)	 		// a --> dihedrals (row)
		{ 						// it fills upper triangular part (including diagonal)
			// Units "outside" (a,b) pair are: (undh[a],undh[b]+1)
			uij_index = undh[a] + ( undh[b] + 1 ) * ( undh[b] + 2 ) / 2;
			// uij_index = undh[a] + undh[b] * ( undh[b] + 1 ) / 2;

			if(debug)
				fprintf(stderr,"a= %d  b= %d  undh[a]= %d  undh[b]+1= %d  uij_index= %d\n", a, b, undh[a], undh[b]+1, uij_index);

			//			fprintf(stderr,"%ld ",uij_index);

			// Given there are two dihedrals per CA (aprox.), up to 4 Uab(ICs) contiguous positions
			// are equal. In "Uij"(units) the same 4 ICs combination share the same Uij, thus
			// the Uij for the 4 ICs combination should be computed only once!
			if( Uij[uij_index] == NULL ) // If "Uij" is not computed yet (compute Uij only once!)
			{
				// "in situ" Uij memory allocation
				if( Uij[uij_index] = (double **) malloc( sizeof(double *) * 6 ) )
				{
					for(int m=0; m<6; m++)
					{
						if( Uij[uij_index][m] = (double *) malloc( sizeof(double) * 6 ) )
						{
							for(int n=0; n<6; n++)
								Uij[uij_index][m][n] = 0.0;
						}
						else
						{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }
					}
				}
				else
				{ printf("Msg(hessianM): Unable to allocate an Uab-matrix element!\n"); exit(1); }

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

				if(debug)
					fprintf(stderr,"\tUij: a=%d b=%d uij1(row=%d col=%d) uij2(row=%d col=%d) uij3(row=%d col=%d) \n",a,b,undh[a],undh[b]+2,undh[a]-1,undh[b]+1,undh[a]-1,undh[b]+2);

				// Computing Uij element
				if( unipa[uij_index] != NULL ) // if "a" and "b" 's units interact
				{				// (whether they have ipas...) --> then, "T" must be computed
					// Computing Tij element "on the fly"
					// (ec.21) Noguti & Go 1983 pp.3685-90
					//	calcTij( T, unipa[uij_index], unipa[uij_index][0], decint, coordCA );


					double r_alpha[3];
					double r_beta[3];
					double v[6];

					for(int i=0; i<6; i++)
						for(int j=0; j<6; j++)
							T[i][j] = 0.0; // T initialization

					for(int i=1; i < unipa[uij_index][0]; i++) // screens inter-unit contacts (k vs. l)
					{
						int in=unipa[uij_index][i];
						int ki = 3*decint[in].k;
						int li = 3*decint[in].l;

						// k-atom position retrieval (alpha)
						r_alpha[0] = coordCA[ki];
						r_alpha[1] = coordCA[ki + 1];
						r_alpha[2] = coordCA[ki + 2];

						// l-atom position retrieval (beta)
						r_beta[0] = coordCA[li];
						r_beta[1] = coordCA[li + 1];
						r_beta[2] = coordCA[li + 2];

						// D-Matrix
						v[0] = r_alpha[1] * r_beta[2] - r_alpha[2] * r_beta[1];
						v[1] = r_alpha[2] * r_beta[0] - r_alpha[0] * r_beta[2];
						v[2] = r_alpha[0] * r_beta[1] - r_alpha[1] * r_beta[0];
						v[3] = r_alpha[0] - r_beta[0];
						v[4] = r_alpha[1] - r_beta[1];
						v[5] = r_alpha[2] - r_beta[2];
						// double d1 = decint[in].C / pow(decint[in].d,2);
						double d1 = decint[in].C / decint[in].d;

						//fprintf(stderr,"%f %f %f %d %d\n", d1, decint[in].C , decint[in].d, ki, li);
						for(int n=0;n<6;n++) // Diadic product
							for(int m=0;m<6;m++)
								T[n][m] += d1 * v[n] * v[m]; // Storing S-matrix
					}




					if(debug)
						fprintf(stderr,"\tTij computed\n");

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

				// at this point, the Uij element is fully computed!

				// Show Uij[]
				if(debug)
				{
					fprintf(stderr,"Uij[%d] a=%d b=%d\n", uij_index, a, b);
					for(int n=0;n<6;n++)
						for(int m=0;m<6;m++)
							fprintf(stderr," %f", Uij[uij_index][n][m]);
					fprintf(stderr,"\n");
				}

			}   	// only if not computed yet!

			// ****************************************************
			// Rab - Matrix Computation
			// Noguti & Go (1983) pp. 3685-90
			// ****************************************************
			// Rab Computation ( the same Single- or Multi- Chain)
			// Table I. Noguti & Go 1983 pp.3685-90
			// backbone vs. backbone --> Rab = Uab

			// (ec.16) Noguti & Go 1983 pp.3685-90
			// fprintf(stderr,"\tRij= ");
			for(int n=0; n<6; n++)
			{
				R[n] = 0.0;
				for(int m=0; m<6; m++)
					R[n] += erx[m][a] * Uij[uij_index][m][n]; // era x Rab
				// fprintf(stderr,"%f ",R[n]);
			}
			// fprintf(stderr,"\n");

			temp = 0.0;
			for(int n=0; n<6; n++)
				temp += R[n] * erx[n][b]; // Rab' x erb = Hessian

			if(debug)
				fprintf(stderr,"\tHab= %f\n",temp);

			// hess_matrix[a+b*(b+1)/2] = temp; // Rab' x erb = Hessian
			hess_matrix[ a * sizex + b ] = temp; // Rab' x erb = Hessian (square matrix upper triangular part including diagonal?)

		}

		/*
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
				}
			}
		}
		 */
	}


	/*
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
				}
			}
		}
	free( Uij );
	 */

	//	for(i=0; i<num_units*(num_units+1)/2; i++)
	//		free( unipa[i] );
	//	free( unipa );

	//	for(i=0;i<6;i++)
	//	{
	//		free( erx[i] ); // erx[6][size]
	//		free( T[i] );
	//	}
	//	free( erx );
	//	free( T );
	//	free( undh );



	// ************* Add entries corresponding to constraints ******************
	// double v[3]; // some vector
	double w[3]; // some vector
	double wcv[3]; // some vector
	double vw[3]; // some vector
	double sum; // some buffers

	// CA (pseudo)atom index of the Ct-anchor
	k = nla + ifpa + 1;

	// v = r_N - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		v[c] = coord[3*(k-1) + c] - coord[3*k + c];
		prod += pow(v[c], 2);
	}
	prod = sqrt(prod); // norm of r_CA - r_N

	for(c=0; c<3; c++)
		v[c] /= prod;

	// w = r_C - r_CA   (of Ct-anchor residue)
	prod = 0.0;
	for(c=0; c<3; c++)
	{
		w[c] = coord[3*(k+1) + c] - coord[3*k + c];
		prod += pow(w[c], 2);
	}
	prod = sqrt(prod); // norm of r_C - r_CA

	for(c=0; c<3; c++)
		w[c] /= prod;

	// cosg = cos(v^w)
	double cosg = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];     // cos(ang(v,w))

	// wcv = w - cosg*v  (difference between w and the projection of w into v; thus wcv is orthogonal to both v and vw, see next...)
	for(c=0; c<3; c++)          // w - cosg*v
		wcv[c] = w[c] - cosg*v[c];

	// vw = v x w
	vw[0] = v[1]*w[2] - v[2]*w[1];
	vw[1] = v[2]*w[0] - v[0]*w[2];
	vw[2] = v[0]*w[1] - v[1]*w[0];    //  v x w


	int isi, jsi;
	for(i=0, isi=0; i<size; i++, isi += sizex) // Mon: Screen "size" dihedrals
	{
		// Mon: CA motion of Ct-anchor is constrained
		hess_matrix[isi+size  ] = der[(nla+1)*size+i].x; // CA pseudo-atom of Ct-anchor
		hess_matrix[isi+size+1] = der[(nla+1)*size+i].y;
		hess_matrix[isi+size+2] = der[(nla+1)*size+i].z;
		// Mon: The effect of all derivatives on changing the wcv vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * wcv[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * wcv[1];
		sum += der[nla*size+i].z * wcv[2];
		hess_matrix[isi+size+3] = sum; // Mon: der * wcv
		// Mon: Constraint the change in the displacement of the N atom of Ct-anchor vw The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[nla*size+i].x * vw[0]; // N pseudo-atom of Ct-anchor
		sum += der[nla*size+i].y * vw[1];
		sum += der[nla*size+i].z * vw[2];
		hess_matrix[isi+size+4] = sum; // Mon: der * vw
		// Mon: The effect of all derivatives on changing the vw vector orientation and/or length must be zero.
		sum  = der[(nla+2)*size+i].x * vw[0]; // CO pseudo-atom of Ct-anchor
		sum += der[(nla+2)*size+i].y * vw[1];
		sum += der[(nla+2)*size+i].z * vw[2];
		hess_matrix[isi+size+5] = sum;
	}
	//***************************************************************

	// fill in lower triangular part
	for(i=0,isi=0; i<sizex; i++,isi+=sizex)
		for(j=0,jsi=0; j<i; j++,jsi+=sizex)
			hess_matrix[isi+j] = hess_matrix[jsi+i];

	return( hess_matrix ); // outputs Hessian Matrix

}




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
		//d = decint[in].C / pow(decint[in].d,2);
		d = decint[in].C / decint[in].d;

		for(n=0;n<6;n++) // Diadic product
			for(m=0;m<6;m++)
				T[n][m] += d * v[n] * v[m]; // Storing S-matrix
	}
}


int diag_dggev(double *eigval, double *eigvect, double *mass_matrix, double *hess_matrix, int size, int nco, int *neig)
{
	bool verb = false;

	/* Lapack subroutine arguments */
	char jobvl, jobvr;
	int  info, lda, ldb, ldvl, ldvr, lwork, nesi, nsix, j, m, n, msi, nsi;
	double *alphar, *alphai, *beta, *vr, *vl, *work;
	int sizex = size + nco;

	/***************************************************/
	/*                                                 */
	/*  now solve the generalized eigenvalue problem   */
	/*                                                 */
	/*  hess_matrix * U  =  lambda * mass_matrix * U   */
	/*                                                 */
	/***************************************************/

	/*
	 * SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
	     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
	 *
	 * In principle, to call a Fortran routine
	 * from C we have to transform the matrix
	 * from row major order to column major order,
	 * but since the matrices are symmetric,
	 * this is not necessary.
	 *
	 *
	 *  -- LAPACK driver routine (version 3.0) --
	 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 *     Courant Institute, Argonne National Lab, and Rice University
	 *     June 30, 1999
	 *
	 *     .. Scalar Arguments ..
	      CHARACTER          JOBVL, JOBVR
	      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
	 *     ..
	 *     .. Array Arguments ..
	      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
	     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
	     $                   VR( LDVR, * ), WORK( * )
	 *     ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
	 *  the generalized eigenvalues, and optionally, the left and/or right
	 *  generalized eigenvectors.
	 *
	 *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
	 *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
	 *  singular. It is usually represented as the pair (alpha,beta), as
	 *  there is a reasonable interpretation for beta=0, and even for both
	 *  being zero.
	 *
	 *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
	 *  of (A,B) satisfies
	 *
	 *                   A * v(j) = lambda(j) * B * v(j).
	 *
	 *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
	 *  of (A,B) satisfies
	 *
	 *                   u(j)**H * A  = lambda(j) * u(j)**H * B .
	 *
	 *  where u(j)**H is the conjugate-transpose of u(j).
	 *
	 *
	 *  Arguments
	 *  =========
	 *
	 *  JOBVL   (input) CHARACTER*1
	 *          = 'N':  do not compute the left generalized eigenvectors;
	 *          = 'V':  compute the left generalized eigenvectors.
	 *
	 *  JOBVR   (input) CHARACTER*1
	 *          = 'N':  do not compute the right generalized eigenvectors;
	 *          = 'V':  compute the right generalized eigenvectors.
	 *
	 *  N       (input) INTEGER
	 *          The order of the matrices A, B, VL, and VR.  N >= 0.
	 *
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
	 *          On entry, the matrix A in the pair (A,B).
	 *          On exit, A has been overwritten.
	 *
	 *  LDA     (input) INTEGER
	 *          The leading dimension of A.  LDA >= max(1,N).
	 *
	 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
	 *          On entry, the matrix B in the pair (A,B).
	 *          On exit, B has been overwritten.
	 *
	 *  LDB     (input) INTEGER
	 *          The leading dimension of B.  LDB >= max(1,N).
	 *
	 *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
	 *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
	 *  BETA    (output) DOUBLE PRECISION array, dimension (N)
	 *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
	 *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
	 *          the j-th eigenvalue is real; if positive, then the j-th and
	 *          (j+1)-st eigenvalues are a complex conjugate pair, with
	 *          ALPHAI(j+1) negative.
	 *
	 *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
	 *          may easily over- or underflow, and BETA(j) may even be zero.
	 *          Thus, the user should avoid naively computing the ratio
	 *          alpha/beta.  However, ALPHAR and ALPHAI will be always less
	 *          than and usually comparable with norm(A) in magnitude, and
	 *          BETA always less than and usually comparable with norm(B).
	 *
	 *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
	 *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
	 *          after another in the columns of VL, in the same order as
	 *          their eigenvalues. If the j-th eigenvalue is real, then
	 *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
	 *          (j+1)-th eigenvalues form a complex conjugate pair, then
	 *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
	 *          Each eigenvector will be scaled so the largest component have
	 *          abs(real part)+abs(imag. part)=1.
	 *          Not referenced if JOBVL = 'N'.
	 *
	 *  LDVL    (input) INTEGER
	 *          The leading dimension of the matrix VL. LDVL >= 1, and
	 *          if JOBVL = 'V', LDVL >= N.
	 *
	 *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
	 *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
	 *          after another in the columns of VR, in the same order as
	 *          their eigenvalues. If the j-th eigenvalue is real, then
	 *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
	 *          (j+1)-th eigenvalues form a complex conjugate pair, then
	 *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
	 *          Each eigenvector will be scaled so the largest component have
	 *          abs(real part)+abs(imag. part)=1.
	 *          Not referenced if JOBVR = 'N'.
	 *
	 *  LDVR    (input) INTEGER
	 *          The leading dimension of the matrix VR. LDVR >= 1, and
	 *          if JOBVR = 'V', LDVR >= N.
	 *
	 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
	 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	 *
	 *  LWORK   (input) INTEGER
	 *          The dimension of the array WORK.  LWORK >= max(1,8*N).
	 *          For good performance, LWORK must generally be larger.
	 *
	 *          If LWORK = -1, then a workspace query is assumed; the routine
	 *          only calculates the optimal size of the WORK array, returns
	 *          this value as the first entry of the WORK array, and no error
	 *          message related to LWORK is issued by XERBLA.
	 *
	 *  INFO    (output) INTEGER
	 *          = 0:  successful exit
	 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
	 *          = 1,...,N:
	 *                The QZ iteration failed.  No eigenvectors have been
	 *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
	 *                should be correct for j=INFO+1,...,N.
	 *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
	 *                =N+2: error return from DTGEVC.
	 *
	 *  =====================================================================
	 */

	jobvl='N';   /* do not compute left eigenvectors */

	jobvr='V';   /* compute the right eigenvectors */

	lda=sizex;
	ldb=sizex;
	ldvl=1;
	ldvr=sizex;

	alphar = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphar);

	alphai = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphai);

	beta = (double *) malloc(sizex * sizeof(double));
	ptr_check(beta);

	vl = (double *) malloc(1*sizex * sizeof(double));
	ptr_check(vl);

	vr = (double *) malloc(sizex*sizex * sizeof(double));   /* eigenvectors */
	ptr_check(vr);

	lwork = 20*sizex;
	work = (double *) malloc(lwork * sizeof(double));
	ptr_check(work);

	dggev_(&jobvl, &jobvr, &sizex, hess_matrix, &lda,
			mass_matrix, &ldb, alphar, alphai, beta,
			vl, &ldvl, vr, &ldvr, work, &lwork, &info);

	// neig will contain the number of "useful" eigenvectors
	*neig = 0;
	for(n=0; n<sizex; n++)
	{
		if(verb>2)
			fprintf(stdout,"ar= %.10e, ai= %.10e, b= %.10e\n",alphar[n], alphai[n], beta[n]);

		if(beta[n]!=0.0 && alphai[n]==0.0) // alphai == 0 means NO imaginary eigenvalue
		{
			if(verb>2)
				fprintf(stdout,"ar= %.10e, ai= %.10e, b= %.10e, lambda= %.10e\n",alphar[n], alphai[n], beta[n], alphar[n]/beta[n]);

			eigval[*neig] = alphar[n] / beta[n];

			//			fprintf(stderr,"Eigenvalue %d: %f, ",n,eigval[*neig]);
			//			if(alphai[n]!=0.0)
			//				fprintf(stderr,"alphai %2d is NOT Zero: %f",n,alphai[n]);

			nesi = (*neig) * size;
			nsix = n*sizex;
			for(j=0; j<size; j++)
				eigvect[nesi+j] = vr[nsix+j];

			//			fprintf(stdout, "Full eigenvector dump:\n");
			//			for(j=0; j<sizex; j++)
			//				fprintf(stdout, " %8.2e", vr[nsix+j]);
			//			fprintf(stdout, " \n");

			(*neig)++;
		}
		//		else
		//		{
		//			// MON
		//			if(beta[n]==0.0)
		//				fprintf(stderr,"beta %2d is Zero: %f, ",n,beta[n]);
		//			if(alphai[n]!=0.0)
		//				fprintf(stderr,"alphai %2d is NOT Zero: %f",n,alphai[n]);
		//			if(alphai[n]==0.0)
		//				fprintf(stderr,"alphai %2d is Zero: %f and alphar is: %f",n,alphai[n],alphar[n]);
		//		}
		//		fprintf(stderr,"\n");
	}

	// sort by increasing eigenvalue
	for(m=0;m<(*neig)-1;m++)
	{
		msi = m*size;
		for(n=m+1;n<(*neig);n++)
			if(eigval[m] > eigval[n])
			{
				SWAPPING(eigval[m], eigval[n], double);
				nsi = n*size;
				for(j=0;j<size;j++)
					SWAPPING(eigvect[msi+j], eigvect[nsi+j], double);
			}
	}

	free(work);
	free(vr);

	free(alphar);
	free(alphai);
	free(beta);
	free(vl);

	return info;
}


inline int diag_dggev(double *eigval, double *eigvect, double *mass_matrix, double *hess_matrix, int size, int nco, int *neig, double *alphar, double *alphai, double *beta, double *vr,  double *vl,  double *work)
{
	bool verb = false;

	/* Lapack subroutine arguments */
	char jobvl, jobvr;
	int  info, lda, ldb, ldvl, ldvr, lwork, nesi, nsix, j, m, n, msi, nsi;

	int sizex = size + nco;

	/***************************************************/
	/*                                                 */
	/*  now solve the generalized eigenvalue problem   */
	/*                                                 */
	/*  hess_matrix * U  =  lambda * mass_matrix * U   */
	/*                                                 */
	/***************************************************/

	/*
	 * SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
	     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
	 *
	 * In principle, to call a Fortran routine
	 * from C we have to transform the matrix
	 * from row major order to column major order,
	 * but since the matrices are symmetric,
	 * this is not necessary.
	 *
	 *
	 *  -- LAPACK driver routine (version 3.0) --
	 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
	 *     Courant Institute, Argonne National Lab, and Rice University
	 *     June 30, 1999
	 *
	 *     .. Scalar Arguments ..
	      CHARACTER          JOBVL, JOBVR
	      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
	 *     ..
	 *     .. Array Arguments ..
	      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
	     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
	     $                   VR( LDVR, * ), WORK( * )
	 *     ..
	 *
	 *  Purpose
	 *  =======
	 *
	 *  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
	 *  the generalized eigenvalues, and optionally, the left and/or right
	 *  generalized eigenvectors.
	 *
	 *  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
	 *  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
	 *  singular. It is usually represented as the pair (alpha,beta), as
	 *  there is a reasonable interpretation for beta=0, and even for both
	 *  being zero.
	 *
	 *  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
	 *  of (A,B) satisfies
	 *
	 *                   A * v(j) = lambda(j) * B * v(j).
	 *
	 *  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
	 *  of (A,B) satisfies
	 *
	 *                   u(j)**H * A  = lambda(j) * u(j)**H * B .
	 *
	 *  where u(j)**H is the conjugate-transpose of u(j).
	 *
	 *
	 *  Arguments
	 *  =========
	 *
	 *  JOBVL   (input) CHARACTER*1
	 *          = 'N':  do not compute the left generalized eigenvectors;
	 *          = 'V':  compute the left generalized eigenvectors.
	 *
	 *  JOBVR   (input) CHARACTER*1
	 *          = 'N':  do not compute the right generalized eigenvectors;
	 *          = 'V':  compute the right generalized eigenvectors.
	 *
	 *  N       (input) INTEGER
	 *          The order of the matrices A, B, VL, and VR.  N >= 0.
	 *
	 *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
	 *          On entry, the matrix A in the pair (A,B).
	 *          On exit, A has been overwritten.
	 *
	 *  LDA     (input) INTEGER
	 *          The leading dimension of A.  LDA >= max(1,N).
	 *
	 *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
	 *          On entry, the matrix B in the pair (A,B).
	 *          On exit, B has been overwritten.
	 *
	 *  LDB     (input) INTEGER
	 *          The leading dimension of B.  LDB >= max(1,N).
	 *
	 *  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
	 *  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
	 *  BETA    (output) DOUBLE PRECISION array, dimension (N)
	 *          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
	 *          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
	 *          the j-th eigenvalue is real; if positive, then the j-th and
	 *          (j+1)-st eigenvalues are a complex conjugate pair, with
	 *          ALPHAI(j+1) negative.
	 *
	 *          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
	 *          may easily over- or underflow, and BETA(j) may even be zero.
	 *          Thus, the user should avoid naively computing the ratio
	 *          alpha/beta.  However, ALPHAR and ALPHAI will be always less
	 *          than and usually comparable with norm(A) in magnitude, and
	 *          BETA always less than and usually comparable with norm(B).
	 *
	 *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
	 *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
	 *          after another in the columns of VL, in the same order as
	 *          their eigenvalues. If the j-th eigenvalue is real, then
	 *          u(j) = VL(:,j), the j-th column of VL. If the j-th and
	 *          (j+1)-th eigenvalues form a complex conjugate pair, then
	 *          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
	 *          Each eigenvector will be scaled so the largest component have
	 *          abs(real part)+abs(imag. part)=1.
	 *          Not referenced if JOBVL = 'N'.
	 *
	 *  LDVL    (input) INTEGER
	 *          The leading dimension of the matrix VL. LDVL >= 1, and
	 *          if JOBVL = 'V', LDVL >= N.
	 *
	 *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
	 *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
	 *          after another in the columns of VR, in the same order as
	 *          their eigenvalues. If the j-th eigenvalue is real, then
	 *          v(j) = VR(:,j), the j-th column of VR. If the j-th and
	 *          (j+1)-th eigenvalues form a complex conjugate pair, then
	 *          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
	 *          Each eigenvector will be scaled so the largest component have
	 *          abs(real part)+abs(imag. part)=1.
	 *          Not referenced if JOBVR = 'N'.
	 *
	 *  LDVR    (input) INTEGER
	 *          The leading dimension of the matrix VR. LDVR >= 1, and
	 *          if JOBVR = 'V', LDVR >= N.
	 *
	 *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
	 *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	 *
	 *  LWORK   (input) INTEGER
	 *          The dimension of the array WORK.  LWORK >= max(1,8*N).
	 *          For good performance, LWORK must generally be larger.
	 *
	 *          If LWORK = -1, then a workspace query is assumed; the routine
	 *          only calculates the optimal size of the WORK array, returns
	 *          this value as the first entry of the WORK array, and no error
	 *          message related to LWORK is issued by XERBLA.
	 *
	 *  INFO    (output) INTEGER
	 *          = 0:  successful exit
	 *          < 0:  if INFO = -i, the i-th argument had an illegal value.
	 *          = 1,...,N:
	 *                The QZ iteration failed.  No eigenvectors have been
	 *                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
	 *                should be correct for j=INFO+1,...,N.
	 *          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
	 *                =N+2: error return from DTGEVC.
	 *
	 *  =====================================================================
	 */

	jobvl='N';   /* do not compute left eigenvectors */

	jobvr='V';   /* compute the right eigenvectors */

	lda=sizex;
	ldb=sizex;
	ldvl=1;
	ldvr=sizex;
	lwork=20*sizex;


	dggev_(&jobvl, &jobvr, &sizex, hess_matrix, &lda,
			mass_matrix, &ldb, alphar, alphai, beta,
			vl, &ldvl, vr, &ldvr, work, &lwork, &info);

	// neig will contain the number of "useful" eigenvectors
	*neig = 0;
	for(n=0; n<sizex; n++)
	{
		if(verb>2)
			fprintf(stdout,"ar= %.10e, ai= %.10e, b= %.10e\n",alphar[n], alphai[n], beta[n]);

		if(beta[n]!=0.0 && alphai[n]==0.0) // alphai == 0 means NO imaginary eigenvalue
		{
			if(verb>2)
				fprintf(stdout,"ar= %.10e, ai= %.10e, b= %.10e, lambda= %.10e\n",alphar[n], alphai[n], beta[n], alphar[n]/beta[n]);

			eigval[*neig] = alphar[n] / beta[n];

			//			fprintf(stderr,"Eigenvalue %d: %f, ",n,eigval[*neig]);
			//			if(alphai[n]!=0.0)
			//				fprintf(stderr,"alphai %2d is NOT Zero: %f",n,alphai[n]);

			nesi = (*neig) * size;
			nsix = n*sizex;
			for(j=0; j<size; j++)
				eigvect[nesi+j] = vr[nsix+j];

			//			fprintf(stdout, "Full eigenvector dump:\n");
			//			for(j=0; j<sizex; j++)
			//				fprintf(stdout, " %8.2e", vr[nsix+j]);
			//			fprintf(stdout, " \n");

			(*neig)++;
		}
		//		else
		//		{
		//			// MON
		//			if(beta[n]==0.0)
		//				fprintf(stderr,"beta %2d is Zero: %f, ",n,beta[n]);
		//			if(alphai[n]!=0.0)
		//				fprintf(stderr,"alphai %2d is NOT Zero: %f",n,alphai[n]);
		//			if(alphai[n]==0.0)
		//				fprintf(stderr,"alphai %2d is Zero: %f and alphar is: %f",n,alphai[n],alphar[n]);
		//		}
		//		fprintf(stderr,"\n");
	}

	// sort by increasing eigenvalue
	for(m=0;m<(*neig)-1;m++)
	{
		msi = m*size;
		for(n=m+1;n<(*neig);n++)
			if(eigval[m] > eigval[n])
			{
				SWAPPING(eigval[m], eigval[n], double);
				nsi = n*size;
				for(j=0;j<size;j++)
					SWAPPING(eigvect[msi+j], eigvect[nsi+j], double);
			}
	}



	return info;
}



//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
void diag_dsygvx(double *hess_matrix,double *mass_matrix,double *eigval, int size, int il, int iu)
{
	/* Lapack input/output */
	char jobz, uplo, range;
	int  itype, lda, ldb, ldz, lwork, info, *iwork, *ifail, nulo, m, evec_size;
	double *work,abstol,*p_evec;

	itype=1;   /* specifies A*x = lambda*B*x */
	jobz='V';  /* Compute eigenvalues and eigenvectors */
	range='I'; // 'I': the IL-th through IU-th eigenvalues will be found.
	uplo='U';  /* Upper triangle of A is stored */
	lda=size;
	ldb=size;
	ldz=size;
	abstol=0;
	lwork = 2*size*size+6*size+1;
	work = (double *) malloc(lwork * sizeof(double));
	for(int i=0; i<lwork; i++) work[i]=0.0;

	iwork = (int *) malloc(5*size * sizeof(int));
	for(int i=0; i<5*size; i++) iwork[i]=0;

	ifail = (int *) malloc(size * sizeof(int));
	for(int i=0; i<size; i++) ifail[i]=0;

	evec_size = (iu-il+1);
	p_evec = (double *) malloc( size * evec_size * sizeof(double) );
	for(int i=0; i<size*evec_size; i++) p_evec[i]=0;

	// Diagonalization
	dsygvx_(&itype, &jobz, &range, &uplo, &size, hess_matrix, &lda,
			mass_matrix, &ldb, &nulo, &nulo, &il, &iu, &abstol, &m,
			eigval, p_evec, &ldz, work, &lwork, iwork,
			ifail, &info);

	fprintf(stdout,"Msg(diag_dsygvx): %d eigenvectors found! info= %d",m,info);

	// Outputting as usual... (eigenvectors inside hess_matrix)
	for(int i=0; i<size*evec_size; i++)
		hess_matrix[i] = p_evec[i]; // copying eigenvectors into hess_matrix

	free(p_evec);
	free(ifail);
	free(work); // <-- "work" could be allocated only once outside!!!
	free(iwork); // <-- "iwork" could be allocated only once outside!!!
}

//SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
//$                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
//$                   LWORK, IWORK, IFAIL, INFO )
//*
//*  -- LAPACK driver routine (version 3.1) --
//*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
//*     November 2006
//*
//*     .. Scalar Arguments ..
//CHARACTER          JOBZ, RANGE, UPLO
//INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
//DOUBLE PRECISION   ABSTOL, VL, VU
//*     ..
//*     .. Array Arguments ..
//INTEGER            IFAIL( * ), IWORK( * )
//DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),
//$                   Z( LDZ, * )
//*     ..
//*
//*  Purpose
//*  =======
//*
//*  DSYGVX computes selected eigenvalues, and optionally, eigenvectors
//*  of a real generalized symmetric-definite eigenproblem, of the form
//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
//*  and B are assumed to be symmetric and B is also positive definite.
//*  Eigenvalues and eigenvectors can be selected by specifying either a
//*  range of values or a range of indices for the desired eigenvalues.
//*
//*  Arguments
//*  =========
//*
//*  ITYPE   (input) INTEGER
//*          Specifies the problem type to be solved:
//*          = 1:  A*x = (lambda)*B*x
//*          = 2:  A*B*x = (lambda)*x
//*          = 3:  B*A*x = (lambda)*x
//*
//*  JOBZ    (input) CHARACTER*1
//*          = 'N':  Compute eigenvalues only;
//*          = 'V':  Compute eigenvalues and eigenvectors.
//*
//*  RANGE   (input) CHARACTER*1
//*          = 'A': all eigenvalues will be found.
//*          = 'V': all eigenvalues in the half-open interval (VL,VU]
//*                 will be found.
//*          = 'I': the IL-th through IU-th eigenvalues will be found.
//*
//*  UPLO    (input) CHARACTER*1
//*          = 'U':  Upper triangle of A and B are stored;
//*          = 'L':  Lower triangle of A and B are stored.
//*
//*  N       (input) INTEGER
//*          The order of the matrix pencil (A,B).  N >= 0.
//*
//*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
//*          On entry, the symmetric matrix A.  If UPLO = 'U', the
//*          leading N-by-N upper triangular part of A contains the
//*          upper triangular part of the matrix A.  If UPLO = 'L',
//*          the leading N-by-N lower triangular part of A contains
//*          the lower triangular part of the matrix A.
//*
//*          On exit, the lower triangle (if UPLO='L') or the upper
//*          triangle (if UPLO='U') of A, including the diagonal, is
//*          destroyed.
//*
//*  LDA     (input) INTEGER
//*          The leading dimension of the array A.  LDA >= max(1,N).
//*
//*  B       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
//*          On entry, the symmetric matrix B.  If UPLO = 'U', the
//*          leading N-by-N upper triangular part of B contains the
//*          upper triangular part of the matrix B.  If UPLO = 'L',
//*          the leading N-by-N lower triangular part of B contains
//*          the lower triangular part of the matrix B.
//*
//*          On exit, if INFO <= N, the part of B containing the matrix is
//*          overwritten by the triangular factor U or L from the Cholesky
//*          factorization B = U**T*U or B = L*L**T.
//*
//*  LDB     (input) INTEGER
//*          The leading dimension of the array B.  LDB >= max(1,N).
//*
//*  VL      (input) DOUBLE PRECISION
//*  VU      (input) DOUBLE PRECISION
//*          If RANGE='V', the lower and upper bounds of the interval to
//*          be searched for eigenvalues. VL < VU.
//*          Not referenced if RANGE = 'A' or 'I'.
//*
//*  IL      (input) INTEGER
//*  IU      (input) INTEGER
//*          If RANGE='I', the indices (in ascending order) of the
//*          smallest and largest eigenvalues to be returned.
//*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
//*          Not referenced if RANGE = 'A' or 'V'.
//*
//*  ABSTOL  (input) DOUBLE PRECISION
//*          The absolute error tolerance for the eigenvalues.
//*          An approximate eigenvalue is accepted as converged
//*          when it is determined to lie in an interval [a,b]
//*          of width less than or equal to
//*
//*                  ABSTOL + EPS *   max( |a|,|b| ) ,
//*
//*          where EPS is the machine precision.  If ABSTOL is less than
//*          or equal to zero, then  EPS*|T|  will be used in its place,
//*          where |T| is the 1-norm of the tridiagonal matrix obtained
//*          by reducing A to tridiagonal form.
//*
//*          Eigenvalues will be computed most accurately when ABSTOL is
//*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
//*          If this routine returns with INFO>0, indicating that some
//*          eigenvectors did not converge, try setting ABSTOL to
//*          2*DLAMCH('S').
//*
//*  M       (output) INTEGER
//*          The total number of eigenvalues found.  0 <= M <= N.
//*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
//*
//*  W       (output) DOUBLE PRECISION array, dimension (N)
//*          On normal exit, the first M elements contain the selected
//*          eigenvalues in ascending order.
//*
//*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
//*          If JOBZ = 'N', then Z is not referenced.
//*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
//*          contain the orthonormal eigenvectors of the matrix A
//*          corresponding to the selected eigenvalues, with the i-th
//*          column of Z holding the eigenvector associated with W(i).
//*          The eigenvectors are normalized as follows:
//*          if ITYPE = 1 or 2, Z**T*B*Z = I;
//*          if ITYPE = 3, Z**T*inv(B)*Z = I.
//*
//*          If an eigenvector fails to converge, then that column of Z
//*          contains the latest approximation to the eigenvector, and the
//*          index of the eigenvector is returned in IFAIL.
//*          Note: the user must ensure that at least max(1,M) columns are
//*          supplied in the array Z; if RANGE = 'V', the exact value of M
//*          is not known in advance and an upper bound must be used.
//*
//*  LDZ     (input) INTEGER
//*          The leading dimension of the array Z.  LDZ >= 1, and if
//*          JOBZ = 'V', LDZ >= max(1,N).
//*
//*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
//*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
//*
//*  LWORK   (input) INTEGER
//*          The length of the array WORK.  LWORK >= max(1,8*N).
//*          For optimal efficiency, LWORK >= (NB+3)*N,
//*          where NB is the blocksize for DSYTRD returned by ILAENV.
//*
//*          If LWORK = -1, then a workspace query is assumed; the routine
//*          only calculates the optimal size of the WORK array, returns
//*          this value as the first entry of the WORK array, and no error
//*          message related to LWORK is issued by XERBLA.
//*
//*  IWORK   (workspace) INTEGER array, dimension (5*N)
//*
//*  IFAIL   (output) INTEGER array, dimension (N)
//*          If JOBZ = 'V', then if INFO = 0, the first M elements of
//*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
//*          indices of the eigenvectors that failed to converge.
//*          If JOBZ = 'N', then IFAIL is not referenced.
//*
//*  INFO    (output) INTEGER
//*          = 0:  successful exit
//*          < 0:  if INFO = -i, the i-th argument had an illegal value
//*          > 0:  DPOTRF or DSYEVX returned an error code:
//*             <= N:  if INFO = i, DSYEVX failed to converge;
//*                    i eigenvectors failed to converge.  Their indices
//*                    are stored in array IFAIL.
//*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
//*                    minor of order i of B is not positive definite.
//*                    The factorization of B could not be completed and
//*                    no eigenvalues or eigenvectors were computed.
//*
//*  Further Details
//*  ===============
//*
//*  Based on contributions by
//*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
//*
//* =====================================================================

// Show square matrix
void show_matrix(double *matrix, int size, char *name, char *fmt)
{
	fprintf(stdout,"%s\n",name);
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)
		{
			fprintf(stdout,fmt,matrix[size*i+j]);
		}
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");
}

void show_cartmode(float *coord, double *Aop, tri *props, int ifr, int nla, char *text, int imod, float factor)
{
	char colortext[110][20] = { {"yellow"}, {"red"}, {"green"}, {"cyan"}, {"white"}, {"orange"}, {"gray"}, {"lime"}, {"pink"}, {"magenta"} };
	float radius = 0.1; // arrow main body radius

	FILE *f_vmd;
	f_vmd = fopen(text, "w");
	if(f_vmd == NULL )
	{
		fprintf(stderr, "> Error: can't open file! %s\n", text);
		exit(1);
	}
	fprintf(f_vmd,"molecule new\n");
	fprintf(f_vmd,"display resetview\n");
	fprintf(f_vmd,"draw color %s\n", colortext[imod]);

	int ifpa = props[ifr].k1; // Index of first (mobile) pseudo-atom of loop

	for(int kp=0; kp<nla+3; kp++)
	{
		int k = kp + ifpa; // Mon: absolute index of "kp" pseudo-atom

		fprintf(f_vmd,"draw cylinder \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20 filled 1\n",
				coord[3*k  ], coord[3*k+1], coord[3*k+2],
				coord[3*k  ] + Aop[imod*3*(nla+3)+3*kp] * factor *0.7,
				coord[3*k+1] + Aop[imod*3*(nla+3)+3*kp+1] * factor *0.7,
				coord[3*k+2] + Aop[imod*3*(nla+3)+3*kp+2] * factor *0.7,
				radius);
		fprintf(f_vmd,"draw cone \"%8.5f %8.5f %8.5f\" \"%8.5f %8.5f %8.5f\" radius %.4f resolution 20\n",
				coord[3*k  ] + Aop[imod*3*(nla+3)+3*kp] *factor *0.7,
				coord[3*k+1] + Aop[imod*3*(nla+3)+3*kp+1] *factor *0.7,
				coord[3*k+2] + Aop[imod*3*(nla+3)+3*kp+2] *factor *0.7,
				coord[3*k  ] + Aop[imod*3*(nla+3)+3*kp] *factor,
				coord[3*k+1] + Aop[imod*3*(nla+3)+3*kp+1] *factor,
				coord[3*k+2] + Aop[imod*3*(nla+3)+3*kp+2] *factor,
				radius*3);
	}
	fclose(f_vmd);
}

// Show a matrix vector-wise...
void show_vectors(FILE *f, double *v, int size, int n, char *name, char *fmt)
{
	fprintf(f,"%s\n",name);
	for(int i=0; i < n;i++)
	{
		fprintf(f,"Vector %d:\n",i+1);
		for(int j=0; j<size; j++) // screen variables
			fprintf(f, fmt, v[i*size + j]);
		fprintf(f,"\n");
	}
}

// Show a single vector
void show_vector(FILE *f, double *v, int size, char *name, char *fmt, bool newline, bool newline2)
{
	//	if(strlen(name) > 0)
	if(newline)
		fprintf(f,"%s\n",name);
	else
		fprintf(f,"%s",name);
	for(int j=0; j<size; j++) // screen variables
		fprintf(f, fmt, v[j]);
	if(newline2)
		fprintf(f,"\n");
}

// Show a single vector (float)
void show_vector(FILE *f, float *v, int size, char *name, char *fmt, bool newline, bool newline2)
{
	//	if(strlen(name) > 0)
	if(newline)
		fprintf(f,"%s\n",name);
	else
		fprintf(f,"%s",name);
	for(int j=0; j<size; j++) // screen variables
		fprintf(f, fmt, v[j]);
	if(newline2)
		fprintf(f,"\n");
}

// Convert IC (dihedral) modes into Cartesian modes
//	masses_sqrt --> Array of the square root of the atomic masses (one per loop atom)
inline double *ic2cart(double *eigvect, int nevec, trd *der, int size, int nla, float *masses_sqrt)
{

	// Cartesian modes
	double *Aop  = (double *) malloc( 3 * nla * nevec * sizeof(double));

	// Compute the matrix of Cartesian normal modes
	double v[3]; // some vector

	if(masses_sqrt == NULL) // Not-weighted Cartesian eigenvectors
	{
		for(int i = 0; i < nevec;i++) // Screen computed modes
		{
			for(int kp=0; kp<nla; kp++) // screen atoms  // Mon added the +3, watch out!
			{
				v[0] = v[1] = v[2] = 0.0;
				for(int j=0; j<size; j++) // screen variables
				{
					v[0] += der[kp*size + j].x * eigvect[i*size + j];
					v[1] += der[kp*size + j].y * eigvect[i*size + j];
					v[2] += der[kp*size + j].z * eigvect[i*size + j];
				}
				Aop[i * 3 * nla + 3*kp ] = v[0];
				Aop[i * 3 * nla + 3*kp + 1] = v[1];
				Aop[i * 3 * nla + 3*kp + 2] = v[2];
			}
		}
	}
	else
	{
		for(int i = 0; i < nevec;i++) // Screen computed modes
		{
			for(int kp=0; kp<nla; kp++) // screen atoms  // Mon added the +3, watch out!
			{
				v[0] = v[1] = v[2] = 0.0;
				for(int j=0; j<size; j++) // screen variables
				{
					v[0] += der[kp*size + j].x * eigvect[i*size + j];
					v[1] += der[kp*size + j].y * eigvect[i*size + j];
					v[2] += der[kp*size + j].z * eigvect[i*size + j];
				}
				Aop[i * 3 * nla + 3*kp ]    =  masses_sqrt[kp] * v[0];
				Aop[i * 3 * nla + 3*kp + 1] = masses_sqrt[kp] * v[1];
				Aop[i * 3 * nla + 3*kp + 2] = masses_sqrt[kp] * v[2];
			}
		}
	}

	return Aop;
}

inline void  *ic2cart(double *eigvect, int nevec, trd *der, int size, int nla, double *Aop, float *masses_sqrt)
{

	// Cartesian modes
	//double *Aop  = (double *) malloc( 3 * nla * nevec * sizeof(double));

	// Compute the matrix of Cartesian normal modes
	double v[3]; // some vector

	if(masses_sqrt == NULL) // Not-weighted Cartesian eigenvectors
	{
		for(int i = 0; i < nevec;i++) // Screen computed modes
		{
			for(int kp=0; kp<nla; kp++) // screen atoms  // Mon added the +3, watch out!
			{
				v[0] = v[1] = v[2] = 0.0;
				for(int j=0; j<size; j++) // screen variables
				{
					v[0] += der[kp*size + j].x * eigvect[i*size + j];
					v[1] += der[kp*size + j].y * eigvect[i*size + j];
					v[2] += der[kp*size + j].z * eigvect[i*size + j];
				}
				Aop[i * 3 * nla + 3*kp ] = v[0];
				Aop[i * 3 * nla + 3*kp + 1] = v[1];
				Aop[i * 3 * nla + 3*kp + 2] = v[2];
			}
		}
	}
	else
	{
		for(int i = 0; i < nevec;i++) // Screen computed modes
		{
			for(int kp=0; kp<nla; kp++) // screen atoms  // Mon added the +3, watch out!
			{
				v[0] = v[1] = v[2] = 0.0;
				for(int j=0; j<size; j++) // screen variables
				{
					v[0] += der[kp*size + j].x * eigvect[i*size + j];
					v[1] += der[kp*size + j].y * eigvect[i*size + j];
					v[2] += der[kp*size + j].z * eigvect[i*size + j];
				}
				Aop[i * 3 * nla + 3*kp ]    =  masses_sqrt[kp] * v[0];
				Aop[i * 3 * nla + 3*kp + 1] = masses_sqrt[kp] * v[1];
				Aop[i * 3 * nla + 3*kp + 2] = masses_sqrt[kp] * v[2];
			}
		}
	}

}

void move_loop_linear(pdbIter *iter, double *cevec, int imod, int ifpa, int nla)

{
	Atom *at;


	// Generate a linear trajectory, including Ct-anchor for checking purposes...
	Tcoor pos;



	int kp = 0;
	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		// Get PDB residue index
		at = (Atom *) iter->get_atom();
		at->getPosition(pos);

		// Apply position increment
		pos[0] += cevec[imod*3*nla+3*kp  ];
		pos[1] += cevec[imod*3*nla+3*kp+1];
		pos[2] += cevec[imod*3*nla+3*kp+2];

		// Set position
		at->setPosition(pos);

		kp++; // Update atom index
	}


}


void move_loop_linear_factor(pdbIter *iter, double *cevec, int imod, int ifpa, int nla, float factor)

{
	Atom *at;


	// Generate a linear trajectory, including Ct-anchor for checking purposes...
	Tcoor pos;



	int kp = 0;
	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		// Get PDB residue index
		at = (Atom *) iter->get_atom();
		at->getPosition(pos);

		// Apply position increment
		pos[0] += cevec[imod*3*nla+3*kp  ] *factor;
		pos[1] += cevec[imod*3*nla+3*kp+1] *factor;
		pos[2] += cevec[imod*3*nla+3*kp+2] *factor;

		// Set position
		at->setPosition(pos);

		kp++; // Update atom index
	}


}

void move_loop_linear_steps(char *file, Macromolecule *mol, double *cevec, int imod, int ifpa, int nla, int steps, float factor)
{
	Atom *at;

	pdbIter *iter = new pdbIter( mol, true, true, true, true ); // iter to screen fragments (residues)

	// Generate a linear trajectory, including Ct-anchor for checking purposes...
	float increment = factor / (float)steps; // 0.05; // 0.05;
	Tcoor pos;
	for(int s = 0; s < steps; s++)
	{
		int kp = 0;
		for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
		{
			// Get PDB residue index
			at = (Atom *) iter->get_atom();
			at->getPosition(pos);

			// Apply position increment
			pos[0] += cevec[imod*3*nla+3*kp  ] * increment;
			pos[1] += cevec[imod*3*nla+3*kp+1] * increment;
			pos[2] += cevec[imod*3*nla+3*kp+2] * increment;

			// Set position
			at->setPosition(pos);

			kp++; // Update atom index
		}

		// Dump current conformation (step) to a Multi-PDB file
		mol->writeMPDB(file,s+1);
	}

	delete iter;
}


float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla)
{

	Tcoor pos,pos2;
	double rmsd = 0.0;
	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		iter2->pos_atom = iter->pos_atom;
		(iter2->get_atom())->getPosition(pos2);

		rmsd += pow(pos2[0] - pos[0], 2) + pow(pos2[1] - pos[1], 2) + pow(pos2[2] - pos[2], 2);
	}

	return sqrt(rmsd / nla);
}

float rmsd_loop_residue(pdbIter *iter, pdbIter *iter2, int ifr, int ifr2, int nla)
{

	Tcoor at, at2;
	double rmsd = 0.0;
	Residue *res, *res2;
	pdbIter *itera, *iterb;

	//fprintf(stderr,"--> %d %d %d\n", ifr,ifr2,nla );

	iter2->pos_fragment = ifr2;

	for( iter->pos_fragment = ifr; iter->pos_fragment <= ifr + nla; iter->next_fragment(),iter2->next_fragment()  ) // just screen mobile loop fragments
	{
		res =  ( Residue * ) iter->get_fragment();
		res2 = ( Residue * ) iter2->get_fragment();

		itera = new pdbIter( res );
		iterb = new pdbIter( res2 );

		for(int i=0; i<3; i++) {
			itera->pos_atom = i;
			( itera->get_atom() )->getPosition( at );
			iterb->pos_atom = i;
			( iterb->get_atom() )->getPosition( at2 );
			//fprintf(stderr,"--> %d %f %f %f %f %f %f   \n", iter->pos_fragment, at[0],at[1],at[2],at2[0],at2[1],at2[2]);

			rmsd += pow(at[0] - at2[0], 2) + pow(at[1] - at2[1], 2) + pow(at[2] - at2[2], 2);

		}
		itera->~pdbIter(); // delete iterator
		iterb->~pdbIter(); // delete iterator


	}

	return sqrt(rmsd / (nla*3));

}




float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, int ifpa2)
{
	Tcoor pos,pos2;
	double rmsd = 0.0;
	int index2 = 0; // loop atom index in 2nd loop
	for( iter->pos_atom = ifpa, iter2->pos_atom = ifpa2 ; iter->pos_atom < ifpa + nla; iter->next_atom(), iter2->next_atom() ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		(iter2->get_atom())->getPosition(pos2);

		rmsd += pow(pos2[0] - pos[0], 2) + pow(pos2[1] - pos[1], 2) + pow(pos2[2] - pos[2], 2);
	}

	return sqrt(rmsd / nla);
}

float rmsd_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta)
{

	Tcoor pos,pos2;
	double rmsd = 0.0;
	int k;
	for( iter->pos_atom = ifpa, k=0; iter->pos_atom < ifpa + nla; iter->next_atom(), k++ ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		iter2->pos_atom = iter->pos_atom;
		(iter2->get_atom())->getPosition(pos2);

		rmsd += pow(pos2[0] - pos[0], 2) + pow(pos2[1] - pos[1], 2) + pow(pos2[2] - pos[2], 2);

		// Compute delta
		delta[ k * 3 ]     = pos2[0] - pos[0];
		delta[ k * 3 + 1 ] = pos2[1] - pos[1];
		delta[ k * 3 + 2 ] = pos2[2] - pos[2];
	}

	return sqrt(rmsd / nla);
}

// Get dihedral angles (an array) from a continuous Macromolecule loop
// WARNING! Dihedral angles order: Psi_NtAnchor, Phi_1st, Psi_1st, Phi_2nd, Psi_2nd, ... , Phi_Nth, Psi_Nth, Phi_CtAnchor
//	iter --> Macromolecule iterator
//	ifr  --> Index (internal) of First Residue
//	nlr  --> Number of (mobile) loop residues
//	p_dihedrals --> Pointer to the array of dihedrals (=NULL forces automatic memory allocation)
void loop_dihedrals(pdbIter *iter, int ifr, int nlr, float **p_dihedrals)
{
	float *dihedrals;
	Residue *res;
	pdbIter *iter3;
	Tcoor atN, atCA, atC, atNn, atCp; // atNn (Next Nt atom), atCp (Previous Ct atom)
	float phi = 9999.0;
	float psi;

	if(*p_dihedrals == NULL) // Not allocated
	{
		*p_dihedrals = (float *) malloc( sizeof( float ) * (nlr+1) * 2 ); // 2 dihedrals (Phi and Psi) per residue (+1 extra dihedral per end)
	}
	// Already allocated
	dihedrals = *p_dihedrals;

	int index = 0; // dihedral index
	for( iter->pos_fragment = ifr-1; iter->pos_fragment <= ifr + nlr; iter->next_fragment() ) // just screen mobile loop fragments
	{
		res = ( Residue * ) iter->get_fragment();

		iter3 = new pdbIter( res );

		if(iter->pos_fragment == ifr-1) // Previous-to-Nt residue (Nt anchor)
		{
			iter3->pos_atom = 0; // N
			( iter3->get_atom() )->getPosition( atN );
			iter3->pos_atom = 1; // CA
			( iter3->get_atom() )->getPosition( atCA );
			iter3->pos_atom = 2; // C
			( iter3->get_atom() )->getPosition( atC );

			// Current C atom will be the previous C atom at next iteration
			atCp[0] = atC[0];
			atCp[1] = atC[1];
			atCp[2] = atC[2];
		}
		else
		{
			iter3->pos_atom = 0; // N (next)
			( iter3->get_atom() )->getPosition( atNn );

			// Compute current Psi angle
			// psi = dihedral_bk( atN, atCA, atC, atNn ) * PDB_PI / 180.0;
			psi = dihedral_bk( atN, atCA, atC, atNn );
			dihedrals[ index ] = psi;
			index++;

			// Update current N
			atN[0] = atNn[0];
			atN[1] = atNn[1];
			atN[2] = atNn[2];

			iter3->pos_atom = 1; // CA
			( iter3->get_atom() )->getPosition( atCA );
			iter3->pos_atom = 2; // C
			( iter3->get_atom() )->getPosition( atC );

			// Compute current Phi angle
			// phi = dihedral_bk( atCp, atN, atCA, atC ) * PDB_PI / 180.0;
			phi = dihedral_bk( atCp, atN, atCA, atC );
			dihedrals[ index ] = phi;
			index++;

			// Current C atom will be previous at next iteration
			atCp[0] = atC[0];
			atCp[1] = atC[1];
			atCp[2] = atC[2];

		}

		iter3->~pdbIter(); // delete iterator
	}

	// Show dihedrals for debugging
	fprintf(stderr,"loop_dihedrals [deg]> ");
	for(int i=0; i < nlr*2; i++)
		fprintf(stderr," %5.1f", dihedrals[i]); // in degrees
	fprintf(stderr,"\n");
}

// Compute the Dihedral angles RMSD between two protein loops
float rmsd_dihedral_loop(pdbIter *iter, pdbIter *iter2, int ifr, int ifr2, int nlr, float *dhs, float *dhs2)
{
	bool remove_array = false;
	bool remove_array2 = false;

	// If "dhs" arrays are NULL, memory must be allocated
	remove_array = dhs == NULL;
	remove_array2 = dhs2 == NULL;

	// Get dihedral angles (an array) from a continuous Macromolecule loop
	// WARNING! Dihedral angles order: Psi_NtAnchor, Phi_1st, Psi_1st, Phi_2nd, Psi_2nd, ... , Phi_Nth, Psi_Nth, Phi_CtAnchor
	loop_dihedrals(iter, ifr, nlr, &dhs);
	loop_dihedrals(iter2, ifr2, nlr, &dhs2);

	// Only remove arrays if not provided by user
	if(remove_array)
		free(dhs);
	if(remove_array2)
		free(dhs2);

	return vector_rmsd(dhs, dhs2, nlr*2 + 2);
}

// Compute the RMSD between two arrays
float vector_rmsd(float *array, float *array2, int size)
{
	double accum = 0.0;
	for(int i=0; i<size; i++)
	{
		accum += powf(array[i] - array2[i], 2);
	}

	return( (float) sqrt(accum / size) );
}

// Compute the RMSD for a given array of increments
float vector_rmsd(double *array, int size)
{
	double accum = 0.0;
	for(int i=0; i<size; i++)
	{
		accum += powf(array[i], 2);
	}

	return( (float) sqrt(accum / size) );
}

void delta_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta)
{
	Tcoor pos,pos2;
	int k;
	for( iter->pos_atom = ifpa, k=0; iter->pos_atom < ifpa + nla; iter->next_atom(), k++ ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		iter2->pos_atom = iter->pos_atom;
		(iter2->get_atom())->getPosition(pos2);

		// Compute delta
		delta[ k * 3 ]     = pos2[0] - pos[0];
		delta[ k * 3 + 1 ] = pos2[1] - pos[1];
		delta[ k * 3 + 2 ] = pos2[2] - pos[2];
	}
}

void delta_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, double *delta, int ifpa2)
{
	Tcoor pos,pos2;
	int k;
	for( iter->pos_atom = ifpa, k=0, iter2->pos_atom = ifpa2; iter->pos_atom < ifpa + nla; iter->next_atom(), k++, iter2->next_atom()) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		(iter2->get_atom())->getPosition(pos2);

		// Compute delta
		delta[ k * 3 ]     = pos2[0] - pos[0];
		delta[ k * 3 + 1 ] = pos2[1] - pos[1];
		delta[ k * 3 + 2 ] = pos2[2] - pos[2];
	}
}

// Compute RMSD from "iffa" to "ifa" and from "ila" to "ilfa"
// iffa2 --> index of flanking first atom in target loop (to consider different number of residues PDBs)
float rmsd_flank(pdbIter *iter, pdbIter *iter2, int iffa, int ifa, int ila, int ilfa, int iffa2)
{
	bool debug = false;
	Tcoor pos,pos2;
	double rmsd = 0.0; // deviations sum
	int index2 = 0; // loop atom index in 2nd loop
	int natoms = 0; // number of atoms

	for( iter->pos_atom = iffa, iter2->pos_atom = iffa2 ; iter->pos_atom < ifa; iter->next_atom(), iter2->next_atom() ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		(iter2->get_atom())->getPosition(pos2);

		if(debug) // dump some debug info
			fprintf(stderr,"Izq A: %10f %10f %10f    B: %10f %10f %10f\n", pos[0], pos[1], pos[2], pos2[0], pos2[1], pos2[2]);

		rmsd += pow(pos2[0] - pos[0], 2) + pow(pos2[1] - pos[1], 2) + pow(pos2[2] - pos[2], 2);
		natoms++; // counting atoms
	}

	for( iter->pos_atom = ila+1, iter2->pos_atom = ila+1 + (iffa2-iffa); iter->pos_atom <= ilfa; iter->next_atom(), iter2->next_atom() ) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);
		(iter2->get_atom())->getPosition(pos2);

		if(debug) // dump some debug info
			fprintf(stderr,"Der A: %10f %10f %10f    B: %10f %10f %10f\n", pos[0], pos[1], pos[2], pos2[0], pos2[1], pos2[2]);

		rmsd += pow(pos2[0] - pos[0], 2) + pow(pos2[1] - pos[1], 2) + pow(pos2[2] - pos[2], 2);
		natoms++; // counting atoms
	}

	return sqrt(rmsd / natoms);
}


// Does some loop atom clash with its environment?
// 	iter,iter2 --> 2 differerent iterators to the same macromolecule
//	ifpa       --> Index of first (pseudo)atom
//	nla        --> Number of loop atoms
//	cut2       --> cutoff distance squared
bool clashed_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, float cut2)
{
	Tcoor pos, pos2; // Position vectors
	float d2; // Distance

	// Loop vs. Environment
	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom()) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);

		for( iter2->pos_atom = 0; iter2->pos_atom < ifpa; iter2->next_atom()) // just screen left-side of environment
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;

			//			if( (d2 = powf(pos2[0] - pos[0], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[1] - pos[1], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[2] - pos[2], 2)) <= cut2)
			//				return true;
		}

		for( iter2->pos_atom = ifpa + nla; !iter2->gend_atom(); iter2->next_atom()) // just screen right-side of environment
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;
			//d2 =+ d2;

			//			d2 = powf(pos2[0] - pos[0], 2);
			//			if( d2 <= cut2 )
			//				return true;
			//
			//			d2 += powf(pos2[1] - pos[1], 2);
			//			if( d2 <= cut2 )
			//				return true;
			//
			//			d2 += powf(pos2[2] - pos[2], 2);
			//			if( d2 <= cut2)
			//				return true;
		}

		for( iter2->pos_atom = iter->pos_atom + 1; iter2->pos_atom < ifpa + nla - 1; iter2->next_atom()) // loop vs. loop (non-redundant)
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;

			//			if( (d2 = powf(pos2[0] - pos[0], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[1] - pos[1], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[2] - pos[2], 2)) <= cut2)
			//				return true;
		}

	}

	return false; // NO clashes detected
}


// Does some loop atom clash with its environment?
// 	iter,iter2 --> 2 differerent iterators to the same macromolecule
//	ifpa       --> Index of first (pseudo)atom in 1st PDB
//	nla        --> Number of loop atoms
//	cut2       --> cutoff distance squared
//  ifpa2      --> Index of first (pseudo)atom in 2nd PDB
bool clashed_loop(pdbIter *iter, pdbIter *iter2, int ifpa, int nla, float cut2, int ifpa2)
{
	Tcoor pos, pos2; // Position vectors
	float d2; // Distance

	// Loop vs. Environment
	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom()) // just screen mobile loop atoms
	{
		(iter->get_atom())->getPosition(pos);

		for( iter2->pos_atom = 0; iter2->pos_atom < ifpa2; iter2->next_atom()) // just screen left-side of environment
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;

			//			if( (d2 = powf(pos2[0] - pos[0], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[1] - pos[1], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[2] - pos[2], 2)) <= cut2)
			//				return true;
		}

		for( iter2->pos_atom = ifpa2 + nla; !iter2->gend_atom(); iter2->next_atom()) // just screen right-side of environment
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;
			//d2 =+ d2;

			//			d2 = powf(pos2[0] - pos[0], 2);
			//			if( d2 <= cut2 )
			//				return true;
			//
			//			d2 += powf(pos2[1] - pos[1], 2);
			//			if( d2 <= cut2 )
			//				return true;
			//
			//			d2 += powf(pos2[2] - pos[2], 2);
			//			if( d2 <= cut2)
			//				return true;
		}

		for( iter2->pos_atom = iter->pos_atom + 1 + ifpa2-ifpa; iter2->pos_atom < ifpa2 + nla - 1; iter2->next_atom()) // loop vs. loop (non-redundant)
		{
			(iter2->get_atom())->getPosition(pos2);

			d2 = powf(pos2[0] - pos[0], 2) + powf(pos2[1] - pos[1], 2) + powf(pos2[2] - pos[2], 2);
			// fprintf(stderr,"d2= %f\n", sqrtf(d2) );

			if(d2 <= cut2)
				return true;

			//			if( (d2 = powf(pos2[0] - pos[0], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[1] - pos[1], 2)) <= cut2)
			//				return true;
			//			if( (d2 += powf(pos2[2] - pos[2], 2)) <= cut2)
			//				return true;
		}

	}

	return false; // NO clashes detected
}


// Scale eigenvectors so that the maximum value of the components is "maxang"
void scale_vectors(double *eigvect, int size, int neig, double maxang)
{
	bool verb = false;
	double maxnorm;
	double temp;

	if(verb)
		fprintf(stdout,"Scaling eigenvectors so that the maximum value is maxang: %f [deg]\n",maxang);

	// maxang *= M_PI / 180.0; // Eigenvectors must be in radians!

	// Normalize all modes
	for(int i=0; i < neig;i++) // Screen all modes
	{
		maxnorm = 0.0;

		for(int j=0;j<size;j++)
			if( ( temp = fabs( eigvect[i*size + j] ) ) > maxnorm )
				maxnorm = temp; // store the maximum

		for(int j=0;j<size;j++)
			eigvect[i*size + j] *= maxang / maxnorm; // the maximum component will have "maxang" value
	}
}

// Compute N-dimensional vector modulus (i.e. vector length)
double vector_modulus(double *v, int size)
{
	double mod = 0.0;
	for(int i=0; i<size; i++)
		mod += pow(v[i],2);
	return sqrt(mod);
}

// Compute N-dimensional vector modulus (i.e. vector length)
float vector_modulus(float *v, int size)
{
	float mod = 0.0;
	for(int i=0; i<size; i++)
		mod += powf(v[i],2);
	return sqrtf(mod);
}

// Element-wise difference between vectors "v" and "w" (d = v - w)
void vector_diff(float *v, float *w, double *d, int size)
{
	for(int i=0; i<size; i++)
		d[i] = v[i] - w[i];
}

// Element-wise difference between dihedral angle [deg] arrays "v" and "w" (d = v - w) taking into account rotation
void dihedrals_diff(float *v, float *w, double *d, int size)
{
	for(int i=0; i<size; i++)
	{
		d[i] = v[i] - w[i];

		//		if( fabs(d[i]) > 180 )
		//		{
		//			if(d[i] >= 180)
		//				d[i] = w[i] - v[i];
		//			else if(d[i] < 180)
		//				d[i] = w[i] - v[i];
		////				d[i] += 180;
		//		}

		//		if( fabs(d[i]) > 180 )
		//		{
		//			if(d[i] >= 180)
		//				d[i] = -(d[i]- 180);
		//			else if(d[i] < 180)
		//				d[i] = -(d[i] + 180);
		//		}

		if( fabs(d[i]) > 180 )
		{
			if(d[i] >= 180)
				d[i] -= 360;
			else
				d[i] += 360;
		}

	}
}

// Compute N-dimensional dot-product between two vectors (i.e. the modulus of the projection between both vectors)
double dotprod(double *v, double *w, int size)
{
	double dot = 0.0;
	for(int i=0; i<size; i++)
		dot += v[i] * w[i];
	return dot;
}

// Compute the normalized (0,1) N-dimensional dot-product between two vectors (i.e. the cosine of the angle between both vectors)
double dotprodnorm(double *v, double *w, int size, double mv, double mw)
{
	double dot = 0.0;
	for(int i=0; i<size; i++)
		dot += v[i] * w[i];
	// fprintf(stderr, "mv= %f  mw= %f\n", mv, mw);
	if(mv == 0.0)
		mv = vector_modulus(v,size); // Compute "v" modulus for normalization, if not already done!
	if(mw == 0.0)
		mw = vector_modulus(w,size); // Compute "w" modulus for normalization, if not already done!
	dot /= mv * mw;
	return dot;
}

// Compute the normalized (0,1) N-dimensional dot-product between two vectors (i.e. the cosine of the angle between both vectors)
float dotprodnorm(float *v, float *w, int size, float mv, float mw)
{
	float dot = 0.0;
	for(int i=0; i<size; i++)
		dot += v[i] * w[i];
	// fprintf(stderr, "mv= %f  mw= %f\n", mv, mw);
	if(mv == 0.0)
		mv = vector_modulus(v,size); // Compute "v" modulus for normalization, if not already done!
	if(mw == 0.0)
		mw = vector_modulus(w,size); // Compute "w" modulus for normalization, if not already done!
	dot /= mv * mw;
	return dot;
}

// Get array of atomic masses from a Macromolecule iterator (num_atoms is for cross-checking purposes)
float *get_masses(pdbIter *iter, int num_atoms)
{
	Atom *at;

	// Allocate the array of masses...
	float *masses; // masses array
	masses = (float *) malloc( sizeof(float) * num_atoms);
	ptr_check(masses);

	// Compute the array of masses...
	for( iter->pos_atom = 0; !iter->gend_atom(); iter->next_atom() ) // screen atoms
	{
		// Get PDB residue index
		at = (Atom *) iter->get_atom();
		masses[iter->pos_atom] = at->getPdbocc(); // Get mass from Macromolecule
	}
	// Some checking
	if(iter->pos_atom != num_atoms)
	{
		fprintf(stderr,"Sorry, there is some mismatch between number of atoms... %d != %d  Forcing exit!", iter->pos_atom, num_atoms);
		exit(1);
	}

	return masses;
}


// Moves atoms given an Internal Coordinates Normal Mode and a Step
// Simple Rotations Scheme (Fast & Multi-Chain/Prot/DNA/RNA/SMOL): 3BB2R & Full-Atom models
// Ref.: Choi, "On Updating Torsion Angles of Molecular Conformations" (2006).
inline void move_loop_dihedral(pdbIter *iter, int ifr, int ilr, tri *props, double *uu, int size, int model, float step)
{
	bool debug = false;
	int j,j2;
	Tcoor e1,nh,ca,co,tv,tvp,T,Ti,Ti2,Ti3,tr;
	Residue *res;
	Atom *atom;
	pdbIter *iter_res_atom;
	double **Ri; // current rotation (Ri-matrix);
	double **Si; // accummulated rotation (Si-matrix)
	double **dummy; // dummy rotation (dummy-matrix)
	double **temp;
	int num_res = 0;
	TMOL fragtype; // enum tmol_protein,tmol_rna,tmol_dna,...
	int indexbase;
	int resn;

	num_res = iter->num_fragment(); // gets number of fragments per current segment (to detect ending)

	// Ri-matrix, Si-matrix, & dummy-matrix allocation
	Ri = (double **) malloc( sizeof(double *) * 3 );
	Si = (double **) malloc( sizeof(double *) * 3 );
	dummy = (double **) malloc( sizeof(double *) * 3 );
	for(int i=0; i<3; i++)
	{
		Ri[i] = (double *) malloc( sizeof(double) * 3 );
		Si[i] = (double *) malloc( sizeof(double) * 3 );
		dummy[i] = (double *) malloc( sizeof(double) * 3 );
	}

	// Ri-matrix, Si-matrix, & dummy-matrix initialization
	for(int i=0; i<3; i++)
	{
		// 1st Ti --> (0,0,0) vector
		Ti[i] = 0.0;
		Ti2[i] = 0.0;

		// 1st Si --> I-matrix
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

	j = 0;
	j2 = 0;
	for( iter->pos_fragment = ifr; iter->pos_fragment <= ilr; iter->next_fragment() ) // screen residues
	{
		res = ( Residue * ) iter->get_fragment();
		iter_res_atom = new pdbIter( res ); // iter residue atoms
		resn = resnum_from_resname( res->getName() );

		if(debug)
			fprintf(stdout,"%4d Residue %s  nan= %d\n",j, res->getName(), props[iter->pos_fragment].nan);

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
				fprintf(stdout,"%4d PHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

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

			// rotating current-PHI affected residue atoms (due to backwards dihedrals and current-Phi)
			for(iter_res_atom->pos_atom = 2; !iter_res_atom->gend_atom(); iter_res_atom->next_atom())
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
		// Note that "resn" is available, you can check numerically Aminoacid identities!!! (faster)
		//		if( type == 2 )
		if( props[iter->pos_fragment].nan == 3 )
			//				if( (props[iter->pos_fragment].nan==3 && iter->pos_fragment != ilr) ||
			//						(props[iter->pos_fragment].nan==2 &&
			//								(iter->pos_fragment==0 ||
			//										(iter->pos_fragment==num_res-1 && model==1) ) ) )
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
				fprintf(stdout,"%4d CHI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

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
				// With 3BB2R, just "O" moves!!!
				iter_res_atom->pos_atom = 4;
				atom = (Atom *) iter_res_atom->get_atom();
				atom->getPosition( tv ) ; // position before rotation
				rmot(Ri, nh, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
				atom->setPosition( tvp ) ; // position after rotation
			}

			j2++;
			j++;
		}

		// PSI
		if( iter->pos_fragment != num_res-1 || model == 2 ) // Segment Last residue never has PSI
		{													// (excepting in Full-Atom!)
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
				fprintf(stdout,"%4d PSI (Updating Ti & Si) angle[%3d]= %f\n",j,j2,uu[j2]*step);

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

			// During PSI-rotation, no current residue atoms are moved !!!

			// Excepting in Full-Atom model!
			// In Full-Atom, the O-atom moves
			if(model == 2)
			{
				// rotating O-atom
				iter_res_atom->pos_atom = 3;
				if( !iter_res_atom->gend_atom() ) // Check whether "O" atom exists
				{
					atom = (Atom *) iter_res_atom->get_atom();
					atom->getPosition( tv ) ; // position before rotation
					rmot(Si, Ti, tv, tvp); // Rigid Motion [ (4x4) x (4x1) matrix multiplication ]
					atom->setPosition( tvp ) ; // position after rotation
				}
			}
			j++;
		}

		delete iter_res_atom;
	}
	// delete iter;

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

// Creates contacts list (IPA) for some loop (intra-loop + loop vs. environment) from two iterators pointing to the same Macromolecule.
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
// MON: This works!!! Either use it instead of the other, or modify it accordingly!
// Creates contacts list (IPA) for some loop (intra-loop + loop vs. environment) from two iterators pointing to the same Macromolecule.
// (The "masses", i.e. occupancies, with zero mass will not be accounted for.)
inline void make_ipas_loop(pdbIter *iterA, pdbIter *iterB, int ifpa, int nla, float cutoff, twid **p_decint, int *p_nipa)
{
	bool debug = false;
	int nipa, index, k, l;
	double d2;
	float cutoff2;
	cutoff2 = cutoff*cutoff; // to speed-up distance evaluations
	twid *decint;
	decint = *p_decint;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	Atom *atA,*atB;
	Tcoor rA,rB;
	// pdbIter *iterA,*iterB;
	// iterA = new pdbIter(mol);
	// iterB = new pdbIter(mol);

	// MON: Bug? (12/4/2022) --> WTF!!! OMG!!!
	// decint = ( twid * ) realloc( decint,  nla*120 * sizeof( twid ) ); //  60 contacts per atom
	decint = ( twid * ) malloc(  nla*500 * sizeof( twid ) ); //  60 contacts per atom

	if(debug)
		fprintf(stdout,"Msg(make_ipas_new): Creating Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

	for(iterA->pos_atom = ifpa; iterA->pos_atom < ifpa + nla; iterA->next_atom() )
	{
		atA = iterA->get_atom();
		if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
		{
			atA->getPosition(rA);
			//  for(iterB->pos_atom = iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom() )
			for(iterB->pos_atom = 0; !iterB->gend_atom(); iterB->next_atom() )
			{
				// MON: Bad Bug? (12/4/2022) --> following lines prevent considering k<l interactions (k=A, l=B)
				// if(iterA->pos_atom >= iterB->pos_atom) // Avoids counting twice intra-loop contacts
				// continue;

				if(iterB->pos_atom >= ifpa && iterB->pos_atom < ifpa+nla && iterA->pos_atom >= iterB->pos_atom) // Avoids counting twice intra-loop contacts
					continue;

				atB = iterB->get_atom();
				if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
				{
					atB->getPosition(rB);
					// d2 = pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2);
					// if (nipa> 400) printf ("nipa ???  %d %d\n", nipa, nla*500);

					if( (d2 = pow(rA[0]-rB[0],2)) > cutoff2) continue;
					else if((d2+= pow(rA[1]-rB[1],2)) > cutoff2)  continue;
					else if((d2+= pow(rA[2]-rB[2],2)) <= cutoff2)
					{
						nipa++; // Counts number of Interacting Pairs of Atoms
						// decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
						decint[nipa - 1].k = iterA->pos_atom; // k-pseudo-atom index (i-atom index)
						decint[nipa - 1].l = iterB->pos_atom; // l-pseudo-atom index (j-atom index)
						decint[nipa - 1].d = d2; // sqrtf(d2); // set distance
						decint[nipa - 1].C = 0.0; // force constant will be set in the future
					}
				}
			}
		}
	}
	//  delete iterA;
	//  delete iterB;

	if(debug)
		printf( "Msg(make_ipas_loop): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d   %f per residue %d \n", nipa, 1.0*nipa/(nla), nla );

	*p_nipa = nipa; // outputs "nipas"
	*p_decint = decint;
}

inline void make_ipas_loop(pdbIter *iterA, pdbIter *iterB, int ifpa, int nla, float cutoff, twid *decint, int *p_nipa)
{
	bool debug = false;
	int nipa, index, k, l;
	double d2;
	float cutoff2;
	cutoff2 = cutoff*cutoff; // to speed-up distance evaluations
	//twid *decint;
	//decint = *p_decint;
	nipa = 0; // counts the number of interacting pseudo-atom pairs
	Atom *atA,*atB;
	Tcoor rA,rB;
	// pdbIter *iterA,*iterB;
	// iterA = new pdbIter(mol);
	// iterB = new pdbIter(mol);
	// decint = ( twid * ) realloc( decint, (ifpa + nla)*120 * sizeof( twid ) ); // resizes contact list-structure

	if(debug)
		fprintf(stdout,"Msg(make_ipas_new): Creating Interacting Pairs of pseudo-Atoms (IPAs) list.\n");

	for(iterA->pos_atom = ifpa; iterA->pos_atom < ifpa + nla; iterA->next_atom() )
	{
		atA = iterA->get_atom();
		if(atA->getPdbocc() != 0.0) // if it's not a virtual atom (modified 24/11/2009)
		{
			atA->getPosition(rA);
			//		  for(iterB->pos_atom = iterA->pos_atom+1; !iterB->gend_atom(); iterB->next_atom() )
			for(iterB->pos_atom = 0; !iterB->gend_atom(); iterB->next_atom() )
			{


				if(iterB->pos_atom >= ifpa && iterB->pos_atom < ifpa+nla && iterA->pos_atom >= iterB->pos_atom) // Avoids counting twice intra-loop contacts
					continue;

				atB = iterB->get_atom();
				if(atB->getPdbocc() != 0.0)  // if it's not a virtual atom
				{
					atB->getPosition(rB);
					//d2 = pow(rA[0]-rB[0],2) + pow(rA[1]-rB[1],2) + pow(rA[2]-rB[2],2);
					//if(d2 <= cutoff2)

					if( (d2 = pow(rA[0]-rB[0],2)) > cutoff2) continue;
					else if((d2+= pow(rA[1]-rB[1],2)) > cutoff2)  continue;
					else if((d2+= pow(rA[2]-rB[2],2)) <= cutoff2)
					{
						nipa++; // Counts number of Interacting Pairs of Atoms
						// if (nipa> 500) printf ("nipa %d\n", nipa);
						//if (d2<1.0) printf ("nipa %d %d %f\n",iterA->pos_atom,iterB->pos_atom, nipa);

						// decint = ( twid * ) realloc( decint, nipa * sizeof( twid ) ); // resizes contact list-structure
						decint[nipa - 1].k = iterA->pos_atom; // k-pseudo-atom index (i-atom index)
						decint[nipa - 1].l = iterB->pos_atom; // l-pseudo-atom index (j-atom index)
						decint[nipa - 1].d = d2; // sqrtf(d2); // set distance // ojo
						decint[nipa - 1].C = 0.0; // force constant will be set in the future
					}
				}
			}
		}
	}
	//  delete iterA;
	//  delete iterB;

	if(debug)
		printf( "Msg(make_ipas_loop): Number of Interacting Pairs of pseudo-Atoms (NIPAs): %d   %f per residue \n", nipa, 1.0*nipa/(ifpa + nla) );


	*p_nipa = nipa; // outputs "nipas"
}



// Get the array of atomic masses for some loop from one macromolecule iterator that points to a Macromolecule.
//	sqrt --> Set "true" to compute the square root of the masses (to mass-weight Cartesian modes), otherwise it just gets the masses
float *get_masses_loop(pdbIter *iter, int ifpa, int nla, bool sqrt)
{
	Atom *at;
	float *masses = (float *) malloc( sizeof( float ) * nla ); // Allocate memory

	int i;
	if(sqrt) // The square root of the masses
		for(iter->pos_atom = ifpa, i=0; iter->pos_atom < ifpa + nla; iter->next_atom(), i++ )
		{
			at = iter->get_atom();
			masses[i] = sqrtf( at->getPdbocc() );
		}
	else // Just the masses
		for(iter->pos_atom = ifpa, i=0; iter->pos_atom < ifpa + nla; iter->next_atom(), i++ )
		{
			at = iter->get_atom();
			masses[i] = at->getPdbocc();
		}

	return masses;
}

// Creates contacts list (IPA) for some loop (intra-loop + loop vs. environment) from two iterators pointing to the same Macromolecule.
//  iter   --> Iterator of first Macromolecule (e.g. reference)
//  iter2  --> Iterator of moving Macromolecule (e.g. target)
//	ilr    --> Index of last residue (internal numeration)
//  p_dist --> OUTPUT distance between C atom of last loop residue and the N atom of Ct-anchor
//  p_ang  --> OUTPUT bond angle between C atom of last loop residue and the N and C atoms of Ct-anchor
void anchor_drift(pdbIter *iter, pdbIter *iter2, tri *props, int ilr, float *p_dist, float *p_ang)
{
	Atom *at;
	Tcoor n, c, ca, c2, c3;

	int ifalr = props[ilr].k1; // Index of the First Atom of the Last Residue
	int ifaa = props[ilr+1].k1; // Index of the First Atom of the Anchor


	// Get coordinates of the C atom from last residue of target macromolecule
	iter2->pos_atom = ifalr + 2; // C-atom index
	(iter2->get_atom())->getPosition(c);

	// Get coordinates of the N atom from Ct-anchor of reference macromolecule
	iter->pos_atom = ifaa; // N-atom index
	(iter->get_atom())->getPosition(n);

	// Get coordinates of the CA atom from Ct-anchor of reference macromolecule
	iter->pos_atom = ifaa + 1; // CA-atom index
	(iter->get_atom())->getPosition(ca);

	// compute drift...
	/*
	iter->pos_atom = ifaa +2 ;
	(iter->get_atom())->getPosition(c2);
    iter2->pos_atom = ifaa + 2; // C-atom index
    (iter2->get_atom())->getPosition(c3);
    fprintf(stdout,"drift %f C atom from last residue\n", sqrtf( powf(c2[0]-c3[0],2) + powf(c2[1]-c3[1],2) + powf(c2[2]-c3[2],2) ) );
    fprintf(stdout,"C1 %f %f %f\n", c2[0], c2[1], c2[2] );
    fprintf(stdout,"C  %f %f %f\n", c3[0], c3[1], c3[2] );
	 */

	// Compute inter-atomic distance
	*p_dist = sqrtf( powf(n[0]-c[0],2) + powf(n[1]-c[1],2) + powf(n[2]-c[2],2) );

	// Compute vector CA-->N
	n[0] -= ca[0];
	n[1] -= ca[1];
	n[2] -= ca[2];

	// Compute vector CA-->C
	c[0] -= ca[0];
	c[1] -= ca[1];
	c[2] -= ca[2];

	// Compute bond angle
	float dot = dotprodnorm(n, c, 3);
	*p_ang = acosf(dot) * 180 / M_PI;
}

// Mode following routine
//	refmode0 --> Reference CC mode to drive motion (mode following)
//	mol --> Input macromolecule
//	model --> Atomic model
//	type --> ICs type (with/without chi)
//	props --> Properties structure
//	masses --> Atomic masses of the whole protein
//	masses_loop --> Atomic masses of the loop
//	ifa --> Index of first atom of the loop
//	ifr --> Index of first residue of the loop
//	ilr --> Index of the last residue of the loop
//	na --> Number of atoms of the loop
//	size --> Number of ICs of the loop
//	nco --> Number of constraints (typically 6)
//	cutoff --> Distance cutoff for elastic network definition
//	nsamples --> Number of samples (frames)
//	eigval --> IC eigenvalues, automatic memory allocation if NULL (OUTPUT)
//	eigvect --> IC eigenvectors, automatic memory allocation if NULL (OUTPUT)
//	rmsd_conv --> (OPTIONAL) Motion amplitude [RMSD] (motion will stop upon "nsamples" or if "rmsd_conv" RMSD is reached)
//	iterini --> (OPTIONAL) Reference pdb iterator to compute RMSD
//	file_movie --> (OPTIONAL) Output Multi-PDB file
//	delta_rmsd --> (OPTIONAL) Delta RMSD in mode-following strategies (high value means disabled)
//	p_fi --> (OPTIONAL) Frame index (model index in Multi-PDB file). Automatically updated inside.
//	chain --> (OPTIONAL) Chain-ID for output Multi-PDB
//	update --> (OPTIONAL) if "true", the "refmode" will be updated by current most overlapping mode
//	RETURN --> Number of non-null eigenpairs
inline int follow_mode(double *refmode0, Macromolecule *mol, int model, int type, tri *props, float *masses, float *masses_loop, int ifa, int ifr, int ilr,
		int na, int size, int nco, double maxang, float cutoff, int nsteps, int max_loops_save, double *eigval, double *eigvect, float rmsd_conv, pdbIter *iterini,
		char *file_movie, float delta_rmsd, int *p_fi, char chain,  bool update)
{
	bool debug = false; // dump debug info
	float *coord;
	int nipa;
	int neig; // Number of non-null eigenvectors
	int ncomps = 3*(na+3); // Number of Cartesian components (x,y,z for each atom of mobile loop + 3 Ct-anchor atoms)

	double *alpha = (double *) malloc( (size+nco) * sizeof(double)); // Alpha overlaps array (allocate memory for the maximum possible)
	double *delta = (double *) malloc( (size+nco) * sizeof(double)); // Delta overlaps array (allocate memory for the maximum possible)
	double *mode = (double *) malloc( size * sizeof(double)); // Merged mode for motion

	double *refmode;

	//	Update the "refmode" and "modref" with the current most overlapping mode
	if(update)
	{
		refmode = (double *) malloc( ncomps * sizeof(double)); // Allocate reference CC vector

		for(int i=0; i<ncomps; i++)
			refmode[i] = refmode0[i]; // Copy the reference mode "refmode0" into "refmode"
	}
	else
		refmode = refmode0; // Use the initial "refmode0"

	// Store the modulus of reference mode "refmode"
	double modref = vector_modulus( refmode, ncomps ); // Reference CC vector modulus


	// Create iterators for RMSD computations
	pdbIter *itermol = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule
	pdbIter *itermol2 = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule to speed up Elastic Network refresh
	int num_atoms = itermol->num_atom();

	if(eigval == NULL) // if not previously allocated
		eigval = (double *) malloc( sizeof(double) * (size + nco) );

	if(eigvect == NULL) // if not previously allocated
		eigvect = (double *) malloc( sizeof(double) * (size + nco) * size );

	float rmsd; // Current RMSD
	float last_rmsd = 0.0; // Last saved RMSD
	mol->coordMatrix( &coord );

	// Allocate derivatives
	trd *der = (trd *) malloc( (na + 3) * size * sizeof(trd) );
	ptr_check(der);

	int sizex = size + nco;
	double *mass_matrix; // kinetic energy matrix
	mass_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(mass_matrix);

	double *rdr = (double *) malloc(size * sizeof(double)); // Auxiliar array
	ptr_check(rdr);

	double *hess_matrix; // Bordered Hessian matrix
	hess_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(hess_matrix);

	double *cevec  = (double *) malloc( 3 * (na+3) * (size - nco) * sizeof(double));

	//decint = ( twid * ) realloc( decint, (ifpa + nla)*12 * sizeof( twid ) ); // resizes contact list-structure

	twid *decint = ( twid * ) malloc(   (na)*500 * sizeof( twid ) ); //  60 contacts per atom
	ptr_check(decint);

	double *alphar, *alphai, *beta, *vr,  *vl, *work;


	alphar = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphar);

	alphai = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphai);

	beta = (double *) malloc(sizex * sizeof(double));
	ptr_check(beta);

	vl = (double *) malloc(1*sizex * sizeof(double));
	ptr_check(vl);

	vr = (double *) malloc(sizex*sizex * sizeof(double));   /* eigenvectors */
	ptr_check(vr);

	work = (double *) malloc(20*sizex * sizeof(double));
	ptr_check(work);

	float ddrift = 0.0; // Distance drift at Ct end
	float adrift = 0.0; // Angle (bond angle) drift at Ct end
	float adist = 0.0; // Anchor distance
	float aang = 0.0; // Anchor angle
	float adist0 = 0.0; // Initial anchor distance
	float aang0 = 0.0; // Initial anchor angle

    anchor_drift(iterini, itermol, props, ilr, &adist0, &aang0);

    double rmsdL, rmsdP;
    rmsdL=100.0;
	for(int f = 0; f < nsteps; f++) // some maximum number of sampling steps
	{

		//if ( fabs(rmsdP-rmsdL) >= 0.02)
		{
		// Initialize Eigenvalues
		memset(eigval, 0, (size+nco)*sizeof(double));


		// Initialize Eigenvectors
		memset(eigvect, 0, (size+nco)*size*sizeof(double));


		// if(debug)
		//	fprintf(stderr, "follow_mode> Getting coordinates single row (pseudo-atom model)\n");

		update_loop_coords(itermol, props[ifr].k1, na, coord);

		//if(debug)
		//	fprintf(stderr, "follow_mode> Computing derivatives...\n");
		//der = drdqC5x(coord, props, ifr, ilr, na, size, model);
		drdqC5x(coord, props, ifr, ilr, na, size, model, der);


		//if(debug)
		//	fprintf(stderr, "follow_mode> Computing Kinetic Energy matrix (masses matrix)...\n");
		//mass_matrix = kineticC5x(der, masses, props, ifr, na, size, nco);
		kineticC5x(der, masses, props, ifr, na, size, nco, mass_matrix );


		//		if(decint != NULL)
		//			free(decint); // Free obsolete contacts list

		// Allocate "decint" to store the contacts list (ipas-list)
		//if( !(decint = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )  // Required for "realloc"
		//{
		//	fprintf(stderr, "Sorry, \"decint\" memory allocation failed!\n");
		//	exit(1);
		//}

		//if(debug)
		//	fprintf(stderr, "follow_mode> Computing Interacting Pair of (non-virtual) Atoms (ipas)\n");

		// INVERSE EXPONENTIAL (power of distance for contact matrix)
		// Making Interacting Pairs of (non-virtual) Atoms (ipas)

	//if ( fabs(rmsdP-rmsdL) >= 0.02)
				{
		make_ipas_loop(itermol, itermol2, ifa, na, cutoff, decint, &nipa); // Updating Elastic network
		// fprintf(stderr, "nipa %d\n", nipa);
		//if(debug)
		//	fprintf(stderr, "follow_mode> Inverse Exponential (%d nipas) cutoff= %.1f, k= %f, x0= %.1f ", nipa, cutoff_k0, cte_k0, x0);
		for(int i=0; i<nipa; i++)
			decint[i].C = Inv_exp( cte_k0, decint[i].d, x0, power); // setting Force Constants
		rmsdL= rmsdP;
				}
		// IPAs checking
		//if(debug) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
		//	for(int i=0; i<nipa; i++)
		//		fprintf(stderr, "ipa %4d: k= %d  l= %d  d= %f  C= %f\n",i,decint[i].k,decint[i].l,decint[i].d,decint[i].C);

		//if(debug)
		//	fprintf(stderr, "follow_mode> Computing Hessian matrix (potential energy matrix)...\n");
		//hessianC5x(coord, der, props, decint, nipa, ifr, na, size, nco, rdr, hess_matrix);
		hessianFast(coord, coord, der, props, decint, nipa, ifr, ilr, size, nco, num_atoms, hess_matrix);

		// free(decint); // free obsolete elastic network

		//if(debug)
		//{
			// Show Hessian matrix
		//	show_matrix(hess_matrix, size + nco, "Hessian:", " %7.2f");
			// Show Kinetic Energy matrix
		//	show_matrix(mass_matrix, size + nco, "Kinetic:", " %7.0f");
		//}


		// COMPUTING THE EIGENVECTORS AND EIGENVALUES
		//int info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig);
		int info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig, alphar, alphai, beta, vr,  vl, work);

		//free(mass_matrix);
		//free(hess_matrix);

		// Some checking...
		if( info ) // if info != 0
		{
			fprintf(stderr, "\nfollow_mode> An error occured in the matrix diagonalization: %d\n", info);
			exit(1);
		}

		// MON: check this, seems unnecessary...
		if(neig < nevec)
		{
			fprintf(stderr," Warning more eigenvectors requested (%d) than available (%d), forcing exit!\n",nevec,neig);
		}

//		if(debug)
//		{
//			fprintf(stderr, "follow_mode> Eigensolver successfully finished!!! (neig=%d)\n", neig);
//			show_vector(stderr,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
//			show_vectors(stderr,eigvect,size,neig,"Dumping Raw Eigenvectors:", " %5.2e");
//		}

		// Scale eigenvectors so that the maximum value of the components is "maxang"
		if(maxang != 0.0)
		{
			scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
//			if(verb > 1)
//				show_vectors(stderr,eigvect,size,neig,"Dumping Scaled Eigenvectors:", " %5.2e");
		}

		// Compute the Cartesian eigenvectors from the Internal Coordinates eigenvectors
		// double *cevec; // Cartesian eigenvectors
		// cevec = ic2cart(eigvect, neig, der, size, na + 3, masses_loop);
		cevec = ic2cart(eigvect, neig, der, size, na + 3, masses_loop);

		//

	}
		//else
	//{
		//fprintf(stdout,"jump\n");

	//}


		//		sprintf(text,"follow_mode> %4d ", f);
//		if(verb > 0)
//			fprintf(stdout,"follow_mode> %4d ", f);

		// Compute "alpha": projecting "refmode" vector into "cevec" modal space to obtain "alpha" projection
		for(int n=0; n<nevec; n++)
			alpha[n] = dotprodnorm( refmode, cevec + n*ncomps, ncomps, modref );
		// show_vector(stdout, alpha, nevec, "", " %7.4f", false, false);

		// Compute "delta" array from "alpha" ("delta" array is the element-wise square of "alpha")
		pow_vector(alpha, delta, nevec, 2);
//		if(verb > 0)
//			show_vector(stdout, delta, nevec, "", " %6.4f", false, false);

		// Compute "delta" value
//		if(verb > 0)
//			fprintf(stdout, "  %7.5f", sqrt( sum_vector(delta, nevec) ));


		// Using Alpha
		for(int k = 0; k < size; k++)
			mode[ k ] = 0.0; // reset
		for(int n = 0; n < nevec; n++)
			for(int k = 0; k < size; k++) // Mon added the +3, watch out!
				mode[ k ] += eigvect[ n * size + k ] * alpha[n]; // Amplitude of n-th mode


		// Generate trajectory (dihedral motion)
		if(maxang != 0.0)
			scale_vectors(mode, size, 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode

		move_loop_dihedral(itermol, ifr, ilr, props, mode, size, model, 1.0);

		// Does some loop atom clash with its environment?
		//		sprintf(dummy, " %d", clashed_loop( itermol, itermol2, ifa, na, 1.0) );
		//		strcat(text, dummy);

		if(iterini != NULL)
		{
			anchor_drift(iterini, itermol, props, ilr, &adist, &aang);
			ddrift = fabs(adist - adist0); // Distance increment wrt. initial distance
			adrift = fabs(aang - aang0); // Angle increment wrt. initial angle

			if (ddrift > 0.03 || adrift> 9)
			{
				    if(verb > 0)
					fprintf(stdout,"follow_mode> Motion tends to overfit = %8f  %8f rmsd %8f \n", ddrift, adrift, rmsd);
			        break;
			}

			rmsdP = rmsd;
			rmsd = rmsd_loop(iterini, itermol, ifa, na);


			//fprintf(stdout, " %7.4f %7.4f %7.4f %7.4f\n", rmsd, rmsdL, rmsdP, (rmsdL - rmsdP) );



			if(verb > 0)
				fprintf(stdout, " %7.4f\n", rmsd);
			//			sprintf(dummy, " %7.4f", rmsd);
			//			strcat(text, dummy);

			// Convergence check
			if(rmsd > rmsd_conv)
			{
				if(verb > 0) fprintf(stdout,"follow_mode> Motion convergence reached! dRMSD = %8f > %8f\n", rmsd, rmsd_conv);
				mol->writeMloop(file_movie, (*p_fi)++, ifr-1, ilr+1, chain);
				if ((*p_fi)>=max_loops_save && max_loops_save!=0) break;
				// exit(0);
				break;
			}

			// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
			if(delta_rmsd < rmsd - last_rmsd)
			{
				if(verb > 0)
					fprintf(stdout,"%s> Structure dumped into Muli-PDB dRMSD = %8f > %8f\n", prog, rmsd - last_rmsd, delta_rmsd);
				//if (ddrift < 0.05 and adrift < 9)
				{
				mol->writeMloop(file_movie, (*p_fi)++, ifr-1, ilr+1, chain);
				last_rmsd = rmsd; // Keep last saved RMSD
				if ((*p_fi)>=max_loops_save && max_loops_save!=0) break;

				}
			}



			//			// Compute anchor drift
			//			anchor_drift(iterini, itermol, props, ilr, &adist, &aang);
			//			ddrift = adist - adist0; // Distance increment wrt. initial distance
			//			adrift = aang - aang0; // Distance increment wrt. initial distance
			//			sprintf(dummy, " %6.3f %5.2f", ddrift, adrift);
			//			strcat(text, dummy);

			// fprintf(stdout,"%s\n", text); // Dump all output for current frame
			// rmsd_old = rmsd;
		}

		// Convergence test
		//				if(!mr_switch && (rmsd - rmsd_old < rmsd_conv || clashed_loop( itermol, itermol2, ifa, na, 1.0)))


		//		fprintf(stdout,"%s\n", text); // Dump all output for current frame
		//		fprintf(stdout, "\n");

		//	Update the "refmode" and "modref" with the current most overlapping mode
		if(update)
		{
			//			int nmax = get_max_index(delta, nevec);

			for(int i=0; i<ncomps; i++)
				refmode[i] = 0.0; // reset "refmode"

			// Update the reference mode "refmode"
			for(int k=0; k<nevec; k++)
			{
				//				if(alpha[k] > 0) // same direction
				//					for(int i=0; i<ncomps; i++)
				//						refmode[i] += cevec[ k * ncomps + i ] * delta[k];
				//				else // reversed direction
				//					for(int i=0; i<ncomps; i++)
				//						refmode[i] -= cevec[ k * ncomps + i ] * delta[k];

				for(int i=0; i<ncomps; i++)
					refmode[i] += cevec[ k * ncomps + i ] * alpha[k];

			}

			// Update the modulus of reference mode "refmode"
			modref = vector_modulus( refmode, ncomps ); // Reference CC vector modulus
		}

	}

	free(coord); // free obsolete raw coordinates
	free(mode); // free mode
	free(alpha);
	free(delta);
	free(mass_matrix);
	free(der); // Free obsolete derivatives
	free(hess_matrix);
	free(rdr);
	free(cevec); // free obsolete modes
	free(decint); // free obsolete elastic network
	free(work);
	free(vr);
	free(alphar);
	free(alphai);
	free(beta);
	free(vl);

	if(update)
		free(refmode);

	delete itermol;
	delete itermol2;

	return neig; // Return the number of non-null eigenpairs
}

inline int MC_mode(double *refmode0, Macromolecule *mol, int model, int type, tri *props, float *masses, float *masses_loop, int ifa, int ifr, int ilr,
		int na, int size, int nco, double maxang, float cutoff, int nsteps, double *eigval, double *eigvect, float rmsd_conv, pdbIter *iterini,
		char *file_movie, float delta_rmsd, int *p_fi, char chain, bool update)
{
	bool debug = false; // dump debug info
	float *coord;
	int nipa;
	int neig; // Number of non-null eigenvectors
	int ncomps = 3*(na+3); // Number of Cartesian components (x,y,z for each atom of mobile loop + 3 Ct-anchor atoms)

	double *alpha = (double *) malloc( (size+nco) * sizeof(double)); // Alpha overlaps array (allocate memory for the maximum possible)
	double *delta = (double *) malloc( (size+nco) * sizeof(double)); // Delta overlaps array (allocate memory for the maximum possible)
	double *mode = (double *) malloc( size * sizeof(double)); // Merged mode for motion

	double *refmode;

	//	Update the "refmode" and "modref" with the current most overlapping mode
	if(update)
	{
		refmode = (double *) malloc( ncomps * sizeof(double)); // Allocate reference CC vector

		for(int i=0; i<ncomps; i++)
			refmode[i] = refmode0[i]; // Copy the reference mode "refmode0" into "refmode"
	}
	else
		refmode = refmode0; // Use the initial "refmode0"

	// Store the modulus of reference mode "refmode"
	double modref = vector_modulus( refmode, ncomps ); // Reference CC vector modulus


	// Create iterators for RMSD computations
	pdbIter *itermol = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule
	pdbIter *itermol2 = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule to speed up Elastic Network refresh
	int num_atoms = itermol->num_atom();

	if(eigval == NULL) // if not previously allocated
		eigval = (double *) malloc( sizeof(double) * (size + nco) );

	if(eigvect == NULL) // if not previously allocated
		eigvect = (double *) malloc( sizeof(double) * (size + nco) * size );

	float rmsd; // Current RMSD
	float last_rmsd = 0.0; // Last saved RMSD
	mol->coordMatrix( &coord );

	// Allocate derivatives
	trd *der = (trd *) malloc( (na + 3) * size * sizeof(trd) );
	ptr_check(der);

	int sizex = size + nco;
	double *mass_matrix; // kinetic energy matrix
	mass_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(mass_matrix);

	double *rdr = (double *) malloc(size * sizeof(double)); // Auxiliar array
	ptr_check(rdr);

	double *hess_matrix; // Bordered Hessian matrix
	hess_matrix = (double *) malloc(sizex*sizex * sizeof(double));
	ptr_check(hess_matrix);

	double *cevec  = (double *) malloc( 3 * (na+3) * (size - nco) * sizeof(double));

	//decint = ( twid * ) realloc( decint, (ifpa + nla)*12 * sizeof( twid ) ); // resizes contact list-structure

	twid *decint = ( twid * ) malloc(   (na)*500 * sizeof( twid ) ); //  60 contacts per atom
	ptr_check(decint);

	double *alphar, *alphai, *beta, *vr,  *vl, *work;


	alphar = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphar);

	alphai = (double *) malloc(sizex * sizeof(double));
	ptr_check(alphai);

	beta = (double *) malloc(sizex * sizeof(double));
	ptr_check(beta);

	vl = (double *) malloc(1*sizex * sizeof(double));
	ptr_check(vl);

	vr = (double *) malloc(sizex*sizex * sizeof(double));   /* eigenvectors */
	ptr_check(vr);

	work = (double *) malloc(20*sizex * sizeof(double));
	ptr_check(work);

	for(int f = 0; f < nsteps; f++) // some maximum number of sampling steps
	{
		// Initialize Eigenvalues
		for(int i=0; i<size+nco; i++)
			eigval[i] = 0.0;

		// Initialize Eigenvectors
		for(int i=0; i<(size+nco)*size; i++)
			eigvect[i] = 0.0;

		if(debug)
			fprintf(stderr, "follow_mode> Getting coordinates single row (pseudo-atom model)\n");

		update_loop_coords(itermol, props[ifr].k1, na, coord);

		if(debug)
			fprintf(stderr, "follow_mode> Computing derivatives...\n");
		//der = drdqC5x(coord, props, ifr, ilr, na, size, model);
		drdqC5x(coord, props, ifr, ilr, na, size, model, der);


		if(debug)
			fprintf(stderr, "follow_mode> Computing Kinetic Energy matrix (masses matrix)...\n");
		//mass_matrix = kineticC5x(der, masses, props, ifr, na, size, nco);
		kineticC5x(der, masses, props, ifr, na, size, nco, mass_matrix );


		//		if(decint != NULL)
		//			free(decint); // Free obsolete contacts list

		// Allocate "decint" to store the contacts list (ipas-list)
		//if( !(decint = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )  // Required for "realloc"
		//{
		//	fprintf(stderr, "Sorry, \"decint\" memory allocation failed!\n");
		//	exit(1);
		//}

		if(debug)
			fprintf(stderr, "follow_mode> Computing Interacting Pair of (non-virtual) Atoms (ipas)\n");

		// INVERSE EXPONENTIAL (power of distance for contact matrix)
		// Making Interacting Pairs of (non-virtual) Atoms (ipas)
		// make_ipas_loop(itermol, itermol2, ifa, na, cutoff, &decint, &nipa); // Updating Elastic network
		make_ipas_loop(itermol, itermol2, ifa, na, cutoff, decint, &nipa); // Updating Elastic network

		// fprintf(stderr, "nipa %d\n", nipa);


		if(debug)
			fprintf(stderr, "follow_mode> Inverse Exponential (%d nipas) cutoff= %.1f, k= %f, x0= %.1f ", nipa, cutoff_k0, cte_k0, x0);

		for(int i=0; i<nipa; i++)
			decint[i].C = Inv_exp( cte_k0, decint[i].d, x0, power); // setting Force Constants

		// IPAs checking
		if(debug) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
			for(int i=0; i<nipa; i++)
				fprintf(stderr, "ipa %4d: k= %d  l= %d  d= %f  C= %f\n",i,decint[i].k,decint[i].l,decint[i].d,decint[i].C);

		if(debug)
			fprintf(stderr, "follow_mode> Computing Hessian matrix (potential energy matrix)...\n");
		hessianC5x(coord, der, props, decint, nipa, ifr, na, size, nco, rdr, hess_matrix);
		//hessianFast(coord, coord, der, props, decint, nipa, ifr, ilr, size, nco, num_atoms, hess_matrix);

		// free(decint); // free obsolete elastic network

		if(debug)
		{
			// Show Hessian matrix
			show_matrix(hess_matrix, size + nco, "Hessian:", " %7.2f");

			// Show Kinetic Energy matrix
			show_matrix(mass_matrix, size + nco, "Kinetic:", " %7.0f");
		}


		// COMPUTING THE EIGENVECTORS AND EIGENVALUES
		//int info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig);
		int info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig, alphar, alphai, beta, vr,  vl, work);

		//free(mass_matrix);
		//free(hess_matrix);

		// Some checking...
		if( info ) // if info != 0
		{
			fprintf(stderr, "\nfollow_mode> An error occured in the matrix diagonalization: %d\n", info);
			exit(1);
		}

		// MON: check this, seems unnecessary...
		if(neig < nevec)
		{
			fprintf(stderr," Warning more eigenvectors requested (%d) than available (%d), forcing exit!\n",nevec,neig);
		}

		if(debug)
		{
			fprintf(stderr, "follow_mode> Eigensolver successfully finished!!! (neig=%d)\n", neig);
			show_vector(stderr,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
			show_vectors(stderr,eigvect,size,neig,"Dumping Raw Eigenvectors:", " %5.2e");
		}

		// Scale eigenvectors so that the maximum value of the components is "maxang"
		if(maxang != 0.0)
		{
			scale_vectors(eigvect,size,neig,maxang * M_PI / 180.0);
			if(verb > 1)
				show_vectors(stderr,eigvect,size,neig,"Dumping Scaled Eigenvectors:", " %5.2e");
		}

		// Compute the Cartesian eigenvectors from the Internal Coordinates eigenvectors
		// double *cevec; // Cartesian eigenvectors
		// cevec = ic2cart(eigvect, neig, der, size, na + 3, masses_loop);
		cevec = ic2cart(eigvect, neig, der, size, na + 3, masses_loop);


		//		sprintf(text,"follow_mode> %4d ", f);
		if(verb > 0)
			fprintf(stdout,"follow_mode> %4d ", f);

		// Compute "alpha": projecting "refmode" vector into "cevec" modal space to obtain "alpha" projection
		for(int n=0; n<nevec; n++)
			alpha[n] = dotprodnorm( refmode, cevec + n*ncomps, ncomps, modref );
		// show_vector(stdout, alpha, nevec, "", " %7.4f", false, false);

		// Compute "delta" array from "alpha" ("delta" array is the element-wise square of "alpha")
		pow_vector(alpha, delta, nevec, 2);
		if(verb > 0)
			show_vector(stdout, delta, nevec, "", " %6.4f", false, false);

		// Compute "delta" value
		if(verb > 0)
			fprintf(stdout, "  %7.5f", sqrt( sum_vector(delta, nevec) ));


		// Using Alpha
		for(int k = 0; k < size; k++)
			mode[ k ] = 0.0; // reset
		for(int n = 0; n < nevec; n++)
			for(int k = 0; k < size; k++) // Mon added the +3, watch out!
				mode[ k ] += eigvect[ n * size + k ] * alpha[n]; // Amplitude of n-th mode


		// Generate trajectory (dihedral motion)
		if(maxang != 0.0)
			scale_vectors(mode, size, 1, fabs(maxang) * M_PI / 180.0); // Scale current merged mode

		move_loop_dihedral(itermol, ifr, ilr, props, mode, size, model, 1.0);

		// Does some loop atom clash with its environment?
		//		sprintf(dummy, " %d", clashed_loop( itermol, itermol2, ifa, na, 1.0) );
		//		strcat(text, dummy);

		if(iterini != NULL)
		{
			rmsd = rmsd_loop(iterini, itermol, ifa, na);
			if(verb > 0)
				fprintf(stdout, " %7.4f\n", rmsd);
			//			sprintf(dummy, " %7.4f", rmsd);
			//			strcat(text, dummy);

			// Convergence check
			if(rmsd > rmsd_conv)
			{
				if(verb > 0) fprintf(stdout,"follow_mode> Motion convergence reached! dRMSD = %8f > %8f\n", rmsd, rmsd_conv);
				mol->writeMloop(file_movie, (*p_fi)++, ifr-1, ilr+1, chain);
				// exit(0);
				break;
			}

			// Write just the indicated loop (from the index of first residue "ifr" to index of last residue "lfr") into a Multi-PDB
			if(delta_rmsd < rmsd - last_rmsd)
			{
				if(verb > 0)
					fprintf(stdout,"%s> Structure dumped into Muli-PDB dRMSD = %8f > %8f\n", prog, rmsd - last_rmsd, delta_rmsd);
				mol->writeMloop(file_movie, (*p_fi)++, ifr-1, ilr+1, chain);
				last_rmsd = rmsd; // Keep last saved RMSD
			}

			//			// Compute anchor drift
			//			anchor_drift(iterini, itermol, props, ilr, &adist, &aang);
			//			ddrift = adist - adist0; // Distance increment wrt. initial distance
			//			adrift = aang - aang0; // Distance increment wrt. initial distance
			//			sprintf(dummy, " %6.3f %5.2f", ddrift, adrift);
			//			strcat(text, dummy);

			// fprintf(stdout,"%s\n", text); // Dump all output for current frame
			// rmsd_old = rmsd;
		}

		// Convergence test
		//				if(!mr_switch && (rmsd - rmsd_old < rmsd_conv || clashed_loop( itermol, itermol2, ifa, na, 1.0)))


		//		fprintf(stdout,"%s\n", text); // Dump all output for current frame
		//		fprintf(stdout, "\n");

		//	Update the "refmode" and "modref" with the current most overlapping mode
		if(update)
		{
			//			int nmax = get_max_index(delta, nevec);

			for(int i=0; i<ncomps; i++)
				refmode[i] = 0.0; // reset "refmode"

			// Update the reference mode "refmode"
			for(int k=0; k<nevec; k++)
			{
				//				if(alpha[k] > 0) // same direction
				//					for(int i=0; i<ncomps; i++)
				//						refmode[i] += cevec[ k * ncomps + i ] * delta[k];
				//				else // reversed direction
				//					for(int i=0; i<ncomps; i++)
				//						refmode[i] -= cevec[ k * ncomps + i ] * delta[k];

				for(int i=0; i<ncomps; i++)
					refmode[i] += cevec[ k * ncomps + i ] * alpha[k];

			}

			// Update the modulus of reference mode "refmode"
			modref = vector_modulus( refmode, ncomps ); // Reference CC vector modulus
		}

	}

	free(coord); // free obsolete raw coordinates
	free(mode); // free mode
	free(alpha);
	free(delta);
	free(mass_matrix);
	free(der); // Free obsolete derivatives
	free(hess_matrix);
	free(rdr);
	free(cevec); // free obsolete modes
	free(decint); // free obsolete elastic network
	free(work);
	free(vr);
	free(alphar);
	free(alphai);
	free(beta);
	free(vl);

	if(update)
		free(refmode);

	delete itermol;
	delete itermol2;

	return neig; // Return the number of non-null eigenpairs
}
// Loop NMA routine to just compute the eigenvectors/values given some macromolecular loop
//	mol --> Input macromolecule
//	model --> Atomic model
//	type --> ICs type (with/without chi)
//	props --> Properties structure
//	masses --> Atomic masses array
//	ifa --> Index of first atom of the loop
//	ifr --> Index of first residue of the loop
//	ilr --> Index of the last residue of the loop
//	na --> Number of atoms of the loop
//	size --> Number of ICs of the loop
//	nco --> Number of constraints (typically 6)
//	cutoff --> Distance cutoff for elastic network definition
//	eigval --> IC eigenvalues, automatic memory allocation if NULL (OUTPUT)
//	eigvect --> IC eigenvectors, automatic memory allocation if NULL (OUTPUT)
//	der --> Derivatives (dr/dq), always allocated here (free elsewhere) (OUTPUT)
//	RETURN --> Number of non-null eigenpairs
int nma_loop(Macromolecule *mol, int model, int type, tri *props, float *masses, int ifa, int ifr, int ilr,
		int na, int size, int nco, float cutoff, double *eigval, double *eigvect, trd **p_der)
{
	bool debug = false; // dump debug info
	double *mass_matrix; // Kinetic energy matrix
	double *hess_matrix; // Hessian matrix
	twid *decint = NULL; // Contacts list structure
	int nipa;

	// Create iterators for RMSD computations
	pdbIter *itermol = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule
	pdbIter *itermol2 = new pdbIter( mol, true, true, true, true ); // Iterator to current (moving) macromolecule to speed up Elastic Network refresh

	if(eigval == NULL) // if not previously allocated
		eigval = (double *) malloc( sizeof(double) * (size + nco) );

	if(eigvect == NULL) // if not previously allocated
		eigvect = (double *) malloc( sizeof(double) * (size + nco) * size );

	// Initialize Eigenvalues
	for(int i=0; i<size+nco; i++)
		eigval[i] = 0.0;

	// Initialize Eigenvectors
	for(int i=0; i<(size+nco)*size; i++)
		eigvect[i] = 0.0;

	float *coord; // (pseudo)atomic coordinates single row vector
	if(debug)
		fprintf(stdout, "%s> Getting coordinates single row (pseudo-atom model)\n", prog);
	mol->coordMatrix( &coord );

	trd *der; // Derivatives
	if(debug) {
		fprintf(stdout, "%s> Computing derivatives...\n", prog);
	}
	der = drdqC5x(coord, props, ifr, ilr, na, size, model); // Derivatives

	if(debug)
		fprintf(stdout, "%s> Computing Kinetic Energy matrix (masses matrix)...\n", prog);
	mass_matrix = kineticC5x(der, masses, props, ifr, na, size, nco);

	// Allocate "decint" to store the contacts list (ipas-list)
	if( !(decint = ( twid * ) malloc( 1 * sizeof( twid ) ) ) )  // Required for "realloc"
	{
		fprintf(stdout,"Sorry, \"decint\" memory allocation failed!\n");
		exit(1);
	}

	// INVERSE EXPONENTIAL (power of distance for contact matrix)
	// Making Interacting Pair of (non-virtual) Atoms (ipas)
	make_ipas_loop(itermol, itermol2, ifa, na, cutoff, &decint, &nipa); // Updating Elastic network
	// fprintf( stdout, "%s> Inverse Exponential (%d nipas) cutoff= %.1f, k= %f, x0= %.1f ", prog, nipa, cutoff_k0, cte_k0, x0);
	for(int i=0; i<nipa; i++)
		decint[i].C = Inv_exp( cte_k0, decint[i].d, x0, power); // setting Force Constants

	// IPAs checking
	if(debug) // If Hessian and Kinetic energy matrices calculation and diagonalization are enabled.
		for(int i=0; i<nipa; i++)
			fprintf(stdout,"ipa %4d: k= %d  l= %d  d= %f  C= %f\n",i,decint[i].k,decint[i].l,decint[i].d,decint[i].C);

	if(debug)
		fprintf(stdout, "%s> Computing Hessian matrix (potential energy matrix)...\n", prog);
	hess_matrix = hessianC5x(coord, der, props, decint, nipa, ifr, na, size, nco);
	free(coord);
	free(decint);
	//	free(der); // Free obsolete derivatives

	// Show Hessian matrix
	if(debug)
		show_matrix(hess_matrix, size + nco, "Hessian:", " %7.2f");

	// Show Kinetic Energy matrix
	if(debug)
		show_matrix(mass_matrix, size + nco, "Kinetic:", " %7.0f");

	// COMPUTING THE EIGENVECTORS AND EIGENVALUES
	int neig = 0; // Number of non-null eigenvectors
	int info = diag_dggev(eigval, eigvect, mass_matrix, hess_matrix, size, nco, &neig);
	free(mass_matrix);
	free(hess_matrix);
	// show_vectors(stdout,eigvect,size,neig,"Dumping Raw Eigenvectors:", " %5.2e");

	//*  DSYGVX computes SELECTED eigenvalues, and optionally, eigenvectors
	//*  of a real generalized symmetric-definite eigenproblem, of the form
	//*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
	//*  and B are assumed to be symmetric and B is also positive definite.
	//*  Eigenvalues and eigenvectors can be selected by specifying either a
	//*  range of values or a range of indices for the desired eigenvalues.
	//	int info = 0;
	//	fprintf(stderr,"size= %d\n",size);
	//	diag_dsygvx(hess_matrix, mass_matrix, eigval, size, nco, 1, 17);
	//	show_vectors(stdout,hess_matrix,size,17,"Dumping Raw Eigenvectors:", " %5.2e");
	//exit(0);



	// Some checking...
	if( info ) // if info != 0
	{
		fprintf(stdout,"\nnma_loop> An error occured in the matrix diagonalization: %d\n", info);
		exit(1);
	}

	if(debug)
	{
		fprintf( stdout, "%s> Eigensolver successfully finished!!! (neig=%d)\n", prog, neig);
		show_vector(stdout,eigval,neig,"Dumping Raw Eigenvalues:", " %5.2e");
		show_vectors(stdout,eigvect,size,neig,"Dumping Raw Eigenvectors:", " %5.2e");
	}

	*p_der = der; // Return derivatives
	return neig; // Return the number of non-null eigenpairs
}

// Get loop coordinates into a pre-allocated "coord" array (intended for copy & paste or coordinates backup)
void get_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord)
{
	Tcoor pos;
	int i = 0;

	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		// Get atom position
		( iter->get_atom() )->getPosition(pos);

		// Apply position increment
		coord[i*3]     = pos[0];
		coord[i*3 + 1] = pos[1];
		coord[i*3 + 2] = pos[2];
		//fprintf(stdout,"%f %f %f\n",pos[0],pos[1],pos[2]);

		i++; // loop atom local index
	}
}




inline void update_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord)
{
	Tcoor pos;
	int i = 0;

	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		// Get atom position
		( iter->get_atom() )->getPosition(pos);
		int index = iter->pos_atom*3;
		//fprintf(stdout,"%s> N %d  %f %f %f \n", prog, index,  coord[index], coord[index + 1], coord[index + 2]);
		//fprintf(stdout,"%s> N %d  %f %f %f \n", prog, index,  pos[0], pos[1], pos[2]);
		coord[index]     = pos[0];
		coord[index + 1] = pos[1];
		coord[index + 2] = pos[2];
	}
}





// Set loop coordinates from a pre-allocated "coord" array (intended for copy & paste or coordinates backup)
void set_loop_coords(pdbIter *iter, int ifpa, int nla, float *coord)
{
	Tcoor pos;
	int i = 0;

	for( iter->pos_atom = ifpa; iter->pos_atom < ifpa + nla; iter->next_atom() ) // just screen mobile loop atoms
	{
		// Apply position increment
		pos[0] = coord[i*3];
		pos[1] = coord[i*3 + 1];
		pos[2] = coord[i*3 + 2];

		// Set atom position
		( iter->get_atom() )->setPosition(pos);

		i++; // loop atom local index
	}
}

// Reads all lines from a text file and returns the number of lines read (automatic memory allocation)
//  file       --> File name
//  p_lines    --> Pointer to the lines array (it will be allocated automatically)
//  linelength --> Length of each line (number of characters)
int readTextLines(char *file, char ***p_lines, int linelength) // Reading all rows from textfile
{
	FILE *f;
	char myline[1024]; // For reading text files
	char **lines; // lines list
	int nlines = 0; // Number of lines
	int nscan; // Number of successful format conversions from sscanf
	lines = (char **) malloc( sizeof(char *) );

	if( !(f = fopen(file,"r")) )
	{
		fprintf(stderr,"readTextLines> ERROR opening %s. Forcing exit!\n",file);
		exit(2);
	}

	while( fgets(myline,1024,f) )
	{
		//fprintf(stdout,"%s\n",myline);
		if(myline[0] != '#')
		{
			// Memory allocations...
			if( !( lines = (char **) realloc( lines, sizeof(char *) * (nlines+1) ) ) )
			{
				fprintf(stderr,"readTextLines> Memory for text list can't be allocated! Forcing exit!\n\n");
				exit(1);
			}

			lines[nlines] = (char *) malloc( sizeof(char) * linelength );

			nscan = sscanf(myline, "%s", lines[nlines]); // Reading line

			// fprintf(stderr,"line %3d: %s\n",nlines+1,lines[nlines]);

			if(nscan != 1)
			{
				fprintf(stderr,"readTextLines> %d successful format conversions from sscanf()! Forcing exit!\n\n",nscan);
				exit(1);
			}

			nlines++;
		}
	}
	fclose(f);

	*p_lines = lines; // output lines
	return nlines;
}

// Compute IC "mode" from "eigvect" IC eigenvectors and "alpha" and "delta" arrays. Requires:
//	"nevec"  --> Number of eigenvectors
//	"size"   --> Number of dihedral coordinates
void make_mode(double *alpha, double *delta, double *eigvect, double *mode, int nevec, int size)
{
	// Initialize current mode
	for(int k = 0; k < size; k++) // Mon added the +3, watch out!
		mode[ k ] = 0.0; // zero initialization

	for(int n=0; n<nevec; n++)
	{
		// Generate "merged" Dihedral-coordinates mode for motion
		if(alpha[n] > 0.0) // Reverse n-th mode
			for(int k = 0; k < size; k++) // Mon added the +3, watch out!
				mode[ k ] += eigvect[ n * size + k ] * delta[n]; // Amplitude of n-th mode
		else
			for(int k = 0; k < size; k++) // Mon added the +3, watch out!
				mode[ k ] -= eigvect[ n * size + k ] * delta[n]; // Reverse current mode to maintain initial direction
	}
}

// Compute element-wise "p" power of "in" vector into "out" vector, both of same length "len"
void pow_vector(double *in, double *out, int len, int p)
{
	for(int n=0; n<len; n++)
		out[n] = pow(in[n], p);
}

// Adds the "len" elements in vector "in"
double sum_vector(double *in, int len)
{
	double sum = 0.0;
	for(int n=0; n<len; n++)
		sum += in[n];
	return sum;
}

// Get the index of the maximum value in the array
int get_max_index(double *v, int n)
{
	double max = -1e99;
	int imax = 0; // index of maximum value

	for(int i=0; i<n; i++)
		if(v[i] > max)
		{
			max = v[i];
			imax = i;
		}

	return imax;
}

