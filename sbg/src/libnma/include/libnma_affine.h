#ifndef LIBNMA_AFFINE_H_
#define LIBNMA_AFFINE_H_

#include "nma.h"

extern "C"
{
	int dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);
	void dpotrs_( char *uplo, int *n, int *nrhs, double *A, int *lda, double *B, int *ldb, int *info );
}

// Normalizes a vector (1 length): c = norm(c)
void normalize(float *v,float *c);
// Returns the vector length (magnitude)
float magnitude(float *v,int dim = 3);
// Returns the vector length (magnitude)
double magnitude(double *v,int dim = 3);
// Returns the squared vector length (magnitude)
float sqrMagnitude(float *v);
// Computes the cross-product: a x b = c
void crossProduct(float *a, float *b, float *c);
// Computes the dot-product: a · b = c
float dotProduct(float *a, float *b);
// Multiply vector by scalar: a *= s
void multByScalar(float *a, float s);
// Multiply vector by scalar: a * s = c
void multByScalar(float *a, float s, float *c);
// Subtract two vectors: a - b = c
void subtract(float *a, float *b, float *c);
// Add two vectors: a - b = c
void add(float *a, float *b, float *c);
// Builds the Rotational Correlation Matrix. See:
// "Automated Illustration of Molecular Flexibility"
// Aaron Bryden, George Phillips Jr., Michael Gleicher
// IEEE Transactions on Visualization and Computer Graphics, Volume 18, Number 1, page 132--145 - Jan 2012
void buildRotationalCorrelationMatrix(Macromolecule *mol, double *vector, float **p_corrMatrix);
// Builds the Dot-Product-based Correlation Matrix.
void buildDotProductCorrelationMatrix(Macromolecule *mol, double *vector, float **p_corrMatrix);
// Solves linear system
void solveCholesky( double *A, int colsA, double *B, int colsB, double *X );
// Assumes bottom row of matrix "m" is [0,0,0,1]
void multMatrix_4x4_x_3(double *m,double *vec, double *ret);
// Assumes bottom row of matrix "m" is [0,0,0,1] (single precision)
void multMatrix_4x4_x_3f(double *m,float *vec, float *ret);
// Generic version... using Cramer's rule
// Note: this assumes the matrix is invertible (det != 0)
void invertMatrix4x4(double *m);
// Transposes a 4x4 matrix "m" into "mat"
void transposeMatrix4x4(double *m,double *mat);
// Multiply two 4x4 matrices ("A" and "B") into "C"
void multMatrix_4x4_x_4x4(double *A, double *B, double *C);

// This function is intended to compute the Affine transformation: v = M·x + t
// (it uses homogeneous coordinates, i.e. a 4 elements vector to define a 3D position)
// coords --> array of macromolecular coordinates
// natoms --> number of atoms
// vector --> pointer to the selected eigenvector
// cluster --> pointer to the array of atom indices belonging to the cluster of atoms
// nclus --> number of cluster atoms
// p_expMat --> pointer to the OUTPUT Affine transformation matrix (=NULL --> automatic memory allocation)
// weights --> array with weighting factors (OPTIONAL)
void getExpMat(float *coords, int natoms, double *vector, int *cluster, int nclus, double **p_expMat, float *weights=NULL);

// Affine model integration.
// A --> Input affine model matrix (LAPACK's standard format)
// X --> Output (integrated) affine model matrix (LAPACK's standard format)
// scale --> Integration range
// iters --> Number of interations (Alexa's method?)
void getExp(double *Af, double *X, double scale, int iters=7);
// Shows a standard format matrix (standard output)
void show_matrix_standard(double *matrix, int size, char *name);
// This calculates the Error (formula XX in Bryden's paper).
// The Error is the distance between the eigenvector component and the affine model evaluation.
// pos --> atom position
// vec --> eigenvector component for such atom
// expMat --> Affine transformation matrix (4x4)
double getError(float *pos, double *vec, double *expMat);
// Projects one point into a plane defined by a point and a normal. Output in "pos2"
void projectPointToPlane(float *pos, float *com, float *normAxis, float *pos2);
// Algorithm to generate a given number of spiral points uniformly distributed on the surface of a sphere.
// double r;		// true radius
// int n;		    // number of points
// float **p_points; // Pointer to the points-array (for dynamic memory allocation)
void spherePoints(double r, int n, float **p_points);
// Compute the squared distance between two points
float sqrDist(float *p1, float *p2);
// Compute the distance between two points
float Dist(float *p1, float *p2);
// Builds a boolean adjacency matrix at atom level.
// mol --> The macromolecule
// p_adj --> Pointer to the matrix (packed storage). Set *adj to NULL for automatic memory allocation.
// cutoff --> Proximity distance cutoff
void buildAdjacencyMatrix(Macromolecule *mol, bool **p_adj, float cutoff);
// Compute adjacency matrix at cluster level. (It also can update, but inefficiently...)
// coords --> array of macromolecular coordinates
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// p_cadj --> pointer to the cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
//void computeClusterAdjMatrix(int **clusters, int nclusters, int *natoms, bool *adj, bool **p_cadj);
// adj --> atomic adjacency matrix
void computeClusterAdjMatrix(float *coords, int **clusters, int nclusters, int *natoms, bool **p_cadj, float cutoff);

// Performs Hierarchical Clustering interleaved with Affine model calculation.
// mol --> Macromolecule (to screen at residue level)
// coords --> array of atomic coordinates
// requested --> Number of clusters requested (in "rclusters")
// evecx --> Eigenvector pointer
// cutoff --> Adjacency cutoff
// p_rclusters --> Pointer to the requested clusters array (OUTPUT)
// p_rnatoms --> Pointer to the requested clusters number of atoms per cluster array (OUTPUT)
// p_rerrors --> Pointer to the requested clusters error array (OUTPUT)
// OPTIONAL INPUT:
//   rigid_bodies --> Number-of-residues sized array with the indices of the rigid bodies detected in a "fixfile"
//   nrigid_bodies --> Number of such rigid bodies (it will be the initial number of clusters)
//   extend_size --> Number of atoms threshold to consider cluster extension (clusters with less atoms will be extended). (Default=5)
//   extend_cutoff --> Distance (in Angstroms) to extend small clusters. (Default=10A)
void doHiearchicalExpMatClustering(Macromolecule *mol, float *coords, int requested, double *evecx, float cutoff, int ****p_rclusters, int ***p_rnatoms, float ***p_rerrors, int *rigid_bodies=NULL, int nrigid_bodies=0, int extend_size=5, float extend_cutoff=10);

// Shows an Adjacency Matrix (boolean)
void showAdjMatrix(bool *adj, int n, char *name);
// Shows all clusters and their atoms
void showClusters(int **clusters, int nclusters, int *natoms, float *errors);
void computeClusterMinDistMatrix(float *coords, int **clusters, int nclusters, int *natoms, float **p_cadj);
// Shows an Adjacency Matrix (boolean)
void showDistMatrix(float *adj, int n, char *name);
// Update adjacency matrix at cluster level upon merging.
// coords --> array of macromolecular coordinates
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// cadj --> cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
// i,j --> merged cluster indices (i<j)
// NOTE: "clusters" and "natoms" should be already updated! (but "nclusters" should not!)
void updateClusterAdjMatrix(float *coords, int **clusters, int nclusters, int *natoms, bool *cadj, float cutoff, int x, int y);
void updateClusterMinDistMatrix(Macromolecule *mol,int **clusters, int nclusters, int *natoms, float *cadj, float cutoff, int x, int y);
// Update adjacency matrix at cluster level upon merging.
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// cadj --> cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
// i,j --> merged cluster indices (i<j)
// NOTE: "clusters", "nclusters" and "natoms" should NOT have been updated yet!
void updateClusterAdjMatrix2(int **clusters, int nclusters, int *natoms, bool *cadj, float cutoff, int x, int y);

#endif /*LIBNMA_AFFINE_H_*/
