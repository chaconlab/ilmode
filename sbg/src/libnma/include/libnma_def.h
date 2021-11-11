#ifndef LIBNMA_DEF_H_
#define LIBNMA_DEF_H_

#include "nma.h"

// Just for back-compatibility with "nmac"
void save_ascii_defmob (int , double *, double *,  double *, char *, char *filename=NULL);
void save_ascii_defmob (int num_atoms, double *deform, double *mob,  double *bf);
void compute_def (int , float *, double *, myfloat *, myfloat *, int, int, double **, double **);
// Mon modified (25/3/2010)
void compute_def(int num_atoms, float *pdb_model, double *dist_matrix,
                  floating *eigval, floating *eigvect, int neigval, int nf1,
                  double **defp, double  **mobp, double dc=15);
void norm_defmob (int , double *);

// They were "extern" before...
void ludcmp(double **, int , int *, double *, double * );
// void hardy(int, double *, double **, double **, int *,double *, double dc=0.001);
void hardy( int num_points, double * dist_matrix, double ** a, double ** b, int * indx, double *vv, double dc=15 );
void grad_hmq(float *, int, double **, double **, int *, double *,
              double *, double *, double *, double * );
double lnorm(double, double, double, double, double, double);
void djacobi3( double [4][4], int , double [4], double [4][4], int * );
void lubksb(double **, int, int *, double []);
void compute_def (int , float *, double *, float *, float *, int, int, double **, double **);
void save_ascii_defmob (int, double *, double *);
void norm_defmob (int , double *);

#endif /*LIBNMA_DEF_H_*/
