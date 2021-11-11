/************************************************************************
*                     LIBRARY: libnma_affine                            *
*************************************************************************
* This library is part of the iMOD package, URL: http://chaconlab.org   *
* (c) Jose Ramon Lopez-Blanco (Mon) and Pablo Chacon.                   *
* Rocasolano's Structural Bioinformatics Group (2011-2013)              *
*************************************************************************
*                                                                       *
*   Library for the Automated Illustration of Molecular Flexibility     *
*   resulting from eigenspace based methods (NMA and ED)                *
*                                                                       *
*   It was Implemented and optimized from:                              *
*      "Automated Illustration of Molecular Flexibility"                *
*      Aaron Bryden, George Phillips Jr., Michael Gleicher              *
*      IEEE Transactions on Visualization and Computer Graphics         *
*      18(1) 132-145 (Jan 2012)                                         *
*                                                                       *
*************************************************************************
* This library is free software; you can redistribute it and/or modify  *
* it under the terms of the GNU General Public License as published by  *
* the Free Software Foundation; either version 2 of the License, or     *
* (at your option) any later version.                                   *
************************************************************************/

#include <libnma_affine.h>


// Builds the Rotational Correlation Matrix. See:
// "Automated Illustration of Molecular Flexibility"
// Aaron Bryden, George Phillips Jr., Michael Gleicher
// IEEE Transactions on Visualization and Computer Graphics, Volume 18, Number 1, page 132--145 - Jan 2012
void buildRotationalCorrelationMatrix(Macromolecule *mol, double *vector, float **p_corrMatrix)
{
	int i,j;
	float *correlationMatrix;
	float *correlationPatternMatrix;
	float *means;
	float *stdevs;
	float a[3],an[3],b[3],bn[3],pa[3],pb[3],axb[3],p1[3],p2[3],dummy1[3],dummy2[3],dummy3[3];
	float closestA[3],closestB[3],angularVelocityA[3],angularVelocityB[3];
	float da,db,ua,ub;
	float currentMean;
	float currentstd;
	float tempval;

	pdbIter *ii,*ij;
	ii = (pdbIter *) new pdbIter(mol);
	ij = (pdbIter *) new pdbIter(mol);

	int dim = mol->get_num_atoms(); //	myENM->consideredAtoms->size();
	correlationMatrix = new float[dim*dim]; // new float[myENM->consideredAtoms->size()*myENM->consideredAtoms->size()];
	means = new float[dim];
	stdevs = new float[dim];
	correlationPatternMatrix = new float[dim*dim];
	*p_corrMatrix = correlationPatternMatrix; // Outputs "correlationPatternMatrix"

	// Initialization
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			correlationMatrix[i+j*dim] = 0.0;

//	for(int eigenVec = 6; eigenVec < 7; eigenVec++)
//	for(int i = 0; i < dim; i++)

	for ( ii->pos_atom = 0; !ii->gend_atom(); ii->next_atom() ) // screen atoms
	{
		i = ii->pos_atom; // it's shorter... (I'm a lazy guy)
//		for(int j = 0; j < dim; j++)
		for ( ij->pos_atom = 0; !ij->gend_atom(); ij->next_atom() ) // screen atoms
		{
			j = ij->pos_atom; // it's shorter... (I'm a lazy guy)
			if(i==j)
			{
				correlationMatrix[i+j*dim] = 1.0;
			}
			else
			{
//				MathVec3D<float> a;// = myENM->firstEigen[i];
//				MathVec3D<float> b;// = myENM->firstEigen[j];

				a[0] = vector[i*3]; // (*(myENM->eigenVectors))(i*3,eigenVec);
				a[1] = vector[i*3+1]; // (*(myENM->eigenVectors))(i*3+1,eigenVec);
				a[2] = vector[i*3+2]; // (*(myENM->eigenVectors))(i*3+2,eigenVec);

				b[0] = vector[j*3]; // (*(myENM->eigenVectors))(j*3,eigenVec);
				b[1] = vector[j*3+1]; // (*(myENM->eigenVectors))(j*3+1,eigenVec);
				b[2] = vector[j*3+2]; //(*(myENM->eigenVectors))(j*3+2,eigenVec);

//				MathVec3D<float> an = a;
//				MathVec3D<float> bn = b;
				normalize(a,an);
				normalize(b,bn);

//				printf("a= %f %f %f  an= %f %f %f\n",a[0],a[1],a[2],an[0],an[1],an[2]);
//				printf("b= %f %f %f  bn= %f %f %f\n",b[0],b[1],b[2],bn[0],bn[1],bn[2]);

//				MathVec3D<float> pa = (*(myENM->consideredAtoms))[i];
//				MathVec3D<float> pb = (*(myENM->consideredAtoms))[j];
				(ii->get_atom())->getPosition(pa);
				(ij->get_atom())->getPosition(pb);
//				printf("pa= %f %f %f  pb= %f %f %f\n",pa[0],pa[1],pa[2],pb[0],pb[1],pb[2]);

//				MathVec3D<float> axb = a^b;
				crossProduct(a,b,axb);

//				if(sqrMagnitude(axb) < 0.000001)
//					continue;

				da = -dotProduct(a,pa);
				db = -dotProduct(b,pb);
//				printf("da= %f  db= %f\n",da,db);

				//calculating intersection line between the 2 planes, p1 and p2 are points on the line
//				MathVec3D<float> p1 = (((db*a) - (da*b))^axb)/axb.sqrMagnitude();
				multByScalar(a,db,dummy1);
				multByScalar(b,da,dummy2);
				subtract(dummy1,dummy2,dummy3);
				crossProduct(dummy3,axb,p1);
				multByScalar(p1,1/magnitude(axb));

//				MathVec3D<float> p2 = p1 + axb;
				add(p1,axb,p2);

//				float ua = (axb.dotProduct(pa-p1))/(axb.dotProduct(p2-p1));
				subtract(pa,p1,dummy1);
				subtract(p2,p1,dummy2);
				ua = dotProduct(axb,dummy1) / dotProduct(axb,dummy2);

//				float ub = (axb.dotProduct(pb-p1))/(axb.dotProduct(p2-p1));
				subtract(pb,p1,dummy1);
				// subtract(p2,p1,dummy2);
				ub = dotProduct(axb,dummy1) / dotProduct(axb,dummy2);

//				MathVec3D<float> closestA = p1 + (axb*ua);
				multByScalar(axb,ua,dummy1);
				add(p1,dummy1,closestA);

//				MathVec3D<float> closestB = p1 + (axb*ub);
				multByScalar(axb,ub,dummy1);
				add(p1,dummy1,closestB);

//				MathVec3D<float> angularVelocityA = ((closestA-pa)^a)/(closestA-pa).sqrMagnitude();
				subtract(closestA,pa,dummy1);
				crossProduct(dummy1,a,dummy2);
				multByScalar(dummy2, 1/sqrMagnitude(dummy1), angularVelocityA);

//				MathVec3D<float> angularVelocityB = ((closestB-pb)^b)/(closestB-pb).sqrMagnitude();
				subtract(closestB,pb,dummy1);
				crossProduct(dummy1,b,dummy2);
				multByScalar(dummy2, 1/sqrMagnitude(dummy1), angularVelocityB);

//				if( angularVelocityA.sqrMagnitude() < angularVelocityB.sqrMagnitude())
				if( sqrMagnitude(angularVelocityA) < sqrMagnitude(angularVelocityB) )
//					correlationMatrix[i+j*dim] += angularVelocityA.dotProduct(angularVelocityB)/sqrt(angularVelocityA.sqrMagnitude()*angularVelocityB.sqrMagnitude()) * angularVelocityA.sqrMagnitude()/angularVelocityB.sqrMagnitude()/1.0;
					correlationMatrix[i+j*dim] += dotProduct(angularVelocityA,angularVelocityB) / sqrt(sqrMagnitude(angularVelocityA)*sqrMagnitude(angularVelocityB)) * sqrMagnitude(angularVelocityA) / sqrMagnitude(angularVelocityB);
				else
//					correlationMatrix[i+j*dim] += angularVelocityA.dotProduct(angularVelocityB)/sqrt(angularVelocityA.sqrMagnitude()*angularVelocityB.sqrMagnitude()) * angularVelocityB.sqrMagnitude()/angularVelocityA.sqrMagnitude()/1.0;
					correlationMatrix[i+j*dim] += dotProduct(angularVelocityA,angularVelocityB) / sqrt(sqrMagnitude(angularVelocityA)*sqrMagnitude(angularVelocityB)) * sqrMagnitude(angularVelocityB) / sqrMagnitude(angularVelocityA);

				//float distancea = ((linePointb - intersectPoint)^(intersectPoint - pa)).magnitude()/(linePointb - interSectPoint).m;

				//float d = -axb.dotProduct(pb);
				//float adistbplane = (pa.dotProduct(axb) +d)/(axb.magnitude());
				//cout << "adistbplane = " << adistbplane << endl;
				//pa = pa - axb*adistbplane;
				//adistbplane = (pa.dotProduct(axb) +d)/(axb.magnitude());
				//cout << "adistbplane after = " << adistbplane << endl;

			}
		}
	}

//	for(int i = 0; i < myENM->consideredAtoms->size(); i++) {
	for(int i = 0; i < dim; i++)
	{
		currentMean = 0;
//		for(int j = 0; j < myENM->consideredAtoms->size(); j++)
		for(int j = 0; j < dim; j++)
			currentMean += correlationMatrix[i+j*dim];

		means[i] = currentMean/dim;
	}
//	for(int i = 0; i < myENM->consideredAtoms->size(); i++)
	for(int i = 0; i < dim; i++)
	{
		currentstd = 0;
//		for(int j = 0; j < myENM->consideredAtoms->size(); j++)
		for(int j = 0; j < dim; j++)
			currentstd += (correlationMatrix[i+j*dim] - means[i])*(correlationMatrix[i+j*dim] - means[i]);

		stdevs[i] = sqrt(currentstd/dim);
	}
//	for(int i = 0; i < myENM->consideredAtoms->size(); i++) {
	for(int i = 0; i < dim; i++)
//		for(int j = 0; j < myENM->consideredAtoms->size(); j++)
		for(int j = 0; j < dim; j++)
		{
			tempval = 0;

			// Mon: This seems to smooth the similarity matrix (correlatinMatrix)...
			for(int k = 0; k < dim; k++)
				tempval += correlationMatrix[i+k*dim]*correlationMatrix[j+k*dim]-means[i]*means[j];

			correlationPatternMatrix[i+j*dim] = (tempval/dim)/(stdevs[i]*stdevs[j]);
		}

//	printf("imodview> Dump some debug output \n");
	FILE *f_affine;
	if( (f_affine = fopen("affine.txt", "w")) == NULL)
//	if( (f_affine = fopen("affineMat.txt", "w")) == NULL)
	{
		printf("imodview> Affine model correlationPatternMatrix file writing error!! \n");
		exit(1);
	}
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			fprintf(f_affine,"%d %d %f\n",i,j,correlationPatternMatrix[i+j*dim]);
//			fprintf(f_affine,"%d %d %f\n",i,j,correlationMatrix[i+j*dim]);
	fclose(f_affine);

//	fstream pgmout, pgmout2;
//	pgmout.open("correlationMatrix.pgm",fstream::out);
//	pgmout2.open("correlationPatternMatrix.pgm",fstream::out);
//	pgmout << "P2" << endl;
//	pgmout << dim << " ";
//	pgmout << dim << endl;
//	pgmout << "255" << endl;
//	for(int i = 0; i < dim; i++) {
//		for(int j = 0; j < dim; j++) {
//			int outval = 255*( correlationMatrix[i+j*dim]*.5 +.5);
//
//			pgmout<< outval << " ";
//		}
//	}
//	pgmout2 << "P2" << endl;
//	pgmout2 << dim << " ";
//	pgmout2 << dim << endl;
//	pgmout2 << "255" << endl;
//	for(int i = 0; i < dim; i++) {
//		for(int j = 0; j < dim; j++) {
//			int outval = 255*( correlationPatternMatrix[i+j*dim]*.5 +.5);
//
//			pgmout2<< outval << " ";
//		}
//	}
//	pgmout.close();
//	pgmout2.close();
}

// Builds the Dot-Product-based Correlation Matrix.
void buildDotProductCorrelationMatrix(Macromolecule *mol, double *vector, float **p_corrMatrix)
{
	int i,j;
	float *correlationMatrix;
	float a[3],b[3],an[3],bn[3];

	pdbIter *ii,*ij;
	ii = (pdbIter *) new pdbIter(mol);
	ij = (pdbIter *) new pdbIter(mol);

	int dim = mol->get_num_atoms(); //	myENM->consideredAtoms->size();
	correlationMatrix = new float[dim*dim]; // new float[myENM->consideredAtoms->size()*myENM->consideredAtoms->size()];
	*p_corrMatrix = correlationMatrix; // Outputs "correlationMatrix"

	// Initialization
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			correlationMatrix[i+j*dim] = 0.0;

	for ( ii->pos_atom = 0; !ii->gend_atom(); ii->next_atom() ) // screen atoms
	{
		i = ii->pos_atom; // it's shorter... (I'm a lazy guy)
		for ( ij->pos_atom = 0; !ij->gend_atom(); ij->next_atom() ) // screen atoms
		{
			j = ij->pos_atom; // it's shorter... (I'm a lazy guy)
			if(i==j)
			{
				correlationMatrix[i+j*dim] = 1.0;
			}
			else
			{
				a[0] = vector[i*3]; // (*(myENM->eigenVectors))(i*3,eigenVec);
				a[1] = vector[i*3+1]; // (*(myENM->eigenVectors))(i*3+1,eigenVec);
				a[2] = vector[i*3+2]; // (*(myENM->eigenVectors))(i*3+2,eigenVec);

				b[0] = vector[j*3]; // (*(myENM->eigenVectors))(j*3,eigenVec);
				b[1] = vector[j*3+1]; // (*(myENM->eigenVectors))(j*3+1,eigenVec);
				b[2] = vector[j*3+2]; //(*(myENM->eigenVectors))(j*3+2,eigenVec);

				normalize(a,an);
				normalize(b,bn);

				correlationMatrix[i+j*dim] += dotProduct(an,bn);
			}
		}
	}

	// Some testing output...
	FILE *f_affine;
	if( (f_affine = fopen("affineDot.txt", "w")) == NULL)
	{
		printf("imodview> Dot product-based correlationMatrix file writing error!! \n");
		exit(1);
	}
	for(int i = 0; i < dim; i++)
		for(int j = 0; j < dim; j++)
			fprintf(f_affine,"%d %d %f\n",i,j,correlationMatrix[i+j*dim]);
	fclose(f_affine);
}

// This function is intended to compute the Affine transformation: v = M·x + t
// (it uses homogeneous coordinates, i.e. a 4 elements vector to define a 3D position)
// coords --> array of macromolecular coordinates
// natoms --> number of atoms
// vector --> pointer to the selected eigenvector
// cluster --> pointer to the array of atom indices belonging to the cluster of atoms
// nclus --> number of cluster atoms
// p_expMat --> pointer to the OUTPUT Affine transformation matrix (=NULL --> automatic memory allocation)
// weights --> array with weighting factors (OPTIONAL)
void getExpMat(float *coords, int natoms, double *vector, int *cluster, int nclus, double **p_expMat, float *weights)
{
	int i,j,k;
	double *expMat;
	double *vec; // current atom eigenvector
	double weight = 1.0;
	float *pos; // atomic position
	double totalWeight = 0.0;
	double jjt[4][4], jb[4][3];
	double jjt_array[16],jb_array[12];
	double x[12]; // 4x3 = 12 elements

	if(*p_expMat == NULL)
	{
		if( !(expMat = (double *) malloc(sizeof(double) * 16)) )
		{
			fprintf(stderr,"getExpMat: Memory allocation failure!\nForcing exit!\n");
			exit(1);
		}
		*p_expMat = expMat; // Output Affine matrix
	}
	else
		expMat = *p_expMat; // Using previously allocated memory...

	jjt[0][0] = jb[0][0] = 0.0;
	jjt[0][1] = jb[0][1] = 0.0;
	jjt[0][2] = jb[0][2] = 0.0;
	jjt[0][3] = 0.0;
	jjt[1][0] = jb[1][0] = 0.0;
	jjt[1][1] = jb[1][1] = 0.0;
	jjt[1][2] = jb[1][2] = 0.0;
	jjt[1][3] = 0.0;
	jjt[2][0] = jb[2][0] = 0.0;
	jjt[2][1] = jb[2][1] = 0.0;
	jjt[2][2] = jb[2][2] = 0.0;
	jjt[2][3] = 0.0;
	jjt[3][0] = jb[3][0] = 0.0;
	jjt[3][1] = jb[3][1] = 0.0;
	jjt[3][2] = jb[3][2] = 0.0;
	jjt[3][3] = 0.0;

	if(weights != NULL) // If there exist weights...
	{
		for(i = 0; i < nclus; i++)
		{
			weight = weights[i];
			totalWeight += weight;
			k = cluster[i]; // cluster's atom index
			vec = vector + 3*k; // the pointer to the eigenvector of the current atom
			pos = coords + 3*k;
			//		fprintf(stderr,"i= %3d   atom= %3d   vec= %f %f %f   pos= %f %f %f\n",i,k,vec[0],vec[1],vec[2],pos[0],pos[1],pos[2]);

			jjt[0][0] += weight*pos[0]*pos[0];
			jjt[0][1] += weight*pos[0]*pos[1];
			jjt[0][2] += weight*pos[0]*pos[2];
			jjt[0][3] += weight*pos[0];
			jjt[1][0] += weight*pos[1]*pos[0];
			jjt[1][1] += weight*pos[1]*pos[1];
			jjt[1][2] += weight*pos[1]*pos[2];
			jjt[1][3] += weight*pos[1];
			jjt[2][0] += weight*pos[2]*pos[0];
			jjt[2][1] += weight*pos[2]*pos[1];
			jjt[2][2] += weight*pos[2]*pos[2];
			jjt[2][3] += weight*pos[2];
			jjt[3][0] += weight*pos[0];
			jjt[3][1] += weight*pos[1];
			jjt[3][2] += weight*pos[2];
			jjt[3][3] += weight;

			jb[0][0] += weight*vec[0]*pos[0];
			jb[1][0] += weight*vec[0]*pos[1];
			jb[2][0] += weight*vec[0]*pos[2];
			jb[3][0] += weight*vec[0];
			jb[0][1] += weight*vec[1]*pos[0];
			jb[1][1] += weight*vec[1]*pos[1];
			jb[2][1] += weight*vec[1]*pos[2];
			jb[3][1] += weight*vec[1];
			jb[0][2] += weight*vec[2]*pos[0];
			jb[1][2] += weight*vec[2]*pos[1];
			jb[2][2] += weight*vec[2]*pos[2];
			jb[3][2] += weight*vec[2];
		}
	}
	else // No weighting...
	{
		for(i = 0; i < nclus; i++)
		{
			k = cluster[i]; // cluster's atom index
			vec = vector + 3*k; // the pointer to the eigenvector of the current atom
			pos = coords + 3*k;
			//		fprintf(stderr,"i= %3d   atom= %3d   vec= %f %f %f   pos= %f %f %f\n",i,k,vec[0],vec[1],vec[2],pos[0],pos[1],pos[2]);

			jjt[0][0] += pos[0]*pos[0];
			jjt[0][1] += pos[0]*pos[1];
			jjt[0][2] += pos[0]*pos[2];
			jjt[0][3] += pos[0];
			jjt[1][0] += pos[1]*pos[0];
			jjt[1][1] += pos[1]*pos[1];
			jjt[1][2] += pos[1]*pos[2];
			jjt[1][3] += pos[1];
			jjt[2][0] += pos[2]*pos[0];
			jjt[2][1] += pos[2]*pos[1];
			jjt[2][2] += pos[2]*pos[2];
			jjt[2][3] += pos[2];
			jjt[3][0] += pos[0];
			jjt[3][1] += pos[1];
			jjt[3][2] += pos[2];
			jjt[3][3] += 1.0;

			jb[0][0] += vec[0]*pos[0];
			jb[1][0] += vec[0]*pos[1];
			jb[2][0] += vec[0]*pos[2];
			jb[3][0] += vec[0];
			jb[0][1] += vec[1]*pos[0];
			jb[1][1] += vec[1]*pos[1];
			jb[2][1] += vec[1]*pos[2];
			jb[3][1] += vec[1];
			jb[0][2] += vec[2]*pos[0];
			jb[1][2] += vec[2]*pos[1];
			jb[2][2] += vec[2]*pos[2];
			jb[3][2] += vec[2];
		}

		totalWeight = (double)nclus;
	}

	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
		{
			jjt_array[4*j + i] = jjt[i][j] / totalWeight;
			if(j<3)
				jb_array[4*j + i] = jb[i][j] / totalWeight;
		}

//	show_matrix_standard(jjt_array,4,"A-matrix before fed to DPTRF");

	solveCholesky(jjt_array,4,jb_array,3,x); // Solves linear system

	for(i=0; i<4; i++)
	{
		expMat[i*4 + 3] = 0.0; // 12 = 3*4
		for(j=0; j<3; j++)
			expMat[i*4 + j] = x[j*4 + i]; // Transposing...
	}

}


// Returns the vector length (magnitude)
float magnitude(float *v,int dim)
{
	double accum = 0.0;
	for(int i=0; i<dim; i++)
		accum += v[i]*v[i];
	return( (float)sqrt(accum) );
}

// Returns the vector length (magnitude)
double magnitude(double *v,int dim)
{
	double accum = 0.0;
	for(int i=0; i<dim; i++)
		accum += v[i]*v[i];
	return( sqrt(accum) );
}

// Returns the squared vector length (magnitude)
float sqrMagnitude(float *v)
{
	return( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
}

// Normalizes a vector (1 length): c = norm(c)
void normalize(float *v,float *c)
{
	float vecLength = magnitude(v);
	if( vecLength != 0 )
	{
		c[0] = v[0] / vecLength;
		c[1] = v[1] / vecLength;
		c[2] = v[2] / vecLength;
	}
	else
	{
		fprintf(stderr,"Warning: called normalize() on the zero vector\n");
	}
}

// Computes the cross-product: a x b = c
void crossProduct(float *a, float *b, float *c)
{
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
}

// Computes the dot-product: a · b = c
float dotProduct(float *a, float *b)
{
	return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

// Multiply vector by scalar: a *= s
void multByScalar(float *a, float s)
{
	a[0] *= s;
	a[1] *= s;
	a[2] *= s;
}

// Multiply vector by scalar: a * s = c
void multByScalar(float *a, float s, float *c)
{
	c[0] = a[0] * s;
	c[1] = a[1] * s;
	c[2] = a[2] * s;
}

// Subtract two vectors: a - b = c
void subtract(float *a, float *b, float *c)
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

// Add two vectors: a - b = c
void add(float *a, float *b, float *c)
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}


// Solves linear system
void solveCholesky( double *A, int colsA, double *B, int colsB, double *X )
{
	int info;

//	show_matrix_standard(A,4,"A-matrix feeded to DPTRF");

	// Compute the Cholesky decomposition of A, leaving A itself unchanged.
	char uplo = 'U';
	dpotrf_( &uplo, &colsA, A, &colsA, &info );

//	show_matrix_standard(A,4,"A-matrix processed from DPTRF");

	if(info != 0)
		fprintf(stderr,"Cholesky Decomposition by LAPACK's dpotrf_ failed! info=%d\n",info);

	// Zero any element below the diagonal.
	for( int i = 0; i < colsA; i++ )
		for( int j = 0; j < i; j++ )
			A[j*colsA + i] = 0.0;

	// Call the LAPACK routine for solving the linear system (dpotrs_).  The function overwrites the right
	// hand side with the solution, so copy B into X and pass in X as the right hand side.
	// X = B;
	for( int i = 0; i < colsA; i++ ) // rows
		for( int j = 0; j < colsB; j++ ) // cols
			X[j*colsA + i] = B[j*colsA + i];
	dpotrs_( &uplo, &colsA, &colsB, A, &colsA, X, &colsA, &info );

	if(info != 0)
		fprintf(stderr,"Solving linear system with LAPACK's dpotrs_ failed! info=%d\n",info);
}

// Assumes bottom row of matrix "m" is [0,0,0,1]
void multMatrix_4x4_x_3(double *m,double *vec, double *ret)
{
	ret[0] = m[0]*vec[0] + m[4]*vec[1] + m[8]*vec[2]  + m[12];
    ret[1] = m[1]*vec[0] + m[5]*vec[1] + m[9]*vec[2]  + m[13];
    ret[2] = m[2]*vec[0] + m[6]*vec[1] + m[10]*vec[2] + m[14];
}

// Assumes bottom row of matrix "m" is [0,0,0,1] (single precision)
void multMatrix_4x4_x_3f(double *m,float *vec, float *ret)
{
	ret[0] = m[0]*vec[0] + m[4]*vec[1] + m[8]*vec[2]  + m[12];
    ret[1] = m[1]*vec[0] + m[5]*vec[1] + m[9]*vec[2]  + m[13];
    ret[2] = m[2]*vec[0] + m[6]*vec[1] + m[10]*vec[2] + m[14];
}

// Generic version... using Cramer's rule
// Note: this assumes the matrix is invertible (det != 0)
void invertMatrix4x4(double *m)
{
//    double *src = getTranspose().m;
    double tmp[12]; /* temp array for pairs */
    double src[16];
    transposeMatrix4x4(m,src);

    /* calculate pairs for first 8 elements (cofactors) */
    tmp[0] = src[10] * src[15];
    tmp[1] = src[11] * src[14];
    tmp[2] = src[9] * src[15];
    tmp[3] = src[11] * src[13];
    tmp[4] = src[9] * src[14];
    tmp[5] = src[10] * src[13];
    tmp[6] = src[8] * src[15];
    tmp[7] = src[11] * src[12];
    tmp[8] = src[8] * src[14];
    tmp[9] = src[10] * src[12];
    tmp[10] = src[8] * src[13];
    tmp[11] = src[9] * src[12];

    /* calculate first 8 elements (cofactors) */
    m[0] =  tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
    m[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
    m[1] =  tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
    m[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
    m[2] =  tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
    m[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
    m[3] =  tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
    m[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
    m[4] =  tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
    m[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
    m[5] =  tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
    m[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
    m[6] =  tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
    m[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
    m[7] =  tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
    m[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];

    /* calculate pairs for second 8 elements (cofactors) */
    tmp[0] = src[2]*src[7];
    tmp[1] = src[3]*src[6];
    tmp[2] = src[1]*src[7];
    tmp[3] = src[3]*src[5];
    tmp[4] = src[1]*src[6];
    tmp[5] = src[2]*src[5];
    tmp[6] = src[0]*src[7];
    tmp[7] = src[3]*src[4];
    tmp[8] = src[0]*src[6];
    tmp[9] = src[2]*src[4];
    tmp[10] = src[0]*src[5];
    tmp[11] = src[1]*src[4];

    /* calculate second 8 elements (cofactors) */
    m[8]  =  tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
    m[8]  -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
    m[9]  =  tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
    m[9]  -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
    m[10] =  tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
    m[10] -= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
    m[11] =  tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
    m[11] -= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
    m[12] =  tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
    m[12] -= tmp[4]*src[11] + tmp[0]*src[9]  + tmp[3]*src[10];
    m[13] =  tmp[8]*src[11] + tmp[0]*src[8]  + tmp[7]*src[10];
    m[13] -= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
    m[14] =  tmp[6]*src[9]  + tmp[11]*src[11] + tmp[3]*src[8];
    m[14] -= tmp[10]*src[11] + tmp[2]*src[8]  + tmp[7]*src[9];
    m[15] =  tmp[10]*src[10] + tmp[4]*src[8]  + tmp[9]*src[9];
    m[15] -= tmp[8]*src[9]   + tmp[11]*src[10] + tmp[5]*src[8];

    /* calculate determinant */
    double det = src[0]*m[0] + src[1]*m[1] + src[2]*m[2] + src[3]*m[3];

    /* calculate matrix inverse */
    for(int j = 0; j < 16; j++)
        m[j] /= det;
}

// Transposes a 4x4 matrix "m" into "mat"
void transposeMatrix4x4(double *m,double *mat)
{
	int ptr = 0;
    for (int row = 0; row < 4; row++)
        for (int col = 0; col < 4; col++)
			mat[col*4 + row] = m[ptr++];
}

// Affine model integration.
// A --> Input affine model matrix (LAPACK's standard format)
// X --> Output (integrated) affine model matrix (LAPACK's standard format)
// scale --> Integration range
// iters --> Number of interations (Alexa's method?)
void getExp(double *Af, double *X, double scale, int iters)
{
	double A[16];
	for(int i=0; i<16; i++)
		A[i] = Af[i];

    double D[16],N[16],cX[16],dummy[16],temp,fact;

    for(int i=0; i<4; i++)
    	for(int j=0; j<4; j++)
    	{
    		if(i==j)
    			temp = 1.0;
    		else
    			temp = 0.0;
    		D[j*4 + i] = temp;
    		N[j*4 + i] = temp;
    		cX[j*4 + i] = temp;
    		X[j*4 + i] = temp;
    	}

	double c = 1.0;

	for(int i = 0; i<16; i++)
		A[i] *= scale;

	double nA = 0.0;
	for(int i = 0; i < 16; i++)
		nA += A[i]*A[i];

	int sqrs = 0;
	double sqv = 1.0;

	while(nA > 1 && (sqrs<50))
	{
		nA /= 2;
		sqrs++;
		sqv /=2;
	}

	for(int i = 0; i<16; i++)
		A[i] *= sqv;

	for (int k = 1; k <= iters; k++)
	{
		c = c*((double)(iters-k+1)) / ((double)(k*(2*iters-k+1)));
		multMatrix_4x4_x_4x4(A, X, dummy);

		for(int i = 0; i<16; i++)
			dummy[i] *= c;

		for(int i = 0; i<16; i++)
			N[i] += dummy[i];

		fact = double((k%2) ? -1 : 1);
		for(int i = 0; i<16; i++)
			dummy[i] *= fact;

		for(int i = 0; i<16; i++)
			D[i] += dummy[i];
	}

	invertMatrix4x4(D);

	multMatrix_4x4_x_4x4(D,N,X);

	while (sqrs-- > 0)
	{
		multMatrix_4x4_x_4x4(X,X,dummy);
		for(int i = 0; i<16; i++)
			X[i] = dummy[i];
    }
}

// Multiply two 4x4 matrices ("A" and "B") into "C"
void multMatrix_4x4_x_4x4(double *A, double *B, double *C)
{
	int col,row,ptr;
    for( col = 0; col < 4; col++ )
        for( row = 0; row < 4; row++ )
        {
    		C[col*4 + row] = 0.0;
			for( ptr = 0; ptr < 4; ptr++ )
				C[col*4 + row] += A[ptr*4 + row] * B[col*4 + ptr];
        }
}

// Shows a standard format matrix (standard output)
void show_matrix_standard(double *matrix, int size, char *name)
{
	int i=0;
	int j=0;
	fprintf(stderr,"%s\n",name);
	for ( i = 0; i < size; i++ )
	{
		for ( j = 0; j < size; j++ )
			fprintf(stderr,"%12.5e ",matrix[j*size + i]);
		fprintf(stderr,"\n");
	}
}

// This calculates the Error (formula XX in Bryden's paper).
// The Error is the distance between the eigenvector component and the affine model evaluation.
// pos --> atom position
// vec --> eigenvector component for said atom
// expMat --> Affine transformation matrix (4x4)
double getError(float *pos, double *vec, double *expMat)
{
	double dist; // distance between two vectors
	float vec2[3]; // Affine matrix evaluated vector

	// Apply affine transformation to given coordinates (pos)
	multMatrix_4x4_x_3f(expMat, pos, vec2); // Assumes bottom row of matrix is [0,0,0,1]

	// The Error is the distance between the two vectors: vec and vec2 (from eigenvector and from affine model, respectively)
	dist = sqrt( (vec2[0]-vec[0])*(vec2[0]-vec[0]) +
				 (vec2[1]-vec[1])*(vec2[1]-vec[1]) +
				 (vec2[2]-vec[2])*(vec2[2]-vec[2])	);

	return(dist);
}

// Projects one point into a plane defined by a point and a normal. Output in "pos2"
void projectPointToPlane(float *pos, float *com, float *normAxis, float *pos2)
{
	float dummy[3];
	// The projection of a point q = (x, y, z) onto a plane given by a point p = (a, b, c) and a normal n = (d, e, f) is
	// q_proj = q - dot(q - p, n) * n                 , where |n|=1.
	// Here: q = pos, p = com, n = normAxis, and q_prog = pos2
	dummy[0] = pos[0] - com[0];
	dummy[1] = pos[1] - com[1];
	dummy[2] = pos[2] - com[2];
	float pnorm = dotProduct(dummy, normAxis);
	// the farthest atom position projected over cluster's middle-plane of the affine transformation.
	pos2[0] = pos[0] - pnorm * normAxis[0];
	pos2[1] = pos[1] - pnorm * normAxis[1];
	pos2[2] = pos[2] - pnorm * normAxis[2];
}

// Algorithm to generate a given number of spiral points uniformly distributed on the surface of a sphere.
// double r;		// true radius
// int n;		    // number of points
// float **p_points; // Pointer to the points-array (for dynamic memory allocation)
void spherePoints(double r, int n, float **p_points)
{
	/*
	Based on O'Rourke's program: "spiral.c", which is based on:
	Reference: E.B. Saff and A.B.J. Kuijlaars, Distributing Many Points on a Sphere,
	The Mathematical Intelligencer, 19(1), Winter (1997);
	-------------------------------------------------------------------------
	This program will generate a given number of spiral points uniformly
	distributed on the surface of a sphere.
	        The idea behind the algorithm is that one can cut the globe with
	N horizontal planes spaced 2/(N-1) units apart, forming N circles of
	latitude on the sphere, each latitude containing one spiral point.  To
	obtain the kth spiral point, one proceeds upward from the (k-1)st point
	(theta(k-1), phi(k-1)) along a great circle to the next latitude and
	travels counterclockwise along ti for a fixed distance to arrive at the
	kth point (theta(k), phi(k)).

	Written by Joseph O'Rourke and Min Xu, June 1997.
	Used in the textbook, "Computational Geometry in C."
	Questions to orourke@cs.smith.edu.
	--------------------------------------------------------------------
	This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
	redistributed in its entirety provided that this copyright notice is
	not removed.
	--------------------------------------------------------------------
	 */
	int k;                /* index */
	double phi1, phi, theta, h, x, y, z;
	double r1, r2, r3;	/* ellipsoid axis lengths */
	float *points;

	if(*p_points == NULL)
		// Allocate memory
		if( !(*p_points = (float *) malloc( sizeof(float)*n*3 ) ) )
		{
			fprintf(stderr,"spherePoints(): Memory allocation failure!\n");
			exit(1);
		}

	points = *p_points;

	r1 = r2 = r3 = r;
	phi1 = 0.0;

	points[0] = 0.0;
	points[1] = 0.0;
	points[2] = -1.0 * r3;

	for ( k = 2; k <= n - 1; k ++ )
	{
		/* Generate a random point on a sphere of radius 1. */
		h = -1.0 + 2.0 * ( k - 1.0 ) / ( double )( n - 1.0 );
		theta = acos ( h );

		if ( theta < 0 || theta > M_PI )
		{
			printf( "spherePoints: Error\n" );
			exit (1);
		}

		phi = phi1 + 3.6 / ( sqrt ( ( double )n * ( 1 - h * h ) ) );
		phi = fmod ( phi, 2 * M_PI );
		phi1 = phi;

		x = cos ( phi ) * sin ( theta );
		y = sin ( phi ) * sin ( theta );
		// z = cos ( theta ); But z==h, so:
		z = h;

		points[(k-1)*3] = r1*x;
		points[(k-1)*3 + 1] = r2*y;
		points[(k-1)*3 + 2] = r3*z;
	}

	points[(n-1)*3] = 0.0;
	points[(n-1)*3 + 1] = 0.0;
	points[(n-1)*3 + 2] = 1.0 * r3;
}

// Compute the squared distance between two points
float sqrDist(float *p1, float *p2)
{
	return( pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2) );
}

// Compute the distance between two points
float Dist(float *p1, float *p2)
{
	return( sqrt(pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2)) );
}

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
//   extend_cutoff --> Distance (in Angstroms) to extend small clusters. (Default=10)
void doHiearchicalExpMatClustering(Macromolecule *mol, float *coords, int requested, double *evecx, float cutoff, int ****p_rclusters, int ***p_rnatoms, float ***p_rerrors, int *rigid_bodies, int nrigid_bodies, int extend_size, float extend_cutoff)
{
	// NOTE: "extend_size" should be greater than "1" to avoid numerical instabilities in "solveCholesky()"...
	// In the paper they use a somewhat high value (extend_size=20), but it slows down clustering process very much.
	// We have found that values as small as 2 work well and that values of 5 or 6 do not slow down the process.
	float extend_cutoff2 = extend_cutoff * extend_cutoff; // Squared distance (Angstroms^2) to extend small clusters.
	bool debug = true;

	// 2. Creating all clusters (each atom/residue belongs to one cluster at start)
	//
	int **clusters; // Clusters array
	int *natoms; // Number of atoms per cluster array
	int num_atoms = mol->get_num_atoms(); // Get number of atoms
	int num_res = mol->get_num_fragments(); // Get number of residues
//	int nclusters = num_atoms; // One atom per cluster at start
//	int nclusters = num_res; // One residue per cluster at start
	int nclusters; // Initial number of clusters

	if(rigid_bodies == NULL || nrigid_bodies == 0)
		nclusters = num_res; // One residue per cluster at start
	else
		nclusters = nrigid_bodies; // One residue per rigid-body at start

	if( !(natoms = (int *) malloc( sizeof(int) * nclusters ) ) ) // Allocating maximum memory for the array of number of atoms...
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}
	float *errors;
	if( !(errors = (float *) malloc( sizeof(float) * nclusters) ) ) // Allocating maximum memory for the array of clusters...
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}
	if( !(clusters = (int **) malloc( sizeof(int *) * nclusters) ) ) // Allocating maximum memory for the array of clusters...
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}

//	// Initializing with one atom per cluster
//	for(int i=0; i<nclusters; i++)
//	{
//		if( !(clusters[i] = (int *) malloc( sizeof(int) * 1 ) ) )
//		{
//			printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
//			exit(1);
//		}
//		clusters[i][0] = i; // one atom per cluster (at start...)
//		natoms[i] = 1;
//		errors[i] = 0.0; // With just one element the affine matrix error should be zero...
//	}

	// CLUSTER INITIALIZATION
	pdbIter *iter = new pdbIter(mol);
	pdbIter *iter2;
	int cont_res = 0;
	int cont_atom = 0;
	int res_atoms;

	if(rigid_bodies == NULL || nrigid_bodies == 0) // Initializing with one residue per cluster
	{
		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() )
		{
			iter2 = new pdbIter( iter->get_fragment() ); // to screen residue atoms
			res_atoms = iter2->num_atom();
			if( !(clusters[cont_res] = (int *) malloc( sizeof(int) * res_atoms ) ) )
			{
				printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
				exit(1);
			}
			natoms[cont_res] = res_atoms;
			errors[cont_res] = 0.0; // With just one element the affine matrix error is set to zero... (who cares anyways!)
			for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
			{
				clusters[cont_res][iter2->pos_atom] = cont_atom; // storing atom indices of current cluster
				cont_atom++; // counting atoms
			}

			delete iter2;
			cont_res++; // counting residues
		}
	}
	else // Initializing with one rigid body per cluster
	{
		int curr_atom; // current cluster atom index
		int rb = -1; // current rigid body index
		int rb_old = -1; // previous rigid body index
		for( iter->pos_fragment = 0; !iter->gend_fragment(); iter->next_fragment() )
		{
			iter2 = new pdbIter( iter->get_fragment() ); // to screen residue atoms
			res_atoms = iter2->num_atom();

			rb_old = rb;
			rb = rigid_bodies[iter->pos_fragment];
			// fprintf(stderr,"res= %d  rb= %d  rb_old= %d  res_atoms= %d\n",iter->pos_fragment,rb,rb_old,res_atoms);
			if(rb != rb_old) // if next rigid-body...
			{   // malloc
				curr_atom = 0; // current cluster atom index
				natoms[rb] = res_atoms;
				errors[rb] = 0.0;
				if( !(clusters[rb] = (int *) malloc( sizeof(int) * res_atoms ) ) )
				{
					printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
					exit(1);
				}
			}
			else // if the same rigid-body...
			{   // realloc
				natoms[rb] += res_atoms;
				if( !(clusters[rb] = (int *) realloc(clusters[rb], sizeof(int) * natoms[rb] ) ) )
				{
					printf("doHiearchicalExpMatClustering: Memory re-allocation error...\nForcing exit!\n");
					exit(1);
				}
			}
			for( iter2->pos_atom = 0; !iter2->gend_atom(); iter2->next_atom() )
			{
				clusters[rb][curr_atom] = cont_atom; // storing atom indices of current cluster
				cont_atom++; // global counting of atoms
				curr_atom++; // count current cluster atoms
			}

			delete iter2;
			cont_res++; // counting residues
		}
	}
	delete iter;
	if(num_res != cont_res || num_atoms != cont_atom)
	{
		printf("ERROR, number of residues or atoms mismatch!\nForcing exit!\n");
		exit(1);
	}

	// Allocate memory for the array of requested clusters arrays.
	int ***rclusters;
	if(*p_rclusters == NULL)
	{
		if( !(rclusters = (int ***) malloc( sizeof(int **) * requested) ) )
		{
			printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_rclusters = rclusters; // output requested clusters
	}
	int **rnatoms;
	if(*p_rnatoms == NULL)
	{
		if( !(rnatoms = (int **) malloc( sizeof(int *) * requested) ) )
		{
			printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_rnatoms = rnatoms; // output requested clusters
	}
	float **rerrors;
	if(*p_rerrors == NULL)
	{
		if( !(rerrors = (float **) malloc( sizeof(float *) * requested) ) )
		{
			printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_rerrors = rerrors; // output requested errors
	}

	// Creating temporary cluster (of maximum attainable size)
	int *temporary; // Temporary cluster
	int ntemporary; // Number of atoms in current temporary cluster
	int nextended; // Number of atoms in current extended cluster
	if( !(temporary = (int *) malloc( sizeof(int) * num_atoms) ) ) // Allocating maximum memory for the temporary cluster...
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}

	float *weights;
	if( !(weights = (float *) malloc( sizeof(float) * num_atoms) ) ) // Allocating maximum memory for weights
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}
	float *raw_weights;
	if( !(raw_weights = (float *) malloc( sizeof(float) * num_atoms) ) ) // Allocating maximum memory for raw_weights
	{
		printf("doHiearchicalExpMatClustering: Memory allocation error...\nForcing exit!\n");
		exit(1);
	}

	// **************************************
	// * THE CLUSTERING ALGORITHM MAIN LOOP *
	// **************************************
	//
	bool progressMade = true;
	int besti, bestj;
    bool *cadj = NULL; // Cluster-wise adjacency matrix
    float *mdist = NULL; // Cluster-wise MIn dist matrix
    float *pos;
    float dist;
	double temp;
	double smallestDiff; // Store the "smallest maximum difference" between clusters
	double maxdiff;
	double diff;
	double avgdiff;
	double *tempMat = NULL;

	// 1. Building the first AdjacencyMatrix (at cluster level)
	computeClusterAdjMatrix(coords,clusters,nclusters,natoms,&cadj,cutoff);

	while(progressMade)
	{
		progressMade = false;

		// 5. Comparing all possible pair of clusters to search for the best candidate (using getError)
		//
		smallestDiff = 1E9; // Store the "smallest maximum difference" between clusters
		besti = -1;
		bestj = -1;
		int compcount = 0;
		int *buffcluster;
		int bufint;

		for(int i = 0; i < nclusters-1; i++)
			for(int j = i+1; j < nclusters; j++)
			{
				// Checking merge criteria 1...?
				if(cadj[i + (j-1)*j/2]) // If clusters "i" and "j" are near...
				{
					// fprintf(stderr,"clusters %d and %d are near\n",i,j);
					maxdiff = 0.0;
					diff = 0.0;
					avgdiff = 0.0;

					// Merging temporary cluster
					ntemporary = 0; // reset temporary cluster number of elements
					for(int k=0; k<natoms[i]; k++,ntemporary++) // screen cluster "i"
						temporary[ntemporary] = clusters[i][k]; // copy atomic indices from cluster "i"
					for(int k=0; k<natoms[j]; k++,ntemporary++) // screen cluster "j"
						temporary[ntemporary] = clusters[j][k]; // copy atomic indices from cluster "j"

					// Extending temporary cluster with neighboring elements "to provide more spatial coherence" (see paper)
					if(ntemporary <= extend_size) // Only checking this when clusters are small enough...
					{
						for(int x=0; x < ntemporary; x++) // current temporary cluster atoms
						{
							pos = coords + 3*temporary[x];
							nextended = ntemporary;
							for(int y = 0; y < num_atoms; y++) // all macromolecule atoms
								raw_weights[y] = sqrDist(pos,coords + 3*y);
						}

						for(int k = 0; k < num_atoms; k++)
							if(raw_weights[k] < extend_cutoff2)
							{
								temporary[nextended] = k; // add "k" atom to extended cluster
								weights[nextended] = (1-sqrt(raw_weights[k])/extend_cutoff) / ntemporary;
								nextended++;
							}

						for(int x=0; x < ntemporary; x++) // current temporary cluster atoms
							weights[x] = 1.0;

						getExpMat(coords,num_atoms,evecx,temporary,nextended,&tempMat,weights); // Computing Affine-matrix for the temporary merged cluster...
					}
					else // if big clusters
					{
//						fprintf(stderr,"i= %d  j= %d  ntemporary= %d\n",i,j,ntemporary);
						getExpMat(coords,num_atoms,evecx,temporary,ntemporary,&tempMat); // Computing Affine-matrix for the temporary merged cluster...
//						show_matrix_standard(tempMat, 4, "Affine matrix");
					}

					// Computing error for current temporary cluster
					for(int k = 0; k < ntemporary; k++)
						avgdiff += getError(coords + 3*temporary[k], evecx + 3*temporary[k], tempMat); // Compute error respect to the Affine-matrix...
					avgdiff /= ntemporary;

					// printf("compcount= %d  maxdiff= %e  avgdiff= %e\n",compcount,maxdiff,avgdiff);
					compcount++;

					// Check if current maximum error (stored in "diff") is smaller than the smallest one found so far...
					if(avgdiff < smallestDiff)
					{
						smallestDiff = avgdiff;
						besti = i;
						bestj = j;
					}
				}
			}

		// 6. If "progress was made", the temporary merge of clusters (i,j) with smallest differences become permanent...
		//
		if(besti > -1)
		{
			progressMade = true;

			// Storing error for current cluster.
			errors[besti] = smallestDiff;

			// 4. Update Cluster-Adjacency-matrix
			//    The "cluster-adjacency-matrix" is a boolean 2D matrix (true= close clusters, false= far clusters)
			//
			updateClusterAdjMatrix2(clusters, nclusters, natoms, cadj, cutoff, besti, bestj);

			ntemporary = 0; // reset temporary cluster number of elements
			for(int k=0; k<natoms[besti]; k++,ntemporary++) // screen cluster "i"
				temporary[ntemporary] = clusters[besti][k]; // copy atomic indices from cluster "i"
			for(int k=0; k<natoms[bestj]; k++,ntemporary++) // screen cluster "j"
				temporary[ntemporary] = clusters[bestj][k]; // copy atomic indices from cluster "j"
			if(debug)
			{
				fprintf(stdout,"\rClusters: %4d --> Merging %4d (%4d) and %4d (%4d) into %4d (%4d), error= %10.4e",nclusters,besti,natoms[besti],bestj,natoms[bestj],besti,ntemporary,smallestDiff);
				fflush(stdout);
			}

			free(clusters[besti]); // delete cluster "besti"
			if( !(clusters[besti] = (int *) malloc(sizeof(int) * ntemporary)) )
			{
				fprintf(stderr,"Failed memory allocation, ntemporary=%d!\nForcing exit!\n",ntemporary);
				exit(1);
			}
			natoms[besti] = ntemporary; // update "i" cluster number of atoms
			for(int k=0; k<ntemporary; k++) // screen cluster "i"
				clusters[besti][k] = temporary[k]; // copy atomic indices of the merged cluster into cluster "besti"

			buffcluster = clusters[bestj];
			clusters[bestj] = clusters[nclusters-1]; // move last cluster into "bestj" cluster
			clusters[nclusters-1] = buffcluster;
			bufint = natoms[bestj];
			natoms[bestj] = natoms[nclusters-1]; // update "j" cluster number of atoms
			natoms[nclusters-1] = bufint;

			// NOTE: "clusters" and "natoms" should be already updated! (but "nclusters" should not!)
			// 4. Update Cluster-Adjacency-matrix
			//    The "cluster-adjacency-matrix" is a boolean 2D matrix (true= close clusters, false= far clusters)
			//
			// updateClusterAdjMatrix(coords, clusters, nclusters, natoms, cadj, cutoff, besti, bestj);

			nclusters--; // decrease number of clusters since two have been merged

			// Leave this commented (do not delete)
			// showClusters(clusters, nclusters, natoms); // Shows all clusters and their atoms
		}

		// STORING OUTPUT
		if(nclusters <= requested)
		{
			if( !(rclusters[nclusters-1] = (int **) malloc(sizeof(int *) * nclusters)) ) // allocating current cluster
			{
				fprintf(stderr,"Failed memory allocation, nclusters=%d! (int**)\nForcing exit!\n",nclusters);
				exit(1);
			}
			if( !(rnatoms[nclusters-1] = (int *) malloc(sizeof(int) * nclusters)) ) // allocating number of atoms array for current cluster
			{
				fprintf(stderr,"Failed memory allocation, nclusters=%d! (int*)\nForcing exit!\n",nclusters);
				exit(1);
			}
			if( !(rerrors[nclusters-1] = (float *) malloc(sizeof(float) * nclusters)) ) // allocating errors array for current cluster
			{
				fprintf(stderr,"Failed memory allocation, nclusters=%d! (float)\nForcing exit!\n",nclusters);
				exit(1);
			}

			for(int i = 0; i < nclusters; i++) // screening clusters
			{
				rnatoms[nclusters-1][i] = natoms[i]; // store current cluster number of atoms
				rerrors[nclusters-1][i] = errors[i]; // store current cluster error
				if( !(rclusters[nclusters-1][i] = (int *) malloc(sizeof(int) * natoms[i])) ) // allocating atom indices array
				{
					fprintf(stderr,"Failed memory allocation, natoms[i]=%d!\nForcing exit!\n",natoms[i]);
					exit(1);
				}
				for(int j = 0; j < natoms[i]; j++) // screening cluster's atoms
					rclusters[nclusters-1][i][j] = clusters[i][j]; // copy atomic indices
			}
		}

		if(nclusters <= 1) // If finished the binary tree...
			progressMade = false;
	}
	if(debug)
		fprintf(stdout,"\n");

	// Free used memory (non-output)
	if(rigid_bodies == NULL || nrigid_bodies == 0)
		for(int i=0; i<num_res; i++) // One residue per cluster at start
			free(clusters[i]);
	else
		for(int i=0; i<nrigid_bodies; i++) // One residue per rigid-body at start
			free(clusters[i]);
	free(clusters);
	free(natoms);
	free(temporary);
	free(weights);
}

// Builds a boolean adjacency matrix at atom level.
// mol --> The macromolecule
// p_adj --> Pointer to the matrix (packed storage). Set *adj to NULL for automatic memory allocation.
// cutoff --> Proximity distance cutoff
void buildAdjacencyMatrix(Macromolecule *mol, bool **p_adj, float cutoff)
{
	bool *adj;
	pdbIter *iti = new pdbIter(mol);
	pdbIter *itj = new pdbIter(mol);
	float posi[3],posj[3];
	float cutoffsqr = cutoff*cutoff;
	int num_atoms = mol->get_num_atoms();

	if(*p_adj == NULL) // Enable memory allocation
	{
		if( !(adj = (bool *) malloc( sizeof(bool) * num_atoms*(num_atoms-1)/2 ) ) ) // Allocating packed storage without diagonal
		{
			printf("buildAdjacencyMatrix: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_adj = adj; // output matrix
	}
	else
		adj = *p_adj; // previously allocated matrix

	// !iter1->gend_atom(); iter1->next_atom()
	for(iti->pos_atom = 0; !iti->gend_atom(); iti->next_atom()) // rows
	{
		(iti->get_atom())->getPosition(posi); // Get i-th atom position
		for(itj->pos_atom = iti->pos_atom+1; !itj->gend_atom(); itj->next_atom()) // cols
		{
			(itj->get_atom())->getPosition(posj); // Get j-th atom position
			if(sqrDist(posi,posj) < cutoffsqr)
				adj[iti->pos_atom + (itj->pos_atom-1)*itj->pos_atom/2] = true;
			else
				adj[iti->pos_atom + (itj->pos_atom-1)*itj->pos_atom/2] = false;
			// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
		}
	}
	delete iti;
	delete itj;

}

// Compute adjacency matrix at cluster level. (It also can update, but inefficiently...)
// coords --> array of macromolecular coordinates
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// p_cadj --> pointer to the cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
void computeClusterAdjMatrix(float *coords, int **clusters, int nclusters, int *natoms, bool **p_cadj, float cutoff)
{
	bool *cadj=NULL;
	bool close = false;
	float *pos;
	float cutoff2 = cutoff*cutoff;

	if(*p_cadj == NULL) // Enable memory allocation
	{
		if( !(cadj = (bool *) malloc( sizeof(bool) * nclusters*(nclusters-1)/2 ) ) ) // Allocating packed storage without diagonal
		{
			printf("computeClusterAdjMatrix: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_cadj = cadj; // output matrix
	}
	else
		cadj = *p_cadj; // previously allocated matrix

	for(int i = 0; i < nclusters-1; i++) // rows
		for(int j = i+1; j < nclusters; j++) // cols
		{
			cadj[i + (j-1)*j/2] = false;
			close = false;
			for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
			{
				pos = coords + 3*clusters[i][a];
				for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
				{
//					fprintf(stderr,"i=%4d  j=%4d  a=%4d  b=%4d\n",i,j,a,b);
					if( sqrDist(pos,coords + 3*clusters[j][b]) < cutoff2 ) // If clusters' atoms are neighboring...
					{
						cadj[i + (j-1)*j/2] = true;
						close = true;
					}
				}
			}
		}
}

// Update adjacency matrix at cluster level upon merging.
// coords --> array of macromolecular coordinates
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// cadj --> cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
// i,j --> merged cluster indices (i<j)
// NOTE: "clusters" and "natoms" should be already updated! (but "nclusters" should not!)
void updateClusterAdjMatrix(float *coords, int **clusters, int nclusters, int *natoms, bool *cadj, float cutoff, int x, int y)
{
	int i,j;
	bool close = false;
	float *pos;
	float cutoff2 = cutoff*cutoff;

	// Copy cluster "nclusters-1" into "j"
	for(i = 0; i < y; i++) // screen rows (vertical)
		cadj[i + (y-1)*y/2] = cadj[i + (nclusters-2)*(nclusters-1)/2];
	for(j = y+1; j < nclusters-1; j++) // screen cols (horizontal)
		cadj[y + (j-1)*j/2] = cadj[j + (nclusters-2)*(nclusters-1)/2];

	// Re-compute: "i" vs. all clusters
	i = x;
	for(j = i+1; j < nclusters-1; j++) // cols
	{
		cadj[i + (j-1)*j/2] = false;
		close = false;
		for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
		{
			pos = coords + 3*clusters[i][a];
			for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
			{
				if( sqrDist(pos, coords+3*clusters[j][b]) < cutoff2 ) // If clusters' atoms are neighboring...
				{
					cadj[i + (j-1)*j/2] = true;
					close = true;
				}
			}
		}
		// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
	}
	j = x;
	for(i = 0; i < x; i++) // cols
	{
		cadj[i + (j-1)*j/2] = false;
		close = false;
		for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
		{
			pos = coords + 3*clusters[i][a];
			for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
			{
				if( sqrDist(pos, coords+3*clusters[j][b]) < cutoff2 ) // If clusters' atoms are neighboring...
				{
					cadj[i + (j-1)*j/2] = true;
					close = true;
				}
			}
		}
		// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
	}
}

// Update adjacency matrix at cluster level upon merging.
// clusters --> array of clusters (with atom indices within each element)
// nclusters --> current number of clusters
// natoms --> current number of atoms per cluster
// cadj --> cluster adjacency matrix
// cutoff --> distance cutoff for adjacency
// i,j --> merged cluster indices (i<j)
// NOTE: "clusters", "nclusters" and "natoms" should NOT have been updated yet!
void updateClusterAdjMatrix2(int **clusters, int nclusters, int *natoms, bool *cadj, float cutoff, int x, int y)
{
	int i,j;
	float *pos;
	bool *buf = (bool *) malloc( sizeof(bool) * nclusters );

	for(j=0; j<nclusters; j++)
		buf[j] = false; // far

	// Checking proximity: "x" vs. all clusters
	for(i = 0; i < x; i++) // rows
		if(cadj[i + (x-1)*x/2])
			buf[i] = true; // close
	for(j = x+1; j < nclusters; j++) // cols
		if(cadj[x + (j-1)*j/2])
			buf[j] = true; // close

	// Checking proximity: "y" vs. all clusters
	for(i = 0; i < y; i++) // rows
		if(cadj[i + (y-1)*y/2])
			buf[i] = true; // close
	for(j = y+1; j < nclusters; j++) // cols
		if(cadj[y + (j-1)*j/2])
			buf[j] = true; // close

	// Copy cluster "nclusters-1" into "y"
	for(i = 0; i < y; i++) // screen rows (vertical)
		cadj[i + (y-1)*y/2] = cadj[i + (nclusters-2)*(nclusters-1)/2];
	for(j = y+1; j < nclusters-1; j++) // screen cols (horizontal)
		cadj[y + (j-1)*j/2] = cadj[j + (nclusters-2)*(nclusters-1)/2];

	// Updating proximity: "x" vs. all clusters
	for(i = 0; i < x; i++) // rows
		cadj[i + (x-1)*x/2] = buf[i]; // close
	for(j = x+1; j < nclusters; j++) // cols
		cadj[x + (j-1)*j/2] = buf[j]; // close

	// Updating proximity: "last cluster" vs. (x+y)
	cadj[x + (y-1)*y/2] = buf[nclusters-1]; // close

	free(buf);
}

void computeClusterMinDistMatrix(float *coords, int **clusters, int nclusters, int *natoms, float **p_cadj)
{
	float *cadj,minDist;
	bool close = false;
	float dist;
	float *pos;

	if(*p_cadj == NULL) // Enable memory allocation
	{
		if( !(cadj = (float *) malloc( sizeof(float) * nclusters*(nclusters-1)/2 ) ) ) // Allocating packed storage without diagonal
		{
			printf("computeClusterAdjMatrix: Memory allocation error...\nForcing exit!\n");
			exit(1);
		}
		*p_cadj = cadj; // output matrix
	}
	else
		cadj = *p_cadj; // previously allocated matrix

	for(int i = 0; i < nclusters-1; i++) // rows
		for(int j = i+1; j < nclusters; j++) // cols
		{
			cadj[i + (j-1)*j/2] = false;
			close = false;
			minDist=1E9;
			for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
			{
				pos = coords + clusters[i][a];
				for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
				{
					dist = Dist(pos,coords + 3*clusters[j][b]);
					if(dist < minDist)
						minDist = dist;
				}
			}
			cadj[i + (j-1)*j/2] = minDist;
			// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
		}
}

void updateClusterMinDistMatrix(Macromolecule *mol,int **clusters, int nclusters, int *natoms, float *cadj, float cutoff, int x, int y)
{
	int i,j;
	bool close = false;
	pdbIter *iter = new pdbIter(mol);
	float pos[3],pos2[3];
	float cutoff2 = cutoff*cutoff;
	float minDist,dist;

	// Copy cluster "nclusters-1" into "j"
	for(i = 0; i < y; i++) // screen rows (vertical)
		cadj[i + (y-1)*y/2] = cadj[i + (nclusters-2)*(nclusters-1)/2];
	for(j = y+1; j < nclusters-1; j++) // screen cols (horizontal)
		cadj[y + (j-1)*j/2] = cadj[j + (nclusters-2)*(nclusters-1)/2];

	// Re-compute: "i" vs. all clusters
	i = x;
	for(j = i+1; j < nclusters-1; j++) // cols
	{
		cadj[i + (j-1)*j/2] = false;
		close = false;
		minDist=1E9;
		for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
		{
			iter->pos_atom = clusters[i][a];
			(iter->get_atom())->getPosition(pos);
			for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
			{
				iter->pos_atom = clusters[j][b];
				(iter->get_atom())->getPosition(pos2);
				dist = Dist(pos,pos2);
				if(dist < minDist)
					minDist = dist;
			}
		}
		cadj[i + (j-1)*j/2] = minDist;
		// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
	}
	j = x;
	for(i = 0; i < x; i++) // rows
	{
		cadj[i + (j-1)*j/2] = false;
		close = false;
		minDist=1E9;
		for(int a = 0; a < natoms[i] && !close; a++) // screening i-th cluster (a-th atoms)
		{
			iter->pos_atom = clusters[i][a];
			(iter->get_atom())->getPosition(pos);
			for(int b = 0; b < natoms[j] && !close; b++) // screening j-th cluster (b-th atoms)
			{
				iter->pos_atom = clusters[j][b];
				(iter->get_atom())->getPosition(pos2);
				dist = Dist(pos,pos2);
				if(dist < minDist)
					minDist = dist;
			}
		}
		cadj[i + (j-1)*j/2] = minDist;
		// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
	}
}

// Shows an Adjacency Matrix (boolean)
void showAdjMatrix(bool *adj, int n, char *name)
{
	printf("%s:\n",name);
	for(int i = 0; i < n; i++) // rows
	{
		for(int j = 0; j < n; j++) // cols
		{
			if( i<n-1 && j >= i+1 )
			{
				if(	adj[i + (j-1)*j/2] )
					printf(" 1");
				else
					printf(" 0");
				// Indexing with diagonal --> hess_matrix[a+b*(b+1)/2]
			}
			else
				printf(" -");
		}
		printf("\n");
	}
}

// Shows an Adjacency Matrix (boolean)
void showDistMatrix(float *adj, int n, char *name)
{
	printf("%s:\n",name);
	for(int i = 0; i < n; i++) // rows
	{
		for(int j = 0; j < n; j++) // cols
		{
			if( i<n-1 && j >= i+1 )
			{
				printf(" %5.2f",adj[i + (j-1)*j/2]);
			}
			else
				printf(" %5s","-");
		}
		printf("\n");
	}
}

// Shows all clusters and their atoms
void showClusters(int **clusters, int nclusters, int *natoms, float *errors)
{
	for(int n=0; n<nclusters; n++)
	{
		printf("Cluster %2d  Error %10.4e  Atoms:",n,errors[n]);
		for(int m=0; m<natoms[n]; m++)
			printf(" %d",clusters[n][m]);
		printf(" (%d)\n",natoms[n]);
	}
}
