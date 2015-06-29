#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "jacobi.h"

#define RND ((double) rand()/RAND_MAX)
#define SIZE 150

int main(void)
{
	int i,j,rot;
	double sum;
	gsl_matrix* A = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* V = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* AV = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* VD = gsl_matrix_alloc(SIZE,SIZE);
	gsl_vector* e = gsl_vector_alloc(SIZE);

	// We construct matrix A and make it symmetric
	for (i = 0; i < SIZE; i++)
	{
		gsl_matrix_set(A,i,i,sin(i)+cos(i*i));
		for (j = i+1; j < SIZE; j++)
		{
			gsl_matrix_set(A, i, j, RND);
			gsl_matrix_set(A, j, i, gsl_matrix_get(A,i,j));
		}
	}

	// perform cyclic jacobi diagonalisation
	rot = jacobi_cyclic(A,e,V,JACOBI_SORT_DESC);

	printf("To evaluate our implementation of Jacobi diagonalisation\n");
	printf("we first consider a symmetric matrix with entries\n");
	printf("\tA_{ij} = sin(i)+cos(j*i) for j >= i\n");
	printf("and size COL=ROW=%d.\n",SIZE);
	printf("A total of %d rotations were neccessary for convergence\n",rot);
	printf("The eigenvalues are\n");
	gsl_vector_fprintf(stdout,e,"\t%g");

	// restore A
	for(i=0; i<SIZE; i++)
	{
		for(j=i+1; j<SIZE; j++)
		{
			gsl_matrix_set(A,i,j,gsl_matrix_get(A,j,i));
		}
	}

	// calculate AV
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,A,V,0,AV);

	// make D = diag(e), saved in A, and calculate VD
	gsl_matrix_set_zero(A);
	for(i=0; i<SIZE; i++)
	{
		gsl_matrix_set(A,i,i,gsl_vector_get(e,i));
	}
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,V,A,0,VD);

	// if AV == VD diagonalisation was succesfull
	// calculate entrywise norm of |AV - VD|
	gsl_matrix_sub(AV,VD);
	sum = 0;
	for(i=0; i<SIZE; i++)
	{
		for(j=0; j<SIZE; j++)
		{
			sum += gsl_matrix_get(AV,i,j)*gsl_matrix_get(AV,i,j);
		}
	}
	printf("Since V^-1AV = D, D = diag(eigenvalues), we expect AV == VD\n");
	printf("The entrywise norm of |AV - DV| is found to be\n");
	printf("\t|AV - DV| = %g\n",sqrt(sum));

	gsl_matrix_view a = gsl_matrix_submatrix(A,0,0,4,4);
	gsl_matrix_view v = gsl_matrix_submatrix(V,0,0,4,4);
	gsl_vector_view l = gsl_vector_subvector(e,0,4);
	// Construct the Hilbert matrix, H_{ij} = 1/(i+j+1)
	for(i=0; i<4; i++)
	{
		gsl_matrix_set(&a.matrix,i,i,1.0/(i+i+1.0));
		for(j=i+1; j<4; j++)
		{
			gsl_matrix_set(&a.matrix,i,j,1.0/(i+j+1.0));
			gsl_matrix_set(&a.matrix,j,i,gsl_matrix_get(&a.matrix,i,j));
		}
	}
	// diagonalise it using Jacobi's algorithm
	rot = jacobi_cyclic(&a.matrix,&l.vector,&v.matrix,JACOBI_SORT_DESC);
	printf("\nWe now consider the Hilbert matrix, H_ij = 1/(i+j+1)\n");
	printf("From GSL's example on eigensystems the eigenvalues for\n");
	printf("the 4th order Hilbert matrix are\n");
	printf("\t9.6702e-05\n\t6.7383e-03\n\t1.6914e-01\n\t1.5002e+00\n");
	printf("Using Jacobi's cyclic algorithm we find them to be\n");
	gsl_vector_fprintf(stdout,&l.vector,"\t%g");
	printf("A total of %d rotations were necessary.\n",rot);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_matrix_free(AV);
	gsl_matrix_free(VD);
	gsl_vector_free(e);
}
