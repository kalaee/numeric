#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "cholesky.h"

#define	SIZE	20

int main(void)
{
	int i,j;
	// define the system
	gsl_matrix* A = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* B = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* AI = gsl_matrix_alloc(SIZE,SIZE);
	gsl_vector* y = gsl_vector_alloc(SIZE);
	gsl_vector* x = gsl_vector_alloc(SIZE);

	// define the entries of A
	for (i = 0; i < SIZE; i++)
	{
		gsl_vector_set(y,i,i*sin(i*i+i));
		for (j = 0; j <= i; j++)
		{
			gsl_matrix_set(A, j, i, 1+i*j+j*sin (i+j) + i*cos (j*i));
			if (i != j)
			{
				gsl_matrix_set(A,j,i,gsl_matrix_get(A,i,j));
			}
		}
	}
	// construct the vector b
	gsl_blas_dgemv(CblasNoTrans,1,A,y,0,x);

	// copy A into B for preservation
	gsl_matrix_memcpy(B,A);

	/* perform Cholesky decomposition on A and solve Ax = Ay for x
	 * determine the determinant and matrix inverse */
	chol_dec(A);
	chol_bak(A,x);
	double d = chol_det(A);
	chol_inv(A,AI);

	// find the norm of the vector x-y, if zero: the solution is correct
	gsl_vector_sub(x,y);
	printf("Solve Ax = b for x, where b = Ay, using Cholesky decomp\n");
	printf("Evaluating the deviation between x and y:\n");
	printf("\t|x-y| =\t%g\n",gsl_blas_dnrm2(x));

	// determinant from GSL using LU decomp
	gsl_matrix_memcpy(A,B);
	gsl_permutation* p = gsl_permutation_alloc(SIZE);
	gsl_linalg_LU_decomp(A,p,&i);
	double dgsl = gsl_linalg_LU_det(A,i);
	printf("\nCompare the algorithm for computation");
	printf(" of the absolute value of the determinant");
	printf(" with that found using LU decomp in GSL\n");
	printf("\tdet(A)/det(A)_gsl - 1 =\t%g\n",d/dgsl-1);

	/* matrix product B*Binv is computed and saved in a.matrix
	 * to evaluate the succes of the computation the entrywise norm
	 * of the matrix B*Binv - I, where I is the identity matrix
	 * is evaluated */
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,B,AI,0,A);
	printf("\nCompute the inverse matrix Binv of B. The entrywise norm of B*Binv - I is\n");
	printf("\t|B*Binv - I| =\t");
	d = 0;
	for(i=0; i < SIZE; i++)
	{
		gsl_matrix_set(A,i,i,gsl_matrix_get(A,i,i)-1);
		for(j=0; j < SIZE; j++)
		{
			d += pow(gsl_matrix_get(A,i,j),2);
		}
	}
	printf("%g\n",sqrt(d));

	gsl_vector_free(y);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(AI);
	gsl_permutation_free(p);
	return 0;
}
