#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "jacobi.h"

int main(int argc, char* argv[])
{
	assert(argc == 2);
	int SIZE = atoi(argv[1]);
	int i,j,rot;
	gsl_matrix* A = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* V = gsl_matrix_alloc(SIZE,SIZE);
	gsl_vector* e = gsl_vector_alloc(SIZE);

	// We construct matrix A and make it symmetric
	for (i = 0; i < SIZE; i++)
	{
		gsl_matrix_set(A,i,i,i*sin(i)+cos(i*i));
		for (j = 0; j < SIZE; j++)
		{
			gsl_matrix_set(A, i, j, 1.0/(i+j+1.0));
		}
	}

	// perform cyclic jacobi diagonalisation
	rot = jacobi_row(A,e,V,JACOBI_SORT_ASC);

	printf("%d\t%d\n",SIZE,rot);

	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(e);
}
