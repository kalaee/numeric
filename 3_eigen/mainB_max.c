#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "jacobi.h"

int main(int argc, char * argv[])
{
	int i,j,rot, SIZE;
	assert(argc == 2);
	SIZE = atoi(argv[1]);
	gsl_matrix* A = gsl_matrix_alloc(SIZE,SIZE);
	gsl_matrix* V = gsl_matrix_alloc(SIZE,SIZE);
	gsl_vector* e = gsl_vector_alloc(SIZE);

	// We construct matrix A and make it symmetric
	for (i = 0; i < SIZE; i++)
	{
		gsl_matrix_set(A,i,i,sin(i)+cos(i*i));
		for (j = i+1; j < SIZE; j++)
		{
			gsl_matrix_set(A, i, j, sin(i) + cos(j*i));
		}
	}

	// perform cyclic jacobi diagonalisation
	rot = jacobi_max_row(A,e,V,JACOBI_SORT_DESC);

	fprintf(stdout,"Size %d,\trot: %d\n",SIZE,rot);
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(e);
}
