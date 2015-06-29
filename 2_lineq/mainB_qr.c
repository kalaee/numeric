#include <gsl/gsl_matrix.h>
#include "qr.h"
#define RND ((double)rand()/RAND_MAX)

int main(int argc, char* argv[])
{
	int i,j, size;
	size = atoi(argv[1]);
	gsl_matrix* A = gsl_matrix_alloc(size,size);
	gsl_matrix* R = gsl_matrix_alloc(size,size);
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			gsl_matrix_set (A, j, i, RND);
		}
	}
	qr_dec(A,R);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	return 0;
}
