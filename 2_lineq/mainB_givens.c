#include <gsl/gsl_matrix.h>
#include "givens.h"

#define RND ((double)rand()/RAND_MAX)

int main(int argc, char* argv[])
{
	int i,j, size;
	size = atoi(argv[1]);
	gsl_matrix* A = gsl_matrix_alloc(size,size);
	for (i = 0; i < size; i++)
	{
		for (j = 0; j < size; j++)
		{
			gsl_matrix_set (A, j, i, RND);
		}
	}
	givens_qr_dec(A);
	gsl_matrix_free(A);
	return 0;
}
