#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "backsub.h"

void chol_dec(gsl_matrix* A)
{
	assert(A->size1 == A->size2);
	int i,j,k;
	double l;

	// generate the lower part, L inplace
	for(i=0; i<A->size1; i++)
	{
		// compute the off-diagonal elements
		for(j=0; j<i; j++)
		{
			l = gsl_matrix_get(A,i,j);
			for(k=0; k<j; k++)
			{
				l -= gsl_matrix_get(A,i,k)*gsl_matrix_get(A,j,k);
			}
			gsl_matrix_set(A,i,j,l/gsl_matrix_get(A,j,j));
		}
		// compute the diagonal elements
		l = gsl_matrix_get(A,i,i);
		for(k=0; k<i; k++)
		{
			l -= pow(gsl_matrix_get(A,i,k),2);
		}
		assert(l >= 0);
		gsl_matrix_set(A,i,i,sqrt(l));
	}

	/* insert the upper part, L^T, in the matrix A
	 * note that the diagonal elements are identical */
	for(i=0; i<A->size1; i++)
	{
		for(j=i+1; j<A->size1; j++)
		{
			gsl_matrix_set(A,i,j,gsl_matrix_get(A,j,i));
		}
	}
	
	return;
}

void chol_bak(gsl_matrix* LL, gsl_vector* b)
{
	/* assert that b is of correct size
	 * it is assumed that LL is square from Cholesky decomp */
	assert(b->size == LL->size1);
	/* backsubstituion from LL^T*x = b to L^T*x = Linv*b
	 * and from L^T*x = Linv*b to x = L^Tinv*Linv*b
	 * solution is stored in b */
	backsub_lower(LL,b);
	backsub_upper(LL,b);
}

double chol_det(gsl_matrix* LL)
{
	assert(LL->size1 == LL->size2);
	int i;
	double d=1;
	for(i=0; i<LL->size1; i++)
	{
		d *= gsl_matrix_get(LL,i,i);
	}
	return d*d;
}

// inverse matrix of A using Cholesky decomp
void chol_inv(gsl_matrix* LL, gsl_matrix* AI)
{
	// assert that AI has the corect dimensions
	assert(AI->size1 == AI->size2 && AI->size1 == LL->size1);
	int i;
	// solve the systems LL^Tx_i = e_i, x_i are columns of the inverse matrix
	for(i = 0; i < AI->size1; i++)
	{
		gsl_vector_view xi = gsl_matrix_column(AI,i);
		gsl_vector_set_basis(&xi.vector,i);
		chol_bak(LL,&xi.vector);
	}
	return;
}
