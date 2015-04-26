#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>
#include "backsub.h"

// QR decomposition using Gram-Schmidt orthogonalization
void qr_dec(gsl_matrix* A, gsl_matrix* R)
{
	// assert that R has the correct dimensions
	assert(R->size1 == R->size2 && R->size1 == A->size2);

	int i,j;
	double r,s;
	for(i = 0; i < A->size2; i++)
	{
		// normalize the i'th column
		gsl_vector_view ai = gsl_matrix_column(A,i);
		r = gsl_blas_dnrm2(&ai.vector);
		gsl_matrix_set(R,i,i,r);
		gsl_vector_scale(&ai.vector,1/r);
		// orthogonalize the remaining columns
		for(j = i+1; j < A->size2; j++)
		{
			gsl_vector_view aj = gsl_matrix_column(A,j);
			gsl_blas_ddot(&ai.vector,&aj.vector,&s);
			gsl_blas_daxpy(-s,&ai.vector,&aj.vector);
			gsl_matrix_set(R,i,j,s);
		}
	}
	return;
}

// Find the solution to the system Ax = b using QR decomp, QRx = b
void qr_bak(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* res)
{
	/* assert that b and res are of correct dimensions
	 * it is assumed that Q and R are correctly dimensioned
	 * from the QR decomposition
	 */
	assert(res->size == R->size1 && b->size == Q->size1);
	// reducing to the system Rx = Q^Tb
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,res);
	// backsubstituion
	backsub_upper(R,res);
}

// compute the absolute value of the determinant
double qr_absdet(gsl_matrix* R)
{
	int i;
	double d = 1;
	for(i=0; i < R->size1; i++)
	{
		d *= gsl_matrix_get(R,i,i);
	}
	return fabs(d);
}

// inverse matrix of A using QR decomp and the workspace vector w of length AI->size1
void qr_inv(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* AI, gsl_vector* w)
{
	/* assert that AI and w has the corect dimensions,
	 * it is assumed that Q and R are correctly dimensioned
	 * from the QR decomposition
	 */
	assert(AI->size1 == AI->size2 && AI->size1 == w->size && AI->size1 == R->size1);
	int i;
	// solve the systems QRx_i = e_i, x_i are columns of the inverse matrix
	for(i = 0; i < AI->size1; i++)
	{
		gsl_vector_view xi = gsl_matrix_column(AI,i);
		gsl_vector_set_basis(w,i);
		qr_bak(Q,R,w,&xi.vector);
	}
	
	return;
}


