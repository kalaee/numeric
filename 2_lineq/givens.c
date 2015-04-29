#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "backsub.h"

// QR decomp using Givens rotations
void givens_qr_dec(gsl_matrix* A)
{
	int i,j,k;
	double t, xi, xj;
	// zero each column at a time
	for(i=0; i < A->size2; i++)
	{
		for(j=i+1; j < A->size1; j++)
		{
			/* only perform rotation if actually needed
			 * will increase speed if matrix is sparse */
			if (gsl_matrix_get(A,j,i)!=0)
			{
				t = atan2(gsl_matrix_get(A,j,i),gsl_matrix_get(A,i,i));
				for(k = i; k < A->size2; k++)
				{
					xi = gsl_matrix_get(A,i,k);
					xj = gsl_matrix_get(A,j,k);
					gsl_matrix_set(A,i,k,xi*cos(t)+xj*sin(t));
					gsl_matrix_set(A,j,k,-xi*sin(t)+xj*cos(t));
				}
				gsl_matrix_set(A,j,i,t);
			}
		}
	}
}

// Solve the linear system QRx=b and save the resulting vector in res
void givens_qr_bak(gsl_matrix* QR, gsl_vector* b, gsl_vector* res)
{
	// Assert that b and res has the correct dimensions
	assert(res->size == QR->size2 && b->size == QR->size1);
	int i,j;
	double t, vi, vj;
	// transform QRx=b -> Rx = Q^T b
	for(i = 0; i < QR->size2; i++)
	{
		for(j=i+1; j<QR->size1; j++)
		{
			t = gsl_matrix_get(QR,j,i);
			vi = gsl_vector_get(b,i);
			vj = gsl_vector_get(b,j);
			gsl_vector_set(b,i,vi*cos(t)+vj*sin(t));
			gsl_vector_set(b,j,-vi*sin(t)+vj*cos(t));
		}
	}
	// solve Rx = Q^T b via backsubstitution and save the result in res
	gsl_vector_view x = gsl_vector_subvector(b,0,QR->size2);
	gsl_matrix_view U = gsl_matrix_submatrix(QR,0,0,QR->size2,QR->size2);
	backsub_upper(&U.matrix,&x.vector);
	gsl_vector_memcpy(res,&x.vector);
}

// determine the determinant of QR using the products of R_{ii}
double givens_qr_det(gsl_matrix* QR)
{
	assert(QR->size1 == QR->size2);
	int i;
	double d=1;
	for(i=0; i < QR->size1; i++)
	{
		d *= gsl_matrix_get(QR,i,i);
	}
	return d;
}

// inverse matrix of A using QR decomp and the workspace vector w of length AI->size1
void givens_qr_inv(gsl_matrix* QR, gsl_matrix* AI, gsl_vector* w)
{
	// assert that AI and w has the corect dimensions
	assert(AI->size1 == AI->size2 && AI->size1 == w->size && AI->size1 == QR->size1);
	int i;
	// solve the systems QRx_i = e_i where x_i are columns of the inverse matrix
	for(i = 0; i < AI->size1; i++)
	{
		gsl_vector_view xi = gsl_matrix_column(AI,i);
		gsl_vector_set_basis(w,i);
		givens_qr_bak(QR,w,&xi.vector);
	}
	return;
}
