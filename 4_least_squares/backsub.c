#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <assert.h>

// solve Ux = b inplace, where U is an upper triangular matrix
void backsub_upper(gsl_matrix* U, gsl_vector* b)
{
	assert(U->size1 == U->size2 && U->size1 == b->size);
	int i,j;
	double s;
	for(i = b->size-1; i >= 0; i--)
	{
		s = gsl_vector_get(b,i);
		for(j=i+1; j < b->size; j++)
		{
			s -= gsl_matrix_get(U,i,j)*gsl_vector_get(b,j);
		}
		gsl_vector_set(b,i,s/gsl_matrix_get(U,i,i));
	}
	return;
}

// solve Lx = b inplace, where L is an lower triangular matrix
void backsub_lower(gsl_matrix* L, gsl_vector* b)
{
	assert(L->size1 == L->size2 && L->size1 == b->size);
	int i,j;
	double s;
	for(i = 0; i<b->size; i++)
	{
		s = gsl_vector_get(b,i);
		for(j=0; j < i; j++)
		{
			s -= gsl_matrix_get(L,i,j)*gsl_vector_get(b,j);
		}
		gsl_vector_set(b,i,s/gsl_matrix_get(L,i,i));
	}
	return;
}

// find the inverse of upper triangular matrix U using backsubstitution
void backsub_upper_inv(gsl_matrix* U, gsl_matrix* Uinv)
{
	// assert that U and Uinv are correctly dimensioned
	assert(U->size1 == U->size2 && Uinv->size1 == Uinv->size2 && U->size1 == Uinv->size1);
	// make columns in Uinv the cartesian basis
	gsl_matrix_set_identity(Uinv);
	// perform backsubstitution for each column
	gsl_vector_view u;
	int i;
	for(i=0; i<U->size1; i++)
	{
		u = gsl_column(Uinv,i);
		backsub_upper(U,&u.vector);
	}
	return;
}
