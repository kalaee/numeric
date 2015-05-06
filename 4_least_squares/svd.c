#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "jacobi.h"

// define the workspace necessary for a linear least square fit
typedef struct
{
	int n, m;
	gsl_matrix* A;
	gsl_matrix* U;
	gsl_matrix* V;
	gsl_vector* b;
} svd_workspace;

// allocate the memory for least square fitting on system with n_data points and nf functions
svd_workspace* svd_workspace_alloc(int n_data, int nf)
{
	svd_workspace* W = (svd_workspace*) malloc(sizeof(svd_workspace));
	W->n = n_data;
	W->m = nf;
	W->A = gsl_matrix_alloc(n_data,nf);
	W->U = gsl_matrix_alloc(n_data,nf);
	W->V = gsl_matrix_alloc(nf,nf);
	W->b = gsl_vector_alloc(n_data);
	return W;
}

void svd_workspace_free(svd_workspace* W)
{
	gsl_matrix_free(W->A);
	gsl_matrix_free(W->U);
	gsl_matrix_free(W->V);
	gsl_vector_free(W->b);
	free(W);
}
// make singular value decomposition of matrix A in preparation of least squares fitc
void svd_prep(gsl_vector* c, gsl_matrix* S, svd_workspace* W)
{
	// calculate A^TA and perform jacobi eigenvalue diagonalisation using cyclic sweeps
	// eigenvalues are stored in c, eigenvectors in W->V
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,W->A,W->A,0,S);
	jacobi_cyclic(S,c,W->V);

	// calculate U = AVD^{1/2} stored in W->U
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,W->A,W->V,0,W->U);
	int i;
	gsl_vector_view v;
	for(i=0; i<W->m; i++)
	{
		v = gsl_matrix_column(W->U,i);
		gsl_vector_scale(&v.vector,pow(gsl_vector_get(c,i),-0.5));
	}

	// calculate the covariance matrix S = VD^{-1}V^T
	// we begin with creating D^{-1}V^T and store it in W->A
	for(i=0; i<W->m; i++)
	{
		v = gsl_matrix_row(W->A,i);
		gsl_matrix_get_col(&v.vector,W->V,i);
		gsl_vector_scale(&v.vector,1/gsl_vector_get(c,i));
	}
	gsl_matrix_view a = gsl_matrix_submatrix(W->A,0,0,W->m,W->m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,W->V,&a.matrix,0,S);
	return;
}

void svd_sol(gsl_vector* c, svd_workspace* W)
{
	// calculate c_fit = VD^{-1/2}U^Tb, where diag of D is saved in vector c
	// after calculation the fitting parameters will be stored in c	
	gsl_vector_view u = gsl_matrix_row(W->A,0);
	gsl_blas_dgemv(CblasTrans,1,W->U,W->b,0,&u.vector);
	int i;
	double val;
	for(i=0; i<W->m; i++)
	{
		val = gsl_vector_get(&u.vector,i)/sqrt(gsl_vector_get(c,i));
		gsl_vector_set(&u.vector,i,val);
	}
	gsl_blas_dgemv(CblasNoTrans,1,W->V,&u.vector,0,c);
	return;
}

// least square fit using singular value decomposition
void svd_fit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double func(int i, double z),
		gsl_vector* c, gsl_matrix* S, svd_workspace* W)
{
	int i,j;
	double xi, dyi;

	// generate vector with weighted y-values
	gsl_vector_memcpy(W->b,y);
	gsl_vector_div(W->b,dy);

	// generate matrix of weighted coefficients for a given xi
	for(i = 0; i < W->n; i++)
	{
		xi = gsl_vector_get(x,i);
		dyi = gsl_vector_get(dy,i);
		for(j = 0; j < W->m; j++)
		{
			gsl_matrix_set(W->A,i,j,func(j,xi)/dyi);
		}
	}

	// perform singular value decomposition
	svd_prep(c,S,W);

	// solve Ac = b from using singular value decomp
	svd_sol(c,W);

	return;
}
