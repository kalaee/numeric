#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "../2_lineq/backsub.h"
#include "../2_lineq/qr.h"

// define the workspace necessary for a linear least square fit
typedef struct
{
	int n, m;
	gsl_matrix* A;
	gsl_matrix* R;
	gsl_matrix* Rinv;
	gsl_vector* b;
} ls_workspace;

// allocate the memory for least square fitting on system with n_data points and nf functions
ls_workspace* ls_workspace_alloc(int n_data, int nf)
{
	ls_workspace* W = (ls_workspace*) malloc(sizeof(ls_workspace));
	W->n = n_data;
	W->m = nf;
	W->A = gsl_matrix_alloc(n_data,nf);
	W->R = gsl_matrix_alloc(nf,nf);
	W->Rinv = gsl_matrix_alloc(nf,nf);
	W->b = gsl_vector_alloc(n_data);
	return W;
}

// free the allocated memory for least square fitting
void ls_workspace_free(ls_workspace* W)
{
	gsl_matrix_free(W->A);
	gsl_matrix_free(W->R);
	gsl_matrix_free(W->Rinv);
	gsl_vector_free(W->b);
	free(W);
	return;
}

// fitting routine
void ls_fit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double func(int i, double z), gsl_vector* c, gsl_matrix* S, ls_workspace* W)
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

	// make Gram-Schmidt QR decomp and solve the linear system
	qr_dec(W->A,W->R);
	qr_bak(W->A,W->R,W->b,c);
	// calculate the inverse of triangular matrix R from which
	// the covariance matrix S is estimated
	backsub_upper_inv(W->R,W->Rinv);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,W->Rinv,W->Rinv,0,S);
	return;
}
