#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "backsub.h"
#include "qr.h"

typedef struct
{
	int n, m;
	gsl_matrix* A, R, Rinv;
	gsl_vector* b;
} ls_workspace;

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

void ls_workspace_free(ls_workspace* W)
{
	gsl_matrix_free(W->A);
	gsl_matrix_free(W->R);
	gsl_matrix_free(W->Rinv);
	gsl_vector_free(W->b);
	free(W);
	return;
}

void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double func(int i, double z),
	gsl_vector* c, gsl_matrix* S, ls_workspace* W)
{
	int i,j;
	double xi, dyi;
	gsl_vector_memcpy(W->b,y);
	gsl_vector_div(W->b,dy);
	for(i = 0; i < W->n; i++)
	{
		xi = gsl_vector_get(x,i);
		dyi = gsl_vector_get(dy,i);
		for(j = 0; j < W->m; j++)
		{
			gsl_matrix_set(A,i,j,f(j,xi)/dyi);
		}
	}
	qr_dec(W->A,W->R);
	backsub_upper_inv(W->R,W->Rinv);
	qr_bak(W->A,W->R,W->b,c);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,W->Rinv,W->Rinv,0,S);
	return;
}
