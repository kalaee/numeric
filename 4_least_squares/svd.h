#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "jacobi.h"

typedef struct
{
	int n, m;
	gsl_matrix* A;
	gsl_matrix* U;
	gsl_matrix* V;
	gsl_vector* b;
} svd_workspace;

svd_workspace* svd_workspace_alloc(int n_data, int nf);

void svd_workspace_free(svd_workspace* W);
void svd_prep(gsl_vector* c, gsl_matrix* S, svd_workspace* W);
void svd_sol(gsl_vector* c, svd_workspace* W);
void svd_fit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double func(int i, double z),gsl_vector* c, gsl_matrix* S, svd_workspace* W);
