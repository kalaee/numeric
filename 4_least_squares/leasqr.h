#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
	int n, m;
	gsl_matrix* A, R, Rinv;
	gsl_vector* b;
} ls_workspace;

ls_workspace* ls_workspace_alloc(int n_data, int nf);

void ls_workspace_free(ls_workspace* W);

void lsfit(gsl_vector* x, gsl_vector* y, gsl_vector* dy, double func(int i, double z),
	gsl_vector* c, gsl_matrix* S, ls_workspace* W);
