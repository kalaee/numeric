#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
	gsl_vector* df;
	gsl_vector* df0;
	gsl_vector* Dx;
	gsl_vector* z;
	gsl_vector* y;
	gsl_matrix* H;
} min_newton_sr1_workspace;
min_newton_sr1_workspace* min_newton_sr1_workspace_alloc(int n);
void min_newton_sr1_workspace_free(min_newton_sr1_workspace* W);
int min_newton_sr1(double f(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), gsl_vector* x, double alpha, double tol, min_newton_sr1_workspace* W);
