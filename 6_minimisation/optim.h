#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
	gsl_vector* df;
	gsl_vector* Dx;
	gsl_vector* z;
	gsl_matrix* H;
} min_newton_workspace;

min_newton_workspace* min_newton_workspace_alloc(int n);
void min_newton_workspace_free(min_newton_workspace* W);
int min_newton(double f(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* x, double alpha, double tol, min_newton_workspace* W);
