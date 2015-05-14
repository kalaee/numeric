#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
	int n;
	gsl_matrix* J;
	gsl_vector* z;
	gsl_vector* Dx;
	gsl_vector* fx;
	gsl_vector* fz;
} roots_newton_workspace;
roots_newton_workspace* roots_newton_workspace_alloc(int n);
void roots_newton_workspace_free(roots_newton_workspace* W);
int roots_newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, roots_newton_workspace* W);
int roots_newton_derivative(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, roots_newton_workspace* W);
int roots_newton_interp(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, roots_newton_workspace* W);
