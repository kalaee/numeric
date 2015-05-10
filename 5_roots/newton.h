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
} newton_workspace;
newton_workspace* newton_workspace_alloc(int n);
void newton_workspace_free(newton_workspace* W);
int newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, newton_workspace* W);
int newton_derivative(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, newton_workspace* W);
int newton_derivative_interp(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, newton_workspace* W);
