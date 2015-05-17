#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#define ODE_SUCCES	0
#define ODE_FAIL	1
#define ODE_FAIL_UPPER	1000

typedef struct
{
	int n;	// number of function
	void* s;	// workspace og step function
	gsl_vector* err;
	gsl_vector* yh;
} ode_workspace;

ode_workspace* ode_workspace_alloc(int n, void* STEPPER_WORKSPACE_ALLOC(int n));
void ode_workspace_free(ode_workspace* W, void STEPPER_WORKSPACE_FREE(void* STEPPER_WORKSPACE));
int ode_driver(void f(double t, gsl_vector* y, gsl_vector* dy),
	void STEPPER(double t, double h, gsl_vector* y,
		void f(double t, gsl_vector* y, gsl_vector* dydt),
		gsl_vector* yh, gsl_vector* err, void* WORKSPACE),
	double* a, double b, double* h0, gsl_vector* ya, double acc, double eps, ode_workspace* W);

