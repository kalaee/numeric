#include <gsl/gsl_vector.h>
#include <math.h>

typedef struct
{
	int n;
	gsl_vector* k0;
	gsl_vector* k1;
	gsl_vector* k2;
	gsl_vector* k3;
	gsl_vector* k4;
	gsl_vector* k5;
	gsl_vector* yh;
} ode_rk5_workspace;

void* ODE_RK5_ALLOC(int n);


void ODE_RK5_FREE(void* WORKSPACE);


void ODE_RK5(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, void* WORKSPACE);

