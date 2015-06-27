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
} ode_rkf45_workspace;

void* ODE_RKF45_ALLOC(int n);
void ODE_RKF45_FREE(void* WORKSPACE);
void ODE_RKF45(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, void* WORKSPACE);

