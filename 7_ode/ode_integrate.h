#include <gsl/gsl_vector.h>
#include "ode_driver.h"
#include "ode_rkf45.h"

double ode_integrate(void f(double t, gsl_vector* y, gsl_vector* dy), double a, double b, double ACC, double EPS, ode_workspace* W);
