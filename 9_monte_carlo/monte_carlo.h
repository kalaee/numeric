#include <math.h>
#include <gsl/gsl_rng.h>

void monte_carlo_sample_vector(double* x, double* a, double* h, int dim, gsl_rng* r);
double monte_carlo_plain(double f(double* t), double* a, double* b, int dim, int N, double *err, gsl_rng* r);
double monte_carlo_tolerance(double f(double* t), double* a, double* b, int dim, double acc, double eps, int THRESHOLD, double *err, gsl_rng* r);
