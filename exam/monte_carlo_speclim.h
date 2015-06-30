double monte_carlo_speclim(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, double ACC, double EPS, double THRESHOLD, double *err, gsl_rng* r);
