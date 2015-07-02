#define ADAPT_NRECUR_MAX 1000
double adapt_nd_speclim_last(double f(double* x), double a, double b, int dim, double* coord, double y2, double y3, double acc, double eps, double *err, int nrecur);
double adapt_nd_speclim_now(double f(double* x), double a, double b, double d(double* x, int p), double u(double* x, int p), int dim, int lvl, double* coord, double y2, double y3, double acc, double eps, double*err, int nrecur);
double adapt_nd_speclim_next(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, int lvl, double* coord, double acc, double eps, double* err);
double adapt_nd_speclim(double f(double* x), double d(double* x, int p), double u (double* x, int p), int dim, double acc, double eps, double *err);

