#include <math.h>
#include <assert.h>
#include <stdio.h>
// Quadratic Adaptive integral using Recursive functions and Open intervals of order 4 with 2 for error estimate
double qaro24(double f(double x), double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur);
double qaro(double f(double), double a, double b, double acc, double eps, double * err);

