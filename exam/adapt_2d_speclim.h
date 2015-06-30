#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ADAPT_NRECUR_MAX 100000

double adapt24_inner_speclim(double f(double x, double y), double x, double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur);
double adapt24_outer_speclim(double f(double x, double y), double ax, double bx, double d(double x), double u(double x), double y2, double y3, double acc, double eps, double *err, int nrecur);
double adapt_2d_speclim(double f(double x, double y), double ax, double bx, double d(double x), double u (double x), double acc, double eps, double* err);
