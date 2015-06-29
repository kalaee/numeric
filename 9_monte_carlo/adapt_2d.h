#include <stdio.h>
#include <math.h>

#define ADAPT_NRECUR_MAX 10000

double adapt24_inner(double f(double x, double y), double x, double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur);
double adapt24_outer(double f(double x, double y), double ax, double ay, double bx, double by, double y2, double y3, double acc, double eps, double *err, int nrecur);
double adapt_2d(double f(double x, double y), double ax, double ay, double bx, double by, double acc, double eps, double* err);
