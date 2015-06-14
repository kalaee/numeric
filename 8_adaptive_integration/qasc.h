#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QASC_MAX_DEPTH	100000

double qasc23(double f(double t), double a, double b, double acc, double eps, double *err);
double qasc(double f(double), double a, double b, double acc, double eps, double * err);

