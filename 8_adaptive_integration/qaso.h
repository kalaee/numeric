#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QASO_MAX_DEPTH	100000

double qaso23(double f(double t), double a, double b, double acc, double eps, double* err);
double qaso(double f(double), double a, double b, double acc, double eps, double * err);

