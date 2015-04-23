#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c, *d;
} cspline;
int cspline_binsearch(cspline * s, double z);
cspline* cspline_alloc(int n, double x[], double y[], double dds0, double ddsn);
double cspline_eval(cspline * s, double z);
void cspline_free(cspline * s);
double cspline_derivative(cspline * s, double z);
double cspline_integral(cspline *s, double z);

