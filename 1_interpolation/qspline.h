#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c, *h;
} qspline;
int qspline_binsearch(qspline * s, double z);
qspline* qspline_alloc(int n, double x[], double y[]);
double qspline_eval(qspline * s, double z);
void qspline_free(qspline * s);
double qspline_derivative(qspline * s, double z);
double qspline_integral(qspline *s, double z);

