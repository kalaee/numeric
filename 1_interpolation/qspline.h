#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c;
} qspline;

// construct the quadratic spline
qspline* qspline_alloc(int n, double x[], double y[])
{
	assert(n>1);

	// allocating memory
	qspline* s = (qspline*) malloc(n*sizeof(qspline));
	s->n = n;
	s->x = (double*) malloc( n*sizeof(double) );
	s->y = (double*) malloc( n*sizeof(double) );
	s->b = (double*) malloc( (n-1)*sizeof(double) );
	s->c = (double*) malloc( (n-1)*sizeof(double) );
	
	// transfering x and y into s
	int i;
	for(i=0; i<n; i++)
	{
		s->x[i] = x[i];
		s->y[i] = y[i];
	}

	// estimating p[i]
	double p[n-1];
	for(i=0; i<n-1; i++)
	{
		p[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
	}

	// forward recursion with c[0] = 0
	s->c[0] = 0;
	for(i = 1; i<n-1; i++)
	{
		s->c[i] = (p[i]-p[i-1]-s->c[i-1]*(x[i]-x[i-1]))/(x[i+1]-x[i]);
	}

	// backward recursion from 1/2*c[n-2]
	s->c[n-2] *= 0.5;
	for(i=n-2; i>0; i--)
	{
		s->c[i-1] = (p[i]-p[i-1]-s->c[i]*(x[i+1]-x[i]))/(x[i]-x[i-1]);
	}

	// estimating b[i]
	for(i=0; i<n-1; i++)
	{
		s->b[i] = p[i] - s->c[i]*(x[i+1]-x[i]);
	}

	return s;
}

// evaluates the quadratic spline
double qspline_eval(qspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// binary search for the correct coordinate
	int i = 0, j = s->n, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > s->x[k])

			i = k;
		else
			j = k;
	}

	return s->y[i] + s->b[i]*(z - s->x[i]) + s->c[i]*pow(z - s->x[i],2);
}

// frees the memory
void qspline_free(qspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

