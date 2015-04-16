#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c;
} qspline;
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
