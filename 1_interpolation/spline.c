#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// Implementation of quadratic spline, "qspline"

// define the qspline structure
typedef struct
{
	int n;
	double x[], y[], b[], c[];
} qspline;

// creating the qspline structure
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
	
	// estimating p[i]
	int i;
	double p[n-1];
	for(i=0; i<n-1; i++)
	{
		p[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
	}

	// forward recursion with c[0] = 0
	for(i = 1, c[0]=0; i<n-1; i++)
	{
		c[i] = (p[i]-p[i-1]-c[i-1]*(x[i]-x[i-1]))/(x[i+1]-x[i]);
	}

	// backward recursion from 1/2*c[n-2]
	c[n-2] *= 0.5;
	for(i=n-2; i>0; i--)
	{
		c[i-1] = (p[i]-p[i-1]-c[i]*(x[i+1]-x[i]))/(x[i]-x[i-1]);
	}

	// estimating b[i]
	for(i=0; i<n-1; i++)
	{
		b[i] = p[i] - c[i]*(x[i+1]-x[i]);
	}

	return s;
}

double qspline_evaluate(qspline * s, double z)
{
	assert(s->x[0] <= z && z <= s->x[s->n-1]);
	
	// binary search for the correct coordinate
	int i = 0, j = n/2, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > x[k])
			i = m;
		else
			j = m;
	}

	return s->y[i] + s->b[i]*(z - s->x[i]) + s->c[i]*pow(z - s->x[i],2);
}

// function for freeing the memory after use
void qspline_free(qspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

// Implementation of linear spline, "lspline"
double lspline(int n, double x[], double y[], double z)
{
	assert(n > 1 && x[0]<= z && z <= x[n-1]);

	// binary search for the correct coordinate
	int i = 0, j = n/2, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > x[k])
			i = m;
		else
			j = m;
	}

	return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}
