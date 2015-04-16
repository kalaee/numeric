#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// implementation of linear spline, "lspline"
double lspline(int n, double x[], double y[], double z)
{
	assert(n > 1 && x[0]<= z && z <= x[n-1]);

	// binary search for the correct coordinate
	int i = 0, j = n, k;
	while(j-i > 1)
	{
		k = (i+j)/2;
		if (z > x[k])
			i = k;
		else
			j = k;
	}

	return y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
}

// Implementation of quadratic spline, "qspline"
// define the qspline structure
typedef struct
{
	int n;
	double *x, *y, *b, *c;
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

//  freeing the memory after use
void qspline_free(qspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

int main(void)
{
	// generating source data
	int length = 10, i;
	double x[10], y[10];
	FILE * f = fopen("A.dat","w");
	for(i=0; i<10; i++)
	{
		x[i] = i+0.5 * sin(i);
		y[i] = i + cos(i*i);
		fprintf(f,"%g\t%g\n",x[i],y[i]);
	}

	// applying linear spline
	int steps = 500;
	double dz = (x[length-1]-x[0])/steps, z[steps];
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		z[i] = x[0]+i*dz;
		fprintf(f,"%g\t%g\n",z[i],lspline(length,x,y,z[i]));
	}
	
	// applying quadratic spline
	qspline * s = qspline_alloc(length, x, y);
	fprintf(f,"\n\n");
	for(i=0; i<steps; i++)
	{
		fprintf(f,"%g\t%g\n",z[i],qspline_eval(s,z[i]));
	}
	fclose(f);
	qspline_free(s);

	return 0;
}



















