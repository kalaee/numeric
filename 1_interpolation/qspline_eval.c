#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c;
} qspline;

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
