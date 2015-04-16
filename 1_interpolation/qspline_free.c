#include <stdio.h>
#include <stdlib.h>

typedef struct
{
	int n;
	double *x, *y, *b, *c;
} qspline;

void qspline_free(qspline * s)
{
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

