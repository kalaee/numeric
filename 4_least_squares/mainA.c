#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "leasqr.h"

#define N_DATA	20
#define XI -1.0
#define XF 9.0
#define NOISE 0.3

double my_func(int i, double x)
{
	switch(i)
	{
		case 0: return cos(x); break;
		case 1: return sin(2*x); break;
		case 2: return x; break;
		case 3: return exp(-x*x); break;
		default: {fprintf(stderr,"my_func: wrong i: %d\n",i); return NAN;}
	}
}

int main(void)
{

	gsl_vector* x = gsl_vector_alloc(N_DATA);
	gsl_vector* y = gsl_vector_calloc(N_DATA);
	gsl_vector* dy = gsl_vector_alloc(N_DATA);
	gsl_vector* c = gsl_vector_alloc(4);
	gsl_matrix* S = gsl_matrix_alloc(4,4);
	int i, j;
	double xi, dx, yi;


	// generate som data with arbitrary "noise" and uncertanties
	dx = (XF - XI)/N_DATA;
	for(i=0; i<N_DATA; i++)
	{
		xi = XI + i*dx;
		gsl_vector_set(x,i,xi);
		yi = NOISE*sin(xi+i);
		for(j=0; j<4; j++)
		{
			yi += my_func(j,xi);
		}
		gsl_vector_set(y,i,yi);
		gsl_vector_set(dy,i,1+0.2*sin(i));
	}

	// do the actual fitting
	ls_workspace w = ls_workspace_alloc(N_DATA,4);
	lsfit(x,y,dy,func,c,S,w)
	ls_workspace_free(w);

	

}
