#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "svd.h"

#define N_DATA	12
#define XI -1.0
#define XF 9.0
#define NOISE 0.3

double my_func(int i, double x)
{
	switch(i)
	{
		case 0: return cos(x); break;
		case 1: return sin(0.5*x); break;
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
	gsl_vector* dc = gsl_vector_alloc(4);
	gsl_matrix* S = gsl_matrix_alloc(4,4);
	int i, j, k;
	double xi, dx, yi, dyi, fjxi;
	double c_true[] = {5, 10, 2, 3};


	// generate some data with arbitrary "noise" and uncertanties
	dx = (XF - XI)/N_DATA;
	for(i=0; i<N_DATA; i++)
	{
		xi = XI + i*dx;
		gsl_vector_set(x,i,xi);
		yi = NOISE*sin(xi+i);
		for(j=0; j<4; j++)
		{
			yi += c_true[j]*my_func(j,xi);
		}
		gsl_vector_set(y,i,yi);
		gsl_vector_set(dy,i,0.6+0.3*sin(i*i));
	}

	// do the actual fitting
	svd_workspace* w = svd_workspace_alloc(N_DATA,4);
	svd_fit(x,y,dy,my_func,c,S,w);
	svd_workspace_free(w);

	// output datapoints for plotting
	for(i=0; i<N_DATA; i++)
	{
		fprintf(stderr,"%g\t%g\t%g\n",gsl_vector_get(x,i),
								gsl_vector_get(y,i),
								gsl_vector_get(dy,i));
	}

	// output fitting parameters to stdout
	printf("The fitting parameters are:\n");
	printf("i\tc_true\tc_fit\tdc\n");
	for(i=0; i<4; i++)
	{
		gsl_vector_set(dc,i,sqrt(gsl_matrix_get(S,i,i)));
		printf("%d\t%.4g\t%.4g\t%.4g\n",i,c_true[i],gsl_vector_get(c,i),gsl_vector_get(dc,i));
	}
	
	// output the fitted function values to file A_fit.dat
	// with an estimate of the uncertainty on y
	dx = (XF - XI)/1000;
	fprintf(stderr,"\n\n");
	for(i=0; i<1000; i++)
	{
		xi = XI + i*dx;
		yi = 0;
		dyi = 0;

		for(j=0; j<4; j++)
		{
			fjxi = my_func(j,xi);
			yi += gsl_vector_get(c,j)*fjxi;
			dyi += gsl_matrix_get(S,j,j)*fjxi*fjxi;
			for(k=0; k<j; k++)
			{
				dyi += 2*gsl_matrix_get(S,j,k)*fjxi*my_func(k,xi);
			}
		}
		dyi = sqrt(dyi);
		fprintf(stderr,"%g\t%g\t%g\t%g\n",xi,yi,yi+dyi,yi-dyi);
	}

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_vector_free(dc);
	gsl_matrix_free(S);
}
