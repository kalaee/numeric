#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "leasqr.h"

double f(int i, double x)
{
	switch(i)
	{
		case 0: return 1; break;
		case 1: return x; break;
		default: {fprintf(stderr,"f: wrong i:%d",i); return NAN;}
	}
}

int main(void)
{
	int i;
	double xi, yi, dyi;
	gsl_vector* x = gsl_vector_alloc(2);
	gsl_vector* y = gsl_vector_alloc(2);
	gsl_vector* c = gsl_vector_alloc(2);
	gsl_vector* dy = gsl_vector_alloc(2);
	gsl_vector* dc = gsl_vector_alloc(2);
	gsl_matrix* S = gsl_matrix_alloc(2,2);
	// set two data points at (0,3) and (1,5) with unc. 0.3 and 0.1
	gsl_vector_set(x,0,0);
	gsl_vector_set(x,1,1);
	gsl_vector_set(y,0,3);
	gsl_vector_set(y,1,5);
	gsl_vector_set(dy,0,0.3);
	gsl_vector_set(dy,1,0.1);
	// do least square fitting
	ls_workspace * W = ls_workspace_alloc(2,2);
	ls_fit(x,y,dy,f,c,S,W);
	ls_workspace_free(W);
	gsl_vector_set(dc,0,gsl_matrix_get(S,0,0));
	gsl_vector_set(dc,1,gsl_matrix_get(S,1,1));
	// print fitting parameters to stdout
	fprintf(stdout,"The fitting parameters to the function y=c0+c1*x is\n");
	fprintf(stdout,"i\tc_i\tdc_i\n");
	fprintf(stdout,"%d\t%g\t%g\n",0,gsl_vector_get(c,0),sqrt(gsl_vector_get(dc,0)));
	fprintf(stdout,"%d\t%g\t%g\n",1,gsl_vector_get(c,1),sqrt(gsl_vector_get(dc,1)));
	// print data to stderr
	fprintf(stderr,"%d\t%d\t%g\n%d\t%d\t%g\n\n\n",0,3,0.3,1,5,0.1);
	for(i=0; i<500; i++)
	{
		xi =-1+i*0.01;
		yi = gsl_vector_get(c,0)+gsl_vector_get(c,1)*xi;
		dyi = sqrt(gsl_vector_get(dc,0)+gsl_vector_get(dc,1)*xi*xi+2*gsl_matrix_get(S,0,1)*xi);
		fprintf(stderr,"%g\t%g\t%g\t%g\n",xi,yi,yi+dyi,yi-dyi);
	}
	// free allocated memory
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(c);
	gsl_vector_free(dc);
	gsl_matrix_free(S);
	return 0;
}
