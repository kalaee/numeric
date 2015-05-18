#include <gsl/gsl_vector.h>
#include <math.h>
#include "../libs/vector.h"

typedef struct
{
	int n;
	gsl_vector* k0;
	gsl_vector* k1;
	gsl_vector* k2;
	gsl_vector* k3;
	gsl_vector* k4;
	gsl_vector* k5;
	gsl_vector* yh;
} ode_rk5_workspace;

void* ODE_RK5_ALLOC(int n)
{
	ode_rk5_workspace* W = (ode_rk5_workspace*) malloc(sizeof(ode_rk5_workspace));
	W->n = n;
	W->k0 = gsl_vector_alloc(n);
	W->k1 = gsl_vector_alloc(n);
	W->k2 = gsl_vector_alloc(n);
	W->k3 = gsl_vector_alloc(n);
	W->k4 = gsl_vector_alloc(n);
	W->k5 = gsl_vector_alloc(n);
	W->yh = gsl_vector_alloc(n);
	return (void*) W;
}

void ODE_RK5_FREE(void* WORKSPACE)
{
	ode_rk5_workspace* W = (ode_rk5_workspace*) WORKSPACE;
	gsl_vector_free(W->k0);
	gsl_vector_free(W->k1);
	gsl_vector_free(W->k2);
	gsl_vector_free(W->k3);
	gsl_vector_free(W->k4);
	gsl_vector_free(W->k5);
	gsl_vector_free(W->yh);
	free(W);
	return;
}

void ODE_RK5(double t, double h, gsl_vector* y, void f(double t, gsl_vector* y, gsl_vector* dydt), gsl_vector* yh, gsl_vector* err, void* WORKSPACE)
{
	ode_rk5_workspace* W = (ode_rk5_workspace*) WORKSPACE;
	double yi, k0i, k1i, k2i, k3i, k4i, k5i;
	int i, k;
	gsl_vector_memcpy(err,y);
	for (k=0; k<3; k++)
	{
		if (k == 1)
		{
			gsl_vector_memcpy(W->yh,yh);
			h /= 2;
		}
		else if (k == 2)
		{
			gsl_vector_memcpy(err,yh);
			t += h;
		}
		// estimate k0
		f(t,err,W->k0);
		gsl_vector_scale(W->k0,h);
		// estimate k1
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi+0.25*k0i);
		}
		f(t+0.25*h,yh,W->k1);
		gsl_vector_scale(W->k1,h);
		// estimate k2
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			k1i = gsl_vector_get(W->k1,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi+3./32.*k0i+9./32.*k1i);
		}
		f(t+3./8.*h,yh,W->k2);
		gsl_vector_scale(W->k2,h);
		// estimate k3
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			k1i = gsl_vector_get(W->k1,i);
			k2i = gsl_vector_get(W->k2,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi+1932./2197.*k0i-7200./2197.*k1i+7296./2197.*k2i);
		}
		f(t+12./13.*h,yh,W->k3);
		gsl_vector_scale(W->k3,h);
		// estimate k4
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			k1i = gsl_vector_get(W->k1,i);
			k2i = gsl_vector_get(W->k2,i);
			k3i = gsl_vector_get(W->k3,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi+439./216.*k0i-8.*k1i+3680./513.*k2i-845./4104.*k3i);
		}
		f(t+h,yh,W->k4);
		gsl_vector_scale(W->k4,h);
		// estimate k5
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			k1i = gsl_vector_get(W->k1,i);
			k2i = gsl_vector_get(W->k2,i);
			k3i = gsl_vector_get(W->k3,i);
			k4i = gsl_vector_get(W->k4,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi-8./27.*k0i+2.*k1i-3544./2565.*k2i+1859./4104.*k3i-11./40.*k4i);
		}
		f(t+0.5*h,yh,W->k5);
		gsl_vector_scale(W->k5,h);
		// estimate y(t+h) and error, i.e. difference between 4th and 5th order
		for(i=0; i<W->n; i++)
		{
			k0i = gsl_vector_get(W->k0,i);
			k1i = gsl_vector_get(W->k1,i);
			k2i = gsl_vector_get(W->k2,i);
			k3i = gsl_vector_get(W->k3,i);
			k4i = gsl_vector_get(W->k4,i);
			k5i = gsl_vector_get(W->k5,i);
			yi = gsl_vector_get(err,i);
			gsl_vector_set(yh,i,yi+16./135.*k0i+6656./12825.*k2i+28561./56430.*k3i-9./50.*k4i+2./55.*k5i);
		}
	}
	vector_sum(1,W->yh,-1,yh,err);
	gsl_vector_scale(err,1./(pow(2,5)-1));
	return;
}
