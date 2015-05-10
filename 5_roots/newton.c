#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "../2_lineq/givens.h"

typedef struct
{
	int n;
	gsl_matrix* J;
	gsl_vector* z;
	gsl_vector* Dx;
	gsl_vector* fx;
	gsl_vector* fz;
} newton_workspace;

newton_workspace* newton_workspace_alloc(int n)
{
	newton_workspace* W = (newton_workspace*) malloc(sizeof(newton_workspace));
	W->n = n;
	W->J = gsl_matrix_alloc(n,n);
	W->z = gsl_vector_alloc(n);
	W->Dx = gsl_vector_alloc(n);
	W->fx = gsl_vector_alloc(n);
	W->fz = gsl_vector_alloc(n);

	return W;
}

void newton_workspace_free(newton_workspace* W)
{
	gsl_matrix_free(W->J);
	gsl_vector_free(W->z);
	gsl_vector_free(W->Dx);
	gsl_vector_free(W->fx);
	gsl_vector_free(W->fz);
	free(W);
}

void vector_add(double alpha, gsl_vector* a, double beta, gsl_vector* b, gsl_vector* c)
{
	int i;
	for(i=0; i<a->size; i++)
	{
		gsl_vector_set(c,i,alpha*gsl_vector_get(a,i)+beta*gsl_vector_get(b,i));
	}
	return;
}

int newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, newton_workspace* W)
{
	int n = W->n;
	int i,j, counter = 0;
	double lambda, normfx;
	do
	{
		counter++;
		f(x,W->fx);
		for(j = 0; j < n; j++)
		{
			gsl_vector_set(x,j,gsl_vector_get(x,j)+dx);
			f(x,W->fz);
			gsl_vector_sub(W->fz,W->fx);
			for(i=0; i<n; i++)
			{
				gsl_matrix_set(W->J,i,j,gsl_vector_get(W->fz,i)/dx);
			}
			gsl_vector_set(x,j,gsl_vector_get(x,j)-dx);
		}
		givens_qr_dec(W->J);
		lambda = 2;
		normfx = gsl_blas_dnrm2(W->fx);
		givens_qr_bak(W->J,W->fx,W->Dx);
		do
		{
			lambda /= 2.;
			vector_add(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
		} while(gsl_blas_dnrm2(W->fz) > normfx && lambda > 0.01);
		gsl_vector_memcpy(x,W->z);
		gsl_vector_memcpy(W->fx,W->fz);
	} while (gsl_blas_dnrm2(W->Dx) > dx && gsl_blas_dnrm2(W->fx) > tol);
	return counter;
}

int newton_derivative(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, newton_workspace* W)
{
	int counter = 0;
	double lambda, normfx;
	do
	{
		counter++;
		f(x,W->fx);
		df(x,W->J);
		givens_qr_dec(W->J);
		lambda = 2;
		normfx = gsl_blas_dnrm2(W->fx);
		givens_qr_bak(W->J,W->fx,W->Dx);
		do
		{
			lambda /= 2.;
			vector_add(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
		} while(gsl_blas_dnrm2(W->fz) > normfx && lambda > 0.01);
		gsl_vector_memcpy(x,W->z);
		gsl_vector_memcpy(W->fx,W->fz);
	} while (gsl_blas_dnrm2(W->fx) > tol);
	return counter;
}

int newton_derivative_interp(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, newton_workspace* W)
{
	int counter = 0;
	double g, gp, gl, c, lambda, normfx, normfz;
	do
	{
		counter++;
		f(x,W->fx);
		df(x,W->J);
		givens_qr_dec(W->J);
		normfx = gsl_blas_dnrm2(W->fx);
		givens_qr_bak(W->J,W->fx,W->Dx);
		lambda = 1;
		vector_add(1,x,-lambda,W->Dx,W->z);
		f(W->z,W->fz);
		normfz = gsl_blas_dnrm2(W->fz);
		while(normfz > normfx && lambda > 0.01)
		{
			g = 0.5*normfx*normfx;
			gp = -normfx*normfx;
			gl = 0.5*normfz*normfz;
			c = (gl-g-gp*lambda)/lambda/lambda;
			lambda = -gp/2/c;
			vector_add(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
			normfz = gsl_blas_dnrm2(W->fz);
		}
		gsl_vector_memcpy(x,W->z);
		gsl_vector_memcpy(W->fx,W->fz);
	} while (normfz > tol);
	return counter;
}
