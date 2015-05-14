#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "../2_lineq/givens.h"
#include "vector.h"
// roots_newton workspace
typedef struct
{
	int n;
	gsl_matrix* J;
	gsl_vector* z;
	gsl_vector* Dx;
	gsl_vector* fx;
	gsl_vector* fz;
} roots_newton_workspace;

// allocate necessary memory for roots_newton root-finding
// it is assumed that a function of n variables returns a vector
// of length n
roots_newton_workspace* roots_newton_workspace_alloc(int n)
{
	roots_newton_workspace* W = (roots_newton_workspace*) malloc(sizeof(roots_newton_workspace));
	W->n = n;
	W->J = gsl_matrix_alloc(n,n);
	W->z = gsl_vector_alloc(n);
	W->Dx = gsl_vector_alloc(n);
	W->fx = gsl_vector_alloc(n);
	W->fz = gsl_vector_alloc(n);

	return W;
}

// free allocated memory for roots_newton root-finding
void roots_newton_workspace_free(roots_newton_workspace* W)
{
	gsl_matrix_free(W->J);
	gsl_vector_free(W->z);
	gsl_vector_free(W->Dx);
	gsl_vector_free(W->fx);
	gsl_vector_free(W->fz);
	free(W);
}

// ordinary roots_newton root-finding for functions with unknown derivative
// f is the function, x is initial guess, on termination solution is stored in x
// dx is length scale for numerical estimate of the derivative and tol the tolerance for convergence
int roots_newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, roots_newton_workspace* W)
{
	int n = W->n;
	int i,j, counter = 0;
	double lambda, normfx;
	do
	{
		counter++;
		f(x,W->fx);
		// numerical estimate of derivative and store in Jacobi matrix, J
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
		// solve system J Dx = f(x)
		givens_qr_dec(W->J);
		lambda = 2;
		normfx = gsl_blas_dnrm2(W->fx);
		givens_qr_bak(W->J,W->fx,W->Dx);
		// find lambda which satisfies f(x-l*Dx) < (1-l/2)*f(x)
		// if not make lambda smaller: lambda = lambda/2
		// if lambda < 0.01 does not satisfy our requirement, move anyway to
		// to proceed from another position
		do
		{
			lambda /= 2.;
			vector_sum(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
		} while(gsl_blas_dnrm2(W->fz) > normfx && lambda > 0.01);
		gsl_vector_memcpy(x,W->z);
	// terminate algorithm if convergence achieved or Dx is smaller than out length scale, dx
	} while (gsl_blas_dnrm2(W->Dx) > dx && gsl_blas_dnrm2(W->fz) > tol);
	// return number of steps
	return counter;
}

// roots_newton method for root-finding with known derivative, df, x is initial guess, on output it will contain solution
// tol is convergence tolerance and W the neccesary workspace
int roots_newton_derivative(void f(gsl_vector* x, gsl_vector* fx), void df(gsl_vector* x, gsl_matrix* J), gsl_vector* x, double tol, roots_newton_workspace* W)
{
	int counter = 0;
	double lambda, normfx;
	do
	{
		counter++;
		f(x,W->fx);
		// calculate Jacobi matrix
		df(x,W->J);
		givens_qr_dec(W->J);
		lambda = 2;
		normfx = gsl_blas_dnrm2(W->fx);
		// solve equation J Dx = f(x)
		givens_qr_bak(W->J,W->fx,W->Dx);
		// find lambda which satisfies f(x-l*Dx) < (1-l/2)*f(x)
		// if not make lambda smaller: lambda = lambda/2
		// if lambda < 0.01 does not satisfy our requirement, move anyway to
		// to proceed from another position
		do
		{
			lambda /= 2.;
			vector_sum(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
		} while(gsl_blas_dnrm2(W->fz) > (1-lambda/2)*normfx && lambda > 0.01);
		gsl_vector_memcpy(x,W->z);
	// if convergence criteria meet, end algorithm
	} while (gsl_blas_dnrm2(W->fz) > tol);
	// return number of steps
	return counter;
}


// roots_newton method for root-finding with unknown derivative and quadratic interpolation
// x is initial guess, on output it will contain the solution of the system
// dx is length for numerical estimate of derivative and tol is tolerance for convergence
int roots_newton_interp(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double tol, roots_newton_workspace* W)
{
	int i,j, n = W->n, counter = 0;
	double g0, gp0, gl, c, lambda, normfx, normfz;
	do
	{
		counter++;
		f(x,W->fx);
		// estimate Jacobi matrix numerically
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
		// solve J Dx = f(x)
		givens_qr_dec(W->J);
		normfx = gsl_blas_dnrm2(W->fx);
		givens_qr_bak(W->J,W->fx,W->Dx);
		// find lambda which satisfies f(x-l*Dx) < (1-l/2)*f(x)
		// if not make lambda smaller using quadratic interpolation for minimization
		// if lambda < 0.01 does not satisfy our requirement, move anyway to
		// to proceed from another position
		// note that quadratic interpolation will only work close to an actual solution
		lambda = 1;
		vector_sum(1,x,-lambda,W->Dx,W->z);
		f(W->z,W->fz);
		g0 = 0.5*normfx*normfx;
		gp0 = -normfx*normfx;
		normfz = gsl_blas_dnrm2(W->fz);
		while(normfz > normfx && lambda > 0.1)
		{
			gl = 0.5*normfz*normfz;
			c = (gl-g0-gp0*lambda)/lambda/lambda;
			lambda = -gp0/2/c;
			vector_sum(1,x,-lambda,W->Dx,W->z);
			f(W->z,W->fz);
			normfz = gsl_blas_dnrm2(W->fz);
		}
		gsl_vector_memcpy(x,W->z);
	// end algorithm if Dx is smaller than length scale dx or convergence is achieved
	} while (gsl_blas_dnrm2(W->Dx) > dx && normfz > tol);
	// return number of steps before solution
	return counter;
}

