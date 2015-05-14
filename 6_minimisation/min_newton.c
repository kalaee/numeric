#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "../2_lineq/givens.h"
#include "../libs/vector.h"

// workspace struceture for minimisation using Newton's method
typedef struct
{
	gsl_vector* df;
	gsl_vector* Dx;
	gsl_vector* z;
	gsl_matrix* H;
} min_newton_workspace;

// allocate memory for newton workspace, n is number of variables
min_newton_workspace* min_newton_workspace_alloc(int n)
{
	min_newton_workspace* W = (min_newton_workspace*) malloc(sizeof(min_newton_workspace));
	W->df = gsl_vector_alloc(n);
	W->Dx = gsl_vector_alloc(n);
	W->z = gsl_vector_alloc(n);
	W->H = gsl_matrix_alloc(n,n);
	return W;
}

// free allocated memory for workspace
void min_newton_workspace_free(min_newton_workspace* W)
{
	gsl_vector_free(W->df);
	gsl_vector_free(W->Dx);
	gsl_vector_free(W->z);
	gsl_matrix_free(W->H);
	free(W);
}

// Minimisation using Newton's method with backtracking, f is function to be minimised,
// gradient a function which calculates gradient of f at x and store it in df,
// hessian calculates the Hessian matrix at x and store it in H
// x is coordinates of initial guess, on output the solution is stored in x
// alpha i the scaling factor for the Armijo condition in backtracking, can be as low as 1e-4
// tol is tolerance for evaluating convergence, if: norm(df) < tol, convergence is achieved
int min_newton(double f(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* x, double alpha, double tol, min_newton_workspace* W)
{
	int steps = 0;
	double lambda, fz, fx, dot;
	// gradient is calculated before do-while loop,
	// for the remainder of the algorithm the gradient for an
	// iteration is calculated in the previous iteration due to
	// evaluation of the while-condition (norm(gradient) > tol)
	gradient(x,W->df);
	do
	{
		steps++;
		fx = f(x);
		hessian(x,W->H);
		// solve the system H Dx = df using QR decomp with Givens rotations
		givens_qr_dec(W->H);
		givens_qr_bak(W->H,W->df,W->Dx);
		// attempt a step x -> z = x+lambda*Dx
		// if f(z) does not satisfy Armijo condition,
		// f(z) < f(x) + alpha*lambda*Dx^T*df,
		// backtrack using lambda -> lambda/2 and attempt again
		// if lambda too small, take step just to continue somewhere else
		gsl_blas_ddot(W->df,W->Dx,&dot);
		lambda = -2;
		do
		{
			lambda /= 2;
			vector_sum(1,x,lambda,W->Dx,W->z);
			fz = f(W->z);
		} while(fz > fx + alpha*lambda*dot && lambda > 0.01);
		gsl_vector_memcpy(x,W->z);
		gradient(x,W->df);
	// if norm of gradient is lower than tol convergence is achieved
	} while(gsl_blas_dnrm2(W->df) > tol);
	// return number of steps taken during minimisation
	return steps;
}

