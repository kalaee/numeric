#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "../2_lineq/givens.h"
#include "../libs/vector.h"

// workspace struceture for minimisation using Newton's method with SR1
typedef struct
{
	gsl_vector* df;
	gsl_vector* df0;
	gsl_vector* Dx;
	gsl_vector* z;
	gsl_vector* y;
	gsl_matrix* H;
} min_newton_sr1_workspace;

// allocate memory for newton SR1 workspace, n is number of variables
min_newton_sr1_workspace* min_newton_sr1_workspace_alloc(int n)
{
	min_newton_sr1_workspace* W = (min_newton_sr1_workspace*) malloc(sizeof(min_newton_sr1_workspace));
	W->df = gsl_vector_alloc(n);
	W->df0 = gsl_vector_alloc(n);
	W->y = gsl_vector_alloc(n);
	W->Dx = gsl_vector_alloc(n);
	W->z = gsl_vector_alloc(n);
	W->H = gsl_matrix_alloc(n,n);
	return W;
}

// free allocated memory for workspace
void min_newton_sr1_workspace_free(min_newton_sr1_workspace* W)
{
	gsl_vector_free(W->df);
	gsl_vector_free(W->df0);
	gsl_vector_free(W->y);
	gsl_vector_free(W->Dx);
	gsl_vector_free(W->z);
	gsl_matrix_free(W->H);
	free(W);
}

// Minimisation using Quasi Newton's method with backtracking and SR1, f is function to be minimised,
// function does not need to calculate Hessian matrix, can be approximated using symmetric rank 1 updates
// gradient is a function which calculates gradient of f at x and store it in df,
// x is coordinates of initial guess, on output the solution is stored in x
// alpha i the scaling factor for the Armijo condition in backtracking, can be as low as 1e-4
// tol is tolerance for evaluating convergence, if: norm(df) < tol, convergence is achieved
// on output the functions returns number of steps taken
int min_newton_sr1(double f(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), gsl_vector* x, double alpha, double tol, min_newton_sr1_workspace* W)
{
	int steps = 0;
	double lambda, fz, fx, dot;
	// gradient is calculated before do-while loop,
	// for the remainder of the algorithm the gradient for an
	// iteration is calculated in the previous iteration due to
	// evaluation of the while-condition (norm(gradient) > tol)
	// The inverse Hessian matrix is initialised as the unity matrix
	gradient(x,W->df);
	gsl_matrix_set_identity(W->H);
	do
	{
		steps++;
		fx = f(x);
		gsl_blas_dsymv(CblasUpper,1,W->H,W->df,0,W->Dx);
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
		} while(fz > fx + alpha*lambda*dot && lambda > 0.5);
		// define y = df - df0 for use in SR1, note that the Hessian matrix is symmetric
		gsl_vector_memcpy(W->df0,W->df);
		gsl_vector_memcpy(x,W->z);
		gradient(x,W->df);
		vector_sum(1,W->df,-1,W->df0,W->y);
		// apply SR1 update to H using knowledge that H is symmetric
		gsl_blas_dsymv(CblasUpper,-1,W->H,W->y,lambda,W->Dx);
		gsl_blas_ddot(W->Dx,W->y,&dot);
		gsl_blas_dsyr(CblasUpper,1./dot,W->Dx,W->H);
	// if norm of gradient smaller than tolerance, convergence is achieved
	} while (gsl_blas_dnrm2(W->df) > tol);
	return steps;
}
