#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include "../2_lineq/givens.h"

typedef struct
{
	int n;
	gsl_vector* df;
	gsl_vector* Dx;
	gsl_vector* z;
	gsl_matrix* H;

} minimisation_newton_workspace;

int minimisation_newton(double f(gsl_vector* x), void gradient(gsl_vector* x, gsl_vector* df), void hessian(gsl_vector* x, gsl_matrix* H), gsl_vector* x, double alpha, double tol, minimisation_newton_workspace* W)
{
	int steps = 0;
	double lambda, fz, fx, dot, normdfx;
	do
	{
		steps++;
		fx = f(x);
		gradient(x,W->df);
		hessian(x,W->H);
		givens_qr_dec(W->H);
		givens_qr_bak(W->H,W->df,W->Dx);
		gsl_blas_ddot(W->df,W->Dx,&dot);
		lambda = -2;
		do
		{
			lambda /= 2;
			vector_sum(1,x,lambda,W->Dx,W->z);
			fz = f(W->z);
		} while(fz > fx + alpha*lambda*dot);
		gsl_vector_memcpy(x,W->z);

	} while(normdfx > tol);
	return steps;
}
