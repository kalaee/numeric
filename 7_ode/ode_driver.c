#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>

#define ODE_SUCCESS	0
#define ODE_FAIL	1
#define ODE_FAIL_UPPER	1000

typedef struct
{
	int n;	// number of function
	void* s;	// workspace og step function
	gsl_vector* err;
	gsl_vector* yh;
} ode_workspace;

ode_workspace* ode_workspace_alloc(int n, void* STEPPER_WORKSPACE_ALLOC(int n))
{
	ode_workspace* W = (ode_workspace*) malloc(sizeof(ode_workspace));
	W->n = n;
	W->err = gsl_vector_alloc(n);
	W->yh = gsl_vector_alloc(n);
	W->s = STEPPER_WORKSPACE_ALLOC(n);
	return W;
}

void ode_workspace_free(ode_workspace* W, void STEPPER_WORKSPACE_FREE(void* STEPPER_WORKSPACE))
{
	gsl_vector_free(W->yh);
	gsl_vector_free(W->err);
	STEPPER_WORKSPACE_FREE(W->s);
	free(W);
	return;
}

int ode_driver(void f(double t, gsl_vector* y, gsl_vector* dy),
	void STEPPER(double t, double h, gsl_vector* y,
		void f(double t, gsl_vector* y, gsl_vector* dydt),
		gsl_vector* yh, gsl_vector* err, void* WORKSPACE),
	double* a, double b, double* h0, gsl_vector* ya, double acc, double eps, ode_workspace* W)
{
	int flag = ODE_FAIL;
	double err, tol;
	// make sure that end-point is actually b
	if(*a == b)
	{
		return ODE_SUCCESS;
	}
	if( *a + *h0 > b)
	{
		*h0 = b - *a;
	}
	
	do
	{
		flag += ODE_FAIL;
		STEPPER(*a,*h0,ya,f,W->yh,W->err,W->s);
		tol = (gsl_blas_dnrm2(W->yh)*eps+acc)*sqrt(*h0/(b-*a));
		err = gsl_blas_dnrm2(W->err);
		if (err < tol)
		{
			flag = ODE_SUCCESS;
			*a += *h0;
			gsl_vector_memcpy(ya,W->yh);
		}
		if (err > 0)
		{
			*h0 *= pow(tol/err,0.25)*0.95;
		}
		else
		{
			*h0 *= 2;
		}
	} while (flag != ODE_SUCCESS && flag < ODE_FAIL_UPPER);
	if (flag != ODE_SUCCESS)
	{
		fprintf(stderr,"ode_driver failed to achieve desired precision for t = %g after %d iterations.\n",*a,flag);
	}
	return flag;
}
