#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
#include <assert.h>

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

// allocate workspace for ode_driver, requires specification of stepper type: ODE_RKF45 or ODE_RK5, n is length of vector
ode_workspace* ode_workspace_alloc(int n, void* STEPPER_WORKSPACE_ALLOC(int n))
{
	ode_workspace* W = (ode_workspace*) malloc(sizeof(ode_workspace));
	W->n = n;
	W->err = gsl_vector_alloc(n);
	W->yh = gsl_vector_alloc(n);
	W->s = STEPPER_WORKSPACE_ALLOC(n);
	return W;
}

// free memory for ode_driver, requires specification of stepper type: ODE_RKF45 or ODE_RK5
void ode_workspace_free(ode_workspace* W, void STEPPER_WORKSPACE_FREE(void* STEPPER_WORKSPACE))
{
	gsl_vector_free(W->yh);
	gsl_vector_free(W->err);
	STEPPER_WORKSPACE_FREE(W->s);
	free(W);
	return;
}

// routine which evolve the differential equaiton exactly one step using adaptive step size
// f is the ode-function, which takes variable t, initial condition in vector y and outputs derivative in dy
// STEPPER designates which step-function is used, options are RKF45 or RK5
// *a is pointer to begining value of t, b is end value of t, *h0 pointer to initial step
// ya vector of initial conditions, acc absolute tolerance, eps relative tolerance and
// W is ode_workspace
// On termination the incremented bottom value is stored in a, next value of y is stored in y
// and adjusted step size is stored in h0
// note that the routine assumes a is smaller than or equal to b
int ode_evolve(void f(double t, gsl_vector* y, gsl_vector* dy),
	void STEPPER(double t, double h, gsl_vector* y,
		void f(double t, gsl_vector* y, gsl_vector* dydt),
		gsl_vector* yh, gsl_vector* err, void* WORKSPACE),
	double* a, double b, double* h0, gsl_vector* ya, double acc, double eps, ode_workspace* W)
{
	int flag = ODE_FAIL;
	double err, tol;
	// do we need to take a step?
	if(*a == b)
	{
		return ODE_SUCCESS;
	}
	// if step bigger than allowed, decrease to maximum
	if( *a + *h0 > b)
	{
		*h0 = b - *a;
	}
	// step-routine. Will continue until step succeded or loop has exhausted allowed number of attempts (ODE_FAIL_UPPER)
	do
	{
		// for every attempt, flag is increase by one ODE_FAIL to counter number of attempts
		flag += ODE_FAIL;
		// attempt the step
		STEPPER(*a,*h0,ya,f,W->yh,W->err,W->s);
		// estimate error
		tol = (gsl_blas_dnrm2(W->yh)*eps+acc)*sqrt(*h0/(b-*a));
		err = gsl_blas_dnrm2(W->err);
		// if error smaller than tolerance, accept result and increase a by the step
		if (err < tol)
		{
			flag = ODE_SUCCESS;
			*a += *h0;
			gsl_vector_memcpy(ya,W->yh);
		}
		// if error is non-zero increase h0 accordingly
		if (err > 0)
		{
			*h0 *= pow(tol/err,0.25)*0.95;
		}
		// if error is zero, attempt to make the next step twice as big
		else
		{
			*h0 *= 2;
		}
	} while (flag != ODE_SUCCESS && flag < ODE_FAIL_UPPER);
	// if terminated loop because of number of attempts, send message to stream: stderr
	if (flag != ODE_SUCCESS)
	{
		fprintf(stderr,"ode_driver failed to achieve desired precision for t = %g after %d iterations.\n",*a,flag);
	}
	return flag;
}


// This routine drives the ode system from a to b using ode_evolve and stores the steps in matrix y
// y must have y0->size+1 columns and sufficiently high number of rows to store every step
// f is the ode-function, which takes variable t, initial condition in vector y and outputs derivative in dy
// STEPPER designates which step-function is used, options are RKF45 or RK5
// a is begining value of t, b is end value of t, h is the guessed initial step
// y0 vector of initial conditions, acc absolute tolerance, eps relative tolerance and
// W is ode_workspace
// note that the routine assumes a is smaller than or equal to b
int ode_driver(void f(double t, gsl_vector* y, gsl_vector* dy),
	void STEPPER(double t, double h, gsl_vector* y,
	void f(double t, gsl_vector* y, gsl_vector* dydt),
	gsl_vector* yh, gsl_vector* err, void* WORKSPACE),
	double a, double b, double h, gsl_vector* y0, gsl_matrix* y, double acc, double eps, ode_workspace* W)
{
	// assert that matrix y and vector y0 have corresponding dimensions
	assert(y->size2 == y0->size + 1);
	int status, i,j;
	// store initial conditions in matrix y
	i = 0;
	gsl_matrix_set(y,0,0,a);
	for(j=0; j<y0->size; j++)
	{
		gsl_matrix_set(y,0,j+1,gsl_vector_get(y0,j));
	}
	// while-loop to drive ode_evolve from a to b
	while(a < b)
	{
		i++;
		// if we have taken too many steps to store them in y-matrix, send message to stderr
		if (!(i<y->size1))
		{
			fprintf(stderr,"Routine ode_driver exceeded capacity of provided matrix.\n");
			break;
		}
		// take next step using ode_evolve
		status = ode_evolve(f,STEPPER,&a,b,&h,y0,acc,eps,W);
		// if step went wrong break ode_drive
		if (status != ODE_SUCCESS)
		{
			break;
		}
		// store new step in matrix
		gsl_matrix_set(y,i,0,a);
		for(j=0; j<y0->size; j++)
		{
			gsl_matrix_set(y,i,j+1,gsl_vector_get(y0,j));
		}
	}
	// return index of last used row in y
	return i;
}
