#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// Quadratic Adaptive integral using Recursive Closed intervals
// and trapezoid evaluations of order 3 with 2 for error estimate

// thus function makes the actual integration of the interval
// f is the integrand function, a and b the first and last coordinates of the interval,
// y1 and y3 are the previously estimated values of the integrand relevant at points
// a and b respectively
// acc allowed absolute error and eps relative allowed error,
// *err is pointer to where the error estimate is stored, nrecur counts the recursion level
// note that qarc23 assumes a > b
double qarc23(double f(double x), double a, double b, double y1, double y3, double acc, double eps, double *err, int nrecur)
{
	// check iteration level is not too deep
	// if too big, send error message and return 0
	if (nrecur > 100000)
	{
		fprintf(stderr,"Too many subdivisions!\n");
		exit(EXIT_FAILURE);
	}
	// make estimate of third and second order
	double h = b-a;
	double y2 = f(a+h/2.);
	double Q = h/4.*(y1+2.*y2+y3);
	double q = h/2.*(y1+y3);
	*err = fabs(Q-q);
	// if error sufficiently small: return integral value and error estimate
	if ( *err < acc + eps*fabs(Q) )
	{
		return Q;
	}
	// if error too big: split interval in two, and integrate each interval
	// recycling estimated point
	else
	{
		acc /= sqrt(2);
		double err1, err2;
		double Q1 = qarc23(f,a,a+h/2,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = qarc23(f,a+h/2,b,y2,y3,acc,eps,&err2,nrecur+1);
		*err = sqrt(err1*err1 + err2*err2);
		return Q1 + Q2;
	}
}

// this function takes the integration arguments and prepares them for qarc23
// f is the integrand, a and b the start and end points of the interval,
// acc and eps are the allowed absolute and relative errors,
// *err is pointer to where the error estimate is stored
// note that due to this algorithm using closed intervals, it cannot apply transformations
// for integrals with infinities in the limits
double qarc(double f(double), double a, double b, double acc, double eps, double * err)
{
	// if a < b, perform integration
	if (a < b)
	{
		double y1 = f(a);
		double y3 = f(b);
		double Q = qarc23(f,a,b,y1,y3,acc,eps,err,0);
		return Q;		
	}
	// if b < a, return the negative of integral from b to a
	else if (b < a)
	{
		return -qarc(f,b,a,acc,eps,err);
	}
	// a == b, integral is by definition zero
	else if (a == b)
	{
		*err = 0;
		return 0;
	}
	// otherwise, return error message about limits
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		exit(EXIT_FAILURE);
	}
}
