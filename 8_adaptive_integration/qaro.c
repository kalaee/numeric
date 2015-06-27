#include <math.h>
#include <assert.h>
#include <stdio.h>

// Quadratic Adaptive integral using Recursive functions and Open intervals
// with trapezoid evaluations of order 4 with 2 for error estimate

// this function makes the actual integration of the interval
// f is the integrand function, a and b the first and last coordinates of the interval,
// y2 and y3 are the previously estimated values of the integrand relevant at points
// a+2*h/6 and b-2*h/6, where h = b - a, respectively
// acc allowed absolute error and eps relative allowed error,
// *err is pointer to where the error estimate is stored, nrecur counts the recursion level
// note that qaro24 assumes a > b
double qaro24(double f(double x), double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	// check iteration level is not too deep
	// if too big, send error message and return 0
	if (nrecur > 100000)
	{
		fprintf(stderr,"Too many subdivisions!\n");
		return 0;
	}
	// make integral estimate of fourth and second order
	double h = b-a;
	double y1 = f(a+h/6);
	double y4 = f(a+h*5/6);
	double Q = h/6*(2*y1+y2+y3+2*y4);
	double q = h/2*(y1+y4);
	*err = fabs(Q-q);
	// if error sufficiently small: return integral value and error estimate
	if ( *err < acc + eps*fabs(Q) )
	{
		return Q;
	}
	// if error too big: split interval in two and integrate each interval
	// recycling estimated points
	else
	{
		acc /= sqrt(2);
		double err1, err2;
		double Q1 = qaro24(f,a,a+h/2,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = qaro24(f,a+h/2,b,y3,y4,acc,eps,&err2,nrecur+1);
		*err  = err1 + err2;
		return Q1 + Q2;
	}
}

// this function takes the integration arguments and prepares them for qaro24
// f is the integrand, a and b the start and end points of the interval,
// acc and eps are the allowed absolute and relative errors,
// *err is pointer to where the error estimate is stored
double qaro(double f(double), double a, double b, double acc, double eps, double * err)
{
	// if a < b, continue to checking for infinities
	if (a < b)
	{
		// the following checks for possible infinities in integration limits
		// and, if found, makes suitable transformation to integrals with finite limits
		if ( !isinf(a) && !isinf(b) )
		{
			double h = b - a;
			double y2 = f(a+h*2/6);
			double y3 = f(a+h*4/6);
			double Q = qaro24(f,a,b,y2,y3,acc,eps,err,0);
			return Q;
		}
		else if ( !isinf(a) && isinf(b) )
		{
			double ft(double t){ return f(a+t/(1-t))/pow(1-t,2); }
			double y2 = ft(2./6.);
			double y3 = ft(4./6.);
			double Q = qaro24(ft,0,1,y2,y3,acc,eps,err,0);
			return Q;
		}
		else if ( isinf(a) && !isinf(b) )
		{
			double ft (double t){ return f(b - (1-t)/t)/t/t; }
			double y2 = ft(2./6.);
			double y3 = ft(4./6.);
			double Q = qaro24(ft,0,1,y2,y3,acc,eps,err,0);
			return Q;
		}
		else if ( isinf(a) && isinf(b) )
		{
			double ft(double t){ return ( f((1-t)/t) + f((t-1)/t))/t/t; }
			double y2 = ft(2./6.);
			double y3 = ft(4./6.);
			double Q = qaro24(ft,0,1,y2,y3,acc,eps,err,0);
			return Q;
		}
		else
		{
				fprintf(stderr,"Something is wrong with the limits!\n");
				return 0;
		}
	}
	// if b < a, make call qaro and return the negative value of integral from b to a
	else if (b < a)
	{
		// reverse order of limits, caro24 and infinity-evalutaion is designed for a < b
		return -qaro(f,b,a,acc,eps,err);
	}
	// if a == b, integral is by definition zero
	else if (a == b )
	{
		*err = 0;
		return 0;
	}
	// otherwise, send error message about limits
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
