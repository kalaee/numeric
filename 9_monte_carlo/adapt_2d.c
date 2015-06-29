#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ADAPT_NRECUR_MAX 100000

// perform integration over the y coordinate for a given x value
// the routine uses adaptive quadratures with open intervals
double adapt24_inner(double f(double x, double y), double x, double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	// check iteration level is not too deep
	// if too big, send error message and return 0
	if (nrecur > ADAPT_NRECUR_MAX)
	{
		fprintf(stderr,"Error: adapt24_inner reached too many subdivisions!\n");
		exit(EXIT_FAILURE);
	}
	// make integral estimate of fourth and second order
	double h = b-a;
	double y1 = f(x,a+h/6.);
	double y4 = f(x,a+h*5./6.);
	double Q = h/6.*(2.*y1+y2+y3+2.*y4);
	double q = h/2.*(y1+y4);
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
		double Q1 = adapt24_inner(f,x,a,a+h/2,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = adapt24_inner(f,x,a+h/2,b,y3,y4,acc,eps,&err2,nrecur+1);
		*err  = sqrt(err1*err1 + err2*err2);
		return Q1 + Q2;
	}
}

// perform integration over x, where for each y1 and y4, integration over y-coordinate is performed via adapt24_inner
double adapt24_outer(double f(double x, double y), double ax, double ay, double bx, double by, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	// check iteration level is not too deep
	// if too big, send error message and return 0
	if (nrecur > ADAPT_NRECUR_MAX)
	{
		fprintf(stderr,"Error: adapt24_outer reached too many subdivisions!\n");
		exit(EXIT_FAILURE);
	}

	// make integral estimate of fourth and second order
	double dummy;	// dummy for error-pointer to adapt24_inner below
	double hx = bx-ax;	// width of interval in x coordinate
	double hy = by-ay;	// width in y-coordinate
	double acc_inner = acc/2.;	// reduce absolute tolerance by sqrt(4) when passing to adapt24_inner
	double x1 = ax+hx/6.;	// coordinate to evaluate y?_x1_inner below
	double x4 = ax+hx*5./6.;	// coordinate to evaluate y?_x4_inner below
	double y2_x1_inner = f(x1,ay+hy*2./6.);
	double y3_x1_inner = f(x1,ay+hy*4./6.);
	double y2_x4_inner = f(x4,ay+hy*2./6.);
	double y3_x4_inner = f(x4,ay+hy*4./6.);
	double y1 = adapt24_inner(f,x1,ay,by,y2_x1_inner,y3_x1_inner,acc_inner,eps,&dummy,0);
	double y4 = adapt24_inner(f,x4,ay,by,y2_x4_inner,y3_x4_inner,acc_inner,eps,&dummy,0);
	double Q = hx/6.*(2.*y1+y2+y3+2*y4); // integral estimate to fourth order
	double q = hx/2.*(y1+y4);	// integral estimate to second order
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
		double Q1 = adapt24_outer(f,ax,ay,ax+hx/2,by,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = adapt24_outer(f,ax+hx/2,ay,bx,by,y3,y4,acc,eps,&err2,nrecur+1);
		*err  = sqrt(err1*err1 + err2*err2);
		return Q1 + Q2;
	}
}
double adapt_2d(double f(double x, double y), double ax, double ay, double bx, double by, double acc, double eps, double* err)
{
	// prepare neccesary arguments for first iteration of integration
	double dummy;
	double hx = bx-ax;
	double hy = by-ay;
	double acc_inner = acc/2.;
	double x2 = ax+hx*2./.6;
	double x3 = ax+hx*4./6.;
	double y2_x2_inner = f(x2,ay+hy*2./6.);
	double y3_x2_inner = f(x2,ay+hy*4./6.);
	double y2_x3_inner = f(x3,ay+hy*2./6.);
	double y3_x3_inner = f(x3,ay+hy*4./6.);
	double y2 = adapt24_inner(f,x2,ay,by,y2_x2_inner,y3_x2_inner,acc_inner,eps,&dummy,0);
	double y3 = adapt24_inner(f,x3,ay,by,y2_x3_inner,y3_x3_inner,acc_inner,eps,&dummy,0);
	// integrate
	double Q = adapt24_outer(f,ax,ay,bx,by,y2,y3,acc,eps,err,0);
	return Q;
}
