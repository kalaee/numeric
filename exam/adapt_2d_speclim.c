#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ADAPT_NRECUR_MAX 100000

// perform integration over the y coordinate for a given x value
// the routine uses adaptive quadratures with open intervals
double adapt24_inner_speclim(double f(double x, double y), double x, double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	if (a == b)
	{
		*err = 0;
		return 0;
	}
	else
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
			double Q1 = adapt24_inner_speclim(f,x,a,a+h/2,y1,y2,acc,eps,&err1,nrecur+1);
			double Q2 = adapt24_inner_speclim(f,x,a+h/2,b,y3,y4,acc,eps,&err2,nrecur+1);
			*err  = sqrt(err1*err1 + err2*err2);
			return Q1 + Q2;
		}
	}
}

// perform integration over x, where for each y1 and y4, integration over y-coordinate is performed via adapt24_inner
double adapt24_outer_speclim(double f(double x, double y), double ax, double bx, double d(double x), double u(double x), double y2, double y3, double acc, double eps, double *err, int nrecur)
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
	double acc_inner = acc/2.;	// reduce absolute tolerance by sqrt(4) when passing to adapt24_inner
	double x1 = ax+hx/6.;	// coordinate to evaluate y?_x1_inner below
	double x4 = ax+hx*5./6.;	// coordinate to evaluate y?_x4_inner below
	double ay_x1 = d(x1);
	double ay_x4 = d(x4);
	double by_x1 = u(x1);
	double by_x4 = u(x4);
	double hy_x1 = by_x1 - ay_x1;
	double hy_x4 = by_x4 - ay_x4;
	double y2_x1_inner = f(x1,ay_x1+hy_x1*2./6.);
	double y3_x1_inner = f(x1,ay_x1+hy_x1*4./6.);
	double y2_x4_inner = f(x4,ay_x4+hy_x4*2./6.);
	double y3_x4_inner = f(x4,ay_x4+hy_x4*4./6.);
	double y1 = adapt24_inner_speclim(f,x1,ay_x1,by_x1,y2_x1_inner,y3_x1_inner,acc_inner,eps,&dummy,0);
	double y4 = adapt24_inner_speclim(f,x4,ay_x4,by_x4,y2_x4_inner,y3_x4_inner,acc_inner,eps,&dummy,0);
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
		double Q1 = adapt24_outer_speclim(f,ax,ax+hx/2,d,u,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = adapt24_outer_speclim(f,ax+hx/2,bx,d,u,y3,y4,acc,eps,&err2,nrecur+1);
		*err  = sqrt(err1*err1 + err2*err2);
		return Q1 + Q2;
	}
}
double adapt_2d_speclim(double f(double x, double y), double ax, double bx, double d(double x), double u (double x), double acc, double eps, double* err)
{
	if (ax == bx)
	{
		*err = 0;
		return 0;
	}
	else
	{
		// prepare neccesary arguments for first iteration of integration
		double dummy;
		double hx = bx-ax;
		double acc_inner = acc/2.;
		double x2 = ax+hx*2./.6;
		double x3 = ax+hx*4./6.;
		double ay_x2 = d(x2);
		double ay_x3 = d(x3);
		double by_x2 = u(x2);
		double by_x3 = u(x3);
		double hy_x2 = by_x2 - ay_x2;
		double hy_x3 = by_x3 - ay_x3;
		double y2_x2_inner = f(x2,ay_x2+hy_x2*2./6.);
		double y3_x2_inner = f(x2,ay_x2+hy_x2*4./6.);
		double y2_x3_inner = f(x3,ay_x3+hy_x3*2./6.);
		double y3_x3_inner = f(x3,ay_x3+hy_x3*4./6.);
		double y2 = adapt24_inner_speclim(f,x2,ay_x2,by_x2,y2_x2_inner,y3_x2_inner,acc_inner,eps,&dummy,0);
		double y3 = adapt24_inner_speclim(f,x3,ay_x3,by_x3,y2_x3_inner,y3_x3_inner,acc_inner,eps,&dummy,0);
		// integrate
		double Q = adapt24_outer_speclim(f,ax,bx,d,u,y2,y3,acc,eps,err,0);
		return Q;
	}
}
