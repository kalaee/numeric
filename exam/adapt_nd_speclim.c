#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define ADAPT_NRECUR_MAX 10000

double adapt_nd_speclim_next(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, int lvl, double* coord, double acc, double eps, double* err);

// integrate over the last coordinate, x[dim-1], using trapezoid adaptive quadrature of order 24
// x1 provides the coordinates x[i], i=0,...,dim-2
// values y2 and y3 are provided from preceding calculations
double adapt_nd_speclim_last(double f(double* x), double a, double b, int dim, double* x1, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	// check iteration level is not too deep
	// if too big, send error message and return 0
	if (nrecur > ADAPT_NRECUR_MAX)
	{
		fprintf(stderr,"adapt_nd_speclim_last: recursion level too deep in last dimension!\n");
		exit(EXIT_FAILURE);
	}
	// prepare integral estimate of fourth and second order
	int i;
	double h, x4[dim], Q, q, y1, y4;
	h = b-a;
	// copy coordinate from x1 to x4
	for (i=0; i<dim-1; i++)
	{
		x4[i] = x1[i];
	}
	// estimate points at x1 and x4 for points y1 and y4 in 24 quadrature
	x1[dim-1] = a + h/6;
	x4[dim-1] = a + h*5/6;
	y1 = f(x1);
	y4 = f(x4);
	// estimate integral values and error
	Q = h/6.*(2.*y1+y2+y3+2.*y4);
	q = h/2.*(y1+y4);
	*err = fabs(Q-q);
	// if error sufficiently small: return integral value and error estimate
	if ( *err < acc + eps*fabs(Q) )
	{
		return Q;
	}
	// if error too big: split interval in two and integrate each interval
	// recycling estimated points and increment recursion number, nrecur
	else
	{
		acc /= sqrt(2);
		double Q1, Q2, err1, err2;
		Q1 = adapt_nd_speclim_last(f,a,a+h/2,dim,x1,y1,y2,acc,eps,&err1,nrecur+1);
		Q2 = adapt_nd_speclim_last(f,a+h/2,b,dim,x4,y3,y4,acc,eps,&err2,nrecur+1);
		*err  = sqrt(err1*err1 + err2*err2);
		return Q1 + Q2;
	}
}

// integrate over coordinate x[lvl], x1 provides the designated coordinates for x[i], i= 0,...,lvl-1
// this routine works by recursively calling itself to integrate over the next coordinate
// if the estimated error is bigger than allowed, the integration interval is split in two
// and the function is recursively called twice, for each half, at the same coordinate lvl
// albeit with an incremented nrecur, for the next recursion level
// needs the values y2 and y3 of this interval from preceding calculations
double adapt_nd_speclim_now(double f(double* x), double a, double b, double d(double* x, int p), double u(double* x, int p), int dim, int lvl, double* x1, double y2, double y3, double acc, double eps, double*err, int nrecur)
{
	// if recursion too deep, assume integration has failed
	if (nrecur > ADAPT_NRECUR_MAX)
	{
		fprintf(stderr,"adapt_nd_speclim_now: recursion level too deep in dim = %d!\n",lvl);
		exit(EXIT_FAILURE);
	}
	int i;
	double h, x4[dim], Q, q, acc_next, err1, err2, y1, y4;
	// find width at given coordinate
	h = b - a;
	// copy coordinate to new array x4
	for(i=0; i<lvl; i++)
	{
		x4[i] = x1[i];
	}
	// integration using trapezoid quadrature of order 24
	x1[lvl] = a + h/6;
	x4[lvl] = a + h*5/6;
	// estimate y1 and y4 from integration over next coordinate level
	// the absolute tolerance is scaled by sqrt(4) for the next level
	// (because of the neccesaery estimation of foure points in order 24)
 	acc_next = acc/2;
	y1 = adapt_nd_speclim_next(f,d,u,dim,lvl+1,x1,acc_next,eps,&err1);
	y4 = adapt_nd_speclim_next(f,d,u,dim,lvl+1,x4,acc_next,eps,&err2);
	Q = h/6*(2*y1+y2+y3+2*y4);
	q = h/2*(y1+y4);
	*err = fabs(Q-q);
	// if error lower than tolerated, return integral value
	if (*err < acc + eps*fabs(Q))
	{
		return Q;
	}
	// if not, divide interval in two and repeat procedure for each
	// scale the absolute tolerance accordingly
	// additionally: increment the recursion number
	else
	{
		acc /= sqrt(2);
		double err1, err2, Q1, Q2;
		Q1 = adapt_nd_speclim_now(f,a,a+h/2,d,u,dim,lvl,x1,y1,y2,acc,eps,&err1,nrecur+1);
		Q2 = adapt_nd_speclim_now(f,a+h/2,b,d,u,dim,lvl,x4,y3,y4,acc,eps,&err2,nrecur+1);
		// add the error quadratically
		*err = sqrt(err1*err1+err2*err2);
		return Q1 + Q2;
	}
}

// this routine initiates integration of next coordinate x[lvl] at the designated coordinates x3[i], i=0,...,lvl-1
double adapt_nd_speclim_next(double f(double* x), double d(double* x, int p), double u(double* x, int p), int dim, int lvl, double* x2, double acc, double eps, double* err)
{

	int i;
	double Q, h, a, b, y2, y3, x3[dim], acc_next, err1, err2;
	// estimate coordinate width at desired lvl,
	// if limits are identical, integral is zero
	a = d(x2,lvl);
	b = u(x2,lvl);
	if (a == b)
	{
		*err = 0;
		return 0;
	}
	h = b - a;
	// 
	for (i=0; i<lvl; i++)
	{
		x3[i] = x2[i];
	}
	// coordinates for point x2 and x3 in order 24 quadrature
	x2[lvl] = a + h*2/6;
	x3[lvl] = a + h*4/6;
	// if x[lvl] is the last coordinate, integrate using 1D routine
	if (lvl == dim-1)
	{
		y2 = f(x2);
		y3 = f(x3);
		Q = adapt_nd_speclim_last(f,a,b,dim,x2,y2,y3,acc,eps,err,0);
		return Q;
	}
	// if not, integrate using adapt_nd_speclim_now routine for arbitrary levels
	// estimate y2 and y3 using this routine for next coordinate
	else
	{
		acc_next = acc/2;
		y2 = adapt_nd_speclim_next(f,d,u,dim,lvl+1,x2,acc_next,eps,&err1);
		y3 = adapt_nd_speclim_next(f,d,u,dim,lvl+1,x3,acc_next,eps,&err2);
		Q = adapt_nd_speclim_now(f,a,b,d,u,dim,lvl,x2,y2,y3,acc,eps,err,0);
		// add errors quadratically
		*err = sqrt((*err)*(*err)+err1*err1+err2*err2);
		return Q;
	}
}

// this initiates integration of multidimensional function f in the iterated integral where the limits of the different
// coordinates are functions of all previous coordinates.
// d(x,p) yields lower limit for the p'th coordinate for a given set of coordinates x[i], i=0,...,p-1
// similarly for u(x,p) which yields the upper limit
// dim is dimension of the integral and acc and eps are absolute and relative tolerance
// err is pointer to where the estimated error is stored
double adapt_nd_speclim(double f(double* x), double d(double* x, int p), double u (double* x, int p), int dim, double acc, double eps, double *err)
{
	double x[dim];
	// if lower and upper limit are equal in first coordinate integal is zero
	if(u(x,0) == d(x,0))
	{
		*err = 0;
		return 0;
	}
	// initiate integration of zeroth coordinate
	else
	{
		double Q;
		Q = adapt_nd_speclim_next(f,d,u,dim,0,x,acc,eps,err);
		return Q;
	}
}
