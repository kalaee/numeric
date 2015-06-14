#include <math.h>
#include <assert.h>
#include <stdio.h>
// Quadratic Adaptive integral using Recursive Closed intervals of order 3 with 2 for error estimate
double qarc23(double f(double x), double a, double b, double y1, double y3, double acc, double eps, double *err, int nrecur)
{
	if (nrecur > 100000)
	{
		fprintf(stderr,"Too many subdivisions!\n");
		return 0;
	}
	double h = b-a;
	double y2 = f(a+h/2.);
	double Q = h/4.*(y1+2.*y2+y3);
	double q = h/2.*(y1+y3);
	*err = fabs(Q-q);
	if ( *err < acc + eps*fabs(Q) )
	{
		return Q;
	}
	else
	{
		acc /= sqrt(2);
		double err1, err2;
		double Q1 = qarc23(f,a,a+h/2,y1,y2,acc,eps,&err1,nrecur+1);
		double Q2 = qarc23(f,a+h/2,b,y2,y3,acc,eps,&err2,nrecur+1);
		*err = err1 + err2;
		return Q1 + Q2;
	}
}

// note that due to this algorithm using closed intervals, it cannot apply transformations
// for integrals with infinities in the limits
double qarc(double f(double), double a, double b, double acc, double eps, double * err)
{
	if (a < b)
	{
		double y1 = f(a);
		double y3 = f(b);
		double Q = qarc23(f,a,b,y1,y3,acc,eps,err,0);
		return Q;		
	}
	else if (b < a)
	{
		return -qarc(f,b,a,acc,eps,err);
	}
	else if (a == b)
	{
		*err = 0;
		return 0;
	}
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
