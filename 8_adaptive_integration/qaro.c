#include <math.h>
#include <assert.h>
#include <stdio.h>
// Quadratic Adaptive integral using Recursive functions and Open intervals of order 4 with 2 for error estimate
double qaro24(double f(double x), double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	double h = b-a;
	double y1 = f(a+h/6);
	double y4 = f(a+h*5/6);
	double Q = h/6*(2*y1+y2+y3+2*y4);
	double q = h/2*(y1+y4);
	*err = fabs(Q-q);
	if ( *err < acc + eps*fabs(Q) )
	{
		return Q;
	}
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

double qaro(double f(double), double a, double b, double acc, double eps, double * err)
{
	if (a == b)
	{
		*err = 0;
		return 0;
	}
	else
	{
		double h = b - a;
		double y2 = f(a+h*2/6);
		double y3 = f(a+h*4/6);
		double Q = qaro24(f,a,b,y2,y3,acc,eps,err,0);
		return Q;
	}
}
