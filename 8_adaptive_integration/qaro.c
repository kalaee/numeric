#include <math.h>
#include <assert.h>
#include <stdio.h>

// Quadratic Adaptive integral using Recursive functions and Open intervals of order 4 with 2 for error estimate
double qaro24(double f(double x), double a, double b, double y2, double y3, double acc, double eps, double *err, int nrecur)
{
	if (nrecur > 100000)
	{
		fprintf(stderr,"Too many subdivisions!\n");
		return 0;
	}
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
	if (a < b)
	{
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
	else if (b < a)
	{
		// reverse order of limits, caro24 and infinity-evalutaion is designed for a < b
		return -qaro(f,b,a,acc,eps,err);
	}
	else if (a == b )
	{
		// integral is by definition zero
		*err = 0;
		return 0;
	}
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
