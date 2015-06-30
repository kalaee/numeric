#include <stdio.h>
#include <math.h>
#include "adapt_2d_speclim.h"

#define ACC	1e-6
#define EPS	1e-6

int counter;

double f1(double x, double y)
{
	counter++;
	return exp(y-x);
}

double d1(double x)
{
	return sin(x)-0.5;
}

double u1(double x)
{
	return cos(0.5*x*x)+0.5;
}


double f2(double x, double y)
{
	counter++;
	return sin(x+y);
}

double d2(double x)
{
	return 1/cos(x)-1;
}

double u2(double x)
{
	if(x == 0)
	{
		return 1;
	}
	else
	{
		return x/tan(x);
	}
}



int main(void)
{
	double exact, estim, err;
	fprintf(stdout,"Tolerances\nACC:\t%g\nEPS:\t%g\n\n",ACC,EPS);


	fprintf(stdout,"Integrand: exp(y-x),\nx from 0 to 2*pi, y from d(x) to u(x),\n");
	fprintf(stdout,"where d(x) = sin(x)-0.5 and u(x) = cos(0.5*x^2)+0.5\n");
	exact = 2.55398;
	fprintf(stdout,"Actual value (WolframAlpha):\t%g\n",exact);
	counter = 0;
	estim = adapt_2d_speclim(f1,0,2*M_PI,d1,u1,ACC,EPS,&err);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);


	fprintf(stdout,"Integrand: sin(x+y),\nx from -pi/4 to pi/4, y from d(x) to u(x),\n");
	fprintf(stdout,"where d(x) = 1/cos(x)-1 and u(x) = x/tan(x)\n");
	exact = 0.55951;
	fprintf(stdout,"Actual value (WolframAlpha):\t%g\n",exact);
	counter = 0;
	estim = adapt_2d_speclim(f2,-M_PI/4.,M_PI/4.,d2,u2,ACC,EPS,&err);
	fprintf(stdout,"Estim:\t%g\n",estim);
	fprintf(stdout,"Estim. err:\t%g\n",err);
	fprintf(stdout,"Actual err:\t%g\n",estim-exact);
	fprintf(stdout,"Calls to function:\t%d\n\n",counter);




	return 0;
}
