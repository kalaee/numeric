#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "qasc.h"
#include "qaso.h"
#include "qaro.h"

#define ACC	1e-12
#define EPS	1e-12

#define ACC1 1e-14
#define EPS1 1e-13
int counter;
// integrand f(x) = log(1+tan(t))
double f1(double t)
{
	counter++;
	return log(1+tan(t));
}

double f2(double t)
{
	counter++;
	return 1/(2+cos(t));
}

double f3(double t)
{
	counter++;
	return cos(2*t)/(1+0.5*cos(t)+0.25*0.25);
}

double f4(double t)
{
	counter++;
	return 4*sqrt(1-pow(1-t,2));
}

double f4gsl(double t, void* param)
{	
	counter++;
	return 4*sqrt(1-pow(1-t,2));
}

double f5(double t)
{
	counter++;
	return exp(-t);
}

double f6(double t)
{
	counter++;
	return t*t*exp(-fabs(t));
}

double f7(double t)
{
	counter++;
	return 1/sqrt(1+pow(t,4));
}

int main(void)
{
	double err, Q, trueval;
	counter = 0;
	trueval = log(2)*M_PI/8;
	fprintf(stdout,"|============================|\n");
	fprintf(stdout,"|   Testing QASC and QASO    |\n");
	fprintf(stdout,"|============================|\n");

	fprintf(stdout,"\nint_0^{pi/4} log(1+tan(t))\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qasc(f1,0,M_PI/4,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QASC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaso(f1,0,M_PI/4,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QASO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	trueval = 2*M_PI/sqrt(3);
	fprintf(stdout,"\n\nint_0^2*PI 1/(2+cos(t))\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	fprintf(stdout,"============================\n");
	Q = qasc(f2,0,M_PI*2,ACC,EPS,&err);
	fprintf(stdout,"QASC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaso(f2,0,M_PI*2,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QASO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	
	counter = 0;
	trueval = M_PI*0.25*0.25/(1-0.25*0.25);
	fprintf(stdout,"\n\nint_0^PI cos(t)/(1+0.5*cos(t)+0.25^2)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	fprintf(stdout,"============================\n");
	Q = qasc(f3,0,M_PI,ACC,EPS,&err);
	fprintf(stdout,"QASC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaso(f3,0,M_PI,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QASO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);


	counter = 0;
	trueval = M_PI;
	fprintf(stdout,"\n\nint_0^1 4*sqrt(1-(1-x)^2)\n");
	fprintf(stdout,"True value:\t%.16g\n",trueval);
	fprintf(stdout,"Tolerance: ACC = %g, EPS = %g\n",ACC1,EPS1);
	fprintf(stdout,"============================\n");
	Q = qasc(f4,0,1,ACC1,EPS1,&err);
	fprintf(stdout,"QASC: %.16g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaso(f4,0,1,ACC1,EPS1,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QASO: %.16g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &f4gsl;
	gsl_integration_qag(&F,0,1,1e-14,1e-13,1000,GSL_INTEG_GAUSS61,w,&Q,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"GSL QAG: %.16g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	fprintf(stdout,"\n|============================|\n");
	fprintf(stdout,"|   Limits with infinities   |\n");
	fprintf(stdout,"|============================|\n");

	counter = 0;
	trueval = 1;
	fprintf(stdout,"\nint_0^INF exp(t)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qaro(f5,0,INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	trueval = -2;
	fprintf(stdout,"\n\nint_0^-INF t^2*exp(t)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qaso(f6,0,-INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QASO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	trueval = 3.708149354602743836867700694390520092435197647043533811171;
	fprintf(stdout,"\n\nint_-INF^INF 1/sqrt(1+t^4)\n");
	fprintf(stdout,"True value:\t%.16g\n",trueval);
	Q = qaro(f7,-INFINITY,INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	gsl_integration_workspace_free(w);

	return 0;
}
