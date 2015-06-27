#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "qarc.h"
#include "qaro.h"

// allowed absolute and relative errors
#define ACC	1e-12
#define EPS	1e-12
#define ACC1 1e-14
#define EPS1 1e-13

// int which counts number of calls to function
int counter;

// the different integrands
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
// function 4 again, because of GSL
double f4gsl(double t, void* param)
{	
	counter++;
	return 4*sqrt(1-pow(1-t,2));
}

int main(void)
{
	double err, Q, trueval;
	// perform integration and output estimate
	// together with error and calls to function

	// f1
	counter = 0;
	trueval = log(2)*M_PI/8;
	fprintf(stdout,"int_0^{pi/4} log(1+tan(t))\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qarc(f1,0,M_PI/4,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaro(f1,0,M_PI/4,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QARO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	// f2
	counter = 0;
	trueval = 2*M_PI/sqrt(3);
	fprintf(stdout,"\n\nint_0^2*PI 1/(2+cos(t))\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	fprintf(stdout,"============================\n");
	Q = qarc(f2,0,M_PI*2,ACC,EPS,&err);
	fprintf(stdout,"QARC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaro(f2,0,M_PI*2,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QARO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	// f3
	counter = 0;
	trueval = M_PI*0.25*0.25/(1-0.25*0.25);
	fprintf(stdout,"\n\nint_0^PI cos(t)/(1+0.5*cos(t)+0.25^2)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	fprintf(stdout,"============================\n");
	Q = qarc(f3,0,M_PI,ACC,EPS,&err);
	fprintf(stdout,"QARC: %g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaro(f3,0,M_PI,ACC,EPS,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QARO: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	// f4
	counter = 0;
	trueval = M_PI;
	fprintf(stdout,"\n\nint_0^1 4*sqrt(1-(1-x)^2)\n");
	fprintf(stdout,"True value:\t%.16g\n",trueval);
	fprintf(stdout,"Tolerance: ACC = %g, EPS = %g\n",ACC1,EPS1);
	fprintf(stdout,"============================\n");
	Q = qarc(f4,0,1,ACC1,EPS1,&err);
	fprintf(stdout,"QARC: %.16g\t\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q-trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	counter = 0;
	Q = qaro(f4,0,1,ACC1,EPS1,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"QARO: %.16g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);
	
	// f4 using GSL
	counter = 0;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_function F;
	F.function = &f4gsl;
	gsl_integration_qag(&F,0,1,1e-14,1e-13,1000,GSL_INTEG_GAUSS61,w,&Q,&err);
	fprintf(stdout,"----------------------------\n");
	fprintf(stdout,"GSL QAG: %.16g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	gsl_integration_workspace_free(w);

	return 0;
}
