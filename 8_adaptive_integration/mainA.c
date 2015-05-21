#include <math.h>
#include <stdio.h>
#include "qarc.h"
#include "qaro.h"

#define ACC	1e-12
#define EPS	1e-12

#define ACC1 1e-19
#define EPS1 1e-19
#define ACC2 1e-11
#define EPS2 1e-11
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

int main(void)
{
	double err, Q, trueval;
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


	counter = 0;
	trueval = M_PI;
	fprintf(stdout,"\n\nint_0^1 4*sqrt(1-(1-x)^2)\n");
	fprintf(stdout,"True value:\t%.16g\n",trueval);
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


	return 0;
}
