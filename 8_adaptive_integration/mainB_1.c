#include <math.h>
#include <stdio.h>
#include "qarc.h"
#include "qaro.h"

#define ACC	1e-12
#define EPS	1e-12

int counter;
double f1(double t)
{
	counter++;
	return exp(-t);
}

double f2(double t)
{
	counter++;
	return t*t*exp(-fabs(t));
}

double f3(double t)
{
	counter++;
	return 1/sqrt(1+pow(t,4));
}

int main(void)
{
	double err, Q, trueval;
	counter = 0;
	trueval = 1;
	fprintf(stdout,"int_0^INF exp(t)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qaro(f1,0,INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	trueval = -2;
	fprintf(stdout,"\n\nint_0^-INF t^2*exp(t)\n");
	fprintf(stdout,"True value:\t%g\n",trueval);
	Q = qaro(f2,0,-INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);

	counter = 0;
	trueval = 3.708149354602743836867700694390520092435197647043533811171;
	fprintf(stdout,"\n\nint_-INF^INF 1/sqrt(1+t^4)\n");
	fprintf(stdout,"True value:\t%.16g\n",trueval);
	Q = qaro(f3,-INFINITY,INFINITY,ACC,EPS,&err);
	fprintf(stdout,"============================\n");
	fprintf(stdout,"QARO estimate: %g\nError estimate: %g\n",Q,err);
	fprintf(stdout,"Actual error:\t%g\n",Q -trueval);
	fprintf(stdout,"Counts:\t%d\n",counter);


	return 0;
}
