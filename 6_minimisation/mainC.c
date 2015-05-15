#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>
#include "min_newton.h"
#include "min_newton_sr1.h"
#include "../libs/roots_newton.h"
#include "simplex.h"

#define	TOL	1e-9
#define ALPHA	1e-4

// Rosenbrock valley function
double f1(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return pow(1-x,2) + 100*pow(y-x*x,2);
}
// Himmelblau function
double f2(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}
// Three hump camel function (THC)
double f3(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return 2*x*x - 1.05*pow(x,4) + pow(x,6)/6+x*y+y*y;
}
int main(void)
{
	int steps;
	simplex_workspace* W = simplex_workspace_alloc(2);
	steps = simplex(f1,-10,10,TOL,W);
	fprintf(stdout,"Minimising Rosenbrock's valley function\n");
	fprintf(stdout,"Domain of initial position: -10 < x,y < 10\n");
	fprintf(stdout,"After %d steps the function is minimised:\n",steps);
	fprintf(stdout,"(x,y) = (%g\t%g)\n",gsl_vector_get(W->ce,0),gsl_vector_get(W->ce,1));
	fprintf(stdout,"f(x,y) = %g\n",f1(W->ce));

	steps = simplex(f2,-10,10,TOL,W);
	fprintf(stdout,"\nMinimising Himmelblau's function\n");
	fprintf(stdout,"Domain of initial position: -10 < x,y < 10\n");
	fprintf(stdout,"After %d steps the function is minimised:\n",steps);
	fprintf(stdout,"(x,y) = (%g\t%g)\n",gsl_vector_get(W->ce,0),gsl_vector_get(W->ce,1));
	fprintf(stdout,"f(x,y) = %g\n",f2(W->ce));
	
	steps = simplex(f3,-5,5,TOL,W);
	fprintf(stdout,"\nMinimising Three-Hump Camel function\n");
	fprintf(stdout,"Domain of initial position: -5 < x,y < 5\n");
	fprintf(stdout,"After %d steps the function is minimised:\n",steps);
	fprintf(stdout,"(x,y) = (%g\t%g)\n",gsl_vector_get(W->ce,0),gsl_vector_get(W->ce,1));
	fprintf(stdout,"f(x,y) = %g\n",f3(W->ce));
	fprintf(stdout,"Note that although THC has a global minimum at (x,y) = (0,0)\n");
	fprintf(stdout,"it also has two local minima close the global one, which might be\n");
	fprintf(stdout,"reached instead.\n");
	simplex_workspace_free(W);
	return 0;
}
