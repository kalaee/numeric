#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>
#include "optim.h"

#define	TOL	1e-9
#define ALPHA	1e-4

// Rosenbrock valley function
double f1(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return pow(1-x,2) + 100*pow(y-x*x,2);
}

// Gradient of Rosenbrock
void Gf1(gsl_vector* x0, gsl_vector* df)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);

	gsl_vector_set(df,0,2*(x-1-200*y*x+200*x*x*x));
	gsl_vector_set(df,1,200*(y-x*x));
	return;
}

// Hessian matrix of Rosenbrock
void Hf1(gsl_vector* x0, gsl_matrix* H)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);

	gsl_matrix_set(H,0,0,-400*y+1200*x*x+2);
	gsl_matrix_set(H,0,1,-400*x);
	gsl_matrix_set(H,1,0,gsl_matrix_get(H,0,1));
	gsl_matrix_set(H,1,1,200);
	return;
}

// Rosenbrock valley function
double f2(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return pow(x*x+y-11,2)+pow(x+y*y-7,2);
}

// Gradient of Rosenbrock
void Gf2(gsl_vector* x0, gsl_vector* df)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);

	gsl_vector_set(df,0,2*(y*y+2*x*(y+x*x-11)+x-7));
	gsl_vector_set(df,1,2*(2*y*(y*y+x-7)+y+x*x-11));
	return;
}

// Hessian matrix of Rosenbrock
void Hf2(gsl_vector* x0, gsl_matrix* H)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);

	gsl_matrix_set(H,0,0,4*y+12*x*x-42);
	gsl_matrix_set(H,0,1,4*(y+x));
	gsl_matrix_set(H,1,0,gsl_matrix_get(H,0,1));
	gsl_matrix_set(H,1,1,12*y*y+4*x-26);
	return;
}


int main(void)
{
	int steps;
	gsl_vector* x = gsl_vector_alloc(2);
	min_newton_workspace* W = min_newton_workspace_alloc(2);

	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	fprintf(stdout,"Rosenbrock's valley function:\n");
	fprintf(stdout,"Initial guess:\n");
	fprintf(stdout,"x =\t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	fprintf(stdout,"f(x) =\t%g\n",f1(x));
	steps = min_newton(f1,Gf1,Hf1,x,ALPHA,TOL,W);
	fprintf(stdout,"After %d steps:\n",steps);
	fprintf(stdout,"x = \t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	fprintf(stdout,"f(x) = \t%g\n",f1(x));

	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	fprintf(stdout,"\nHimmelblau's function:\n");
	fprintf(stdout,"Initial guess:\n");
	fprintf(stdout,"x =\t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	fprintf(stdout,"f(x) =\t%g\n",f2(x));
	steps = min_newton(f2,Gf2,Hf2,x,ALPHA,TOL,W);
	fprintf(stdout,"After %d steps:\n",steps);
	fprintf(stdout,"x = \t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	fprintf(stdout,"f(x) = \t%g\n",f2(x));

	gsl_vector_free(x);
	min_newton_workspace_free(W);
	return 0;
}
