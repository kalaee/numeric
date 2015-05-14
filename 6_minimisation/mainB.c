#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <stdio.h>
#include "min_newton.h"
#include "min_newton_sr1.h"
#include "../libs/roots_newton.h"

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

// Three hump camel function (THC)
double f3(gsl_vector* x0)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	return 2*x*x - 1.05*pow(x,4) + pow(x,6)/6+x*y+y*y;
}

// Gradient of THC
void Gf3(gsl_vector* x0, gsl_vector* df)
{
	double	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_vector_set(df,0,y + pow(x,5)-4.2*pow(x,3)+4*x);
	gsl_vector_set(df,1,2*y+x);
	return;
}

// Hessian matrix of THC
void Hf3(gsl_vector* x0, gsl_matrix* H)
{
	double	x = gsl_vector_get(x0,0);
	gsl_matrix_set(H,0,0,5*pow(x,4)-12.6*pow(x,2)+4);
	gsl_matrix_set(H,0,1,1);
	gsl_matrix_set(H,1,0,1);
	gsl_matrix_set(H,1,1,2);
	return;
}

int main(void)
{
	int steps;
	gsl_vector* x = gsl_vector_alloc(2);
	min_newton_sr1_workspace* WSR1 = min_newton_sr1_workspace_alloc(2);
	min_newton_workspace* W = min_newton_workspace_alloc(2);
	roots_newton_workspace* Wroot = roots_newton_workspace_alloc(2);

	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	fprintf(stdout,"\nRosenbrock's function:\n");
	fprintf(stdout,"Initial guess:\n");
	fprintf(stdout,"x =\t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	steps = min_newton(f1,Gf1,Hf1,x,ALPHA,TOL,W);
	fprintf(stdout,"True Solution\nx =\t%g,\t%g\n",1.,1.);
	fprintf(stdout,"Newton's minimisation method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	steps = min_newton_sr1(f1,Gf1,x,ALPHA,TOL,WSR1);
	fprintf(stdout,"Quasi-Newton's with SR1:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	steps = roots_newton(Gf1,x,1e-6,TOL,Wroot);
	fprintf(stdout,"Newton's root-finding method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));


	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	fprintf(stdout,"\nHimmelblau's function:\n");
	fprintf(stdout,"Initial guess:\n");
	fprintf(stdout,"x =\t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	steps = min_newton(f2,Gf2,Hf2,x,ALPHA,TOL,W);
	fprintf(stdout,"True Solution\nx =\t%g,\t%g\n",-2.805118,3.131312);
	fprintf(stdout,"Note that the function has more global minima,\n");
	fprintf(stdout,"yet, the above solution appears to be the one relevant from our initial guess\n");
	fprintf(stdout,"Newton's minimisation method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	steps = min_newton_sr1(f2,Gf2,x,ALPHA,TOL,WSR1);
	fprintf(stdout,"Quasi-Newton's with SR1:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-3);
	gsl_vector_set(x,1,10);
	steps = roots_newton(Gf2,x,1e-6,TOL,Wroot);
	fprintf(stdout,"Newton's root-finding method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));


	gsl_vector_set(x,0,-0.5);
	gsl_vector_set(x,1,10);
	fprintf(stdout,"\nThree-Hump Camel function (THC):\n");
	fprintf(stdout,"Initial guess:\n");
	fprintf(stdout,"x =\t%g,\t%g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));
	steps = min_newton(f3,Gf3,Hf3,x,ALPHA,TOL,W);
	fprintf(stdout,"True Solution\nx =\t%g,\t%g\n",0.,0.);
	fprintf(stdout,"Note that although THC has one global minimum it also has two other local minima\n");
	fprintf(stdout,"close to the global minimum. Depending on the initial guess the algorithms might\n");
	fprintf(stdout,"find these instead and terminate.\n");
	fprintf(stdout,"Newton's minimisation method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-0.5);
	gsl_vector_set(x,1,10);
	steps = min_newton_sr1(f3,Gf3,x,ALPHA,TOL,WSR1);
	fprintf(stdout,"Quasi-Newton's with SR1:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));
	gsl_vector_set(x,0,-0.5);
	gsl_vector_set(x,1,10);
	steps = roots_newton(Gf3,x,1e-6,TOL,Wroot);
	fprintf(stdout,"Newton's root-finding method:\n\tsteps=%d\n\tx= %g\t%g\n",steps,gsl_vector_get(x,0),gsl_vector_get(x,1));

	gsl_vector_free(x);
	min_newton_sr1_workspace_free(WSR1);
	min_newton_workspace_free(W);
	roots_newton_workspace_free(Wroot);
	return 0;
}
