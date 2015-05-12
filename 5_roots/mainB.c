#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include "newton.h"

#define TOL	1e-9


void f1(gsl_vector* x0, gsl_vector* fx)
{
	double 	A = 1.0e4,
			x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_vector_set(fx,0,A*x*y-1);
	gsl_vector_set(fx,1,exp(-x)+exp(-y)-1-1/A);
	return;
}

void df1(gsl_vector* x0, gsl_matrix* J)
{
	double	A = 10000,
			x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_matrix_set(J,0,0,A*y);
	gsl_matrix_set(J,0,1,A*x);
	gsl_matrix_set(J,1,0,-exp(-x));
	gsl_matrix_set(J,1,1,-exp(-y));
	return;
}

void f2(gsl_vector* x0, gsl_vector* fx)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_vector_set(fx,0,2*(x-1)-400*x*(y-x*x));
	gsl_vector_set(fx,1,200*(y-x*x));
	return;
}

void df2(gsl_vector* x0, gsl_matrix* J)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_matrix_set(J,0,0,2-400*y+1200*x*x);
	gsl_matrix_set(J,0,1,-400*x);
	gsl_matrix_set(J,1,0,-400*x);
	gsl_matrix_set(J,1,1,200);
}

void f3(gsl_vector* x0, gsl_vector* fx)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_vector_set(fx,0,4*(x*x+y-11)*x+2*(x+y*y-7));
	gsl_vector_set(fx,1,2*(x*x+y-11)+4*(x+y*y-7)*y);
	return;
}

void df3(gsl_vector* x0, gsl_matrix* J)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1);
	gsl_matrix_set(J,0,0,12*x*x+4*y-42);
	gsl_matrix_set(J,0,1,4*x+4*y);
	gsl_matrix_set(J,1,0,4*x+4*y);
	gsl_matrix_set(J,1,1,-26+4*x+12*y*y);
}

void f4(gsl_vector* x0, gsl_vector* fx)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1),
			z = gsl_vector_get(x0,2);
	gsl_vector_set(fx,0,exp(x)+exp(y/z)-1);
	gsl_vector_set(fx,1,sin(y)-cos(y)*z),
	gsl_vector_set(fx,2,x*y*z-100);
	return;
}

void df4(gsl_vector* x0, gsl_matrix* J)
{
	double 	x = gsl_vector_get(x0,0),
			y = gsl_vector_get(x0,1),
			z = gsl_vector_get(x0,2);
	gsl_matrix_set(J,0,0,exp(x));
	gsl_matrix_set(J,0,1,1/z*exp(y/z));
	gsl_matrix_set(J,0,2,-y*exp(y/z)/z/z);
	gsl_matrix_set(J,1,0,0);
	gsl_matrix_set(J,1,1,cos(y)+sin(y)*z);
	gsl_matrix_set(J,1,2,-cos(y));
	gsl_matrix_set(J,2,0,y*z);
	gsl_matrix_set(J,2,1,x*z);
	gsl_matrix_set(J,2,2,x*y);
}

int main(void)
{
	int n = 2, calls;
	gsl_vector* x0 = gsl_vector_alloc(n);
	gsl_vector* fx = gsl_vector_alloc(n);
	gsl_vector_set(x0,0,-2);
	gsl_vector_set(x0,1,8);
	newton_workspace* W = newton_workspace_alloc(n);

	calls = newton_derivative(f1,df1,x0,TOL,W);
	f1(x0,fx);
	fprintf(stdout,"Solving the system of equations\n");
	fprintf(stdout,"\tAxy = 1\n\texp(-x)+exp(-y)=1+1/A\n");
	fprintf(stdout,"where A = 10000\n");
	fprintf(stdout,"calls:\t%d\n",calls);
	fprintf(stdout,"(x,y) = (%g,\t%g),\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1));
	fprintf(stdout,"f(x,y) = (%g,\t%g)\n",gsl_vector_get(fx,0),gsl_vector_get(fx,1));

	gsl_vector_set(x0,0,-2);
	gsl_vector_set(x0,1,8);
	calls = newton_derivative(f2,df2,x0,TOL,W);
	f2(x0,fx);
	fprintf(stdout,"\nFinding the minimum of Rosenbrock's valley function\n");
	fprintf(stdout,"\tf(x,y) = (1-x)^2+100*(y-x^2)^2\n");
	fprintf(stdout,"calls:\t%d\n",calls);
	fprintf(stdout,"(x,y) = (%g,\t%g),\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1));
	fprintf(stdout,"f(x,y) = (%g,\t%g)\n",gsl_vector_get(fx,0),gsl_vector_get(fx,1));

	gsl_vector_set(x0,0,-2);
	gsl_vector_set(x0,1,8);
	calls = newton_derivative(f3,df3,x0,TOL,W);
	f3(x0,fx);
	fprintf(stdout,"\nFinding the minimum of Himmelblau's function\n");
	fprintf(stdout,"\tf(x,y) = (x^2+y-11)^2+(x+y^2-7)^2\n");
	fprintf(stdout,"calls:\t%d\n",calls);
	fprintf(stdout,"(x,y) = (%g,\t%g),\n",gsl_vector_get(x0,0),gsl_vector_get(x0,1));
	fprintf(stdout,"f(x,y) = (%g,\t%g)\n",gsl_vector_get(fx,0),gsl_vector_get(fx,1));

	gsl_vector* x1 = gsl_vector_alloc(3);
	gsl_vector* fx1 = gsl_vector_alloc(3);
	newton_workspace* W1 = newton_workspace_alloc(3);
	gsl_vector_set(x1,0,-2);
	gsl_vector_set(x1,1,8);
	gsl_vector_set(x1,2,2);
	fprintf(stdout,"\nFinding the solution to the system of equation\n");
	fprintf(stdout,"\texp(x)+exp(y/z) = 1\n");
	fprintf(stdout,"\tsin(y) = cos(x)*z\n");
	fprintf(stdout,"\tx*y*z = 100\n");
	fprintf(stdout,"Initial guess\n");
	f4(x1,fx1);
	fprintf(stdout,"(x,y) = (%g,\t%g,\t%g),\n",gsl_vector_get(x1,0),gsl_vector_get(x1,1),gsl_vector_get(x1,2));
	fprintf(stdout,"f(x,y) = (%g,\t%g,\t%g)\n",gsl_vector_get(fx1,0),gsl_vector_get(fx1,1),gsl_vector_get(fx1,2));
	calls = newton_derivative(f4,df4,x1,TOL,W1);
	f4(x1,fx1);fprintf(stdout,"calls:\t%d\n",calls);
	fprintf(stdout,"(x,y) = (%g,\t%g,\t%g),\n",gsl_vector_get(x1,0),gsl_vector_get(x1,1),gsl_vector_get(x1,2));
	fprintf(stdout,"f(x,y) = (%g,\t%g,\t%g)\n",gsl_vector_get(fx1,0),gsl_vector_get(fx1,1),gsl_vector_get(fx1,2));

	gsl_vector_free(x0);
	gsl_vector_free(fx);
	gsl_vector_free(x1);
	gsl_vector_free(fx1);
	newton_workspace_free(W1);
	newton_workspace_free(W);

	return 0;
}
