#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// maximal iteration depth allowed for the nodes
#define QASC_MAX_DEPTH	100000

// the node structure represent data about an interval
// lvl is recursion level of the node
// a is first coordinate, b last coordinate,
// y1 and y3 are known values of the integrand at a and b
// acc is allowed absolute error for the interval
// *next is pointer to the node representing next
// part of the whole interval of integration
typedef struct
{
	int lvl;
	double a, b, y1, y3, acc;
	void* next;
} node;

// qasc23 performs integration using closed interval and trapezoid
// evaluations of third and second order
// f is the integrand, a and b the start and end points of the interval,
// acc and eps are the allowed absolute and relative errors,
// *err is pointer to where the error estimate is stored
double qasc23(double f(double t), double a, double b, double acc, double eps, double *err)
{
	double h, y2, Q, q, diff, sqrt2 = sqrt(2);
	// long double is neccesary due to the nonrecursive algorithm
	// adding every contribution directly to I, thus loosing precision
	// from the machine epsilon
	long double I = 0;
	*err = 0;

	// allocate the memory for the node representing at first the whole interval
	// and prepare the node
	node* now = (node*) malloc(sizeof(node));
	node* next;
	node* done;
	now->lvl = 0;
	now->a = a;
	now->b = b;
	now->y1 = f(a);
	now->y3 = f(b);
	now->acc = acc;
	// the last node will always point to NULL
	now->next = NULL;

	// the do-while loop continues until the integral value
	// for the node with NULL pointer has been evaluated
	do {
		// check that the node level is not too deep
		if (now->lvl < QASC_MAX_DEPTH)
		{
			// estimate integral and error, saved in diff
			h = now->b - now->a;
			y2 = f(now->a+0.5*h);
			Q = h*0.25*(now->y1+2*y2+now->y3);
			q = h*0.5*(now->y1+now->y3);
			diff = fabs(Q-q);

			// if diff not too big, add integral value to I,
			// free memory for evaluated node and continue to
			// next node in line
			if (diff < now->acc + eps*fabs(Q))
			{
				I += Q;
				*err += diff;
				done = now;
				now = (node*) done->next;
				free(done);
			}
			// if diff too big
			else
			{
				// allocte memory for new node
				next = (node*) malloc(sizeof(node));
				// copy information from old node to new
				memcpy(next,now,sizeof(node));
				// adjust the information of the now to represent the
				// first half of the interval
				now->lvl++;
				now->b = now->a+h*0.5;
				now->y3 = y2;
				now->acc /= sqrt2;
				now->next = (void*) next;
				// adjust the information of the next to represent second
				// half of the interval
				next->lvl = now->lvl;
				next->a = now->b;
				next->y1 = y2;
				next->acc = now->acc;
			}
		}
		// if node level too deep, assume integration is not possible,
		// send error message and free memory of unevaluated nodes
		else
		{
			fprintf(stderr,"Too many subdivisions during integration!\nTerminating integrator\t...\t");
			// free memory until null pointer is reached
			while (now != NULL)
			{
				done = now;
				now = (node*) done->next;
				free(done);
			}
			fprintf(stderr,"done!\n");
		}
	// terminate when null pointer is reached
	} while ( now != NULL);
	return I;
}

// this function takes the integration arguments and prepares them for qasc23
// f is the integrand, a and b the start and end points of the interval,
// acc and eps are the allowed absolute and relative errors,
// *err is pointer to where the error estimate is stored
// note that due to this algorithm using closed intervals, it cannot apply transformations
// for integrals with infinities in the limits
double qasc(double f(double), double a, double b, double acc, double eps, double * err)
{
	// check integral limits, if a < b, perform integration
	if (a < b)
	{
		double Q = qasc23(f,a,b,acc,eps,err);
		return Q;
	}
	// if b < a, return negative of integral from b to a
	else if (b < a)
	{
		return -qasc(f,b,a,acc,eps,err);
	}
	// if a == b, integral is by definition zero
	else if (a == b)
	{
		*err = 0;
		return 0;
	}
	// if nothing fits, send error message about intervals
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
