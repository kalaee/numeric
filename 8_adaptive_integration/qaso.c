#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// maximal iteration depth allowed for the nodes
#define QASO_MAX_DEPTH	100000

// the node structure represents data about an interval
// lvl is recursion level of the node
// a and b are frist and last coordinate of the interval
// y2 and y3 are integrand values at a+2h/6 and b-2h/6,
// where h = b - a
// acc is allowed absolute error and next is pointer to
// the node representing next part of the whole integraion interval
typedef struct
{
	int lvl;
	double a, b, y2, y3, acc;
	void* next;
} stack;


// qaso24 perform integration using open intervals and trapezoid
// evaluations of fourth and second order
// f is the integrand, a and b the start and end of interval
// acc and eps are allowed absolute and relative deviations
// *err is pointer to where estimated error is stored
// note that qaso24 asssumes a < b
double qaso24(double f(double t), double a, double b, double acc, double eps, double* err)
{
	double h, y1, y4, Q, q, diff, sqrt2 = sqrt(2);
	// long double is neccesary due to the nonrecursive algorithm
	// adding every contribution directly to I, thus loosing precision
	// from the machine epsilon
	long double I = 0;

	// allocate memory for the node representing the whole
	// interval and prepare the node
	*err = 0;
	stack* now = (stack*) malloc(sizeof(stack));
	stack* next;
	stack* done;
	now->lvl = 0;
	now->a = a;
	now->b = b;
	h = b-a;
	now->y2 = f(a+h*2/6);
	now->y3 = f(a+h*4/6);
	now->acc = acc;
	// the last node will always point to null
	now->next = NULL;
	
	// the do-while loop continues until the integral value
	// for the node with null pointer has been evaluated
	do
	{
		// check that the node level is not too deep
		if (now->lvl < QASO_MAX_DEPTH)
		{
			// estimate integral and error, stored in diff
			h = (now->b - now->a);
			y1 = f(now->a + h/6);
			y4 = f(now->a + h*5/6);
			Q = h/6*(2*y1+now->y2+now->y3+2*y4);
			q = h/2*(y1+y4);
			diff = fabs(Q-q);
			// if diff not too big, add contribution to I
			// free memory of current node and continue to next node
			if (diff < now->acc + eps*fabs(Q))
			{
				I += Q;
				*err += diff*diff;
				done = now;
				now = (stack*) done->next;
				free(done);
			}
			// if diff too big
			else
			{
				// allocate memory for new node
				next = (stack*) malloc(sizeof(stack));
				// copy information from current node to new
				memcpy(next,now,sizeof(stack));
				// let "now" node represent first half of interval
				now->lvl++;
				now->b = now->a+h/2;
				now->y3 = now->y2;
				now->y2 = y1;
				now->acc /= sqrt2;
				now->next = (void*) next;
				// let next node represent second half of interval
				next->lvl = now->lvl;
				next->a = now->b;
				next->y2 = next->y3;
				next->y3 = y4;
				next->acc = now->acc;
			}
		}
		// if node level too deep, assume integration is not possible
		// send error message to stream stderr and free memory
		else
		{
			fprintf(stderr,"Too many subdivisions during integration!\nTerminating integrator\t...\t");
			// free memory until null pointer is reached
			while(now != NULL)
			{
				done = now;
				now = (stack*) done->next;
				free(done);
			}
			fprintf(stderr,"done!\n");
			exit(EXIT_FAILURE);
		}
	// terminate loop when null pointer is reached
	} while (now != NULL);
	*err = sqrt(*err);
	return I;
}

// qaso takes the integration argument and prepares them for qaso24
// f is the integrand, a and b the start and end of interval
// acc and eps are allowed absolute and relative deviations
// *err is pointer to where estimated error is stored
double qaso(double f(double), double a, double b, double acc, double eps, double * err)
{
	// check that a < b, qaso42 assumes this to be the case
	if (a < b)
	{
		// determine whether the limits contain infinities
		// if yes, assign the appropriate transformation of the integral
		if ( !isinf(a) && !isinf(b) )
		{
			return qaso24(f,a,b,acc,eps,err);
		}
		else if ( !isinf(a) && isinf(b) )
		{
			double ft(double t){ return f(a+t/(1-t))/pow(1-t,2); }
			return qaso24(ft,0,1,acc,eps,err);
		}
		else if ( isinf(a) && !isinf(b) )
		{
			double ft (double t){ return f(b - (1-t)/t)/t/t; }
			return qaso24(ft,0,1,acc,eps,err);
		}
		else if ( isinf(a) && isinf(b) )
		{
			double ft(double t){ return ( f((1-t)/t) + f((t-1)/t))/t/t; }
			return qaso24(ft,0,1,acc,eps,err);
		}
		else
		{
			fprintf(stderr,"Something is wrong with the limits!\n");
			exit(EXIT_FAILURE);
		}
	}
	// if a > b, return negative of integral with limits switched
	else if (b < a)
	{
		return -qaso(f,b,a,acc,eps,err);
	}
	// if a == b, integral is by definition zero
	else if (a == b )
	{
		*err = 0;
		return 0;
	}
	// if nothing fits, send error message about the limits
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
