#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QASC_MAX_DEPTH	100000

typedef struct
{
	int lvl;
	double a, b, y1, y3, acc;
	void* next;
} stack;

double qasc23(double f(double t), double a, double b, double acc, double eps, double *err)
{
	double h, y2, Q, q, diff, sqrt2 = sqrt(2);
	// long double is neccesary due to the nonrecursive algorithm
	// adding every contribution directly to I, thus loosing precision
	// from the machine epsilon
	long double I = 0;
	*err = 0;
	stack* now = (stack*) malloc(sizeof(stack));
	stack* next;
	stack* done;
	now->lvl = 0;
	now->a = a;
	now->b = b;
	now->y1 = f(a);
	now->y3 = f(b);
	now->acc = acc;
	now->next = NULL;
	do {
		if (now->lvl < QASC_MAX_DEPTH)
		{
			h = now->b - now->a;
			y2 = f(now->a+0.5*h);
			Q = h*0.25*(now->y1+2*y2+now->y3);
			q = h*0.5*(now->y1+now->y3);
			diff = fabs(Q-q);

			if (diff < now->acc + eps*fabs(Q))
			{
				I += Q;
				*err += diff;
				done = now;
				now = (stack*) done->next;
				free(done);
			}
			else
			{
				next = (stack*) malloc(sizeof(stack));
				memcpy(next,now,sizeof(stack));

				now->lvl++;
				now->b = now->a+h*0.5;
				now->y3 = y2;
				now->acc /= sqrt2;
				now->next = (void*) next;

				next->lvl = now->lvl;
				next->a = now->b;
				next->y1 = y2;
				next->acc = now->acc;
			}
		}
		else
		{
			fprintf(stderr,"Too many subdivisions during integration!\nTerminating integrator\t...\t");
			while (now != NULL)
			{
				done = now;
				now = (stack*) done->next;
				free(done);
			}
			fprintf(stderr,"done!\n");
		}
	} while ( now != NULL);
	return I;
}

// note that due to this algorithm using closed intervals, it cannot apply transformations
// for integrals with infinities in the limits
double qasc(double f(double), double a, double b, double acc, double eps, double * err)
{
	if (a < b)
	{
		double Q = qasc23(f,a,b,acc,eps,err);
		return Q;
	}
	else if (b < a)
	{
		return -qasc(f,b,a,acc,eps,err);
	}
	else if (a == b)
	{
		*err = 0;
		return 0;
	}
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
