#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define QASO_MAX_DEPTH	100000

typedef struct
{
	int lvl;
	double a, b, y2, y3, acc;
	void* next;
} stack;

double qaso24(double f(double t), double a, double b, double acc, double eps, double* err)
{
	double h, y1, y4, Q, q, diff, sqrt2 = sqrt(2);
	// long double is neccesary due to the nonrecursive algorithm
	// adding every contribution directly to I, thus loosing precision
	// from the machine epsilonlong double	I = 0;
	long double I = 0;
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
	now->next = NULL;
	do
	{
		if (now->lvl < QASO_MAX_DEPTH)
		{
			h = (now->b - now->a);
			y1 = f(now->a + h/6);
			y4 = f(now->a + h*5/6);
			Q = h/6*(2*y1+now->y2+now->y3+2*y4);
			q = h/2*(y1+y4);
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
				now->b = now->a+h/2;
				now->y3 = now->y2;
				now->y2 = y1;
				now->acc /= sqrt2;
				now->next = (void*) next;

				next->lvl = now->lvl;
				next->a = now->b;
				next->y2 = next->y3;
				next->y3 = y4;
				next->acc = now->acc;
			}
		}
		else
		{
			fprintf(stderr,"Too many subdivisions during integration!\nTerminating integrator\t...\t");
			while(now != NULL)
			{
				done = now;
				now = (stack*) done->next;
				free(done);
			}
			fprintf(stderr,"done!\n");
		}
	} while (now != NULL);
	return I;
}

double qaso(double f(double), double a, double b, double acc, double eps, double * err)
{
	if (a < b)
	{
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
			return 0;
		}
	}
	else if (b < a)
	{
		// reverse order of limits, caro24 and infinity-evalutaion is designed for a < b
		return -qaso(f,b,a,acc,eps,err);
	}
	else if (a == b )
	{
		// integral is by definition zero
		*err = 0;
		return 0;
	}
	else
	{
		fprintf(stderr,"Something is wrong with the limits!\n");
		return 0;
	}
}
