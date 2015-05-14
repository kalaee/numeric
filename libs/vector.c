#include <gsl/gsl_vector.h>

void vector_sum(double alpha, gsl_vector* a, double beta, gsl_vector* b, gsl_vector* c)
{
	int i;
	for(i=0; i<a->size; i++)
	{
		gsl_vector_set(c,i,alpha*gsl_vector_get(a,i)+beta*gsl_vector_get(b,i));
	}
	return;
}
