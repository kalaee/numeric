#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#define JACOBI_SORT_ASC 1
#define JACOBI_SORT_DESC -1

int jacobi_cyclic(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V);
int jacobi_row(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int SORT);
int jacobi_max_row(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int SORT);
