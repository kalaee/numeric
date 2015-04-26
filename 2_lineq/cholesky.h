#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

void chol_dec(gsl_matrix* A);
void chol_bak(gsl_matrix* LL, gsl_vector* b);
double chol_det(gsl_matrix* LL);
void chol_inv(gsl_matrix* LL, gsl_matrix* AI);
