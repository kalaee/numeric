#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void backsub_upper(gsl_matrix* U, gsl_vector* b);
void backsub_lower(gsl_matrix* L, gsl_vector* b);
