#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void qr_dec(gsl_matrix* A, gsl_matrix* R);
void qr_bak(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* res);
double qr_absdet(gsl_matrix* R);
void qr_inv(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* AI, gsl_vector* w);
