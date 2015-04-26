#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

void givens_qr_dec(gsl_matrix* A);
void givens_qr_bak(gsl_matrix* QR, gsl_vector* b, gsl_vector* res);
double givens_qr_det(gsl_matrix* QR);
void givens_qr_inv(gsl_matrix* QR, gsl_matrix* AI, gsl_vector* w);
