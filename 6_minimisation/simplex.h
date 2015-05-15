#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct
{
	int lo, hi, n; // dimension of the system
	gsl_matrix* simplex; // polytope
	gsl_vector* p1; // reflected
	gsl_vector* p2; // expanded/contracted
	gsl_vector* ce; // centroid
	gsl_vector* fp; // f(p_i) values
} simplex_workspace;

simplex_workspace* simplex_workspace_alloc(int n);
void simplex_workspace_free(simplex_workspace* W);
int simplex(double f(gsl_vector* x),double lower, double upper, double simplex_goal_size, simplex_workspace* W);
