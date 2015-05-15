#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>

// pseduo-random numbers from 0 to 1, sufficient for this routine
#define RND ((double)rand()/RAND_MAX)

// structure definition for simplex workspace
typedef struct
{
	int lo, hi, n; // dimension of the system
	gsl_matrix* simplex; // polytope
	gsl_vector* p1; // reflected
	gsl_vector* p2; // expanded/contracted
	gsl_vector* ce; // centroid
	gsl_vector* fp; // f(p_i) values
} simplex_workspace;

// allocate memory for a system in n dimensions
simplex_workspace* simplex_workspace_alloc(int n)
{
	simplex_workspace* W = (simplex_workspace*) malloc(sizeof(simplex_workspace));
	W->n = n;
	W->simplex = gsl_matrix_alloc(n+1,n);
	W->p1 = gsl_vector_alloc(n);
	W->p2 = gsl_vector_alloc(n);
	W->ce = gsl_vector_alloc(n);
	W->fp = gsl_vector_alloc(n+1);
	return W;
}

// free allocated workspace memory
void simplex_workspace_free(simplex_workspace* W)
{
	gsl_matrix_free(W->simplex);
	gsl_vector_free(W->p1);
	gsl_vector_free(W->p2);
	gsl_vector_free(W->ce);
	gsl_vector_free(W->fp);
	free(W);
	return;
}

// generate coordinates to vertices of simplex amoeba pseudo-randomly
// in interval lower <= x < upper for all coordinates.
void simplex_generate(double lower, double upper, simplex_workspace* W)
{
	int i,j, m = W->n+1;
	double width = upper - lower;
	for(i=0; i< m; i++)
	{
		for(j=0; j< W->n; j++)
		{
			gsl_matrix_set(W->simplex,i,j,lower+width*RND);
		}
	}
	return;
}

// simpelx reflection
void reflection(simplex_workspace* W)
{
	int i;
	double hii, cei;
	for(i = 0; i < W->n; i++)
	{
		hii = gsl_matrix_get(W->simplex,W->hi,i);
		cei = gsl_vector_get(W->ce,i);
		gsl_vector_set(W->p1,i,2*cei-hii);
	}
	return;
}

// simplex expansion
void expansion(simplex_workspace* W)
{
	int i;
	double hii, cei;
	for(i=0; i < W->n; i++)
	{
		hii = gsl_matrix_get(W->simplex,W->hi,i);
		cei = gsl_vector_get(W->ce,i);
		gsl_vector_set(W->p2,i,3*cei-2*hii);
	}
	return;
}

// simplex contraction
void contraction(simplex_workspace* W)
{
	int i;
	double hii, cei;
	for(i=0; i<W->n; i++)
	{
		hii = gsl_matrix_get(W->simplex,W->hi,i);
		cei = gsl_vector_get(W->ce,i);
		gsl_vector_set(W->p2,i,0.5*(hii+cei));
	}
	return;
}

// simplex reduction
void reduction(double f(gsl_vector* x), simplex_workspace* W)
{
	int i,k, m = W->n+1;
	double ki, loi;
	gsl_vector_view v;
	// reduce vertices
	for(i=0; i < W->n; i++)
	{
		loi = gsl_matrix_get(W->simplex,W->lo,i);
		for(k = 0; k < m; k++)
		{
			if(k != W->lo)
			{
				ki = gsl_matrix_get(W->simplex,k,i);
				gsl_matrix_set(W->simplex,k,i,0.5*(loi+ki));
			}
		}
	}
	// recalculate function values in affected vertices
	for(i=0; i < m; i++)
	{
		if (i != W->lo)
		{
			v = gsl_matrix_row(W->simplex,i);
			gsl_vector_set(W->fp,i,f(&v.vector));
		}
	}
	return;
}

// update simplex for higher, lower and centroid
void simplex_update(simplex_workspace* W)
{
	int i, m = W->n+1;
	W->hi = 0;
	W->lo = 0;
	double 	highest = gsl_vector_get(W->fp,0),
			lowest = gsl_vector_get(W->fp,0),
			fpi;
	gsl_vector_view v;
	// find higher and lower
	for(i=1; i<m; i++)
	{
		fpi = gsl_vector_get(W->fp,i);
		if(fpi > highest)
		{
			highest = fpi;
			W->hi = i;
		}
		if(fpi < lowest)
		{
			lowest = fpi;
			W->lo = i;
		}
	}
	// find centroid
	gsl_vector_set_zero(W->ce);
	for(i=0; i<m; i++)
	{
		if(i != W->hi)
		{
			v = gsl_matrix_row(W->simplex,i);
			gsl_vector_add(W->ce,&v.vector);
		}
	}
	gsl_vector_scale(W->ce,1./W->n);
	return;
}

// calculate distance between vertices a and b as the euclidean norm
// of the vector from a to b
double simplex_dist(simplex_workspace* W, int a, int b)
{
	double s = 0;
	int i;
	for(i=0; i<W->n; i++)
	{
		s += pow(gsl_matrix_get(W->simplex,a,i)-gsl_matrix_get(W->simplex,b,i),2);
	}
	return sqrt(s);
}

// estimate the size of the simplex as the largest distance between two vertices
double simplex_size(simplex_workspace* W)
{
	double dist, s = 0;
	int i, m = W->n+1;
	for(i=1; i<m; i++)
	{
		dist = simplex_dist(W,0,i);
		if(dist > s)
		{
			s = dist;
		}
	}
	return s;
}

// initial calculation of function values at all vertices
void simplex_initialize(double f(gsl_vector* x), simplex_workspace* W)
{
	int i, m= W->n+1;
	gsl_vector_view v;
	for(i=0; i < m; i++)
	{
		v = gsl_matrix_row(W->simplex,i);
		gsl_vector_set(W->fp,i,f(&v.vector));
	}
}


// simplex routine
// f is n dimensional function to be minimised, lower and upper are bounds of coordinates for initial vertices
// simplex_goal_size is tolerance for convergene and W is workspace for simplex routine of dimension n
// out termination the vector W->ce will contain coordinates for lowest vertex
int simplex(double f(gsl_vector* x),double lower, double upper, double simplex_goal_size, simplex_workspace* W)
{
	int steps = 0;
	double fp1, fp2, flo, fhi;
	// initialize system by generating vertices and
	// finding their function values
	simplex_generate(lower,upper,W);
	simplex_initialize(f,W);
	do
	{
		// make an update for higher, lower and centroid
		simplex_update(W);
		fhi = gsl_vector_get(W->fp,W->hi);
		flo = gsl_vector_get(W->fp,W->lo);
		// make reflection
		reflection(W);
		fp1 = f(W->p1);
		if (fp1 < flo)
		{
			// if f(reflected) < f(lower) attempt expansion
			expansion(W);
			fp2 = f(W->p2);
			if (fp2 < fp1)
			{
				// if f(expanded) < f(reflecred) accept expansion
				gsl_matrix_set_row(W->simplex,W->hi,W->p2);
				gsl_vector_set(W->fp,W->hi,fp2);
			}
			else
			{
				// if not, accept reflection
				gsl_matrix_set_row(W->simplex,W->hi,W->p1);
				gsl_vector_set(W->fp,W->hi,fp1);
			}
		}
		else
		{
			if (fp1 < fhi)
			{
				// if f(reflected) < f(higher) accept reflection
				gsl_matrix_set_row(W->simplex,W->hi,W->p1);
				gsl_vector_set(W->fp,W->hi,fp1);			
			}
			else
			{
				// if not, attempt contraction
				contraction(W);
				fp2 = f(W->p2);
				if (fp2 < fhi)
				{
					// if f(contracted) < f(higher), accept contraction
					gsl_matrix_set_row(W->simplex,W->hi,W->p2);
					gsl_vector_set(W->fp,W->hi,fp2);
				}
				else
				{
					// if not, we must be in a valley, perform reduction
					reduction(f,W);
				}
			}
		}
	steps++;
	// if simplex has reduced sufficiently, that is size(simplex) < simplex_goal_size
	// convergence is achieved. Return number of steps before convergence.
	} while (simplex_size(W) > simplex_goal_size);
	// copy lowest vertex to centroid vector
	gsl_matrix_get_row(W->ce,W->simplex,W->lo);
	return steps;
}
