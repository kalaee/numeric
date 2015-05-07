#include <assert.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#define PI 3.14159265358979323846264338327950288419716939937510

/* perform rotations in row-wise cyclic sweeps on symmetric
 * matrix A to diagonalize it, using the Jacobi eigenvalue algorithm.
 * Upper triagnular part is destroyed but lower is conserved
 * The eigenvalues and corresponding vectors are
 * stored in eig and V respectively.
 * The function output is number of rotations before convergence */
int jacobi_cyclic(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int SORT)
{
	assert(A->size1 == A->size2 && V->size1 == V->size2 && A->size1 == V->size1 && eig->size == V->size1);
	int flag, rot, p, q, i, n = A->size1;
	double app, aqq, apq, phi, s, c, app1, aqq1;
	double aip, aiq, api, aqi, vip, viq;
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}
	gsl_matrix_set_identity(V);
	rot = 0;
	do {
		flag = 0;
		for(p=0; p<n; p++)
		{
			for(q=p+1; q<n; q++)
			{
				app = gsl_vector_get(eig,p);
				aqq = gsl_vector_get(eig,q);
				apq = gsl_matrix_get(A,p,q);
				phi = (atan2(2*apq,aqq-app)+SORT*PI)/2;
				s = sin(phi);
				c = cos(phi);
				app1 = c*c*app-2*s*c*apq+s*s*aqq;
				aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
				if (app1 != app || aqq1 != aqq)
				{
					flag = 1; rot++;
					gsl_vector_set(eig,p,app1);
					gsl_vector_set(eig,q,aqq1);
					gsl_matrix_set(A,p,q,0);
					for(i=0; i<p; i++)
					{
						aip = gsl_matrix_get(A,i,p);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,i,p,c*aip - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*aip);
					}
					for(i=p+1; i<q; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,p,i,c*api - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*api);
					}
					for(i=q+1; i<n; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aqi = gsl_matrix_get(A,q,i);
						gsl_matrix_set(A,p,i,c*api - s*aqi);
						gsl_matrix_set(A,q,i,c*aqi + s*api);
					}
					for(i=0; i<n; i++)
					{
						vip = gsl_matrix_get(V,i,p);
						viq = gsl_matrix_get(V,i,q);
						gsl_matrix_set(V,i,p,c*vip - s*viq);
						gsl_matrix_set(V,i,q,c*viq + s*vip);
					}
				}
			}
		}	
	} while (flag != 0);

	return rot;
}

/* perform rotations in cyclic sweeps, minimising all off-diagonal elements
 * in a single row before continuing to the next. The integer SORT determines
 * whether to write the algorithm so as to calculate the eigenvalues in
 * ascending or descending manner. The integers are defined in jacobi.h and
 * are signified JACOBI_SORT_ASC and JACOBI_SORT_DESC. Mathematically they
 * assign the quadrant in which to find the angle from tangent.
 * THe used algorithm is Jacobi Eigenvalue Algorithm.
 * Upper triagnular part is destroyed but lower is conserved
 * The eigenvalues and corresponding vectors are stored in eig and V respectively.
 * The function output is number of rotations before convergence */
int jacobi_row(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int SORT)
{
	assert(A->size1 == A->size2 && V->size1 == V->size2 && A->size1 == V->size1 && eig->size == V->size1);
	int flag, rot, p, q, i, n = A->size1;
	double app, aqq, apq, phi, s, c, app1, aqq1;
	double aip, aiq, api, aqi, vip, viq;
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}
	gsl_matrix_set_identity(V);
	rot = 0;
	for(p=0; p<n; p++)
	{
		do
		{
			flag = 0;
			for(q=p+1; q<n; q++)
			{
				app = gsl_vector_get(eig,p);
				aqq = gsl_vector_get(eig,q);
				apq = gsl_matrix_get(A,p,q);
				phi = (atan2(2*apq,aqq-app)+SORT*PI)/2;
				s = sin(phi);
				c = cos(phi);
				app1 = c*c*app-2*s*c*apq+s*s*aqq;
				aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
				if (app1 != app)
				{
					flag = 1; rot++;
					gsl_vector_set(eig,p,app1);
					gsl_vector_set(eig,q,aqq1);
					gsl_matrix_set(A,p,q,0);
					for(i=0; i<p; i++)
					{
						aip = gsl_matrix_get(A,i,p);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,i,p,c*aip - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*aip);
					}
					for(i=p+1; i<q; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,p,i,c*api - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*api);
					}
					for(i=q+1; i<n; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aqi = gsl_matrix_get(A,q,i);
						gsl_matrix_set(A,p,i,c*api - s*aqi);
						gsl_matrix_set(A,q,i,c*aqi + s*api);
					}
					for(i=0; i<n; i++)
					{
						vip = gsl_matrix_get(V,i,p);
						viq = gsl_matrix_get(V,i,q);
						gsl_matrix_set(V,i,p,c*vip - s*viq);
						gsl_matrix_set(V,i,q,c*viq + s*vip);
					}
				}
			}
		} while (flag != 0);
	}

	return rot;
}

/* perform rotations in cyclic sweeps, minimising the biggest off-diagonal element
 * in a single row before continuing to the next. The integer SORT determines
 * whether to write the algorithm so as to calculate the eigenvalues in
 * ascending or descending manner. The integers are defined in jacobi.h and
 * are signified JACOBI_SORT_ASC and JACOBI_SORT_DESC. Mathematically they
 * assign the quadrant in which to find the angle from tangent.
 * The used algorithm is Jacobi Eigenvalue Algorithm.
 * Upper triagnular part is destroyed but lower is conserved
 * The eigenvalues and corresponding vectors are stored in eig and V respectively.
 * The function output is number of rotations before convergence */
int jacobi_max_row(gsl_matrix* A, gsl_vector* eig, gsl_matrix* V, int SORT)
{
	assert(A->size1 == A->size2 && V->size1 == V->size2 && A->size1 == V->size1 && eig->size == V->size1);
	int flag, rot, p, q, i, n = A->size1;
	double app, aqq, apq, phi, s, c, app1, aqq1;
	double aip, aiq, api, aqi, vip, viq, max;
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}

	gsl_matrix_set_identity(V);

	rot = 0;
	for(p=0; p<n-1; p++)
	{
		do
		{
			flag = 0;
			q = p+1;
			max = gsl_matrix_get(A,p,q);
			for(i=p+2; i<n; i++)
			{
				if(max < fabs(gsl_matrix_get(A,p,i)))
				{
					max = fabs(gsl_matrix_get(A,p,i));
					q = i;
				}
			}
			if(q < n)
			{
				app = gsl_vector_get(eig,p);
				aqq = gsl_vector_get(eig,q);
				apq = gsl_matrix_get(A,p,q);
				phi = (atan2(2*apq,aqq-app)+SORT*PI)/2;
				s = sin(phi);
				c = cos(phi);
				app1 = c*c*app-2*s*c*apq+s*s*aqq;
				aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
				if (app1 != app)
				{
					flag = 1; rot++;
					gsl_vector_set(eig,p,app1);
					gsl_vector_set(eig,q,aqq1);
					gsl_matrix_set(A,p,q,0);
					for(i=0; i<p; i++)
					{
						aip = gsl_matrix_get(A,i,p);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,i,p,c*aip - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*aip);
					}
					for(i=p+1; i<q; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aiq = gsl_matrix_get(A,i,q);
						gsl_matrix_set(A,p,i,c*api - s*aiq);
						gsl_matrix_set(A,i,q,c*aiq + s*api);
					}
					for(i=q+1; i<n; i++)
					{
						api = gsl_matrix_get(A,p,i);
						aqi = gsl_matrix_get(A,q,i);
						gsl_matrix_set(A,p,i,c*api - s*aqi);
						gsl_matrix_set(A,q,i,c*aqi + s*api);
					}
					for(i=0; i<n; i++)
					{
						vip = gsl_matrix_get(V,i,p);
						viq = gsl_matrix_get(V,i,q);
						gsl_matrix_set(V,i,p,c*vip - s*viq);
						gsl_matrix_set(V,i,q,c*viq + s*vip);
					}
				}
			}
		} while (flag != 0);
	}

	return rot;
}
