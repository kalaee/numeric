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
	// store diagonal values
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}
	gsl_matrix_set_identity(V);
	rot = 0;
	// perform jacobi rotation
	do {
		flag = 0;
		// consider each element in order
		for(p=0; p<n; p++)
		{
			for(q=p+1; q<n; q++)
			{
				// find angle to rotate A(p,q) to zero
				app = gsl_vector_get(eig,p);
				aqq = gsl_vector_get(eig,q);
				apq = gsl_matrix_get(A,p,q);
				phi = (atan2(2*apq,aqq-app)+SORT*PI)/2;
				s = sin(phi);
				c = cos(phi);
				app1 = c*c*app-2*s*c*apq+s*s*aqq;
				aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
				// if rotation affected the entries, apply rotation to whole matrix
				// change flag and increment rot by one
				if (app1 != app || aqq1 != aqq)
				{
					flag = 1; rot++;
					gsl_vector_set(eig,p,app1);
					gsl_vector_set(eig,q,aqq1);
					gsl_matrix_set(A,p,q,0);
					// first three performs rotation
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
					// calculate eigenvectors
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
	// continure until no entries can be reduced anymore
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
	// store diagonal entries
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}
	gsl_matrix_set_identity(V);
	rot = 0;
	// reduce completely one row at a time
	for(p=0; p<n; p++)
	{
		do
		{
			flag = 0;
			for(q=p+1; q<n; q++)
			{
				// estimate rotation angle
				app = gsl_vector_get(eig,p);
				aqq = gsl_vector_get(eig,q);
				apq = gsl_matrix_get(A,p,q);
				phi = 0.5*(atan2(2*apq,aqq-app)+SORT*PI);
				s = sin(phi);
				c = cos(phi);
				app1 = c*c*app-2*s*c*apq+s*s*aqq;
				aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
				// if rotation does anything, apply
				// change flag and increment rot
				if (app1 != app || aqq1 != aqq)
				{
					flag = 1; rot++;
					gsl_vector_set(eig,p,app1);
					gsl_vector_set(eig,q,aqq1);
					gsl_matrix_set(A,p,q,0);
					// first three loops apply rotation
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
					// last one calculates eigenvectors
					for(i=0; i<n; i++)
					{
						vip = gsl_matrix_get(V,i,p);
						viq = gsl_matrix_get(V,i,q);
						gsl_matrix_set(V,i,p,c*vip - s*viq);
						gsl_matrix_set(V,i,q,c*viq + s*vip);
					}
				}
			}
		// continue until row completely reduced
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
	// store diagonal values
	for(i=0; i<n; i++)
	{
		gsl_vector_set(eig,i,gsl_matrix_get(A,i,i));
	}
	gsl_matrix_set_identity(V);

	rot = 0;
	// consider each row until completely reduced
	for(p=0; p<n-1; p++)
	{
		do
		{
			// find the largest element of the row
			flag = 0;
			q = p+1;
			max = fabs(gsl_matrix_get(A,p,q));
			for(i=p+2; i<n; i++)
			{
				if(max < fabs(gsl_matrix_get(A,p,i)))
				{
					max = fabs(gsl_matrix_get(A,p,i));
					q = i;
				}
			}
			// estimate neccesary rotation angle
			app = gsl_vector_get(eig,p);
			aqq = gsl_vector_get(eig,q);
			apq = gsl_matrix_get(A,p,q);
			phi = 0.5*(atan2(2*apq,aqq-app)+SORT*PI);
			s = sin(phi);
			c = cos(phi);
			app1 = c*c*app-2*s*c*apq+s*s*aqq;
			aqq1 = s*s*app+2*s*c*apq+c*c*aqq;
			// if rotation changes values, apply
			if (app1 != app || aqq1 != aqq)
			{
				flag = 1; rot++;
				gsl_vector_set(eig,p,app1);
				gsl_vector_set(eig,q,aqq1);
				gsl_matrix_set(A,p,q,0);
				// first three loops rotates matrix
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
				// last loop estimates eigenvectors
				for(i=0; i<n; i++)
				{
					vip = gsl_matrix_get(V,i,p);
					viq = gsl_matrix_get(V,i,q);
					gsl_matrix_set(V,i,p,c*vip - s*viq);
					gsl_matrix_set(V,i,q,c*viq + s*vip);
				}
			}
		} while (flag != 0);
	}

	return rot;
}
