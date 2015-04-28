#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include "qr.h"

#define	ROW	100
#define COL	90

int main(void)
{
	int i,j;
	// define the system
	gsl_matrix* A = gsl_matrix_alloc(ROW,COL);
	gsl_matrix* B = gsl_matrix_alloc(COL,COL);
	gsl_matrix* BI = gsl_matrix_alloc(COL,COL);
	gsl_matrix* R = gsl_matrix_alloc(COL,COL);
	gsl_vector* b = gsl_vector_alloc(ROW);
	gsl_vector* y = gsl_vector_alloc(COL);
	gsl_vector* x = gsl_vector_alloc(COL);

	// define the entries of matrix A and vector y
	for (i = 0; i < COL; i++)
	{
		gsl_vector_set(y,i,sin(i)+cos(i*i));
		for (j = 0; j < ROW; j++)
		{
			gsl_matrix_set (A, j, i, j*sin (i) + cos (j*i));
		}
	}

	// copy COLxCOL submatrix of A into B for later
	gsl_matrix_view a = gsl_matrix_submatrix(A,0,0,COL,COL);
	gsl_matrix_memcpy(B,&a.matrix);

	// construct the vector b
	gsl_blas_dgemv(CblasNoTrans,1,A,y,0,b);

	// perform QR decomposition on A and solve Ax = b for x
	qr_dec(A,R);
	qr_bak(A,R,b,x);

	// find the norm of the vector x-y, if zero: the solution is correct
	gsl_vector_sub(x,y);
	
	printf("Solve Ax = b for x, where b = Ay, using QR decomp\n");
	printf("Evaluating the deviation between x and y:\n");
	printf("\t|x-y| =\t%g\n",gsl_blas_dnrm2(x));

	// the abs. val of determinant from QR decomp and the inverse
	gsl_matrix_memcpy(&a.matrix,B);
	qr_dec(B,R);
	double d = qr_absdet(R);
	qr_inv(B,R,BI,x);

	// to evaluate Binv, we measure the entrywise norm of B*Binv - I
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&a.matrix,BI,0,R);
	printf("\nCompute the inverse matrix Binv of B. The entrywise norm of B*Binv - I is\n");
	printf("\t|B*Binv - I| =\t");
	double sum = 0;
	for(i=0; i < COL; i++)
	{
		gsl_matrix_set(R,i,i,gsl_matrix_get(R,i,i)-1);
		for(j=0; j < COL; j++)
		{
			sum += pow(gsl_matrix_get(R,i,j),2);
		}
	}
	printf("%g\n",sqrt(sum));

	// determinant from GSL using LU decomp
	gsl_permutation* p = gsl_permutation_alloc(COL);
	gsl_linalg_LU_decomp(&a.matrix,p,&i);
	double dgsl = fabs(gsl_linalg_LU_det(&a.matrix,i));
	printf("\nCompare the algorithm for computation");
	printf(" of the absolute value of the determinant\n");
	printf("\t|det(A)|/|det(A)_gsl| - 1 =\t%g\n",d/dgsl-1);

	gsl_vector_free(y);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(B);
	gsl_matrix_free(BI);
	gsl_permutation_free(p);
	return 0;
}
