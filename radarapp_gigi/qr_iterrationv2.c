#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "qr_iterrationv2.h"
#include "extrae_user_events.h"

#define VECT 0
/* Decomposes the matrix A into QR */
void QRdecompose(float complex A[16][16], float complex Q[16][16], float complex R[16][16])
{
    float complex T[16][1], S[16][1];

    for (int i = 0; i < 16; i++)
    {

        //Qi = Ui
        matrix_copy_column(16, 16, A, i, 16, 16, Q, i);

        for (int j = 0; j < i; j++)
        {
	    R[i][j] = 0.0 + 0.0 * I;
            //r[j,i] = Qj^T * Ui

            matrix_copy_column(16,16, Q, j, 16, 1, T, 0);
            matrix_copy_column(16,16, A, i, 16, 1, S, 0);
            float complex r = 0.0 + 0.0 * I;

	    #if VECT
	    #pragma clang loop vectorize(enable)
	    #endif
            for (int k = 0; k < 16; k++)
            {
                r += T[k][0] * S[k][0];
	    }

	    R[j][i] = r;
	    matrix_column_multiply(16, 1, T,0,r);
	    matrix_column_subtract(16, 16, Q, i, 16, 1, T, 0);

        }

	//r[i,i] = ||Qi||
        R[i][i] = vector_length(Q, i);
        //Qi = Qi/r[i,i]
        matrix_column_divide(Q, i, R[i][i]);
    }

}

/* Copies a matrix column from msrc at column col1 to mdst at column col2 */
void matrix_copy_column(int m1, int n1, float complex msrc[m1][n1], int col1, int m2, int n2, float complex mdst[m2][n2], int col2)
{
    for (int i = 0; i < 16; i++)
    {
        mdst[i][col2] = msrc[i][col1];
    }
}

/* Subtracts m2's column c2 from m1'si column c1 */
void matrix_column_subtract(int l1, int n1, float complex m1[l1][n1], int c1, int l2, int n2, float complex m2[l2][n2], int c2)
{
    const int rows = 16;	
    #if VECT
    #pragma clang loop vectorize(assume_safety)
    #endif	
    for (int i = 0; i < rows; i++)
    {
        m1[i][c1] -= m2[i][c2];
    }

}

/* Returns the length of the vector column in m */
float complex vector_length(float complex m[16][16], int column)
{
    float complex length = 0 + 0 * I;
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (int row = 0; row < 16; row++)
    {
        length += m[row][column] * m[row][column];
    }
    return csqrt(length);
}

/* Divides the matrix column c in m by k */
void matrix_column_divide(float complex m[16][16], int c, float complex k)
{
    const int rows = 16;	
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif	
    for (int i = 0; i < rows; i++)
    {
        m[i][c] /= k;
    }
}

/* Multiplies the matrix column c in m by k */
void matrix_column_multiply(int l, int n, float complex m[l][n], int c, float complex k)
{
    const int rows = 16;		
    #if VECT
    #pragma clang loop vectorize(assume_safety)	
    #endif
    for (int i = 0; i < rows; i++)
    {
        m[i][c] *= k;
    }

}

/* Debugging purposes only */

void print_matrix_test(int n, float complex m[n][n], bool isMatrixQ)
{
    if(isMatrixQ)
	    printf("Matrix Q: \n\n");
    else
	    printf("Matrix R: \n\n");

    for (int row = 0; row < 16; row++)
    {
        printf("[");
        for (int col = 0; col < 16; col++)
        {
            printf("(%.5f + %.5f) ", creal(m[row][col]), cimag(m[row][col]));
        }
        printf("]\n");
    }
    printf("\n");
}
 
bool test_QR(int n, float complex Q[n][n], float complex R[n][n], float complex A[n][n])
{
    float complex res[n][n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            res[i][j] = 0.0;
      	    #if VECT
	    #pragma clang loop vectorize(enable)
	    #endif
      	    for (int k = 0; k < n; k++)
            {
                res[i][j] += Q[i][k] * R[k][j];
            }
        }
    }
    //matrix *O = create_matrix_from_array(n, n, res);
    float precision = 0.00001;
    for(int i=0; i < n; i++) {
    	for(int j=0; j < n; j++) {
		if(!(((creal(res[i][j]) - precision) < creal(A[i][j])) && ((cimag(res[i][j]) - precision) < cimag(A[i][j])) && ((creal(res[i][j]) + precision) > creal(A[i][j])) && ((cimag(res[i][j]) + precision) > cimag(A[i][j]))))
			return false;
	}
    } 
    return true;
}

void compute_qr(int n, float complex U[n][n]) {
    float complex Q[n][n], R[n][n]; 

    QRdecompose(U, Q, R);
    print_matrix_test(n, Q, true);
    print_matrix_test(n, R, false);
/*
    bool isQRValid = test_QR(n, Q, R, U);
    if(isQRValid)
	    printf("QR Iterration is valid! Input matrix is equal to Q*R \n");
    else
	    printf("QR Iterration not valid! \n");
*/
}

