#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include "covariance.h"
#include "covariance_v2.h"
#include "hessenberg.h"
#include "unitary.h"
#include "eigen.h"
#include "eigen_v2.h"
#include "qr_iterrationv2.h"
#include "MUSIC_algorithm.h"
#include "MUSIC_algorithm_v3.h"
#include "MUSIC_algorithm_v2.h"
//#include "cycles.h"
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include "extrae_user_events.h"

Complex* read_snapshotMatrix(int m, int n, int ld_snapshot);
void print_matrix(char *desc, int m, int n, Complex *C, int lda);
uint64_t start_cycles, start_instret, end_cycles, end_instret;

int main(int argc, char *argv[])
{
	int m, n, ld_snapshot;
	
	
	m = 16;
	n = 8;
	ld_snapshot = 16; 
	
	Complex *snapshotMatrix;
	snapshotMatrix = (Complex *)malloc(sizeof(Complex) * m * n);
	snapshotMatrix = read_snapshotMatrix(m, n, ld_snapshot);

	float covarianceReal[m][m], covarianceImag[m][m];
	float complex hessenbergMatrix[m][m];
	float complex tau[m-1];
	float D[16], E[16];

	//Step 1: covariance matrix
	compute_covariance(m, n, snapshotMatrix, ld_snapshot, covarianceReal, covarianceImag);

	//Step 2: hessenberg matrix   
	compute_hessenberg(m, covarianceReal, covarianceImag, hessenbergMatrix, tau, D, E);

	//Step 3: unitary matrix
	compute_unitary(m, hessenbergMatrix, tau);
	//Step 4: eigen
	compute_eigenv2(m, &covarianceReal[0][0], &covarianceImag[0][0]);

	//Step 5: QR iterration
	compute_qr(m, hessenbergMatrix);	
	
	//THIS WILL BE CHANGED AFTER STEP 4 IS FIXED (EIGENVECTORS FROM STEP 4 WILL BE THE INPUT FOR STEP 6)
	//THIS MATRIX IS JUST FOR TESTING PURPOSES BUT THE LOGIC FOR STEP 6 WILL REMAIN THE SAME
	float complex dummyMatrix[m][m];
	for(int i=0; i<m; i++) {
		for(int j=0; j<m; j++) {
			dummyMatrix[i][j] = covarianceReal[i][j] + covarianceImag[i][j] * I;
		}	
	}

	//step 6: MUSIC algorithm/
//	double time_spent = 0.0;
//	clock_t begin = clock();
	compute_MUSIC(m, dummyMatrix);
//	clock_t end = clock();
//	time_spent += (double)(end-begin) / CLOCKS_PER_SEC;
	free(snapshotMatrix);
	return 1;
}

Complex* read_snapshotMatrix(int m, int n, int ld_snapshot)
{
	Complex *snapshotMatrix;
	int flag = 0;
	int i, k, num;
	float real, imag;
	char sign;

	FILE *inputFile;
	inputFile = fopen("matrices/Scenario_1.dat", "r");

	if (inputFile == NULL)
	{
		fprintf(stderr, "Can't open the file !\n");
		exit(1);
	}

	real = 0;
	imag = 0;

	snapshotMatrix = (Complex *)malloc(sizeof(Complex) * m * n);

	if (flag == EOF)
	{
		fprintf(stderr, "Error has been encountered while reading values from twidder_factor.txt!\n");
		exit(1);
	}

	for (i = 0; i < m; i++)
	{
		for (k = 0; k < n; k++)
		{
			if ((num = fscanf(inputFile, "%f %c %fi", &real, &sign, &imag)) > 0)
			{
				if (sign == '-')
					imag *= -1;
				snapshotMatrix[i * n + k].re = real;
				snapshotMatrix[i * n + k].im = imag;
			}
			else
			{
				k--;
			}
		}
		printf("\n");
	}

	print_matrix("Snapshot matrix", m, n, snapshotMatrix, n);


	return snapshotMatrix;
}

void print_matrix(char *desc, int m, int n, Complex *C, int lda)
{
	int i, j;
	printf("\n %s\n", desc);
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if ((C[i * lda + j].im) < 0)
			{
				printf("%.8g-%.8gi  ", (C[i * lda + j].re), (C[i * lda + j].im) * -1);
			}
			else
				printf("%.8g+%.8gi  ", (C[i * lda + j].re), (C[i * lda + j].im));
		}
		printf("\n");
	}
}
