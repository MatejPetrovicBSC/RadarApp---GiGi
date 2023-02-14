#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include "qr_iterration.h"
#include "unitary.h"

//#include "/apps/riscv/fpga-sdv/ila2prv/include/trigger.h"

#define MAX(x,y) (((x) >= (y)) ? (x) : (y))
#define MIN(x,y) (((x) <= (y)) ? (x) : (y))
#define VECT 0

int ILACLC(int M, int N, float complex A[M][N], int LDA) {
    float complex ZERO = 0.0 + 0 * I;
    int i;
	int out;
    if (N == 0) {
        out = N;
	} else if ((A[0][N-1] != ZERO || A[M-1][N-1] != ZERO)) {
        out = N;
    } else {
        for (out = N-1; out >= 0; out--) {
            for (i = 0; i < M; i++) {
                if (A[i][out] != ZERO) return out;
            }
        }
    }
    return out;
}

bool LSAME(char CA, char CB) {
    int INTA, INTB;
	char charA, charB;
    bool out;
	out = (CA == CB);
    if (out)
		return out;
    
    INTA = CA;
    INTB = CB;
	charA = CA;
	charB = CB;
    if (INTA >= 97 && INTA <= 122) charA = INTA - 32;
    if (INTB >= 97 && INTB <= 122) charB = INTB - 32;
    out = (charA == charB);
	return out;
}

void CGEMVU(char TRANS, int M, int N, float complex ALPHA, int LDA, float complex A[LDA][LDA], float complex X[M], int INCX, float complex BETA, float complex Y[N], int INCY) {
    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;
    float complex TEMP;
    int i, INFO, j, JX, JY, KX, KY, LENX, LENY;
    bool NOCONJ;
    INFO = 0;
	if ((M == 0) || (N == 0) || ((ALPHA == ZERO) && (BETA == ONE)))
		return;	
	
    LENX = M;
    LENY = N;
    KX = 0;
    KY = 0;
    if (BETA != ONE) {
        for (i = 0; i < LENY; i++)
            Y[i] = ZERO;
    }
    if (ALPHA == ZERO)
		return;
    JY = KY;	
    for (j = 0; j < N; j++) {
        TEMP = ZERO;
	for (i = 0; i < M; i++) {	
            TEMP = TEMP + conj(A[i][j])*X[i];
	}	
        Y[JY] = Y[JY] + ALPHA * TEMP;
        JY = JY + INCY;
    }
	
    return;
}

void CGERCU(int M, int N, float complex ALPHA, int LDA, float complex X[LDA], int INCX, float complex Y[LDA], int INCY, float complex A[LDA][LDA]) {
    float complex ZERO = 0.0 + 0 * I;
    float complex TEMP;
    int i, INFO, j, JY;
    INFO = 0;
    
    
    if ((M == 0) || (N == 0) || (ALPHA == ZERO))
		return;
    JY = 0;

    for (j = 0; j < N; j++) {
        if (Y[JY] != ZERO) {
            TEMP = ALPHA * conj(Y[JY]);
		#if VECT
		#pragma clang loop vectorize(assume_safety)
		#endif
	    	for (i = 0; i < M; i++) {
                	A[i][j] = A[i][j] + X[i] * TEMP;	
	    	}		
        }
        JY = JY + INCY;
    }	
		
    return;
}

void CLARFU(char SIDE, int M, int N, int LDC, float complex V[M], int INCV, float complex TAU, float complex C[LDC][LDC], float complex WORK[LDC]) {

    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;	
	
    int i, LASTV, LASTC;
    LASTV = 0;
    LASTC = 0;
  
    
    if (TAU != ZERO) {
        LASTV = M;
        i = 1 + (LASTV - 1) * INCV;
        while ((LASTV > 0) && (V[i - 1] == ZERO)) {
            LASTV = LASTV - 1;
            i = i - INCV;
        }
        LASTC = ILACLC(LASTV, N, C, LDC);		
    }
	
    if (LASTV > 0) {
        CGEMVU('C', LASTV, LASTC, ONE, LDC, C, V, INCV, ZERO, WORK, 1);
        CGERCU(LASTV, LASTC, -TAU, LDC, V, INCV, WORK, 1, C);
    }
	
    return;
}

void CSCAL(int N, float complex CA, float complex CX[N], int INCX) {
    int i;
    if (N <= 0 || INCX <= 0) 
	    return;	

    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (i = 0; i < N; i++) {
	CX[i] = CA * CX[i];
    }
    return;
}

void CUNG2L(int M, int N, int K, int LDA, float complex A[LDA][LDA], float complex TAU[K], float complex WORK[LDA], int INFO) {
	float complex ONE = 1.0 + 0 * I;
	float complex ZERO = 0.0 + 0 * I;
    int i, ii, j, l;
    INFO = 0;
    for (j=0; j < N-K; j++) {
		for (l=0; l < M; l++)
            		A[l][j] = ZERO;
        	A[M - N + j][j] = ONE;
    }


    for (i=0; i < K; i++) {
        ii = N - K + i;
        A[M - N + ii][ii] = ONE;
		
		float complex tmpArray[LDA];
		for (int tt=0; tt < LDA; tt++)
			tmpArray[tt] = A[tt][ii];
        	
		CLARFU('L', M - N + ii + 1, ii, LDA, tmpArray, 1, TAU[i], A, WORK);

		for (int tt=0; tt < LDA; tt++) {
			A[tt][ii] = tmpArray[tt];
		}

		float complex tmpArray2[M - N + ii];
		for (int tt=0; tt < M - N + ii; tt++) {
			tmpArray2[tt] = A[tt][ii];
		}

		CSCAL(M - N + ii, -TAU[i], tmpArray2, 1);
		for (int tt=0; tt < M - N + ii; tt++) {
			A[tt][ii] = tmpArray2[tt];
		}
        A[M - N + ii][ii] = ONE - TAU[i];
		
        for (l=M-N+ii+1; l < M; l++) {		
            A[l][ii] = ZERO;
	}
    }	
    return;
}

void CUNGQL(int M, int N, int K, int LDA, float complex A[LDA][LDA], float complex TAU[K], int LWORK, float complex WORK[LWORK], int INFO) {
    float complex ZERO = 0.0 + 0 * I;
    bool LQUERY;
    int i, IB, IINFO, IWS, j, KK, L, LDWORK, LWKOPT, NB, NBMIN, NX;	
	
    INFO = 0;
    
    if (INFO == 0) {
		NB = 32;
		LWKOPT = N*NB;
    }
	
    NBMIN = 2;
    NX = 0;
    IWS = N;
    if (NB > 1 && NB < K) {
        NX = 128;
        if (NX < K) {
            LDWORK = N;
            IWS = LDWORK * NB;
            if (LWORK < IWS) {
                NB = LWORK / LDWORK;
                NBMIN = 2;
            }
        }
    }
    if (NB >= NBMIN && NB < K && NX < K) {
        KK = MIN(K, ((K - NX + NB - 1) / NB) * NB);
        for (j=0; j < N-KK; j++)
			for (i=M-KK+1; i < M; i++)
                		A[i][j] = ZERO;
    } else {
        KK = 0;
    }	
    
    CUNG2L(M - KK, N - KK, K - KK, LDA, A, TAU, WORK, IINFO);
     
    WORK[0] = IWS;
    return;
}
			
void CUNGTR(char UPLO, int N, float complex A[N][N], int LDA, float complex TAU[N], float complex WORK[N], int LWORK, int INFO) {
    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;
    bool LQUERY, UPPER;
    int i, IINFO, j, LWKOPT, NB;
    INFO = 0;
    LQUERY = (LWORK == -1);
    
    UPPER = LSAME(UPLO, 'U');
    
    if (INFO == 0) {
		NB = 32;
        LWKOPT = MAX(1, N - 1)*NB;
        WORK[0] = LWKOPT;
    }
    	
    if (UPPER) {

        for (j=0; j < N-1; j++) {
			for (i=0; i < j; i++)
                		A[i][j] = A[i][j+1];
            A[N-1][j] = ZERO;
        }
        for (i=0; i < N-1; i++)
            A[i][N-1] = ZERO;
        A[N-1][N-1] = ONE;
		
        CUNGQL(N - 1, N - 1, N - 1, LDA, A, TAU, LWORK, WORK, IINFO);
    } 
    WORK[0] = LWKOPT;
    return;
}


void compute_unitary(int N, complex float A[N][N], complex float TAU[N-1]) {
	char UPLO = 'U';
	int LDA = N;
	float complex WORK[N];
	int LWORK = N;
	int INFO;
   

	CUNGTR(UPLO, N, A, LDA, TAU, WORK, LWORK, INFO);

	printf("\nUnitary Matrix\n\n");

	for(int i=0; i < LDA; i++)
	{
		for(int j=0; j < LDA; j++) {
			printf(" (%f, %f) ", creal(A[i][j]), cimag(A[i][j]));
		}
		printf("\n");
	}

	printf(" \n ");

	//	qr_iterration(N, A);
}
