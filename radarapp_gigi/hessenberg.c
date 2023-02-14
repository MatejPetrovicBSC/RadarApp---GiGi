#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "hessenberg.h"
#include "eigen.h"
#include "unitary.h"
//#include <vehave-control.h>
#define VECT 0
#include "cycles.h"


void compute_hessenberg(int n, float resultReal[n][n], float resultImag[n][n], float complex a[n][n], float complex tau[n-1], float D[n], float Er[n])
{

    float E[16];
	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i][j] = resultReal[i][j] + resultImag[i][j] * I;
		}
	}
	/*
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf(" %f %f ", creal(a[i][j]), cimag(a[i][j]));
		}
	}*/

	//Step 2, diagonal and off-diagonal elements of Hessenberg
	ZHETRD('u', n, a, resultReal, resultImag, D, E, tau);

//    Extrae_fini();
/*
	printf("D values:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t\t%.8g\n", D[i]);
	}
	printf("\n\n");
*/
		
	for(int i = n; i> 0; i--) {
		Er[i] = E[i-1];
	}
	Er[0] = 0.0;
/*	
	
	printf("E values:\n");
	for (int i = 0; i < n; i++)
	{
		printf("\t\t%.8g\n", Er[i]);
	}
*/
	printf("\nHessenberg matrix \n\n");
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf(" %g + %gi ", creal(a[i][j]), cimag(a[i][j]));
		}
		printf("\n");
	}

}

float complex conjunction(float xreal, float ximag)
{
    float Ireal = creal(I);
    return xreal + (ximag * -1) * I;
}

void DLADIV(float A, float B, float C, float D, float *P, float *Q)
{
    float E, F;

    if (fabs(D) < fabs(C))
    {
        E = D / C;
        F = C + D * E;
        *P = (A + B * E) / F;
        *Q = (B - A * E) / F;
    }
}
float complex ZLADIV(float complex x, float complex y)
{
    float zi, zr;
    DLADIV(creal(x), cimag(x), creal(y), cimag(y), &zr, &zi);

    return zr + zi * I;
}

float DZNRM2(int n, float complex x[n], float xreal[n], float ximag[n])
{
    float out;
    out = 0.0;
    float temp[n], scale, ssq;

    if (n < 1)
    {
        return 0;
    }
    else
    {
        scale = 0;
        ssq = 1;
        for (int ix = 0; ix < n; ix++)
        {
            if (creal(x[ix]) != 0.0)
            {
                temp[ix] = fabs(xreal[ix]);
                if (scale < temp[ix])
                {
                    ssq = 1.0 + ssq * (scale / temp[ix]) * (scale / temp[ix]);
                    scale = temp[ix];
                }
                else
                {
                    ssq = ssq + (temp[ix] / scale) * (temp[ix] / scale);
                }
            }
            if (cimag(x[ix]) != 0.0)
            {
                temp[ix] = fabs(cimag(x[ix]));
                if (scale < temp[ix])
                {
                    ssq = 1.0 + ssq * (scale / temp[ix]) * (scale / temp[ix]);
                    scale = temp[ix];
                }
                else
                {
                    ssq = ssq + (temp[ix] / scale) * (temp[ix] / scale);
                }
            }
        }
        return (scale * sqrt(ssq));
    }
}

float DLAPY3(float X, float Y, float Z)
{
    return sqrt(X * X + Y * Y + Z * Z);
}

void ZAXPY(int n, float complex *za, float complex zx[n], float zxreal[n], float zximag[n], float complex zy[n])
{
    float zareal, zaimag;
    float zyreal[n], zyimag[n];

    zareal = za[0];
    zaimag = za[1];

//    Extrae_eventandcounters(1000,1);
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (int i = 0; i < n; i++)
    {

        zyreal[i] = creal(zy[i]);
        zyimag[i] = cimag(zy[i]);

        float real = zareal * zxreal[i] - zaimag * zximag[i];
        float imag = zareal * zximag[i] + zaimag * zxreal[i];

        zy[i] = (zyreal[i] + real) + (zyimag[i] + imag) * I;

        //zy[i] = zy[i] + za * zx[i];
    }

//    Extrae_eventandcounters(1000,0);
}

float complex ZDOTC(int n, float complex zx[n], float complex zy[n])
{
    float complex ztemp;
    float xconjreal[n], xconjimag[n];
    float zxreal[n], zximag[n], zyreal[n], zyimag[n];


    ztemp = 0.0 + 0.0 * I;

//    Extrae_eventandcounters(1000,2);
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (int i = 0; i < n; i++)
    {
        zxreal[i] = creal(zx[i]);
        zximag[i] = cimag(zx[i]);

        zyreal[i] = creal(zy[i]);
        zyimag[i] = cimag(zy[i]);

        xconjreal[i] = creal(conjunction(zxreal[i], zximag[i]));
        xconjimag[i] = cimag(conjunction(zxreal[i], zximag[i]));

        float real = xconjreal[i] * zyreal[i] - xconjimag[i] * zyimag[i];
        float imag = xconjreal[i] * zyimag[i] + xconjimag[i] * zyreal[i];

        ztemp = (creal(ztemp) + real) + (cimag(ztemp) + imag) * I;
        //ztemp = ztemp + xconj[i] * zy[i];
    }
//    Extrae_eventandcounters(1000,0);
    return ztemp;
}

void ZSCAL(int n, float complex za, float complex zx[n], float zxreal[n], float zximag[n])
{
    float zareal, zaimag;
    zareal = creal(za);
    zaimag = cimag(za);

//    Extrae_eventandcounters(1000,3);
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (int i = 0; i < n; i++)
    {
        zx[i] = za * zx[i];
        //zx[i] = real + imag * I;
    }
//    Extrae_eventandcounters(1000,0);
}

void ZHEMV(char uplo, int n, float complex alpha, float complex A[n][n], float complex x[n], float complex beta, float complex y[n])
{
    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;
    float complex temp1[n], temp2;

    float xreal[n], ximag[n];
    float Areal[n][n], Aimag[n][n];
    float real[n], imag[n];

    float tmp1real[n], tmp1imag[n];

    float Aconjreal, Aconjimag;
    float alphareal = creal(alpha);
    float alphaimag = cimag(alpha);

    for (int i = 0; i < n; i++)
    {
        y[i] = ZERO;
    }

//    Extrae_eventandcounters(1000,4);
    #if VECT
    #pragma clang loop vectorize(enable)
    #endif
    for (int j = 0; j < n; j++)
    {
        xreal[j] = creal(x[j]);
        ximag[j] = cimag(x[j]);

        tmp1real[j] = alphareal * xreal[j] - alphaimag * ximag[j];
        tmp1imag[j] = alphareal * ximag[j] + alphaimag * xreal[j];

        temp1[j] = tmp1real[j] + tmp1imag[j] * I;
    }
//    Extrae_eventandcounters(1000,0);

 //   Extrae_eventandcounters(1000,5);
    for (int j = 0; j < n; j++)
    {
	temp2 = ZERO;
//        Extrae_eventandcounters(1006,j+1);
	#if VECT
	#pragma clang loop vectorize(enable)
	#endif
	for (int i = 0; i <= j - 1; i++)
        {
            Areal[i][j] = creal(A[i][j]);
            Aimag[i][j] = cimag(A[i][j]);

            real[i] = tmp1real[j] * Areal[i][j] - tmp1imag[j] * Aimag[i][j];
            imag[i] = tmp1real[j] * Aimag[i][j] + tmp1imag[j] * Areal[i][j];

            y[i] = (creal(y[i]) + real[i]) + (cimag(y[i]) + imag[i]) * I;
            //y[i] = y[i] + temp1 * A[i][j];

            Aconjreal = creal(conjunction(Areal[i][j], Aimag[i][j]));
            Aconjimag = cimag(conjunction(Areal[i][j], Aimag[i][j]));

            float lhs = Aconjreal * xreal[i] - Aconjimag * ximag[i];
            float rhs = Aconjreal * ximag[i] + Aconjimag * xreal[i];

            temp2 = (creal(temp2) + lhs) + (cimag(temp2) + rhs) * I;
            //temp2 = temp2 + conjunction(A[i][j]) * x[i];
        }
//    	Extrae_eventandcounters(1006,0);
    }
//    Extrae_eventandcounters(1000,0);

}

void ZHER2(char uplo, int n, float complex alpha, float complex x[n], float complex y[n], float complex A[n][n])
{

    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;
    float complex temp1, temp2;

    float alphareal = creal(alpha);
    float alphaimag = cimag(alpha);

//    Extrae_eventandcounters(1000,6);
    for (int j = 0; j < n; j++)
    {
        if (x[j] != ZERO || y[j] != ZERO)
	{
            temp1 = alpha * conjunction(creal(y[j]), cimag(y[j]));
            temp2 = conjunction(creal(alpha * x[j]), cimag(alpha * x[j]));
 //   	    Extrae_eventandcounters(1008,j+1);
	    #if VECT
	    #pragma clang loop vectorize(assume_safety)
	    #endif
	    for (int i = 0; i <= j - 1; i++)
            {
                A[i][j] = A[i][j] + x[i] * temp1 + y[i] * temp2;
            }
//    	    Extrae_eventandcounters(1008,0);
            
	    A[j][j] = creal(A[j][j]) + creal(x[j] * temp1 + y[j] * temp2);
        }
        else
        {
            A[j][j] = creal(A[j][j]);
        }
    }
//    Extrae_eventandcounters(1000,0);
}

void ZLARFG(int n, float complex *alpha, float complex x[n - 1], float *xreal, float ximag[n - 1], float complex *tau)
{
    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;

    float xnorm = DZNRM2(n - 1, x, xreal, ximag);

    float alphr = creal(*alpha);
    float alphi = cimag(*alpha);
    float beta;


    if (xnorm == ZERO && alphi == ZERO)
    {
        *tau = ZERO;
    }
    else
    {
        if (alphr >= 0)
        {
            beta = -fabs(DLAPY3(alphr, alphi, xnorm));
        }
        else
        {
            beta = fabs(DLAPY3(alphr, alphi, xnorm));
        }
        *tau = ((beta - alphr) / beta) + (-alphi / beta) * I;
        *alpha = ZLADIV(ONE, *alpha - beta);
        ZSCAL(n - 1, *alpha, x, xreal, ximag);

        *alpha = beta;
    }
}

void ZHETD2(char uplo, int n, float complex A[n][n], float Areal[n][n], float Aimag[n][n], float D[n], float E[n - 1], float complex tau[n - 1])
{
    float complex ONE = 1.0 + 0 * I;
    float complex ZERO = 0.0 + 0 * I;
    float complex HALF = 0.5 + 0 * I;
    bool upper = (uplo == 'u' || uplo == 'U');
    float complex taui;


    if (uplo == 'u' || uplo == 'U')
    {
        A[n - 1][n - 1] = creal(A[n - 1][n - 1]);
        for (int i = n - 2; i >= 0; i--)
        {
            float complex alpha = A[i][i + 1];
            float complex tempArray[16] = {0.0 + 0.0*I};
	    float tempArrayreal[16] = {0.0};
	    float tempArrayimag[16] = {0.0};
	    //float complex tempArray[i];
            //float tempArrayreal[i];
            //float tempArrayimag[i];

            for (int j = 0; j < i + 1; j++)
            {

                tempArray[j] = A[j][i + 1];

                tempArrayreal[j] = creal(A[j][i + 1]);
                tempArrayimag[j] = cimag(A[j][i + 1]);
            }

            ZLARFG(i + 1, &alpha, tempArray, tempArrayreal, tempArrayimag, &taui);

            for (int j = 0; j < i + 1; j++)
            {
                A[j][i + 1] = tempArray[j];
            }

            E[i] = creal(alpha);
            if (taui != ZERO)
            {
                A[i][i + 1] = ONE;

		float complex tempArray[16] = {0.0 + 0.0*I};
		float tempArrayreal[16] = {0.0};
		float tempArrayimag[16] = {0.0};
                //float complex tempArray[i];
                //float tempArrayreal[i];
                //float tempArrayimag[i];

                for (int j = 0; j < i + 1; j++)
                {
                    tempArray[j] = A[j][i + 1];
                    tempArrayreal[j] = creal(A[j][i + 1]);
                    tempArrayimag[j] = cimag(A[j][i + 1]);
                }

                float complex temp2Array[i + 1][i + 1];
                float complex tempArray2[i + 1];

                for (int j = 0; j < i + 1; j++)
                {
                    for (int k = 0; k < i + 1; k++)
                    {
                        temp2Array[j][k] = A[j][k];
                    }
                    tempArray2[j] = tau[j];
                }

                ZHEMV(uplo, i + 1, taui, temp2Array, tempArray, ZERO, tempArray2);

                alpha = -HALF * taui * ZDOTC(i + 1, tempArray2, tempArray);

                ZAXPY(i + 1, &alpha, tempArray, tempArrayreal, tempArrayimag, tempArray2);

                ZHER2(uplo, i + 1, -ONE, tempArray, tempArray2, temp2Array);

		//#pragma clang loop vectorize(enable)
                for (int j = 0; j < i + 1; j++)
                {
                    for (int k = 0; k < i + 1; k++)
                    {
                        A[j][k] = temp2Array[j][k];
                    }
                    A[j][i + 1] = tempArray[j];
                    tau[j] = tempArray2[j];
                }
            }
            else
            {
                for (int j = 0; j < i + 1; j++)
                {
                    for (int k = 0; k < i + 1; k++)
                    {
                        A[j][k] = creal(A[j][k]);
                    }
                }
            }

            A[i][i + 1] = E[i] + 0 * I;

            D[i + 1] = creal(A[i + 1][i + 1]);
            tau[i] = taui;
        }
        D[0] = creal(A[0][0]);
    }
    else
    {
        A[0][0] = creal(A[0][0]);
        for (int i = 0; i < n - 1; i++)
        {
            float complex alpha = A[i + 1][i];
            float complex tempArray[n - i - 1];
            float tempArrayreal[i];
            float tempArrayimag[i];

            for (int j = 0; j < n - i - 1; j++)
            {
                tempArray[j] = A[i + 2 + j][i];
                tempArrayreal[j] = creal(A[i + 2 + j][i]);
                tempArrayimag[j] = cimag(A[i + 2 + j][i]);
            }
            ZLARFG(n - i - 1, &alpha, tempArray, tempArrayreal, tempArrayimag, &taui);

            for (int j = 0; j < n - i - 1; j++)
            {
                A[i + 2 + j][i] = tempArray[j];
            }

            E[i] = creal(alpha);

            if (taui != ZERO)
            {
                A[i + 1][i] = ONE;

                float complex temp2Array[n - i - 1][n - i - 1];
                float complex tempArray[n - i - 1];
                float tempArrayreal[n - i - 1];
                float tempArrayimag[n - i - 1];
                float complex tempArray2[n - i - 1];

                for (int j = 0; j < n - i - 1; j++)
                {
		    #if VECT
		    #pragma clang loop vectorize(enable)	
		    #endif
		    for (int k = 0; k < n - i - 1; k++)
                    {
                        temp2Array[j][k] = A[j + i + 1][k + i + 1];
                    }
                    tempArray[j] = A[i + j + 1][i];
                    tempArrayreal[j] = creal(A[i + j + 1][i]);
                    tempArrayimag[j] = cimag(A[i + j + 1][i]);
                    tempArray2[j] = tau[i + j + 1];
                }

                ZHEMV(uplo, n - i - 1, taui, temp2Array, tempArray, ZERO, tempArray2);
                alpha = -HALF * taui * ZDOTC(n - i - 1, tempArray2, tempArray);
                ZAXPY(n - i - 1, &alpha, tempArray, tempArrayreal, tempArrayimag, tempArray2);
                ZHER2(uplo, n - i - 1, -ONE, tempArray, tempArray2, temp2Array);

                for (int j = 0; j < n - i - 1; j++)
                {
                    for (int k = 0; k < n - i - 1; k++)
                    {
                        A[j + i + 1][k + i + 1] = temp2Array[j][k];
                    }
                    A[i + j + 1][i] = tempArray[j];
                    tau[i + j] = tempArray2[j];
                }
            }
            else
            {
                for (int j = 0; j < n - i; j++)
                {
                    for (int k = 0; k < n - i; k++)
                    {
                        A[j + i][k + i] = creal(A[j + i][k + i]);
                    }
                }
            }

            A[i + 1][i] = E[i] + 0 * I;

            D[i] = creal(A[i][i]);
            tau[i] = taui;
        }
        D[n - 1] = creal(A[n - 1][n - 1]);
    }
}

void ZHETRD(char uplo, int n, float complex A[n][n], float Areal[n][n], float Aimag[n][n], float D[n], float E[n - 1], float complex tau[n - 1])
{
    ZHETD2(uplo, n, A, &Areal[0], &Aimag[0], D, E, tau);
}
