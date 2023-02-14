#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "MUSIC_algorithm.h"
#include "covariance.h"
//#include <vehave-control.h>
#ifndef M_PI
    #define M_PI 3.141592
#endif

#define VECT 0
void conjugate_transpose(int n, int m, Complex D[n][m], Complex Dtrans[m][n]);
float largest(int n, float arr[n]);
Complex transpose_SS(Complex SSvalue, Complex SStransvalue);

void compute_MUSIC(int m, float complex eigenVectors[m][m])
{
    int N = 200; //snapshot

    float doa[2] = {20.0 / 180 * M_PI, 60.0 / 180 * M_PI}; //direction of arrival
    float w[2] = {M_PI / 4, M_PI / 3};                     //frequency

    int M = 16.0;                     //number of array elements
    int P = sizeof(w) / sizeof(w[0]); //the number of signal
    float lambda = 150.0;             //wave length
    float d = lambda / 2;             //element spacing
    float snr = 20.0;                 //SNA

    Complex D[P][M];
    float Dreal;
    float Dimag;
    float z_im;

    for (int i = 0; i < P; i++)
    {
        z_im = -1 * 2 * M_PI * d * sin(doa[i]);

        if (z_im == 0)
        {
            z_im = 0.0;
        }
        else
        {
            z_im /= lambda;
        }

        for (int j = 0; j < M; j++)
        {
            Dimag = (float)j * z_im;
            if (Dimag == 0.0)
            {
                Dreal = 1.0;
                Dimag = 0.0;
            }
            else
            {
                Dreal = cos(Dimag);
                Dimag = sin(Dimag);
            }
            D[i][j].re = Dreal;
            D[i][j].im = Dimag;
        }
    }

    //copnjugate transpose matrix
    Complex Dtrans[M][P];
    conjugate_transpose(P, M, D, Dtrans);


    Complex NN[M][M - P]; //estimate noise subspace
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M - P; j++)
        {
            NN[i][j].re = creal(eigenVectors[i][j]);
	        NN[i][j].im = cimag(eigenVectors[i][j]);
        }
    }

    float angle = -90.0;
    float theta[361]; //peak search (-90:0.5:90)

    for (int i = 0; i < 361; i++)
    {
        theta[i] = angle;
        angle += 0.5;
    }

    Complex NNtrans[M - P][M];
    conjugate_transpose(M, M - P, NN, NNtrans);

    Complex SS[M];
    Complex SStrans[M];

    Complex NN_trans_res[M][M];
    Complex SS_NN_NNtrans[M];

    Complex PP;

    float Pmusic[sizeof(theta) / sizeof(theta[0])];

    float SSreal;
    float SSimag;
   
    for(int j=0; j < M; j++) {
	    for(int m=0; m < M; m++) {
		    NN_trans_res[j][m].re = 0.0;
		    NN_trans_res[j][m].im = 0.0;
	}
    for(int l=0; l < M-P; l++) {
        #if VECT
	   	#pragma clang loop vectorize(assume_safety)
	   	#endif
		for(int k=0; k < M; k++) {      		
			NN_trans_res[j][k].re += NN[j][l].re * NNtrans[l][k].re - NN[j][l].im * NNtrans[l][k].im;
			NN_trans_res[j][k].im += NN[j][l].re * NNtrans[l][k].im + NN[j][l].im * NNtrans[l][k].re;
      	   	}
   	}

   }	
   
    for (int i = 0; i < sizeof(theta) / sizeof(theta[0]); i++)
    {

        for (int j = 0; j < M; j++)
        {
            SSimag = -1.0 * 2.0 * j * M_PI * d * sin(theta[i] / 180 * M_PI);

            if (SSimag == 0.0)
            {
                SS[j].im = 0.0;
            }
            else
            {
                SS[j].im = SSimag / lambda;
            }

            if (SS[j].im == 0.0)
            {
                SS[j].re = 1.0;
                SS[j].im = 0.0;
            }
            else
            {
                SSimag = SS[j].im;
                SSreal = SS[j].im;
                SS[j].re = cos(SSimag);
                SS[j].im = sin(SSreal);
            }

            //SS'
            SStrans[j] = transpose_SS(SS[j], SStrans[j]);
	}

    	for (int t = 0; t < M; t++)
    	{
            SS_NN_NNtrans[t].re = 0.0;
            SS_NN_NNtrans[t].im = 0.0;
	    	#if VECT  
	    	#pragma clang loop vectorize(enable)
	    	#endif
	    	for (int g = 0; g < M; g++)
            	{
                	SS_NN_NNtrans[t].re += SS[g].re * NN_trans_res[g][t].re - SS[g].im * NN_trans_res[g][t].im;
                	SS_NN_NNtrans[t].im += SS[g].re * NN_trans_res[g][t].im + SS[g].im * NN_trans_res[g][t].re;
            	}
    	}	
   
  	PP.re=0.0;
    PP.im=0.0;

	#if VECT
	#pragma clang loop vectorize(enable)
   	#endif
    for (int u = 0; u < M; u++)
    {
        PP.re += SS_NN_NNtrans[u].re * SStrans[u].re - SS_NN_NNtrans[u].im * SStrans[u].im;
        PP.im += SS_NN_NNtrans[u].re * SStrans[u].im + SS_NN_NNtrans[u].im * SStrans[u].re;
    }

	Pmusic[i] = 1 / PP.re;
        if (Pmusic[i] < 0)
           Pmusic[i] *= -1;

    } 
  
    float max = largest(361, Pmusic);

    for (int i = 0; i < 361; i++)
    {
        Pmusic[i] = 10.0 * log10(Pmusic[i] / max);
    }
     
    printf("PMusic: \n");
    for (int i = 0; i < 361; i++)
    {
        printf(" %f ", Pmusic[i]);
    }

}

Complex transpose_SS(Complex SSvalue, Complex SStransvalue)
{

    SStransvalue.re = SSvalue.re;
    if (SSvalue.im != 0)
        SStransvalue.im = SSvalue.im * -1;
    return SStransvalue;
}

float largest(int n, float arr[n])
{
    int i;

    float max = arr[0];
    for (i = 1; i < n; i++)
        if (arr[i] > max)
            max = arr[i];

    return max;
}

void conjugate_transpose(int n, int m, Complex D[n][m], Complex Dtrans[m][n])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            Dtrans[j][i].re = D[i][j].re;
            if (D[i][j].im != 0)
                Dtrans[j][i].im = D[i][j].im * -1;
        }
    }
}
