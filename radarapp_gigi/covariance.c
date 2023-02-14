#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "covariance.h"
#include "hessenberg.h"
//#include <vehave-control.h>

#define VECT 0

void compute_covariance(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float covarianceReal[m][m], float covarianceImag[m][m])
{
	float meanReal[m], meanImag[m];
	compute_Means(m, n, snapshotMatrix, ld_snapshot, meanReal, meanImag);
	
	float diffReal[m][n], diffImag[m][n];
	compute_diffMatrix(m, n, snapshotMatrix, ld_snapshot, meanReal, meanImag, diffReal, diffImag);

//    	#if VECT
//		compute_resultMatrixVect(m, n, diffReal, diffImag, covarianceReal, covarianceImag);
  //  	#else
		compute_resultMatrix(m, n, diffReal, diffImag, covarianceReal, covarianceImag);
	
//	#endif
	

}
void compute_Means(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float meanReal[m], float meanImag[m])
{
    float meanRealValue;
    float meanImagValue;

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            meanRealValue += (snapshotMatrix[i * n + j]).re;
            meanImagValue += (snapshotMatrix[i * n + j]).im;
        }
        meanReal[i] = meanRealValue / n;
        meanImag[i] = meanImagValue / n;

        meanRealValue = 0.0;
        meanImagValue = 0.0;
    }
}

void compute_diffMatrix(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float meanReal[m], float meanImag[m], float diffReal[m][n], float diffImag[m][n])
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            diffReal[i][j] = (snapshotMatrix[i * n + j]).re;
            diffImag[i][j] = (snapshotMatrix[i * n + j]).im;
        }
    }
}


void compute_resultMatrix(int m, int n, float diffReal[m][n], float diffImag[m][n], float resultReal[m][m], float resultImag[m][m])
{
    //accumulation that in these values and then returned to resultReal and Imag, so the vectorization can work! otherwise, vectorization fails
    float x = 0;
    float y = 0;
 
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            resultReal[i][j] = 0;
            resultImag[i][j] = 0;
	    x=0;
	    y=0;
	    
//	    #pragma clang loop vectorize(enable)
            for (int k = 0; k < n; k++)
            {

                x += diffReal[i][k] * diffReal[j][k] + diffImag[i][k] * diffImag[j][k];
                y += -diffReal[i][k] * diffImag[j][k] + diffImag[i][k] * diffReal[j][k];
            }
	    
	    
	    x /= n;
	    y /= n;

            resultReal[i][j] = x;
            resultImag[i][j] = y;
	        	
        }
    }

    printf("\nCovariance Matrix:\n\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            printf(" ( %f, %f )", resultReal[i][j], resultImag[i][j]);
        }
        printf("\n");
    }
}

/*
void compute_resultMatrixVect(int m, int n, float diffReal[m][n], float diffImag[m][n], float resultReal[m][m], float resultImag[m][m])
{
//	vehave_trace(1000, 1);
	long gvl = __builtin_epi_vsetvl(128, __epi_e32, __epi_m1);
	long gvl2 = __builtin_epi_vsetvl(256, __epi_e32, __epi_m1);

	__epi_2xf32 vdiffReal = __builtin_epi_vload_2xf32(&diffReal[0][0], gvl);
	__epi_2xf32 vdiffImag = __builtin_epi_vload_2xf32(&diffImag[0][0], gvl);

	#if __riscv_vector_version == 700
		#define __builtin_epi_vmv_v_f_2xf32 __builtin_epi_vbroadcast_2xf32
	#endif

	__epi_2xf32 vresultReal = __builtin_epi_vfmv_v_f_2xf32(0.0f, gvl2);
	__epi_2xf32 vresultImag = __builtin_epi_vfmv_v_f_2xf32(0.0f, gvl2);

	for(int i = m-1; i >= 0; i--) {
		__epi_2xf32 vReal8elemsI = __builtin_epi_vslidedown_2xf32(vdiffReal, i*8, gvl);
		__epi_2xf32 vImag8elemsI = __builtin_epi_vslidedown_2xf32(vdiffImag, i*8, gvl);

		for(int j = m-1; j >=0; j--) {
			__epi_2xf32 vReal8elemsJ = __builtin_epi_vslidedown_2xf32(vdiffReal, j*8, gvl);
			__epi_2xf32 vImag8elemsJ = __builtin_epi_vslidedown_2xf32(vdiffImag, j*8, gvl);

			__epi_2xf32 tmp1 = __builtin_epi_vfmul_2xf32(vReal8elemsI, vReal8elemsJ, 8);

			tmp1 = __builtin_epi_vfmacc_2xf32(tmp1, vImag8elemsI, vImag8elemsJ, 8);

			__epi_2xf32 tmpaccum = __builtin_epi_vfmv_v_f_2xf32(0.0f, 8);

			tmp1 = __builtin_epi_vfredsum_2xf32(tmp1, tmpaccum, 8);

			float resultRealIJ = __builtin_epi_vfmv_f_s_2xf32(tmp1);
			resultReal[i][j] = resultRealIJ;

			__epi_2xf32 tmp2 = __builtin_epi_vfmul_2xf32(vReal8elemsI, vImag8elemsJ, 8);

			tmp2 = __builtin_epi_vfmsac_2xf32(tmp2, vImag8elemsI, vReal8elemsJ, 8);

			tmp2 = __builtin_epi_vfredsum_2xf32(tmp2, tmpaccum, 8);

			float resultImagIJ = __builtin_epi_vfmv_f_s_2xf32(tmp2);
			resultImag[i][j] = resultImagIJ;
		}
	}

	vresultReal = __builtin_epi_vload_2xf32(&resultReal[0][0], gvl2);
	vresultImag = __builtin_epi_vload_2xf32(&resultImag[0][0], gvl2);

	__epi_2xf32 vdivisorN = __builtin_epi_vfmv_v_f_2xf32(8.0f, gvl2);

	vresultReal = __builtin_epi_vfdiv_2xf32(vresultReal, vdivisorN, gvl2);
	vresultImag = __builtin_epi_vfdiv_2xf32(vresultImag, vdivisorN, gvl2);

	__builtin_epi_vstore_2xf32(&resultReal[0][0], vresultReal, gvl2);
	__builtin_epi_vstore_2xf32(&resultImag[0][0], vresultImag, gvl2);

//	vehave_trace(1000, 0);

	printf("\nCovariance matrix: \n\n");
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < m; j++) {
			printf("( %f, %f )", resultReal[i][j], resultImag[i][j]);
		}
		printf("\n");
	}
}
*/
/*
void compute_resultMatrixVect(int m, int n, float diffReal[m][n], float diffImag[m][n], float resultReal[m][m], float resultImag[m][m])
{
//	vehave_trace(1000, 1);
	long gvl = __builtin_epi_vsetvl(128, __epi_e32, __epi_m1);
	long gvl2 = __builtin_epi_vsetvl(256, __epi_e32, __epi_m1);

	__epi_2xf32 vdiffReal = __builtin_epi_vload_2xf32(&diffReal[0][0], gvl);
	__epi_2xf32 vdiffImag = __builtin_epi_vload_2xf32(&diffImag[0][0], gvl);

	#if __riscv_vector_version == 700
		#define __builtin_epi_vmv_v_f_2xf32 __builtin_epi_vbroadcast_2xf32
	#endif

	__epi_2xf32 vresultReal = __builtin_epi_vbroadcast_2xf32(0.0f, gvl2);
	__epi_2xf32 vresultImag = __builtin_epi_vbroadcast_2xf32(0.0f, gvl2);

	for(int i = m-1; i >= 0; i--) {
		__epi_2xf32 vReal8elemsI = __builtin_epi_vslidedown_2xf32(vdiffReal, i*8, gvl);
		__epi_2xf32 vImag8elemsI = __builtin_epi_vslidedown_2xf32(vdiffImag, i*8, gvl);

		for(int j = m-1; j >=0; j--) {
			__epi_2xf32 vReal8elemsJ = __builtin_epi_vslidedown_2xf32(vdiffReal, j*8, gvl);
			__epi_2xf32 vImag8elemsJ = __builtin_epi_vslidedown_2xf32(vdiffImag, j*8, gvl);

			__epi_2xf32 tmp1 = __builtin_epi_vfmul_2xf32(vReal8elemsI, vReal8elemsJ, 8);

			tmp1 = __builtin_epi_vfmacc_2xf32(tmp1, vImag8elemsI, vImag8elemsJ, 8);

			__epi_2xf32 tmpaccum = __builtin_epi_vbroadcast_2xf32(0.0f, 8);

			tmp1 = __builtin_epi_vfredsum_2xf32(tmp1, tmpaccum, 8);

			float resultRealIJ = __builtin_epi_vgetfirst_2xf32(tmp1);
			resultReal[i][j] = resultRealIJ;

			__epi_2xf32 tmp2 = __builtin_epi_vfmul_2xf32(vReal8elemsI, vImag8elemsJ, 8);

			tmp2 = __builtin_epi_vfmsac_2xf32(tmp2, vImag8elemsI, vReal8elemsJ, 8);

			tmp2 = __builtin_epi_vfredsum_2xf32(tmp2, tmpaccum, 8);

			float resultImagIJ = __builtin_epi_vgetfirst_2xf32(tmp2);
			resultImag[i][j] = resultImagIJ;
		}
	}

	vresultReal = __builtin_epi_vload_2xf32(&resultReal[0][0], gvl2);
	vresultImag = __builtin_epi_vload_2xf32(&resultImag[0][0], gvl2);

	__epi_2xf32 vdivisorN = __builtin_epi_vbroadcast_2xf32(8.0f, gvl2);

	vresultReal = __builtin_epi_vfdiv_2xf32(vresultReal, vdivisorN, gvl2);
	vresultImag = __builtin_epi_vfdiv_2xf32(vresultImag, vdivisorN, gvl2);

	__builtin_epi_vstore_2xf32(&resultReal[0][0], vresultReal, gvl2);
	__builtin_epi_vstore_2xf32(&resultImag[0][0], vresultImag, gvl2);

//	vehave_trace(1000, 0);

	printf("\nCovariance matrix: \n\n");
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < m; j++) {
			printf("( %f, %f )", resultReal[i][j], resultImag[i][j]);
		}
		printf("\n");
	}
	
}
*/

