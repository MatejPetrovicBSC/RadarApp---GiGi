typedef struct Complex
{
	float re;
	float im;
} Complex;

void compute_covariance(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float covarianceReal[m][m], float covarianceImag[m][m]);
void compute_Means(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float meanReal[m], float meanImag[m]);
void compute_diffMatrix(int m, int n, Complex *snapshotMatrix, int ld_snapshot, float meanReal[m], float meanImag[m], float diffReal[m][n], float diffImag[m][n]);
void compute_resultMatrixVect(int m, int n, float diffReal[m][n], float diffImag[m][n], float covarianceReal[m][m], float covarianceImag[m][m]);
void compute_resultMatrix(int m, int n, float diffReal[m][n], float diffImag[m][n], float covarianceReal[m][m], float covarianceImag[m][m]);
