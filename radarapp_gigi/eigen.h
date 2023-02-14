float Sign(float a, float b);
void TQL2(int N, float *D, float *E, float Z[N+1][N+1], int *IER);
void multiply_mat_vec(int n, float mat[n][n], float *vec, float *res);
void multiply_scalar(float scalar, float *array, float *res);
void transpose_matrix(int N, float tran[N+1][N+1], float Z[N+1][N+1]);
void compute_eigen(int n, float *D, float E[n]);
void shift_mat(float **Z);
