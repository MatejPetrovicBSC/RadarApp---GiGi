typedef struct
{
    int rows;              /* Number of rows */
    int cols;              /* Number of columns */
    float complex *array; /* Pointer to an array of type TYPE */
} matrix;
/*
matrix *create_matrix(int rows, int cols);
matrix *matrix_column_subtract(matrix *m1, int c1, matrix *m2, int c2);
matrix *create_matrix_from_array(int rows, int cols, float complex m[rows][cols]);
float complex vector_length(matrix *m, int column);
matrix *matrix_column_divide(matrix *m, int c, float complex k);
matrix *matrix_column_multiply(matrix *m, int c, float complex k);
void matrix_copy_column(matrix *msrc, int col1, matrix *mdst, int col2);
void print_matrix_test(matrix *m, bool isMatrixQ);
bool test_QR(int n, matrix *Q, matrix *R, matrix *A);
void compute_qr(int n, float complex m[n][n]);
*/
float complex *create_matrix(int rows, int cols);
void matrix_column_subtract(int l1, int n1, float complex m1[l1][n1], int c1, int l2, int n2, float complex m2[l2][n2], int c2);
float complex *create_matrix_from_array(int rows, int cols, float complex m[rows][cols]);
float complex vector_length(float complex m[16][16], int column);
void matrix_column_divide(float complex m[16][16], int c, float complex k);
void matrix_column_multiply(int l, int n, float complex m[l][n], int c, float complex k);
void matrix_copy_column(int m1, int n1, float complex msrc[m1][n1], int col1, int m2, int n2, float complex mdst[m2][n2], int col2);
void print_matrix_test(int n, float complex m[n][n], bool isMatrixQ);
bool test_QR(int n, float complex Q[n][n], float complex R[n][n], float complex A[n][n]);
void compute_qr(int n, float complex m[n][n]);
