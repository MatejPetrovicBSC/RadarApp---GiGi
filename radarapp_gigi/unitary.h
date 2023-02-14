void XERBLA(char SRNAME[], int INFO);
int ILACLC(int M, int N, float complex A[M][N], int LDA);
int ILACLR(int M, int N, float complex A[M][N], int LDA);
bool LSAME(char CA, char CB);
void CGEMVU(char TRANS, int M, int N, float complex ALPHA, int LDA, float complex A[LDA][LDA], float complex X[M], int INCX, float complex BETA, float complex Y[N], int INCY);
void CGERCU(int M, int N, float complex ALPHA, int LDA, float complex X[LDA], int INCX, float complex Y[LDA], int INCY, float complex A[LDA][LDA]);
void CLARFU(char SIDE, int M, int N, int LDC, float complex V[M], int INCV, float complex TAU, float complex C[LDC][LDC], float complex WORK[LDC]);
void CGEMVD(char TRANS, int M, int N, float complex ALPHA, int LDA, float complex A[M][N], float complex X[M], int INCX, float complex BETA, float complex Y[N], int INCY);
void CGERCD(int M, int N, float complex ALPHA, int LDA, float complex X[M], int INCX, float complex Y[N], int INCY, float complex A[M][N]);
void CLARFD(char SIDE, int M, int N, int LDC, float complex V[M], int INCV, float complex TAU, float complex C[M][N], float complex WORK[LDC]);
void CSCAL(int N, float complex CA, float complex CX[N], int INCX);
void CUNG2R(int M, int N, int K, float complex A[N][N], int LDA, float complex TAU[K], float complex WORK[N], int INFO);
void CTRMM(char SIDE, char UPLO, char TRANSA, char DIAG, int M, int N, float complex ALPHA, float complex A[N][N], int LDA, float complex B[M][N], int LDB);
void CLACGV(int N, float complex X[N], int INCX);
void CCOPY(int N, float complex CX[N], int INCX, float complex CY[N], int INCY);
void CGEMM(char TRANSA, char TRANSB, int M, int N, int K, float complex ALPHA, int LDA, float complex A[LDA][LDA], float complex B[LDA][LDA], int LDB, float complex BETA, float complex C[LDA][LDA], int LDC);
void CLARFB(char SIDE, char TRANS, char DIRECT, char STOREV, int M, int N, int K, float complex V[K][K], int LDV, float complex T[K][K], int LDT, float complex C[N][N], int LDC, float complex WORK[K][K], int LDWORK);
void CTRMV(char UPLO, char TRANS, char DIAG, int N, float complex A[N][N], int LDA, float complex X[N], int INCX);
void CLARFT(char UPLO, char DIRECT, char STOREV, int N, int K, float complex V[N][N], int LDV, float complex TAU[N], float complex T[N][N], int LDT);
void CUNG2L(int M, int N, int K, int LDA, float complex A[LDA][LDA], float complex TAU[K], float complex WORK[LDA], int INFO);
void CUNGQL(int M, int N, int K, int LDA, float complex A[LDA][LDA], float complex TAU[K], int LWORK, float complex WORK[LWORK], int INFO);
void CUNGQR(int M, int N, int K, int LDA, float complex A[LDA][LDA], float complex TAU[K], int LWORK, float complex WORK[N], int INFO);
void CUNGTR(char UPLO, int N, float complex A[N][N], int LDA, float complex TAU[N], float complex WORK[N], int LWORK, int INFO);
void compute_unitary(int N, float complex A[N][N], float complex TAU[N-1]);