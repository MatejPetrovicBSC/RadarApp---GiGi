#include <stddef.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct
{
  float re, im;
} Complexnum;

/* Function Declarations */
static void b_xgemv(int m, int n, const Complexnum A_data[], int ia0, int lda,
                    const Complexnum x_data[], int ix0, Complexnum y_data[]);
static void b_xgerc(int m, int n, const Complexnum alpha1, int ix0, const Complexnum y_data[], Complexnum A_data[], int ia0, int lda);
static void b_xzlarf(int m, int n, int iv0, const Complexnum tau, Complexnum C_data[],
                     int ic0, int ldc, Complexnum work_data[]);
static Complexnum b_xzlarfg(Complexnum *alpha1, Complexnum *x);
static void c_sqrt(Complexnum *x);
static void diagDiagUpperHessNoImag(Complexnum D_data[], int D_size[2]);
static int eml_zlahqr(Complexnum h_data[], int h_size[2], Complexnum z_data[], int z_size[2]);
static int ilazlc(int m, int n, const Complexnum A_data[], int ia0, int lda);
static int ilazlr(int m, int n, const Complexnum A_data[], int ia0, int lda);
static bool ishermitian(const Complexnum A_data[], const int A_size[2]);
static Complexnum recip(const Complexnum y);
static float rt_hypotd_snf(float u0, float u1);
static void schur(const Complexnum A_data[], const int A_size[2], Complexnum V_data[],
                  int V_size[2], Complexnum T_data[], int T_size[2]);
static float xdlapy3(float x1, float x2, float x3);
static void xgehrd(Complexnum a_data[], int a_size[2], Complexnum tau_data[], int tau_size[1]);
static void xgemv(int m, int n, const Complexnum A_data[], int ia0, int lda, const Complexnum x_data[], int ix0, Complexnum y_data[]);
static void xgerc(int m, int n, const Complexnum alpha1, const Complexnum x_data[],
                  int iy0, Complexnum A_data[], int ia0, int lda);
static int xhseqr(Complexnum h_data[], int h_size[2], Complexnum z_data[], int z_size[2]);
static float xnrm2(int n, const Complexnum x_data[], int ix0);
static void xscal(int n, const Complexnum a, Complexnum x_data[], int ix0);
static void xungorghr(int n, int ihi, Complexnum A_data[], int A_size[2], int lda,
                      const Complexnum tau_data[]);
static void xzlarf(int m, int n, int iv0, const Complexnum tau, Complexnum C_data[],
                   int ic0, int ldc, Complexnum work_data[]);
static Complexnum xzlarfg(int n, Complexnum *alpha1, Complexnum x_data[], int ix0);


static void xzungqr(int m, int n, int k, Complexnum A_data[], int A_size[2], int ia0, int lda, const Complexnum tau_data[]);
extern void eigenValue(Complexnum *A_data, const int A_size[2], Complexnum *V_data, int V_size[2], Complexnum *D_data, int D_size[2]);
void compute_eigenv2(int m, float *covarianceReal, float *covarianceImag);
