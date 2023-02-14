#include <math.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "eigen_v2.h"
#include <stdbool.h>
//#include <vehave-control.h>
/*
 * Arguments    : int m
 *                int n
 *                const Complexnum A_data[]
 *                int ia0
 *                int lda
 *                const Complexnum x_data[]
 *                int ix0
 *                Complexnum y_data[]
 * Return Type  : void
 */
#define VECT 0

static void b_xgemv(int m, int n, const Complexnum A_data[], int ia0, int lda,
		const Complexnum x_data[], int ix0, Complexnum y_data[])
{
	int iy;
	int i11;
	int iac;
	int ix;
	float c_re;
	float c_im;
	int i12;
	int ia;
	if (n != 0)
	{
		for (iy = 1; iy <= n; iy++)
		{
			y_data[iy - 1].re = 0.0;
			y_data[iy - 1].im = 0.0;
		}

		iy = 0;
		i11 = ia0 + lda * (n - 1);
		iac = ia0;
		while ((lda > 0) && (iac <= i11))
		{
			ix = ix0 - 1;
			c_re = 0.0;
			c_im = 0.0;
			i12 = (iac + m) - 1;
		    
			#if VECT
			#pragma clang loop vectorize(enable)
			#endif
			for (ia = iac - 1; ia < i12; ia++)
			{
				c_re += A_data[ia].re * x_data[ix].re + A_data[ia].im * x_data[ix].im;
				c_im += A_data[ia].re * x_data[ix].im - A_data[ia].im * x_data[ix].re;
				ix++;
			}
			
			y_data[iy].re += c_re - 0.0 * c_im;
			y_data[iy].im += c_im + 0.0 * c_re;
			iy++;
			iac += lda;
		}
	}
}

/*
 * Arguments    : int m
 *                int n
 *                const Complexnum alpha1
 *                int ix0
 *                const Complexnum y_data[]
 *                Complexnum A_data[]
 *                int ia0
 *                int lda
 * Return Type  : void
 */
static void b_xgerc(int m, int n, const Complexnum alpha1, int ix0, const Complexnum y_data[], Complexnum A_data[], int ia0, int lda)
{

	int jA;
	int jy;
	int j;
	float temp_re;
	float temp_im;
	int ix;
	int i13;
	int ijA;
	float A_data_im;
	if (!((alpha1.re == 0.0) && (alpha1.im == 0.0)))
	{
		jA = ia0 - 1;
		jy = 0;
	        
		for (j = 1; j <= n; j++)
		{
			if ((y_data[jy].re != 0.0) || (y_data[jy].im != 0.0))
			{
				temp_re = y_data[jy].re * alpha1.re + y_data[jy].im * alpha1.im;
				temp_im = y_data[jy].re * alpha1.im - y_data[jy].im * alpha1.re;
				ix = ix0;
				i13 = m + jA;
				#if VECT
				#pragma clang loop vectorize(enable)
				#endif
				for (ijA = jA; ijA < i13; ijA++)
				{
					A_data_im = A_data[ix - 1].re * temp_im + A_data[ix - 1].im * temp_re;
					A_data[ijA].re += A_data[ix - 1].re * temp_re - A_data[ix - 1].im *
						temp_im;
					A_data[ijA].im += A_data_im;
					ix++;
				}

			}

			jy++;
			jA += lda;
		}
	}
}

/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                const Complexnum tau
 *                Complexnum C_data[]
 *                int ic0
 *                int ldc
 *                Complexnum work_data[]
 * Return Type  : void
 */
static void b_xzlarf(int m, int n, int iv0, const Complexnum tau, Complexnum C_data[],
		int ic0, int ldc, Complexnum work_data[])
{
	int lastv;
	int lastc;
	Complexnum b_tau;
	if ((tau.re != 0.0) || (tau.im != 0.0))
	{
		lastv = m;
		lastc = iv0 + m;
		while ((lastv > 0) && ((C_data[lastc - 2].re == 0.0) && (C_data[lastc - 2].im == 0.0)))
		{
			lastv--;
			lastc--;
		}

		lastc = ilazlc(lastv, n, C_data, ic0, ldc);
	}
	else
	{
		lastv = 0;
		lastc = 0;
	}

	if (lastv > 0)
	{
		b_xgemv(lastv, lastc, C_data, ic0, ldc, C_data, iv0, work_data);
		b_tau.re = -tau.re;
		b_tau.im = -tau.im;
		b_xgerc(lastv, lastc, b_tau, iv0, work_data, C_data, ic0, ldc);
	}
}

/*
 * Arguments    : Complexnum *alpha1
 *                Complexnum *x
 * Return Type  : Complexnum
 */
static Complexnum b_xzlarfg(Complexnum *alpha1, Complexnum *x)
{
	Complexnum tau;
	float xnorm;
	float beta1;
	int knt;
	float ai;
	Complexnum b_alpha1;
	float x_re;
	float x_im;
	int k;
	tau.re = 0.0;
	tau.im = 0.0;
	xnorm = rt_hypotd_snf(x->re, x->im);
	if ((xnorm != 0.0) || (alpha1->im != 0.0))
	{
		beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
		if (alpha1->re >= 0.0)
		{
			beta1 = -beta1;
		}

		if (fabs(beta1) < 1.0020841800044864E-292)
		{
			knt = 0;
			do
			{
				knt++;
				x->re *= 9.9792015476736E+291;
				x->im *= 9.9792015476736E+291;
				beta1 *= 9.9792015476736E+291;
				alpha1->re *= 9.9792015476736E+291;
				alpha1->im *= 9.9792015476736E+291;
			} while (!(fabs(beta1) >= 1.0020841800044864E-292));

			beta1 = xdlapy3(alpha1->re, alpha1->im, rt_hypotd_snf(x->re, x->im));
			if (alpha1->re >= 0.0)
			{
				beta1 = -beta1;
			}

			xnorm = beta1 - alpha1->re;
			ai = 0.0 - alpha1->im;
			if (ai == 0.0)
			{
				tau.re = xnorm / beta1;
				tau.im = 0.0;
			}
			else if (xnorm == 0.0)
			{
				tau.re = 0.0;
				tau.im = ai / beta1;
			}
			else
			{
				tau.re = xnorm / beta1;
				tau.im = ai / beta1;
			}

			b_alpha1.re = alpha1->re - beta1;
			b_alpha1.im = alpha1->im;
			*alpha1 = recip(b_alpha1);
			xnorm = alpha1->re;
			ai = alpha1->im;
			x_re = x->re;
			x_im = x->im;
			x->re = xnorm * x_re - ai * x_im;
			x->im = xnorm * x_im + ai * x_re;
	
			//float temp = 0;
			//looop not vectorized
			//#if VECT
			//#pragma clang loop vectorize(enable)
			//#endif
			for (k = 1; k <= knt; k++)
			{
				beta1 *= 1.0020841800044864E-292;
			}
			//beta1 = temp;
			alpha1->re = beta1;
			alpha1->im = 0.0;
		}
		else
		{
			xnorm = beta1 - alpha1->re;
			ai = 0.0 - alpha1->im;
			if (ai == 0.0)
			{
				tau.re = xnorm / beta1;
				tau.im = 0.0;
			}
			else if (xnorm == 0.0)
			{
				tau.re = 0.0;
				tau.im = ai / beta1;
			}
			else
			{
				tau.re = xnorm / beta1;
				tau.im = ai / beta1;
			}

			b_alpha1.re = alpha1->re - beta1;
			b_alpha1.im = alpha1->im;
			*alpha1 = recip(b_alpha1);
			xnorm = alpha1->re;
			ai = alpha1->im;
			x_re = x->re;
			x_im = x->im;
			x->re = xnorm * x_re - ai * x_im;
			x->im = xnorm * x_im + ai * x_re;
			alpha1->re = beta1;
			alpha1->im = 0.0;
		}
	}

	return tau;
}

/*
 * Arguments    : Complexnum *x
 * Return Type  : void
 */
static void c_sqrt(Complexnum *x)
{
	float xr;
	float xi;
	float absxi;
	float absxr;
	xr = x->re;
	xi = x->im;
	if (xi == 0.0)
	{
		if (xr < 0.0)
		{
			absxi = 0.0;
			xr = sqrt(-xr);
		}
		else
		{
			absxi = sqrt(xr);
			xr = 0.0;
		}
	}
	else if (xr == 0.0)
	{
		if (xi < 0.0)
		{
			absxi = sqrt(-xi / 2.0);
			xr = -absxi;
		}
		else
		{
			absxi = sqrt(xi / 2.0);
			xr = absxi;
		}
	}
	else
	{
		absxr = fabs(xr);
		absxi = fabs(xi);
		if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
		{
			absxr *= 0.5;
			absxi *= 0.5;
			absxi = rt_hypotd_snf(absxr, absxi);
			if (absxi > absxr)
			{
				absxi = sqrt(absxi) * sqrt(1.0 + absxr / absxi);
			}
			else
			{
				absxi = sqrt(absxi) * 1.4142135623730951;
			}
		}
		else
		{
			absxi = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
		}

		if (xr > 0.0)
		{
			xr = 0.5 * (xi / absxi);
		}
		else
		{
			if (xi < 0.0)
			{
				xr = -absxi;
			}
			else
			{
				xr = absxi;
			}

			absxi = 0.5 * (xi / xr);
		}
	}

	x->re = absxi;
	x->im = xr;
}

/*
 * Arguments    : Complexnum D_data[]
 *                int D_size[2]
 * Return Type  : void
 */
static void diagDiagUpperHessNoImag(Complexnum D_data[], int D_size[2])
{
	int n;
	int j;
	int i;
	n = D_size[1];
	D_data[0].im = 0.0;
	for (j = 1; j < n; j++)
	{
		D_data[j + D_size[1] * j].im = 0.0;
		D_data[j + D_size[1] * (j - 1)].re = 0.0;
		D_data[j + D_size[1] * (j - 1)].im = 0.0;
		for (i = 1; i <= j; i++)
		{
			D_data[(i + D_size[1] * j) - 1].re = 0.0;
			D_data[(i + D_size[1] * j) - 1].im = 0.0;
		}
	}
}

/*
 * Arguments    : Complexnum h_data[]
 *                int h_size[2]
 *                Complexnum z_data[]
 *                int z_size[2]
 * Return Type  : int
 */
static int eml_zlahqr(Complexnum h_data[], int h_size[2], Complexnum z_data[], int z_size[2])
{
	int info;
	int n;
	int L;
	int itmax;
	int ldh;
	int ldz;
	int j;
	int i;
	float SMLNUM;
	float tst;
	bool exitg1;
	float htmp1;
	float ba;
	Complexnum u2;
	bool goto140;
	bool exitg2;
	int k;
	bool exitg3;
	Complexnum y;
	float t_re;
	float ab;
	bool goto70;
	int m;
	float u_re;
	float u_im;
	float aa;
	float s;
	int b_k;
	float b_SMLNUM;
	Complexnum v[2];
	int c;
	int i14;
	float b_u_re;
	n = h_size[1];
	L = h_size[1];
	if (10 > L)
	{
		L = 10;
	}

	itmax = 30 * L;
	ldh = h_size[1];
	ldz = z_size[1];
	info = 0;
	if ((h_size[1] != 0) && (1 != h_size[1]))
	{

		for (j = 0; j < n - 3; j++)
		{
			h_data[(j + h_size[1] * j) + 2].re = 0.0;
			h_data[(j + h_size[1] * j) + 2].im = 0.0;
			h_data[(j + h_size[1] * j) + 3].re = 0.0;
			h_data[(j + h_size[1] * j) + 3].im = 0.0;
		}

		if (1 <= n - 2)
		{
			h_data[(n + h_size[1] * (n - 3)) - 1].re = 0.0;
			h_data[(n + h_size[1] * (n - 3)) - 1].im = 0.0;
		}

		SMLNUM = 2.2250738585072014E-308 * ((float)n / 2.2204460492503131E-16);
		i = n - 1;
		exitg1 = false;
		while ((!exitg1) && (i + 1 >= 1))
		{
			L = -1;
			goto140 = false;
			ldz = 0;
			exitg2 = false;
			while ((!exitg2) && (ldz <= itmax))
			{
				k = i;
				exitg3 = false;
				while ((!exitg3) && ((k + 1 > L + 2) && (!(fabs(h_data[k + h_size[1] *
										(k - 1)]
										.re) +
									fabs(h_data[k + h_size[1] * (k - 1)].im) <=
									SMLNUM))))
				{
					tst = (fabs(h_data[(k + h_size[1] * (k - 1)) - 1].re) + fabs(h_data[(k + h_size[1] * (k - 1)) - 1].im)) + (fabs(h_data[k + h_size[1] * k].re) + fabs(h_data[k + h_size[1] * k].im));
					if (tst == 0.0)
					{
						if (k - 1 >= 1)
						{
							tst = fabs(h_data[(k + h_size[1] * (k - 2)) - 1].re);
						}

						if (k + 2 <= n)
						{
							tst += fabs(h_data[(k + h_size[1] * k) + 1].re);
						}
					}

					if (fabs(h_data[k + h_size[1] * (k - 1)].re) <= 2.2204460492503131E-16 * tst)
					{
						htmp1 = fabs(h_data[k + h_size[1] * (k - 1)].re) + fabs(h_data[k +
								h_size[1] * (k - 1)]
								.im);
						tst = fabs(h_data[(k + h_size[1] * k) - 1].re) + fabs(h_data[(k +
									h_size[1] * k) -
								1]
								.im);
						if (htmp1 > tst)
						{
							ab = htmp1;
							ba = tst;
						}
						else
						{
							ab = tst;
							ba = htmp1;
						}

						htmp1 = fabs(h_data[k + h_size[1] * k].re) + fabs(h_data[k + h_size
								[1] *
								k]
								.im);
						tst = fabs(h_data[(k + h_size[1] * (k - 1)) - 1].re - h_data[k +
								h_size[1] * k]
								.re) +
							fabs(h_data[(k + h_size[1] * (k - 1)) - 1].im - h_data[k + h_size[1] * k].im);
						if (htmp1 > tst)
						{
							aa = htmp1;
							htmp1 = tst;
						}
						else
						{
							aa = tst;
						}

						s = aa + ab;
						tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
						if ((SMLNUM > tst))
						{
							b_SMLNUM = SMLNUM;
						}
						else
						{
							b_SMLNUM = tst;
						}

						if (ba * (ab / s) <= b_SMLNUM)
						{
							exitg3 = true;
						}
						else
						{
							k--;
						}
					}
					else
					{
						k--;
					}
				}

				L = k - 1;
				if (k + 1 > 1)
				{
					h_data[k + h_size[1] * (k - 1)].re = 0.0;
					h_data[k + h_size[1] * (k - 1)].im = 0.0;
				}

				if (k + 1 >= i + 1)
				{
					goto140 = true;
					exitg2 = true;
				}
				else
				{
					if (ldz == 10)
					{
						t_re = 0.75 * fabs(h_data[(k + h_size[1] * k) + 1].re) + h_data[k +
							h_size[1] * k]
							.re;
						ba = h_data[k + h_size[1] * k].im;
					}
					else if (ldz == 20)
					{
						t_re = 0.75 * fabs(h_data[i + h_size[1] * (i - 1)].re) + h_data[i +
							h_size[1] * i]
							.re;
						ba = h_data[i + h_size[1] * i].im;
					}
					else
					{
						t_re = h_data[i + h_size[1] * i].re;
						ba = h_data[i + h_size[1] * i].im;
						y = h_data[(i + h_size[1] * i) - 1];
						c_sqrt(&y);
						u2 = h_data[i + h_size[1] * (i - 1)];
						c_sqrt(&u2);
						u_re = y.re * u2.re - y.im * u2.im;
						u_im = y.re * u2.im + y.im * u2.re;
						s = fabs(u_re) + fabs(u_im);
						if (s != 0.0)
						{
							htmp1 = 0.5 * (h_data[(i + h_size[1] * (i - 1)) - 1].re - h_data[i + h_size[1] * i].re);
							ab = 0.5 * (h_data[(i + h_size[1] * (i - 1)) - 1].im - h_data[i +
									h_size[1] * i]
									.im);
							aa = fabs(htmp1) + fabs(ab);
							tst = fabs(htmp1) + fabs(ab);
							if (!((s > tst)))
							{
								s = tst;
							}

							if (ab == 0.0)
							{
								t_re = htmp1 / s;
								ba = 0.0;
							}
							else if (htmp1 == 0.0)
							{
								t_re = 0.0;
								ba = ab / s;
							}
							else
							{
								t_re = htmp1 / s;
								ba = ab / s;
							}

							tst = t_re;
							t_re = t_re * t_re - ba * ba;
							ba = tst * ba + ba * tst;
							if (u_im == 0.0)
							{
								u2.re = u_re / s;
								u2.im = 0.0;
							}
							else if (u_re == 0.0)
							{
								u2.re = 0.0;
								u2.im = u_im / s;
							}
							else
							{
								u2.re = u_re / s;
								u2.im = u_im / s;
							}

							y.re = t_re + (u2.re * u2.re - u2.im * u2.im);
							y.im = ba + (u2.re * u2.im + u2.im * u2.re);
							c_sqrt(&y);
							y.re *= s;
							y.im *= s;
							if (aa > 0.0)
							{
								if (ab == 0.0)
								{
									t_re = htmp1 / aa;
									ba = 0.0;
								}
								else if (htmp1 == 0.0)
								{
									t_re = 0.0;
									ba = ab / aa;
								}
								else
								{
									t_re = htmp1 / aa;
									ba = ab / aa;
								}

								if (t_re * y.re + ba * y.im < 0.0)
								{
									y.re = -y.re;
									y.im = -y.im;
								}
							}

							ba = htmp1 + y.re;
							aa = ab + y.im;
							if (aa == 0.0)
							{
								if (u_im == 0.0)
								{
									b_u_re = u_re / ba;
									tst = 0.0;
								}
								else if (u_re == 0.0)
								{
									b_u_re = 0.0;
									tst = u_im / ba;
								}
								else
								{
									b_u_re = u_re / ba;
									tst = u_im / ba;
								}
							}
							else if (ba == 0.0)
							{
								if (u_re == 0.0)
								{
									b_u_re = u_im / aa;
									tst = 0.0;
								}
								else if (u_im == 0.0)
								{
									b_u_re = 0.0;
									tst = -(u_re / aa);
								}
								else
								{
									b_u_re = u_im / aa;
									tst = -(u_re / aa);
								}
							}
							else
							{
								ab = fabs(ba);
								tst = fabs(aa);
								if (ab > tst)
								{
									s = aa / ba;
									tst = ba + s * aa;
									b_u_re = (u_re + s * u_im) / tst;
									tst = (u_im - s * u_re) / tst;
								}
								else if (tst == ab)
								{
									if (ba > 0.0)
									{
										htmp1 = 0.5;
									}
									else
									{
										htmp1 = -0.5;
									}

									if (aa > 0.0)
									{
										tst = 0.5;
									}
									else
									{
										tst = -0.5;
									}

									b_u_re = (u_re * htmp1 + u_im * tst) / ab;
									tst = (u_im * htmp1 - u_re * tst) / ab;
								}
								else
								{
									s = ba / aa;
									tst = aa + s * ba;
									b_u_re = (s * u_re + u_im) / tst;
									tst = (s * u_im - u_re) / tst;
								}
							}

							t_re = h_data[i + h_size[1] * i].re - (u_re * b_u_re - u_im * tst);
							ba = h_data[i + h_size[1] * i].im - (u_re * tst + u_im * b_u_re);
						}
					}

					goto70 = false;
					m = i;
					exitg3 = false;
					while ((!exitg3) && (m > k + 1))
					{
						u2.re = h_data[(m + h_size[1] * (m - 1)) - 1].re - t_re;
						u2.im = h_data[(m + h_size[1] * (m - 1)) - 1].im - ba;
						tst = h_data[m + h_size[1] * (m - 1)].re;
						s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
						if (u2.im == 0.0)
						{
							u2.re /= s;
							u2.im = 0.0;
						}
						else if (u2.re == 0.0)
						{
							u2.re = 0.0;
							u2.im /= s;
						}
						else
						{
							u2.re /= s;
							u2.im /= s;
						}

						tst /= s;
						v[0] = u2;
						v[1].re = tst;
						v[1].im = 0.0;
						if (fabs(h_data[(m + h_size[1] * (m - 2)) - 1].re) * fabs(tst) <=
								2.2204460492503131E-16 * ((fabs(u2.re) + fabs(u2.im)) * ((fabs(h_data[(m + h_size[1] * (m - 1)) - 1].re) + fabs(h_data[(m +
													h_size[1] * (m - 1)) -
												1]
												.im)) +
										(fabs(h_data[m + h_size[1] * m].re) + fabs(h_data[m + h_size[1] * m].im)))))
						{
							goto70 = true;
							exitg3 = true;
						}
						else
						{
							m--;
						}
					}

					if (!goto70)
					{
						u2.re = h_data[k + h_size[1] * k].re - t_re;
						u2.im = h_data[k + h_size[1] * k].im - ba;
						tst = h_data[(k + h_size[1] * k) + 1].re;
						s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
						if (u2.im == 0.0)
						{
							u2.re /= s;
							u2.im = 0.0;
						}
						else if (u2.re == 0.0)
						{
							u2.re = 0.0;
							u2.im /= s;
						}
						else
						{
							u2.re /= s;
							u2.im /= s;
						}

						tst /= s;
						v[0] = u2;
						v[1].re = tst;
						v[1].im = 0.0;
					}

		    		//	Extrae_eventandcounters(1000,3);
					for (b_k = m; b_k <= i; b_k++)
					{
						if (b_k > m)
						{
							v[0] = h_data[(b_k + h_size[1] * (b_k - 2)) - 1];
							v[1] = h_data[b_k + h_size[1] * (b_k - 2)];
						}

						u2 = b_xzlarfg(&v[0], &v[1]);
						if (b_k > m)
						{
							h_data[(b_k + h_size[1] * (b_k - 2)) - 1] = v[0];
							h_data[b_k + h_size[1] * (b_k - 2)].re = 0.0;
							h_data[b_k + h_size[1] * (b_k - 2)].im = 0.0;
						}

						tst = u2.re * v[1].re - u2.im * v[1].im;
//		    				Extrae_eventandcounters(1004,1);
						#if VECT
						#pragma clang loop vectorize(assume_safety)
						#endif
						for (j = b_k - 1; j < n; j++)
						{
							t_re = (u2.re * h_data[(b_k + h_size[1] * j) - 1].re - -u2.im *
									h_data[(b_k + h_size[1] * j) - 1].im) +
								tst * h_data[b_k +
								h_size[1] * j]
								.re;
							ba = (u2.re * h_data[(b_k + h_size[1] * j) - 1].im + -u2.im *
									h_data[(b_k + h_size[1] * j) - 1].re) +
								tst * h_data[b_k +
								h_size[1] * j]
								.im;
							h_data[(b_k + h_size[1] * j) - 1].re -= t_re;
							h_data[(b_k + h_size[1] * j) - 1].im -= ba;
							h_data[b_k + h_size[1] * j].re -= t_re * v[1].re - ba * v[1].im;
							h_data[b_k + h_size[1] * j].im -= t_re * v[1].im + ba * v[1].re;
						}
			
//		    				Extrae_eventandcounters(1004,0);
						if (b_k + 2 < i + 1)
						{
							i14 = b_k;
						}
						else
						{
							i14 = i - 1;
						}

//		    				Extrae_eventandcounters(1005,1);
						#if VECT
						#pragma clang loop vectorize(assume_safety)
						#endif
						for (j = 0; j < i14 + 2; j++)
						{
							t_re = (u2.re * h_data[j + h_size[1] * (b_k - 1)].re - u2.im *
									h_data[j + h_size[1] * (b_k - 1)].im) +
								tst * h_data[j +
								h_size[1] * b_k]
								.re;
							ba = (u2.re * h_data[j + h_size[1] * (b_k - 1)].im + u2.im *
									h_data[j + h_size[1] * (b_k - 1)].re) +
								tst * h_data[j +
								h_size[1] * b_k]
								.im;
							h_data[j + h_size[1] * (b_k - 1)].re -= t_re;
							h_data[j + h_size[1] * (b_k - 1)].im -= ba;
							h_data[j + h_size[1] * b_k].re -= t_re * v[1].re - ba * -v[1].im;
							h_data[j + h_size[1] * b_k].im -= t_re * -v[1].im + ba * v[1].re;
						}

//		    				Extrae_eventandcounters(1005,0);
//		    				Extrae_eventandcounters(1006,1);
						#if VECT
						#pragma clang loop vectorize(assume_safety)
						#endif
						for (j = 0; j < n; j++)
						{
							t_re = (u2.re * z_data[j + z_size[1] * (b_k - 1)].re - u2.im *
									z_data[j + z_size[1] * (b_k - 1)].im) +
								tst * z_data[j +
								z_size[1] * b_k]
								.re;
							ba = (u2.re * z_data[j + z_size[1] * (b_k - 1)].im + u2.im *
									z_data[j + z_size[1] * (b_k - 1)].re) +
								tst * z_data[j +
								z_size[1] * b_k]
								.im;
							z_data[j + z_size[1] * (b_k - 1)].re -= t_re;
							z_data[j + z_size[1] * (b_k - 1)].im -= ba;
							z_data[j + z_size[1] * b_k].re -= t_re * v[1].re - ba * -v[1].im;
							z_data[j + z_size[1] * b_k].im -= t_re * -v[1].im + ba * v[1].re;
						}

//		    				Extrae_eventandcounters(1006,0);
						if ((b_k == m) && (m > k + 1))
						{
							u2.re = 1.0 - u2.re;
							u2.im = 0.0 - u2.im;
							ba = rt_hypotd_snf(u2.re, u2.im);
							if (u2.im == 0.0)
							{
								u2.re /= ba;
								u2.im = 0.0;
							}
							else if (u2.re == 0.0)
							{
								u2.re = 0.0;
								u2.im /= ba;
							}
							else
							{
								u2.re /= ba;
								u2.im /= ba;
							}

							tst = h_data[m + h_size[1] * (m - 1)].re;
							htmp1 = h_data[m + h_size[1] * (m - 1)].im;
							h_data[m + h_size[1] * (m - 1)].re = tst * u2.re - htmp1 * -u2.im;
							h_data[m + h_size[1] * (m - 1)].im = tst * -u2.im + htmp1 * u2.re;
							if (m + 2 <= i + 1)
							{
								tst = h_data[(m + h_size[1] * m) + 1].re;
								htmp1 = h_data[(m + h_size[1] * m) + 1].im;
								h_data[(m + h_size[1] * m) + 1].re = tst * u2.re - htmp1 * u2.im;
								h_data[(m + h_size[1] * m) + 1].im = tst * u2.im + htmp1 * u2.re;
							}

							for (j = m; j <= i + 1; j++)
							{
								if (j != m + 1)
								{
									if (n > j)
									{
										c = j + j * ldh;
										if (!(ldh < 1))
										{
											i14 = c + ldh * ((n - j) - 1);
										        const int cc = i14;	
											//loop not vectorized
											//#if VECT
											//#pragma clang loop vectorize(assume_safety)
											//#endif
											while (c <= cc)
											{
												tst = h_data[c - 1].re;
												htmp1 = h_data[c - 1].im;
												h_data[c - 1].re = u2.re * tst - u2.im * htmp1;
												h_data[c - 1].im = u2.re * htmp1 + u2.im * tst;
												c += ldh;
											}
										}
									}

									c = (j - 1) * ldh;
									y.re = u2.re;
									y.im = -u2.im;
									i14 = (c + j) - 1;
//		    							Extrae_eventandcounters(1007,1);
									#if VECT
									#pragma clang loop vectorize(assume_safety)
									#endif
									while (c + 1 <= i14)
									{
										tst = h_data[c].re;
										htmp1 = h_data[c].im;
										h_data[c].re = y.re * tst - y.im * htmp1;
										h_data[c].im = y.re * htmp1 + y.im * tst;
										c++;
									}
//		    							Extrae_eventandcounters(1007,0);

									c = (j - 1) * ldh;
									y.re = u2.re;
									y.im = -u2.im;
									i14 = c + n;
//		    							Extrae_eventandcounters(1008,1);
									#if VECT
									#pragma clang loop vectorize(assume_safety)
									#endif
									while (c + 1 <= i14)
									{
										tst = z_data[c].re;
										htmp1 = z_data[c].im;
										z_data[c].re = y.re * tst - y.im * htmp1;
										z_data[c].im = y.re * htmp1 + y.im * tst;
										c++;
									}
//		    							Extrae_eventandcounters(1008,0);
								}
							}
						}
					}

	//	    			Extrae_eventandcounters(1000,0);

					u2 = h_data[i + h_size[1] * (i - 1)];
					if (h_data[i + h_size[1] * (i - 1)].im != 0.0)
					{
						tst = rt_hypotd_snf(h_data[i + h_size[1] * (i - 1)].re, h_data[i +
								h_size[1] * (i - 1)]
								.im);
						h_data[i + h_size[1] * (i - 1)].re = tst;
						h_data[i + h_size[1] * (i - 1)].im = 0.0;
						if (u2.im == 0.0)
						{
							u2.re /= tst;
							u2.im = 0.0;
						}
						else if (u2.re == 0.0)
						{
							u2.re = 0.0;
							u2.im /= tst;
						}
						else
						{
							u2.re /= tst;
							u2.im /= tst;
						}

						if (n > i + 1)
						{
							c = i + (i + 1) * ldh;
							y.re = u2.re;
							y.im = -u2.im;
							if (!(ldh < 1))
							{
								i14 = (c + ldh * ((n - i) - 2)) + 1;
								//loop not vectorized
								//#if VECT
								//#pragma clang loop vectorize(assume_safety)
								//#endif
								while (c + 1 <= i14)
								{
									tst = h_data[c].re;
									htmp1 = h_data[c].im;
									h_data[c].re = y.re * tst - y.im * htmp1;
									h_data[c].im = y.re * htmp1 + y.im * tst;
									c += ldh;
								}
							}
						}

						c = i * ldh;
						i14 = c + i;
	//	    				Extrae_eventandcounters(1000,4);
						#if VECT
						#pragma clang loop vectorize(assume_safety)
						#endif
						while (c + 1 <= i14)
						{
							tst = h_data[c].re;
							htmp1 = h_data[c].im;
							h_data[c].re = u2.re * tst - u2.im * htmp1;
							h_data[c].im = u2.re * htmp1 + u2.im * tst;
							c++;
						}

	//	    				Extrae_eventandcounters(1000,0);
						c = i * ldh;
						i14 = c + n;
	//	    				Extrae_eventandcounters(1000,5);
						#if VECT
						#pragma clang loop vectorize(assume_safety)
						#endif
						while (c + 1 <= i14)
						{
							tst = z_data[c].re;
							htmp1 = z_data[c].im;
							z_data[c].re = u2.re * tst - u2.im * htmp1;
							z_data[c].im = u2.re * htmp1 + u2.im * tst;
							c++;
						}
	//	    				Extrae_eventandcounters(1000,0);
					}

					ldz++;
				}
			}

			if (!goto140)
			{
				info = i + 1;
				exitg1 = true;
			}
			else
			{
				i = L;
			}
		}
	}

	return info;
}

/*
 * Arguments    : int m
 *                int n
 *                const Complexnum A_data[]
 *                int ia0
 *                int lda
 * Return Type  : int
 */
static int ilazlc(int m, int n, const Complexnum A_data[], int ia0, int lda)
{
	int j;
	bool exitg2;
	int coltop;
	int ia;
	int exitg1;
	j = n;
	exitg2 = false;
	while ((!exitg2) && (j > 0))
	{
		coltop = ia0 + (j - 1) * lda;
		ia = coltop;
		do
		{
			exitg1 = 0;
			if (ia <= (coltop + m) - 1)
			{
				if ((A_data[ia - 1].re != 0.0) || (A_data[ia - 1].im != 0.0))
				{
					exitg1 = 1;
				}
				else
				{
					ia++;
				}
			}
			else
			{
				j--;
				exitg1 = 2;
			}
		} while (exitg1 == 0);

		if (exitg1 == 1)
		{
			exitg2 = true;
		}
	}

	return j;
}

/*
 * Arguments    : int m
 *                int n
 *                const Complexnum A_data[]
 *                int ia0
 *                int lda
 * Return Type  : int
 */
static int ilazlr(int m, int n, const Complexnum A_data[], int ia0, int lda)
{
	int i;
	bool exitg2;
	int rowleft;
	int ia;
	int exitg1;
	i = m;
	exitg2 = false;
	while ((!exitg2) && (i > 0))
	{
		rowleft = (ia0 + i) - 1;
		ia = rowleft;
		do
		{
			exitg1 = 0;
			if ((lda > 0) && (ia <= rowleft + (n - 1) * lda))
			{
				if ((A_data[ia - 1].re != 0.0) || (A_data[ia - 1].im != 0.0))
				{
					exitg1 = 1;
				}
				else
				{
					ia += lda;
				}
			}
			else
			{
				i--;
				exitg1 = 2;
			}
		} while (exitg1 == 0);

		if (exitg1 == 1)
		{
			exitg2 = true;
		}
	}

	return i;
}

/*
 * Arguments    : const Complexnum A_data[]
 *                const int A_size[2]
 * Return Type  : bool
 */
static bool ishermitian(const Complexnum A_data[], const int A_size[2])
{
	bool p;
	int j;
	bool exitg2;
	int i;
	int exitg1;
	p = (A_size[1] == A_size[0]);
	if (p)
	{
		j = 0;
		exitg2 = false;
		while ((!exitg2) && (j <= A_size[0] - 1))
		{
			i = 0;
			do
			{
				exitg1 = 0;
				if (i <= j)
				{
					if (!((A_data[i + A_size[1] * j].re == A_data[j + A_size[1] * i].re) &&
								(A_data[i + A_size[1] * j].im == -A_data[j + A_size[1] * i].im)))
					{
						p = false;
						exitg1 = 1;
					}
					else
					{
						i++;
					}
				}
				else
				{
					j++;
					exitg1 = 2;
				}
			} while (exitg1 == 0);

			if (exitg1 == 1)
			{
				exitg2 = true;
			}
		}
	}

	return p;
}

/*
 * Arguments    : const Complexnum y
 * Return Type  : Complexnum
 */
static Complexnum recip(const Complexnum y)
{
	Complexnum z;
	float brm;
	float bim;
	float d;
	brm = fabs(y.re);
	bim = fabs(y.im);
	if (y.im == 0.0)
	{
		z.re = 1.0 / y.re;
		z.im = 0.0;
	}
	else if (y.re == 0.0)
	{
		z.re = 0.0;
		z.im = -1.0 / y.im;
	}
	else if (brm > bim)
	{
		bim = y.im / y.re;
		d = y.re + bim * y.im;
		z.re = 1.0 / d;
		z.im = -bim / d;
	}
	else if (brm == bim)
	{
		bim = 0.5;
		if (y.re < 0.0)
		{
			bim = -0.5;
		}

		d = 0.5;
		if (y.im < 0.0)
		{
			d = -0.5;
		}

		z.re = bim / brm;
		z.im = -d / brm;
	}
	else
	{
		bim = y.re / y.im;
		d = y.im + bim * y.re;
		z.re = bim / d;
		z.im = -1.0 / d;
	}

	return z;
}

/*
 * Arguments    : float u0
 *                float u1
 * Return Type  : float
 */
static float rt_hypotd_snf(float u0, float u1)
{
	float y;
	float a;
	float b;
	a = fabs(u0);
	b = fabs(u1);
	if (a < b)
	{
		a /= b;
		y = b * sqrt(a * a + 1.0);
	}
	else if (a > b)
	{
		b /= a;
		y = a * sqrt(b * b + 1.0);
	}
	else
	{
		y = a * 1.4142135623730951;
	}

	return y;
}

/*
 * Arguments    : const Complexnum A_data[]
 *                const int A_size[2]
 *                Complexnum V_data[]
 *                int V_size[2]
 *                Complexnum T_data[]
 *                int T_size[2]
 * Return Type  : void
 */
static void schur(const Complexnum A_data[], const int A_size[2], Complexnum V_data[],
		int V_size[2], Complexnum T_data[], int T_size[2])
{
	
	int i2;
	int loop_ub;
	signed char iv1[2];
	Complexnum tau_data[15];
	int tau_size[1];

	T_size[1] = A_size[1];
	T_size[0] = A_size[0];
	loop_ub = A_size[1] * A_size[0];
	if (0 <= loop_ub - 1)
	{
		memcpy(&T_data[0], &A_data[0], (unsigned int)(loop_ub * (int)sizeof(Complexnum)));
	}

	xgehrd(T_data, T_size, tau_data, tau_size);
	V_size[1] = T_size[1];
	V_size[0] = T_size[0];
	loop_ub = T_size[1] * T_size[0];
	if (0 <= loop_ub - 1)
	{
		memcpy(&V_data[0], &T_data[0], (unsigned int)(loop_ub * (int)sizeof(Complexnum)));
	}

	xungorghr(A_size[1], A_size[1], V_data, V_size, A_size[1], tau_data);
	xhseqr(T_data, T_size, V_data, V_size);
}

/*
 * Arguments    : float x1
 *                float x2
 *                float x3
 * Return Type  : float
 */
static float xdlapy3(float x1, float x2, float x3)
{
	float y;
	float a;
	float b;
	float c;
	a = fabs(x1);
	b = fabs(x2);
	c = fabs(x3);
	if ((a > b))
	{
		y = a;
	}
	else
	{
		y = b;
	}

	if (c > y)
	{
		y = c;
	}

	if ((y > 0.0))
	{
		a /= y;
		b /= y;
		c /= y;
		y *= sqrt((a * a + c * c) + b * b);
	}
	else
	{
		y = (a + b) + c;
	}

	return y;
}

/*
 * Arguments    : Complexnum a_data[]
 *                int a_size[2]
 *                Complexnum tau_data[]
 *                int tau_size[1]
 * Return Type  : void
 */
static void xgehrd(Complexnum a_data[], int a_size[2], Complexnum tau_data[], int tau_size[1])
{
	int n;
	int loop_ub;
	int im1n;
	Complexnum work_data[16];
	int in;
	Complexnum alpha1;
	int u0;
	Complexnum b_tau_data;
	n = a_size[1];
	if (a_size[1] < 1)
	{
		tau_size[0] = 0;
	}
	else
	{
		tau_size[0] = (signed char)(a_size[1] - 1);
	}

	loop_ub = (signed char)a_size[1];
	for (im1n = 0; im1n < loop_ub; im1n++)
	{
		work_data[im1n].re = 0.0;
		work_data[im1n].im = 0.0;
	}

	//loop not vectorized
	//#if VECT
	//#pragma clang loop vectorize(assume_safety)
 	//#endif
	for (loop_ub = 1; loop_ub < n; loop_ub++)
	{
		im1n = (loop_ub - 1) * n + 1;
		in = loop_ub * n + 1;
		alpha1 = a_data[loop_ub + a_size[1] * (loop_ub - 1)];
		u0 = loop_ub + 2;
		if (!(u0 < n))
		{
			u0 = n;
		}

		tau_data[loop_ub - 1] = xzlarfg(n - loop_ub, &alpha1, a_data, u0 + (loop_ub - 1) * n);
		a_data[loop_ub + a_size[1] * (loop_ub - 1)].re = 1.0;
		a_data[loop_ub + a_size[1] * (loop_ub - 1)].im = 0.0;
		xzlarf(n, n - loop_ub, loop_ub + im1n, tau_data[loop_ub - 1], a_data, in, n,
				work_data);
		b_tau_data.re = tau_data[loop_ub - 1].re;
		b_tau_data.im = -tau_data[loop_ub - 1].im;
		b_xzlarf(n - loop_ub, n - loop_ub, loop_ub + im1n, b_tau_data, a_data,
				loop_ub + in, n, work_data);
		a_data[loop_ub + a_size[1] * (loop_ub - 1)] = alpha1;
	}
}

/*
 * Arguments    : int m
 *                int n
 *                const Complexnum A_data[]
 *                int ia0
 *                int lda
 *                const Complexnum x_data[]
 *                int ix0
 *                Complexnum y_data[]
 * Return Type  : void
 */
static void xgemv(int m, int n, const Complexnum A_data[], int ia0, int lda, const Complexnum x_data[], int ix0, Complexnum y_data[])
{
	int iy;
	int ix;
	int i8;
	int iac;
	float c_re;
	float c_im;
	int i9;
	int ia;
	if (m != 0)
	{
		for (iy = 1; iy <= m; iy++)
		{
			y_data[iy - 1].re = 0.0;
			y_data[iy - 1].im = 0.0;
		}

		ix = ix0;
		i8 = ia0 + lda * (n - 1);
		iac = ia0;
		
	//	Extrae_eventandcounters(1000,6);
		while ((lda > 0) && (iac <= i8))
		{
			c_re = x_data[ix - 1].re - 0.0 * x_data[ix - 1].im;
			c_im = x_data[ix - 1].im + 0.0 * x_data[ix - 1].re;
			iy = 0;
			i9 = (iac + m) - 1;

			#if VECT
			#pragma clang loop vectorize(enable)
			#endif
			for (ia = iac; ia <= i9; ia++)
			{
				y_data[iy].re += A_data[ia - 1].re * c_re - A_data[ia - 1].im * c_im;
				y_data[iy].im += A_data[ia - 1].re * c_im + A_data[ia - 1].im * c_re;
				iy++;
			}


			ix++;
			iac += lda;
		}
	//	Extrae_eventandcounters(1000,0);
	}
}

/*
 * Arguments    : int m
 *                int n
 *                const Complexnum alpha1
 *                const Complexnum x_data[]
 *                int iy0
 *                Complexnum A_data[]
 *                int ia0
 *                int lda
 * Return Type  : void
 */
static void xgerc(int m, int n, const Complexnum alpha1, const Complexnum x_data[],
		int iy0, Complexnum A_data[], int ia0, int lda)
{
	int jA;
	int jy;
	int j;
	float temp_re;
	float temp_im;
	int ix;
	int i10;
	int ijA;
	if (!((alpha1.re == 0.0) && (alpha1.im == 0.0)))
	{
		jA = ia0 - 1;
		jy = iy0 - 1;
		
	//	Extrae_eventandcounters(1000,7);
		for (j = 1; j <= n; j++)
		{
			if ((A_data[jy].re != 0.0) || (A_data[jy].im != 0.0))
			{
				temp_re = A_data[jy].re * alpha1.re + A_data[jy].im * alpha1.im;
				temp_im = A_data[jy].re * alpha1.im - A_data[jy].im * alpha1.re;
				ix = 0;
				i10 = m + jA;
//		    		Extrae_eventandcounters(1012,1);
				#if VECT
				#pragma clang loop vectorize(enable)
				#endif
				for (ijA = jA; ijA < i10; ijA++)
				{
					A_data[ijA].re += x_data[ix].re * temp_re - x_data[ix].im * temp_im;
					A_data[ijA].im += x_data[ix].re * temp_im + x_data[ix].im * temp_re;
					ix++;
				}
//		    		Extrae_eventandcounters(1012,0);

			}

			jy++;
			jA += lda;
		}
	//	Extrae_eventandcounters(1000,0);
	}
}

/*
 * Arguments    : Complexnum h_data[]
 *                int h_size[2]
 *                Complexnum z_data[]
 *                int z_size[2]
 * Return Type  : int
 */
static int xhseqr(Complexnum h_data[], int h_size[2], Complexnum z_data[], int z_size[2])
{
	int info;
	int m;
	int istart;
	int jend;
	int j;
	int i;
	info = eml_zlahqr(h_data, h_size, z_data, z_size);
	m = h_size[1];
	if ((h_size[1] == 0) || (h_size[0] == 0) || (3 >= h_size[1]))
	{
	}
	else
	{
		istart = 4;
		if (h_size[1] - 4 < h_size[0] - 1)
		{
			jend = h_size[1] - 3;
		}
		else
		{
			jend = h_size[0];
		}

		for (j = 1; j <= jend; j++)
		{
			for (i = istart; i <= m; i++)
			{
				h_data[(i + h_size[1] * (j - 1)) - 1].re = 0.0;
				h_data[(i + h_size[1] * (j - 1)) - 1].im = 0.0;
			}

			istart++;
		}
	}

	return info;
}

/*
 * Arguments    : int n
 *                const Complexnum x_data[]
 *                int ix0
 * Return Type  : float
 */
static float xnrm2(int n, const Complexnum x_data[], int ix0)
{
	float y;
	float scale;
	int kend;
	int k;
	float absxk;
	float t;
	y = 0.0;
	if (!(n < 1))
	{
		if (n == 1)
		{
			y = rt_hypotd_snf(x_data[ix0 - 1].re, x_data[ix0 - 1].im);
		}
		else
		{
			scale = 3.3121686421112381E-170;
			kend = (ix0 + n) - 1;
			for (k = ix0; k <= kend; k++)
			{
				absxk = fabs(x_data[k - 1].re);
				if (absxk > scale)
				{
					t = scale / absxk;
					y = 1.0 + y * t * t;
					scale = absxk;
				}
				else
				{
					t = absxk / scale;
					y += t * t;
				}

				absxk = fabs(x_data[k - 1].im);
				if (absxk > scale)
				{
					t = scale / absxk;
					y = 1.0 + y * t * t;
					scale = absxk;
				}
				else
				{
					t = absxk / scale;
					y += t * t;
				}
			}

			y = scale * sqrt(y);
		}
	}

	return y;
}

/*
 * Arguments    : int n
 *                const Complexnum a
 *                Complexnum x_data[]
 *                int ix0
 * Return Type  : void
 */
static void xscal(int n, const Complexnum a, Complexnum x_data[], int ix0)
{
	int i7;
	int k;
	float x_data_re;
	float x_data_im;
	i7 = (ix0 + n) - 1;


//	Extrae_eventandcounters(1000,8);
	#if VECT
	#pragma clang loop vectorize(enable)
	#endif
	for (k = ix0; k <= i7; k++)
	{
		x_data_re = x_data[k - 1].re;
		x_data_im = x_data[k - 1].im;
		x_data[k - 1].re = a.re * x_data_re - a.im * x_data_im;
		x_data[k - 1].im = a.re * x_data_im + a.im * x_data_re;
	}
//	Extrae_eventandcounters(1000,0);
}

/*
 * Arguments    : int n
 *                int ihi
 *                Complexnum A_data[]
 *                int A_size[2]
 *                int lda
 *                const Complexnum tau_data[]
 * Return Type  : void
 */
static void xungorghr(int n, int ihi, Complexnum A_data[], int A_size[2], int lda,
		const Complexnum tau_data[])
{
	int j;
	int i;
	int c;
	int iajm1;
	if (n != 0)
	{
		for (j = ihi; j > 1; j--)
		{
			c = (j - 1) * lda;
			for (i = 1; i < j; i++)
			{
				A_data[(c + i) - 1].re = 0.0;
				A_data[(c + i) - 1].im = 0.0;
			}

			iajm1 = c - lda;
			for (i = j; i < ihi; i++)
			{
				A_data[c + i] = A_data[iajm1 + i];
			}

			for (i = ihi; i < n; i++)
			{
				A_data[c + i].re = 0.0;
				A_data[c + i].im = 0.0;
			}
		}

		for (i = 1; i <= n; i++)
		{
			A_data[i - 1].re = 0.0;
			A_data[i - 1].im = 0.0;
		}

		A_data[0].re = 1.0;
		A_data[0].im = 0.0;
		for (j = ihi; j < n; j++)
		{
			c = j * lda;
			for (i = 1; i <= n; i++)
			{
				A_data[(c + i) - 1].re = 0.0;
				A_data[(c + i) - 1].im = 0.0;
			}

			A_data[c + j].re = 1.0;
			A_data[c + j].im = 0.0;
		}

		xzungqr(ihi - 1, ihi - 1, ihi - 1, A_data, A_size, 2 + lda, lda, tau_data);
	}
}
/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                const Complexnum tau
 *                Complexnum C_data[]
 *                int ic0
 *                int ldc
 *                Complexnum work_data[]
 * Return Type  : void
 */
static void xzlarf(int m, int n, int iv0, const Complexnum tau, Complexnum C_data[],
		int ic0, int ldc, Complexnum work_data[])
{
	int lastv;
	int lastc;
	Complexnum b_tau;
	if ((tau.re != 0.0) || (tau.im != 0.0))
	{
		lastv = n;
		lastc = iv0 + n;
		while ((lastv > 0) && ((C_data[lastc - 2].re == 0.0) && (C_data[lastc - 2].im == 0.0)))
		{
			lastv--;
			lastc--;
		}

		lastc = ilazlr(m, lastv, C_data, ic0, ldc);
	}
	else
	{
		lastv = 0;
		lastc = 0;
	}

	if (lastv > 0)
	{
		xgemv(lastc, lastv, C_data, ic0, ldc, C_data, iv0, work_data);
		b_tau.re = -tau.re;
		b_tau.im = -tau.im;
		xgerc(lastc, lastv, b_tau, work_data, iv0, C_data, ic0, ldc);
	}
}

/*
 * Arguments    : int n
 *                Complexnum *alpha1
 *                Complexnum x_data[]
 *                int ix0
 * Return Type  : Complexnum
 */
static Complexnum xzlarfg(int n, Complexnum *alpha1, Complexnum x_data[], int ix0)
{
	Complexnum tau;
	float xnorm;
	float beta1;
	int knt;
	float ai;
	int i6;
	int k;
	Complexnum b_alpha1;
	float x_data_re;
	float x_data_im;
	tau.re = 0.0;
	tau.im = 0.0;
	if (!(n <= 0))
	{
		xnorm = xnrm2(n - 1, x_data, ix0);
		if ((xnorm != 0.0) || (alpha1->im != 0.0))
		{
			beta1 = xdlapy3(alpha1->re, alpha1->im, xnorm);
			if (alpha1->re >= 0.0)
			{
				beta1 = -beta1;
			}

			if (fabs(beta1) < 1.0020841800044864E-292)
			{
				knt = 0;
				i6 = (ix0 + n) - 2;
				do
				{
					knt++;
					//loop vectorized but trace file not created; floating point warning when compiled
					#if VECT
					#pragma clang loop vectorize(enable)
					#endif
					for (k = ix0; k <= i6; k++)
					{
						x_data_re = x_data[k - 1].re;
						x_data_im = x_data[k - 1].im;
						x_data[k - 1].re = 9.9792015476736E+291 * x_data_re - 0.0 *
							x_data_im;
						x_data[k - 1].im = 9.9792015476736E+291 * x_data_im + 0.0 *
							x_data_re;
					}

					beta1 *= 9.9792015476736E+291;
					alpha1->re *= 9.9792015476736E+291;
					alpha1->im *= 9.9792015476736E+291;
				} while (!(fabs(beta1) >= 1.0020841800044864E-292));

				beta1 = xdlapy3(alpha1->re, alpha1->im, xnrm2(n - 1, x_data, ix0));
				if (alpha1->re >= 0.0)
				{
					beta1 = -beta1;
				}

				xnorm = beta1 - alpha1->re;
				ai = 0.0 - alpha1->im;
				if (ai == 0.0)
				{
					tau.re = xnorm / beta1;
					tau.im = 0.0;
				}
				else if (xnorm == 0.0)
				{
					tau.re = 0.0;
					tau.im = ai / beta1;
				}
				else
				{
					tau.re = xnorm / beta1;
					tau.im = ai / beta1;
				}

				b_alpha1.re = alpha1->re - beta1;
				b_alpha1.im = alpha1->im;
				*alpha1 = recip(b_alpha1);
				i6 = (ix0 + n) - 2;
				//trace not created!
				//#if VECT
				//#pragma clang loop vectorize(enable)
				//#endif
				for (k = ix0; k <= i6; k++)
				{
					xnorm = alpha1->re;
					ai = alpha1->im;
					x_data_re = x_data[k - 1].re;
					x_data_im = x_data[k - 1].im;
					x_data[k - 1].re = xnorm * x_data_re - ai * x_data_im;
					x_data[k - 1].im = xnorm * x_data_im + ai * x_data_re;
				}

				for (k = 1; k <= knt; k++)
				{
					beta1 *= 1.0020841800044864E-292;
				}

				alpha1->re = beta1;
				alpha1->im = 0.0;
			}
			else
			{
				xnorm = beta1 - alpha1->re;
				ai = 0.0 - alpha1->im;
				if (ai == 0.0)
				{
					tau.re = xnorm / beta1;
					tau.im = 0.0;
				}
				else if (xnorm == 0.0)
				{
					tau.re = 0.0;
					tau.im = ai / beta1;
				}
				else
				{
					tau.re = xnorm / beta1;
					tau.im = ai / beta1;
				}

				b_alpha1.re = alpha1->re - beta1;
				b_alpha1.im = alpha1->im;
				*alpha1 = recip(b_alpha1);
				i6 = (ix0 + n) - 2;
//				Extrae_eventandcounters(1000,9);
				#if VECT
				#pragma clang loop vectorize(enable)
				#endif
				for (k = ix0; k <= i6; k++)
				{
					xnorm = alpha1->re;
					ai = alpha1->im;
					x_data_re = x_data[k - 1].re;
					x_data_im = x_data[k - 1].im;
					x_data[k - 1].re = xnorm * x_data_re - ai * x_data_im;
					x_data[k - 1].im = xnorm * x_data_im + ai * x_data_re;
				}
//				Extrae_eventandcounters(1000,0);
				alpha1->re = beta1;
				alpha1->im = 0.0;
			}
		}
	}

	return tau;
}

/*
 * Arguments    : int m
 *                int n
 *                int k
 *                Complexnum A_data[]
 *                int A_size[2]
 *                int ia0
 *                int lda
 *                const Complexnum tau_data[]
 * Return Type  : void
 */
static void xzungqr(int m, int n, int k, Complexnum A_data[], int A_size[2], int ia0, int lda, const Complexnum tau_data[])
{
	int j;
	int itau;
	int ia;
	int i;
	int iaii;
	Complexnum work_data[16];
	Complexnum b_tau_data;
	if (!(n < 1))
	{
		for (j = k; j < n; j++)
		{
			ia = ia0 + j * lda;
			for (i = 0; i < m; i++)
			{
				A_data[(ia + i) - 1].re = 0.0;
				A_data[(ia + i) - 1].im = 0.0;
			}

			A_data[(ia + j) - 1].re = 1.0;
			A_data[(ia + j) - 1].im = 0.0;
		}

		itau = k - 1;
		ia = (signed char)A_size[0];
		for (iaii = 0; iaii < ia; iaii++)
		{
			work_data[iaii].re = 0.0;
			work_data[iaii].im = 0.0;
		}

		for (i = k; i >= 1; i--)
		{
			iaii = ((ia0 + i) + (i - 1) * lda) - 2;
			if (i < n)
			{
				A_data[iaii].re = 1.0;
				A_data[iaii].im = 0.0;
				b_xzlarf((m - i) + 1, n - i, iaii + 1, tau_data[itau], A_data, (iaii + lda) + 1, lda, work_data);
			}

			if (i < m)
			{
				b_tau_data.re = -tau_data[itau].re;
				b_tau_data.im = -tau_data[itau].im;
				xscal(m - i, b_tau_data, A_data, iaii + 2);
			}

			A_data[iaii].re = 1.0 - tau_data[itau].re;
			A_data[iaii].im = 0.0 - tau_data[itau].im;
			for (j = 1; j < i; j++)
			{
				A_data[iaii - j].re = 0.0;
				A_data[iaii - j].im = 0.0;
			}

			itau--;
		}
	}
}

/*
 * Arguments    : const Complexnum A_data[]
 *                const int A_size[2]
 *                Complexnum V_data[]
 *                int V_size[2]
 *                Complexnum D_data[]
 *                int D_size[2]
 * Return Type  : void
 */
void eigenValue(Complexnum *A_data, const int A_size[2], Complexnum V_data[],
		int V_size[2], Complexnum D_data[], int D_size[2])
{
	static Complexnum dc0 = {
		0.0, /* re */
		0.0  /* im */
	};

	int b_A_size[2];
	int info;
	int i0;
	int loop_ub;
	Complexnum b_A_data[256];
	int b_V_size[2];
	int i1;
	Complexnum alpha1_data[16];
	int alpha1_size[1];
	Complexnum beta1_data[16];
	int beta1_size[1];
	Complexnum b_V_data[256];
	signed char iv0[2];
	b_A_size[1] = A_size[0];
	b_A_size[0] = A_size[1];
	info = A_size[1];

	for (i0 = 0; i0 < info; i0++)
	{
		loop_ub = A_size[0];
		for (i1 = 0; i1 < loop_ub; i1++)
		{
			b_A_data[i1 + b_A_size[1] * i0] = A_data[i0 + A_size[1] * i1];
		}
	}
	if (ishermitian(b_A_data, b_A_size))
	{
		schur(b_A_data, b_A_size, b_V_data, b_V_size, D_data, D_size);
		b_A_size[1] = D_size[1];
		b_A_size[0] = D_size[0];
		info = D_size[1] * D_size[0];
		if (0 <= info - 1)
		{
			memcpy(&b_A_data[0], &D_data[0], (unsigned int)(info * (int)sizeof(Complexnum)));
		}

		diagDiagUpperHessNoImag(b_A_data, b_A_size);
	}

	D_size[1] = b_A_size[0];
	D_size[0] = b_A_size[1];
	info = b_A_size[1];
	for (i0 = 0; i0 < info; i0++)
	{
		loop_ub = b_A_size[0];
		for (i1 = 0; i1 < loop_ub; i1++)
		{
			D_data[i1 + D_size[1] * i0] = b_A_data[i0 + b_A_size[1] * i1];
		}
	}

	V_size[1] = b_V_size[0];
	V_size[0] = b_V_size[1];
	info = b_V_size[1];
	for (i0 = 0; i0 < info; i0++)
	{
		loop_ub = b_V_size[0];
		for (i1 = 0; i1 < loop_ub; i1++)
		{
			//printf("(%f, %f ) ", b_V_data[i0 + b_V_size[1] * i1].re, b_V_data[i0 + b_V_size[1] * i1].im);
			V_data[i1 + V_size[1] * i0] = b_V_data[i0 + b_V_size[1] * i1];
		}
	}
}

Complexnum* read_snapshotMatrixzz(int m, int n, int ld_snapshot)
{
	Complexnum *snapshotMatrix;
	int flag = 0;
	int i, k, num;
	float real, imag;
	char sign;

	FILE *inputFile;
	inputFile = fopen("matrices/Covariance_1.dat", "r");

	snapshotMatrix = (Complexnum *)malloc(sizeof(Complexnum)*m*n);

	if (inputFile == NULL)
	{
		fprintf(stderr, "Can't open the file !\n");
		exit(1);
	}

	real = 0;
	imag = 0;


	if (flag == EOF)
	{
		fprintf(stderr, "Error has been encountered while reading values from twidder_factor.txt!\n");
		exit(1);
	}

	for (i = 0; i < m; i++)
	{
		for (k = 0; k < n; k++)
		{
			if ((num = fscanf(inputFile, "%f %c %fj", &real, &sign, &imag)) > 0)
			{
				if (sign == '-')
					imag *= -1;
				snapshotMatrix[i * ld_snapshot + k].re = real;
				snapshotMatrix[i * ld_snapshot + k].im = imag;
			}
			else
			{
				k--;
			}
		}
		printf("\n");
	}
	return snapshotMatrix;
}

void compute_eigenv2(int m, float *covarianceReal, float *covarianceImag)
{
	Complexnum *A_data;
	int A_size[2];
	Complexnum V_data[m*m];
	int V_size[2];
	Complexnum D_data[m*m];
	int D_size[2];

	int i, k;

	A_data = (Complexnum *)malloc(sizeof(Complexnum)*m*m);
//	V_data = (Complexnum *)malloc(sizeof(Complexnum)*16*16);
//	D_data = (Complexnum *)malloc(sizeof(Complexnum)*16*16);

	//works with reading covariance matrix but not with passed; fix!!
	Complexnum *A_data1;
	A_data1 = read_snapshotMatrixzz(16,16,16);
	
	for(i=0; i < m; i++) 
	{
		for(k=0; k < m; k++)
		{
//			eigenMatrix[i][k] = covarianceReal[i][k] + covarianceImag[i][k] * I;
//			A_data[i*m + k].re = creal(eigenMatrix[i][k]); 
//			A_data[i*m + k].im = cimag(eigenMatrix[i][k]); 
					
			A_data[i*m + k].re = covarianceReal[i*m + k]; 
			A_data[i*m + k].im = covarianceImag[i*m + k];
		}
	}	
/*
	printf("\n");

	for(i=0; i < m; i++) 
	{
		for(k=0; k < m; k++)
		{
			printf("(%f %f)", A_data1[i*m + k].re, A_data1[i*m + k].im);
		}
		printf("\n");
	}	


	printf("\n");
	for(i=0; i < m; i++) 
	{
		for(k=0; k < m; k++)
		{
			printf("(%f %f)", A_data[i*m + k].re, A_data[i*m + k].im);
		}
		printf("\n");
	}
*/
	A_size[0] = m; //need fixing
	A_size[1] = m; //need fixing

	/* Call the entry-point 'eigenValue'. */
	
	
	eigenValue(A_data1, A_size, V_data, V_size, D_data, D_size);
    	
	printf("Output for Eigen Vector...\n");
	for (i = m-1; i >= 0; i--)
	{ //need fixing 16 as size
		for (k = m-1; k >= 0; k--)
		{ //need fixing 16 as size
			printf("%lf ", V_data[i * m + k].re);
			if (!(V_data[i * m + k].im < 0))
				printf("+");
			else {
				printf("-");
				V_data[i * m + k].im *= -1;
			}
			printf("% lfi", V_data[i * m + k].im);
			printf("\t");
		}
		printf("\n");
	}

	printf("\nOutput for Eigen Value...\n");
	for (i = 0; i < m; i++)
	{ //need fixing 16 as size
		for (k = 0; k < m; k++)
		{ //need fixing 16 as size
			if (i == k)
			{
				printf("%f ", D_data[i * m + k].re);
				if (!(D_data[i * m + k].im < 0))
					printf("+");
				else {
					printf("-");
					D_data[i * m + k].im *= -1;
				}
				printf(" %fi", D_data[i * m + k].im);
				printf("\t");
			}
		}
		printf("\n");
	}
	
}

