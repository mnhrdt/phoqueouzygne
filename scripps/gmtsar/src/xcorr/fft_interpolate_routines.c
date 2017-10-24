#include "gmtsar.h"
/*--------------------------------------------------------------------------------------*/
/* some fft-based interpolation routines - rjm July 2009				*/
/*											*/	
/* fft_interpolate_1d(vin, N, vout, factor)	uses fft to interpolate a signal	*/
/*				vin - complex vector length N				*/
/*				vout - complex vector length M = factor*N		*/
/*				factor - interpolation factor				*/
/*				N must be even						*/
/*				adjust amplitude					*/
/*											*/
/* fft_interpolate_2d(min, N, mout, factor)	2D but uses 1D fft		*/
/*				min - complex matrix length NxN				*/
/*				mout - complex matrix length MxM = factor*N		*/
/*											*/
/* fft_arrange_interpolate(in, N, out, factor)	rearrange terms of transformed signal	*/
/*											*/
/*	all vectors and matrices must be pre-allocated					*/
/*	all use the SIO fcomplex and call cfft1d or cfft2d				*/
/* 	fcomplex - struct {float r; float i} where r is real and i imag			*/
/*--------------------------------------------------------------------------------------*/

void cfft1d_(int *,struct FCOMPLEX *, int *);
void cfft1d(int *,struct FCOMPLEX *, int *);
void fft_interpolate_1d(struct FCOMPLEX *, int, struct FCOMPLEX *, int);
void fft_arrange_interpolate(struct FCOMPLEX *, int, struct FCOMPLEX *, int, int);
void fft_interpolate_2d(struct FCOMPLEX *, int, int, struct FCOMPLEX *, int, int, int);

/*------------------------------------------------------------------------*/
void fft_interpolate_1d(struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int ifactor)
{
int	i, dir, M;

	M = ifactor * N;

	/* forward fft */
	dir = -1;
	cfft1d_(&N, in, &dir);

	/* re-arrange values in fft */
	fft_arrange_interpolate(in, N, out, M, ifactor);

	/* backward fft */
	dir = 1;
	cfft1d_(&M, out, &dir);

	/* scale amplitude */
	for (i=0; i<M; i++) {
		out[i].r = ((float) ifactor ) * out[i].r;
		out[i].i = ((float) ifactor ) * out[i].i;
		}
}
/*------------------------------------------------------------------------*/
void fft_arrange_interpolate(struct FCOMPLEX *in, int N, struct FCOMPLEX *out, int M, int ifactor)
{
int	i, N2;

	N2 = N / 2;

	for (i=0; i<M; i++) out[i].r = out[i].i = 0.0;

	for (i=0; i<N2; i++) {
		out[i].r = in[i].r;
		out[i].i = in[i].i;
		out[i + (2*ifactor-1)*N2].r = in[i + N2].r;
		out[i + (2*ifactor-1)*N2].i = in[i + N2].i;
		}
}
/*--------------------------------------------------------------------------------------*/
void fft_interpolate_2d(struct FCOMPLEX *in, int N1, int M1, struct FCOMPLEX *out, int N,  int M, int ifactor)
{
int	i, j, error_flag;
struct FCOMPLEX *tmp1,*tmp2,*tmp3;

	/* sanity checks */
	if (N != (N1 * ifactor)) error_flag = 1;
	if (M != (M1 * ifactor)) error_flag = 1;

	tmp1 = (struct FCOMPLEX *) malloc(N1 * M * sizeof(struct FCOMPLEX));
	tmp2 = (struct FCOMPLEX *) malloc(N1 * sizeof(struct FCOMPLEX));
	tmp3 = (struct FCOMPLEX *) malloc(N * sizeof(struct FCOMPLEX));

	for (i=0; i<N1*M;i++) tmp1[i].i = tmp1[i].r = 0.0;

	if (debug) print_complex(in, N1, M1, 0);

	for (i=0; i<N1; i++) fft_interpolate_1d(&in[i*M1], M1, &tmp1[i*M], ifactor);

	if (debug) print_complex(tmp1, N1, M, 0);

	/* now do columns - need to transpose */
	for (i=0; i<M; i++) {
		for (j=0; j<N1; j++) tmp2[j] = tmp1[j*M+i];
		fft_interpolate_1d(tmp2, N1, tmp3, ifactor);
		for (j=0; j<N; j++) out[j*M+i] = tmp3[j];
		}

	if (debug) print_complex(out, N, M, 0);

	free((char *) tmp1);
	free((char *) tmp2);
	free((char *) tmp3);
}
/*--------------------------------------------------------------------------------------*/
