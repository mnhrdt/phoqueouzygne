#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "s1a.h"
#include "xmalloc.c"
#include "xfopen.c"

static void asc_write(char *fname, float *x, int w, int h)
{
	FILE *f = xfopen(fname, "w");
	fprintf(f, "%d %d 1 1\n", w, h);
	for (int i = 0; i < w*h; i++)
		fprintf(f, "%g\n", x[i]);
	xfclose(f);
}

static void pgm_write(char *fname, uint8_t *x, int w, int h)
{
	FILE *f = xfopen(fname, "w");
	fprintf(f, "P5\n%d %d\n255\n", w, h);
	fwrite(x, h, w, f);
	xfclose(f);
}

static void pgm_write_wcrop(char *s, uint8_t *x,int w, int h, int x0, int xf)
{
	FILE *f = xfopen(s, "w");
	fprintf(f, "P5\n%d %d\n255\n", 1+xf-x0, h);
	for (int j = 0; j < h; j++)
	for (int i = 0; i <= xf-x0; i++)
		fputc(x[j*w+i+x0], f);
	xfclose(f);
}

#include "iio.h"
static void float_write_wcrop(char *s, float *x, int w, int h, int pd,
		int x0, int xf)
{
	if (x0 < 0) exit(fprintf(stderr,"ERROR: bad x0=%d\n",x0));
	if (xf >=w) exit(fprintf(stderr,"ERROR: bad xf=%d\n",xf));

	int w2 = 1 + xf - x0;
	int h2 = h;
	float *x2 = xmalloc(w2 * h2 * pd * sizeof*x2);
	for (int j = 0; j < h2; j++)
	for (int i = 0; i < w2; i++)
	for (int l = 0; l < pd; l++)
		x2[(j*w2 + i)*pd+l] = x[(j*w + i + x0)*pd+l];
	iio_write_image_float_vec(s, x2, w2, h2, 2);
	free(x2);
}



#include "pickopt.c"
int main(int c, char *v[])
{
	// process named arguments
	// x: raw image
	// y: range-focused image
	// z: azimuth-focused image
	char *filename_x = pick_option(&c, &v, "x", "");
	char *filename_y = pick_option(&c, &v, "y", "");
	char *filename_z = pick_option(&c, &v, "z", "");

	// process positional parameters
	if (c != 4 && c != 6) return fprintf(stderr, "usage:\n\t"
			"%s raw.dat lin0 linf [col0 colf]\n", *v);
	//                0 1       2    3     4    5
	char *filename_in  = v[1];
	int n_first = atoi(v[2]);
	int n_last  = atoi(v[3]);
	int col_first  = c > 5 ? atoi(v[4]) : 0 ;
	int col_last   = c > 5 ? atoi(v[5]) : 0 ;
	if (n_first > n_last)
		return fprintf(stderr, "first(%d) > last(%d)\n",n_first,n_last);

	struct s1a_file f[1];
	s1a_load_whole_datafile_trunc(f, filename_in, 1+n_last);

	int max_nq = 0;
	for (int i = n_first; i <= n_last; i++)
	{
		int nq = f->t[i].secondary_header.field.number_of_quads;
		if (nq > max_nq)
			max_nq = nq;
	}
	fprintf(stderr, "max_nq = %d\n", max_nq);
	int w = 2 * max_nq;
	int h = 1 + n_last - n_first;
	complex float *x = xmalloc(w * h * sizeof*x);
	complex float *y = xmalloc(w * h * sizeof*y);
	complex float *Y = xmalloc(w * h * sizeof*y);
	complex float *z = xmalloc(w * h * sizeof*z);
	for (int i = 0; i < w*h; i++) x[i] = y[i] = Y[i] = z[i] = 0;


	// debug focusing parameters
	{
	double param_TXPSF = s1a_extract_datum_TXPSF (f->t + n_first);
	double param_TXPRR = s1a_extract_datum_TXPRR (f->t + n_first);
	double param_TXPL  = s1a_extract_datum_TXPL  (f->t + n_first);
	int    param_TXPL3 = s1a_extract_datum_TXPL3 (f->t + n_first);
	int    param_NF    = s1a_extract_datum_NF    (f->t + n_first);
	double param_SWL   = s1a_extract_datum_SWL  (f->t + n_first);
	int    param_SWL3  = s1a_extract_datum_SWL3 (f->t + n_first);
	fprintf(stderr, "TXPRR = %.18e\n", param_TXPRR);
	fprintf(stderr, "TXPL = %.18e\n", param_TXPL);
	fprintf(stderr, "TXPSF = %.18e\n", param_TXPSF);
	fprintf(stderr, "TXPL3 = %d\n", param_TXPL3);
	fprintf(stderr, "SWL = %.18e\n", param_SWL);
	fprintf(stderr, "SWL3 = %d\n", param_SWL3);
	fprintf(stderr, "NF = %d\n", param_NF);
	}

	// decode and focus lines (RANGE FOCUSING)
	fprintf(stderr, "w = %d\n", w);
	for (int i = 0; i < h; i++)
	{
		s1a_decode_line_fancy(x + w*i,
				//x_block + w*i, x_brc   + w*i, x_thidx + w*i,
				NULL, NULL, NULL,
				f->t + n_first + i);
		s1a_focus_decoded_line(y + w*i, x + w*i,
				f->t + n_first + i);
	}

	// select columns of concern (AZIMUTH BLOCK)
	int wmin = 0;
	int wmax = w;
	if (col_first || col_last) {
		wmin = col_first;
		wmax = col_last+1;
	}

	// compute DC frequency polynomial, needed later
	// formula 5-18
	complex double ACCC[wmax];
	for (int i = 0; i < wmax; i++) ACCC[i] = 0;
	for (int i = wmin; i < wmax; i++)
	{
		complex double c = 0;
		for (int j = 0; j < h-1; j++)
			c += y[w*j+i] * conj(y[w*(j+1)+i]);
		fprintf(stderr, "c[i] = %g %g, angle = %g\n",
				creal(c), cimag(c), carg(c));
		ACCC[i] = c;
	}
	complex long double ACCC_avg = 0;
	for (int i = wmin; i < wmax; i++)
		ACCC_avg += ACCC[i];
	ACCC_avg /= (wmax - wmin);
	fprintf(stderr, "ACCC_avg = %g %g, angle = %g\n",
			creal(ACCC_avg), cimag(ACCC_avg), carg(ACCC_avg));
	double phi_accc = carg(ACCC_avg); // TODO: use a more robust average

	double TXPSF = s1a_extract_datum_TXPSF(f->t + n_first);
	double TXPRR = s1a_extract_datum_TXPRR(f->t + n_first);
	double TXPL  = s1a_extract_datum_TXPL (f->t + n_first);
	double phi_1 = TXPSF + TXPRR*TXPL / 2; // TODO: check!


	// azimuth FFT
	for (int i = wmin; i < wmax; i++)
	{
		complex float t[h], T[h];
		for (int j = 0; j < h; j++) t[j] = y[w*j + i];
		fft(T, t, h);
		for (int j = 0; j < h; j++) Y[w*j + i] = T[j];
	}


	// azimuth zero-padding


	// secondary range compression SRC (TODO)
	;

	// range cell migration correction RCMC (TODO)
	;

	// formula 5-19
	double f_eta_c = -s1a_extract_datum_PRF(f->t+n_first)*phi_accc/(2*M_PI);
	fprintf(stderr, "f_eta_c = %g\n", f_eta_c);

	// page 6-22 section 6.3.4.1.2
	//
	// XXX TODO FIXME CHECK wether SWST is indeed the fast time
	// (and that it is always the same inside this azimuth block)
	//
	double tau0 = s1a_extract_datum_SWST(f->t + n_first);
	double Fr = s1a_extract_datum_Fr(f->t + n_first);
	double V_s = 7500;
	double V_g = 0.88 * V_s;
	double V_r = sqrt(V_s*V_g); // TODO: obtain from metadata
	//for (int j = wmin; j < wmax; j++)
	for (int j = 0; j < h; j++)
	{
		double tau = tau0 + (j + n_first) / Fr; // 6-37 (contains ERROR)
		double R = tau * SPEED_OF_LIGHT / 2.0;

		// Compute the azimuth frequency vector AFV
		// page 6-23, first line
		int Mfft = h;//100; // should be h, or maybe smaller
		double AFV[Mfft];
		double Dvec[h];
		for (int i = 0; i < Mfft; i++)
		{
			double Fa = s1a_extract_datum_PRF(f->t + n_first);
			AFV[i] = f_eta_c - Fa/2 + i*Fa/Mfft;
		}

		complex float H[h]; // azimuth matched filter 6-35 page 6-22
		for (int i = 0; i < Mfft; i++)
		{
			// formula 6-29
			double f_eta = AFV[i];
			double f_0 = phi_1; // ERROR FIXME UNKNOWN XXX TODO
			double q = SPEED_OF_LIGHT * f_eta / (2 * V_r * f_0);
			if (q*q > 1) {
				fprintf(stderr, "f_eta = %g\n", f_eta);
				fprintf(stderr, "V_r = %g\n", V_r);
				fprintf(stderr, "f_0 = %g\n", f_0);
				fprintf(stderr, "q = %g\n", q);
				fflush(stderr);
			}
			assert(q*q <= 1);
			double D_f_theta_V_r = sqrt(1 - q*q);
			double DD = pow(D_f_theta_V_r, 2);
			Dvec[i] = DD;
			double W = 1;
			H[i] = W * cexp(I * 4*M_PI*R*DD*f_0 / SPEED_OF_LIGHT);
		}
		assert(h == Mfft);

		// normalize the match filter, formula 6-39 page 6-23
		long double N = 0;
		for (int i = 0; i < Mfft; i++)
			N += H[i] * conj(H[i]);
		double E = sqrt(N / Mfft);
		for (int i = 0; i < Mfft; i++)
			H[i] /= E;

		static int iHvals_saved = 0;
		if (!iHvals_saved)
		{
			iHvals_saved = 1;
			complex float iH[h];
			ifft(iH, H, h);
			FILE *ff = fopen("/tmp/iHvals.txt", "w");
			for (int i = 0; i < h; i++)
				fprintf(ff, "%g %g\n", creal(iH[i]),
						cimag(iH[i]));
			fclose(ff);
			ff = fopen("/tmp/DDvals.txt", "w");
			for (int i = 0; i < h; i++)
				fprintf(ff, "%g\n", Dvec[i]);
			fclose(ff);
		}

		// apply the filter, step 2 of page 6-21 section 6.3.4
		for (int i = 0; i < Mfft; i++)
			Y[w*j + i] *= conj(H[i]);
	}

	// azimuth iFFT
	for (int i = wmin; i < wmax; i++)
	{
		complex float t[h], T[h];
		for (int j = 0; j < h; j++) t[j] = Y[w*j + i];
		ifft(T, t, h);
		for (int j = 0; j < h; j++) z[w*j + i] = T[j];
	}


	//// focus columns (AZIMUTH FOCUSING)
	//for (int i = wmin; i < wmax; i++)
	//{
	//	fprintf(stderr, "focusing column %d\n", i);
	//	complex float t1[h], t2[h];
	//	for (int j = 0; j < h; j++)
	//		t1[j] = y[w*j + i];
	//	s1a_focus_column(t2, t1, h, 0, 0, 0, 50);
	//	for (int j = 0; j < h; j++)
	//		z[w*j + i] = t2[j];
	//}

	fprintf(stderr, "not going to free header mem\n");
	//s1a_file_free_memory(f);
	//fprintf(stderr, "i freed the mem!\n");

	int x0 = 0;
	int xf = w - 1;
	if (col_last)
	{
		x0 = col_first;
		xf = col_last;
	}
	if (*filename_x) float_write_wcrop(filename_x, (float*)x, w,h,2, x0,xf);
	if (*filename_y) float_write_wcrop(filename_y, (float*)y, w,h,2, x0,xf);
	if (*filename_z) float_write_wcrop(filename_z, (float*)z, w,h,2, x0,xf);


	free(x);
	free(y);
	free(z);

	return 0;
}
