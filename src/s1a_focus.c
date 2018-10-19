#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>

#include "s1a.h"

static void fft(complex float *X, complex float *x, int n)
{
	fftwf_plan p = fftwf_plan_dft_1d(n, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
}

static void ifft(complex float *x, complex float *X, int n)
{
	fftwf_plan p = fftwf_plan_dft_1d(n, X, x, FFTW_BACKWARD,FFTW_ESTIMATE);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
}

#include "smapa.h"
SMART_PARAMETER(FHACK,37);
static void fill_chirp(complex float *c_out, int n,
		double TXPRR,
		double TXPSF,
		double TXPL,
		int NF
		)
{
	assert(NF < n);
	assert(NF >= 0);
	assert(isfinite(TXPRR));
	assert(isfinite(TXPSF));
	assert(isfinite(TXPL));
	complex float c[NF];

	// SENT-TN-52-7445 sec.4.2.1.1.1 (p.4-6)
	double phi1 = TXPSF + TXPRR*TXPL / 2; // TODO: check!
	double phi2 = TXPRR / 2;
	for (int i = 0; i < NF; i++)
	{
		double t = TXPL*(-0.5 + i/(float)NF);
		double s = phi1 * t + phi2 * t * t;
		//fprintf(stderr, "t,s,s/t[%d] = %g %g\n", i, t, phi1+t*phi2);

		c[i] = cexp(2*M_PI*I*s)/NF;
	}

	for (int i = 0; i < n; i++)
		c_out[i] = 0;
	for (int i = 0; i < NF; i++)
	{
		int j = i < NF/2 ? n-NF/2+i : i-NF/2;
		assert(j < n);
		assert(j >= 0);
		c_out[j] = c[i];
	}

	static int chirp_saved = 0;
	if (!chirp_saved)
	{
		chirp_saved = 1;
		FILE *ff = fopen("/tmp/chirpvals.txt", "w");
		for (int i = 0; i < n; i++)
			fprintf(ff, "%g %g\n", creal(c_out[i]),
					cimag(c_out[i]));
		fclose(ff);
	}
}

static void focus_one_line(
		complex float *y, // output (focused line)
		complex float *x, // input (raw data for a line)
		int n,            // number of complex samples
		double TXPRR,
		double TXPSF,
		double TXPL,
		int NF
		)
{
	complex float c[n];
	fill_chirp(c, n, TXPRR, TXPSF, TXPL, NF);

	complex float C[n], X[n], Y[n];
	fft(C, c, n);
	fft(X, x, n);

	for (int i = 0; i < n; i++)
		Y[i] = X[i] * C[i];

	ifft(y, Y, n);
}


int s1a_focus_decoded_line(complex float *out, complex float *in,
		struct s1a_isp *x)
{
	// length of the line
	int n = 2 * x->secondary_header.field.number_of_quads;

	// chirp parameters
	double param_TXPSF = s1a_extract_datum_TXPSF(x);
	double param_TXPRR = s1a_extract_datum_TXPRR(x);
	double param_TXPL  = s1a_extract_datum_TXPL(x);
	int    param_NF    = s1a_extract_datum_NF(x);

	// focus
	focus_one_line(out, in, n,
			param_TXPRR, param_TXPSF, param_TXPL, param_NF);

	return 1;
}

//int s1a_focus_column(complex float *out, complex float *in, int n,
//		float k, int l, float f)
//{
//	focus_one_line(out, in, n, k, l, f);
//}

