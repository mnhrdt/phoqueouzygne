#include <complex.h>
#include <math.h>
#include <fftw3.h>

static void fft(complex float *X, complex float *x, int n)
{
	fftwf_plan p = fftwf_plan_dft_1d(n, x, X, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
}

static void ifft(complex float *x, complex float *X, int n)
{
	fftwf_plan p = fftwf_plan_dft_1d(n, X, x, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(p);
	fftwf_destroy_plan(p);
}

static void fill_chirp(complex float *c, int n, float k)
{
	for (int i = 0; i < n; i++)
		c[i] = cos(k*i*i);
}

static void focus_one_line(
		complex float *y, // output (focused line)
		complex float *x, // input (raw data for a line)
		int n,            // number of complex samples
		float k           // pulse parameter K
		)
{
	complex float c[n];
	fill_chirp(c, n, k);

	complex float C[n], X[n], Y[n];
	fft(C, c, n);
	fft(X, x, n);

	for (int i = 0; i < n; i++)
		Y[i] = X[i] * C[i];

	ifft(y, Y, n);
}

#include "s1a.h"

int s1a_focus_decoded_line(complex float *out, complex float *in,
		struct s1a_isp *x)
{
	// length of the line
	int n = 2 * x->secondary_header.field.number_of_quads;

	// chirp parameters
	float k = 17;

	// focus
	focus_one_line(out, in, n, k);

	return 1;
}
