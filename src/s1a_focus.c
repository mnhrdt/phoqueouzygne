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
static void fill_chirp(complex float *c, int n, float k, int l, float f)
{
	for (int i = 0; i < n; i++)
	{
		int j = i < n/2 ? i : n - i;
		double t = j / f;//FHACK();//FILTER_REF_FREQ;
		if (abs(j) < l/2)
			c[i] = c[n-i] = cos(k*t*t);
		else
			c[i] = 0;
	}

	//FILE *ff = fopen("/tmp/chirpvals.txt", "w");
	//for (int i = 0; i < n; i++)
	//	fprintf(ff, "%g\n", creal(c[i]));
	//fclose(ff);
}

static void focus_one_line(
		complex float *y, // output (focused line)
		complex float *x, // input (raw data for a line)
		int n,            // number of complex samples
		float k,          // pulse parameter K
		float l,          // filter length TXPL3
		float f           // frequency scaling
		)
{
	complex float c[n];
	fill_chirp(c, n, k, l, f);

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
	float k = s1a_extract_datum_TXPRR(x);
	int   l = s1a_extract_datum_TXPL3(x);
	float f = FHACK();

	static int printed = 0;
	if (!printed) {
		printed = 1;
		fprintf(stderr, "k = %g\n", k);
		fprintf(stderr, "l = %d\n", l);
	}

	// focus
	focus_one_line(out, in, n, k, l, f);

	return 1;
}

int s1a_focus_column(complex float *out, complex float *in, int n,
		float k, int l, float f)
{
	focus_one_line(out, in, n, k, l, f);
}

