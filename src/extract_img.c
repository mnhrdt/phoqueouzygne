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


#include "smapa.h"
SMART_PARAMETER(WMIN,0)
SMART_PARAMETER(WMAX,-1)
SMART_PARAMETER(FHACK2,INFINITY)


#include "iio.h"
#include "pickopt.c"
int main(int c, char *v[])
{
	char *filename_ynorm = pick_option(&c, &v, "y", "y_norm.asc");
	if (c != 4) return fprintf(stderr, "usage:\n\t%s raw.dat l0 lf\n", *v);
	//                                             0 1       2  3
	char *filename_x = v[1];
	int n_first = atoi(v[2]);
	int n_last  = atoi(v[3]);

	if (n_first > n_last)
		return fprintf(stderr, "first(%d) > last(%d)\n",n_first,n_last);

	struct s1a_file f[1];

	s1a_load_whole_datafile_trunc(f, filename_x, 1+n_last);

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
	for (int i = 0; i < w*h; i++) x[i] = 0;
	float *x_norm    = xmalloc(w * h * sizeof*x_norm);
	float *x_real    = xmalloc(w * h * sizeof*x_real);
	float *x_imag    = xmalloc(w * h * sizeof*x_imag);
	float *y_norm    = xmalloc(w * h * sizeof*y_norm);
	uint8_t *x_block = xmalloc(w*h);
	uint8_t *x_brc   = xmalloc(w*h);
	uint8_t *x_thidx = xmalloc(w*h);
	for (int i = 0; i < w*h; i++) x_block[i]=x_brc[i]=x_thidx[i]=42;


	// read focusing parameters
	{
	double param_TXPSF = s1a_extract_datum_TXPSF(f->t + n_first);
	double param_TXPRR = s1a_extract_datum_TXPRR(f->t + n_first);
	double param_TXPL = s1a_extract_datum_TXPL(f->t + n_first);
	int    param_NF = s1a_extract_datum_NF(f->t + n_first);
	fprintf(stderr, "TXPRR = %.18e\n", param_TXPRR);
	fprintf(stderr, "TXPL = %.18e\n", param_TXPL);
	fprintf(stderr, "TXPSF = %.18e\n", param_TXPSF);
	fprintf(stderr, "NF = %d\n", param_NF);
	}

	// decode and focus lines
	fprintf(stderr, "w = %d\n", w);
	fprintf(stderr, "x = %p\n", (void*)x);
	for (int i = 0; i < h; i++)
	{
		s1a_decode_line_fancy(x + w*i,
				x_block + w*i,
				x_brc   + w*i,
				x_thidx + w*i,
				f->t + n_first + i);
		s1a_focus_decoded_line(y + w*i, x + w*i,
				f->t + n_first + i);
	}

	//// focus columns
	//int wmin = WMIN();
	//int wmax = WMAX();
	//for (int i = wmin; i < wmax; i++)
	//{
	//	fprintf(stderr, "focusing column %d\n", i);
	//	complex float t1[h], t2[h];
	//	for (int j = 0; j < h; j++)
	//		t1[j] = y[w*j + i];
	//	s1a_focus_column(t2, t1, h,
	//			param_TXPRR, param_TXPSF, param_TXPL, param_NF);
	//	for (int j = 0; j < h; j++)
	//		y[w*j + i] = t2[j];
	//}

	fprintf(stderr, "going to free mem\n");
	s1a_file_free_memory(f);
	fprintf(stderr, "i freed the mem!\n");

	for (int i = 0; i < w*h; i++)
	{
		x_norm[i] = cabs(x[i]);
		x_real[i] = creal(x[i]);
		x_imag[i] = cimag(x[i]);
		y_norm[i] = cabs(y[i]);
	}
	free(x);
	fprintf(stderr, "i freed more mem!\n");
	//if (isfinite(FHACK2()))
	//	iio_write_image_float(filename_ynorm, y_norm, w, h);

	//pgm_write("x_block.asc", x_block, w, h);
	//pgm_write("x_brc.asc"  , x_brc  , w, h);
	//pgm_write("x_thidx.asc", x_thidx, w, h);
	fprintf(stderr, "going to write stuff\n");
	iio_write_image_float("x_norm.tif", x_norm, w, h);
	iio_write_image_float("x_real.tif", x_real, w, h);
	iio_write_image_float("x_imag.tif", x_imag, w, h);
	iio_write_image_float("y_norm.tif", y_norm, w, h);
	//asc_write("x_real.asc", x_real, w, h);
	//asc_write("x_imag.asc", x_imag, w, h);


	free(x_norm);
	free(x_real);
	free(x_imag);
	free(x_block);
	free(x_brc);
	free(x_thidx);

	return 0;
}
