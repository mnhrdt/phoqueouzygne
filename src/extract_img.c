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
	// x: raw image
	// y: range-focused image
	// z: azimuth-focused image
	char *filename_x = pick_option(&c, &v, "x", "");
	char *filename_y = pick_option(&c, &v, "y", "");
	char *filename_z = pick_option(&c, &v, "z", "");
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
	complex float *z = xmalloc(w * h * sizeof*z);
	for (int i = 0; i < w*h; i++) x[i] = y[i] = z[i] = 0;
	//uint8_t *x_block = xmalloc(w*h);
	//uint8_t *x_brc   = xmalloc(w*h);
	//uint8_t *x_thidx = xmalloc(w*h);
	//for (int i = 0; i < w*h; i++) x_block[i]=x_brc[i]=x_thidx[i]=42;


	// read focusing parameters
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

	// azimuth zero-padding
	;

	// azimuth FFT
	;

	// focus columns (AZIMUTH FOCUSING)
	int wmin = 0;
	int wmax = w;
	if (col_first || col_last) {
		wmin = col_first;
		wmax = col_last+1;
	}
	for (int i = wmin; i < wmax; i++)
	{
		fprintf(stderr, "focusing column %d\n", i);
		complex float t1[h], t2[h];
		for (int j = 0; j < h; j++)
			t1[j] = y[w*j + i];
		s1a_focus_column(t2, t1, h, 0, 0, 0, 30);
		for (int j = 0; j < h; j++)
			z[w*j + i] = t2[j];
	}

	fprintf(stderr, "going to free header mem\n");
	s1a_file_free_memory(f);
	fprintf(stderr, "i freed the mem!\n");

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


//	if (!col_first && !col_last) {
//		iio_write_image_float_vec(filename_out, (float*)x, w, h, 2);
//		iio_write_image_float_vec(filename_ynorm, (float*)y, w, h, 2);
//	}
//	else {
//		int w2 = 1 + col_last - col_first;
//		int h2 = h;
//		if (w2 < 0 || w2 > 50000) return fprintf(stderr, "ERROR:bad w2=%d\n", w2);
//		if (h2 < 0 || h2 > 50000) return fprintf(stderr, "ERROR:bad h2=%d\n", h2);
//		fprintf(stderr, "output file of size %dx%d\n", w2, h2);
//		complex float *x2 = xmalloc(w2 * h2 * sizeof*x2);
//		complex float *y2 = xmalloc(w2 * h2 * sizeof*y2);
//		fprintf(stderr, "a\n");
//		for (int j = 0; j < h2; j++)
//		for (int i = 0; i < w2; i++) {
//			x2[j*w2 + i] = x[j*w + i + col_first];
//			y2[j*w2 + i] = y[j*w + i + col_first];
//		}
//		fprintf(stderr, "a\n");
//		iio_write_image_float_vec(filename_out, (float*)x2, w2, h2, 2);
//		iio_write_image_float_vec(filename_y, (float*)y2, w2, h2, 2);
//		iio_write_image_float_vec(filename_z, (float*)x2, w2, h2, 2);
//		fprintf(stderr, "a\n");
//		free(y2);
//		free(x2);
//		fprintf(stderr, "a\n");
//	}

	//for (int i = 0; i < w*h; i++)
	//{
	//	x_norm[i] = cabs(x[i]);
	//	x_real[i] = creal(x[i]);
	//	x_imag[i] = cimag(x[i]);
	//	y_norm[i] = cabs(y[i]);
	//}
	//free(x);
	//fprintf(stderr, "i freed more mem!\n");
	////if (isfinite(FHACK2()))
	////	iio_write_image_float(filename_ynorm, y_norm, w, h);

	//int x0 = 0;
	//int xf = w - 1;
	//if (col_last)
	//{
	//	x0 = col_first;
	//	xf = col_last;
	//}
	//pgm_write_wcrop("x_block.pgm", x_block, w, h, x0, xf);
	//pgm_write_wcrop("x_brc.pgm"  , x_brc  , w, h, x0, xf);
	//pgm_write_wcrop("x_thidx.pgm", x_thidx, w, h, x0, xf);
	//fprintf(stderr, "going to write stuff\n");
	//iio_write_image_float("x_norm.tif", x_norm, w, h);
	//iio_write_image_float("x_real.tif", x_real, w, h);
	//iio_write_image_float("x_imag.tif", x_imag, w, h);
	//iio_write_image_float("y_norm.tif", y_norm, w, h);
	////asc_write("x_real.asc", x_real, w, h);
	////asc_write("x_imag.asc", x_imag, w, h);


	free(x);
	free(y);
	free(z);
	//free(x_block);
	//free(x_brc);
	//free(x_thidx);

	return 0;
}
