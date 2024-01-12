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
	// process positional parameters
	if (c != 5) return fprintf(stderr, "usage:\n\t"
			"%s lin0 linf in.dat out.npy]\n", *v);
	//                0 1    2    3      4
	int n_first = atoi(v[1]);
	int n_last  = atoi(v[2]);
	char *filename_in  = v[3];
	char *filename_out  = v[4];
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
	//fprintf(stderr, "max_nq = %d\n", max_nq);
	int w = 2 * max_nq;
	int h = 1 + n_last - n_first;
	complex float *x = xmalloc(w * h * sizeof*x);
	complex float *y = xmalloc(w * h * sizeof*y);
	complex float *Y = xmalloc(w * h * sizeof*y);
	complex float *z = xmalloc(w * h * sizeof*z);
	for (int i = 0; i < w*h; i++) x[i] = y[i] = Y[i] = z[i] = 0;


	// debug focusing parameters
	if (0) {
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

	// decode lines (RANGE FOCUSING)
	//fprintf(stderr, "w = %d\n", w);
	for (int i = 0; i < h; i++)
		s1a_decode_line_fancy(x + w*i,
				//x_block + w*i, x_brc   + w*i, x_thidx + w*i,
				NULL, NULL, NULL,
				f->t + n_first + i);

	//fprintf(stderr, "not going to free header mem\n");
	//s1a_file_free_memory(f);
	//fprintf(stderr, "i freed the mem!\n");

	if (*filename_out)
		iio_write_image_float_vec(filename_out, (float*)x, w, h, 2);

	//free(x);

	return 0;
}
