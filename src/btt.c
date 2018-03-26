#include <stdlib.h>
#include "iio.h"
#include "xmalloc.c"
int main(void)
{
	int w = 24076;
	int h = 50000;
	float *x = xmalloc(w * h * sizeof*x);
	for (int i = 0; i < 1000*1000; i++)
		x[i] = rand() / (RAND_MAX + 1.0);
	iio_write_image_float("big.tiff", x, w, h);
	free(x);
	return 0;
}
