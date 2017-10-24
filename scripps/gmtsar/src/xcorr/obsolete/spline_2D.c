#include<stdio.h>
#include<stdlib.h>
#include<math.h>
main()
{
/*
 * nn	size of array
 * x 	x coordinates
 * y	coordinate at which value is desired
 * u	sample values
 * s	derivatives	[output]
 * a	working space [size: nn]
 * eval	value of u at point y
 */

int	i, j, k, ifc;
int	istart, nn, nn2;
double	*x, *u, *s, *a, y, eval;
float	*data, *interp_data;

nn = 4;
ifc = 4;
nn2 = ifc * nn;

data = (float *) malloc (nn * nn * sizeof(double));
interp_data = (float *) malloc (nn2 * nn2 * sizeof(double));

x = (double *) malloc (nn * sizeof(double));
u = (double *) malloc (nn * sizeof(double));
s = (double *) malloc (nn * sizeof(double));
a = (double *) malloc (nn * sizeof(double));

for (i=0; i<nn2*nn2; i++) interp_data[i] = 0.0;
for (i=0; i<nn; i++){
	for (j=0; j<nn; j++) data[i*nn+j] = (float) (i + j);
	}

for (i=0; i<nn; i++){
	for (j=0; j<nn; j++) fprintf(stderr, " %f ", data[i*nn+j]);
	fprintf(stderr, "\n");
	}

fprintf(stderr, "\n");
/* interpolate along rows */
for (i=0; i<nn; i++){
	for (j=0; j<nn; j++){
		x[j] = (double) j;
		u[j] = (double) data[i*nn + j];
		interp_data[ifc*nn2*i + ifc*j] = data[i*nn + j];
		}

	/* calculate s */
	spline_(&istart, &nn, x, u, s, a);

	/* evaluate at point y along row */
	for (j=0; j<nn; j++) {
		for (k=0; k<ifc; k++) {
			y = x[j] + (double) (1.0 * k / ifc);
			evals_(&istart, &y, &nn, x, u, s, &eval);
			interp_data[ifc*nn2*i + ifc*j + k] = eval;
			}
		}
	}

/* interpolate along columns */
for (j=0; j<nn2; j++){
	for (i=0; i<nn2; i+=ifc){
		x[i/ifc] = (double) i;
		u[i/ifc] = (double) interp_data[i*nn2 + j];
		}

	spline_(&istart, &nn, x, u, s, a);

	for (i=0; i<nn; i++) {
		for (k=0; k<ifc; k++) {
			y = x[i] + (double) (1.0 * k);
			evals_(&istart, &y, &nn, x, u, s, &eval);
			interp_data[(ifc*i+k)*nn2 + j] = eval;
			}
		}
	}

/* print */
for (i=0; i<nn2; i++){
	for (j=0; j<nn2; j++) fprintf(stderr, " %4.2f ", interp_data[i*nn2+j]);
	fprintf(stderr, "\n");
	}
}
