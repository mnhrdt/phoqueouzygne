#include<stdio.h>
#include<math.h>
#include "gmtsar.h"
/*------------------------------------------------------------------------*/
/*	calculates 2D fft by doing 1D over rows, transposing, and then 	  */
/*	columns								  */
/*------------------------------------------------------------------------*/
void transpose_complex_NM(struct FCOMPLEX *, int, int);
int cfft2d(int *N, int *M, struct FCOMPLEX *cin, int *dir)
{
void cfft1d_(int *, struct FCOMPLEX *, int *);
int	i;
static	int flag = 1;

	if (flag == 0){
		fprintf(stderr,"using fftpack \n");
		flag = 1;
		}

	if (debug) print_complex(cin, *N, *M, 1);

	/* forward 2D */
	if (*dir == -1) {
		/* forward rows */
		for (i=0; i<(*N); i++) cfft1d_(M, &cin[i*(*M)], dir); 

		/* transpose */
		transpose_complex_NM(cin, *N, *M);

		/* forward columns */
		for (i=0; i<(*M); i++) cfft1d_(N, &cin[i*(*N)], dir); 

		/* transpose */
		transpose_complex_NM(cin, *M, *N);
		}

	/* inverse 2D */
	if (*dir == 1) {
		/* forward rows */
		for (i=0; i<(*N); i++) cfft1d_(M, &cin[i*(*M)], dir); 

		/* transpose */
		transpose_complex_NM(cin, *N, *M);

		/* forward columns */
		for (i=0; i<(*M); i++) cfft1d_(N, &cin[i*(*N)], dir); 

		/* transpose */
		transpose_complex_NM(cin, *M, *N);
		}

	if (debug) print_complex(cin, *N, *M, 1);

return 0;

}
/*------------------------------------------------------------------------*/
void transpose_complex_NM(struct FCOMPLEX *in, int n, int m)
{
int	i, j;
struct FCOMPLEX	*tmp;

	tmp = (struct FCOMPLEX *) malloc(n * m * sizeof(struct FCOMPLEX));

	for (i=0; i<n; i++) {
		for (j=0; j<m; j++) {
			tmp[j*n+i] = in[i*m+j];
			}
		}

	for (i=0; i<(n*m); i++) in[i] = tmp[i];

	free((char *) tmp);
}
/*------------------------------------------------------------------------*/
