#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmtsar.h"
#include "gmt.h"
/*-------------------------------------------------------------------------------------	*/
int print_cpatch(int nx, int ny, struct FCOMPLEX *m)
{
int i;
	for (i=0; i<nx*ny; i++){
		fprintf(stderr," (%4.2f %4.2f) ", m[i].r, m[i].i);
		if ((i+1)/nx  == (i+1)/(float) nx) fprintf(stderr,"\n");
		}
	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int print_patch(int nx, int ny, float *m)
{
int i;
	for (i=0; i<nx*ny; i++){
		fprintf(stderr," %4.2f", m[i]);
		if ((i+1)/nx  == (i+1)/(float) nx) fprintf(stderr,"\n");
		}
	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
