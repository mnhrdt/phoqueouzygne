#include "gmtsar.h"
#include "lib_functions.h"
/*--------------------------------------------------------------*/
/* read i2 SLC complex data					*/					
/* SLCfile - SLC file						*/
/* name - name of SLC_file					*/
/* sdata - input i2 data					*/
/* cdata - output complex data (float) 				*/
/* xdim - length of row						*/
/*								*/
int read_SLC_short2float(FILE *SLCfile, char *name, short *sdata, fcomplex *cdata, int xdim, int psize, double dfact)
{
long     k;

	if (fread(sdata, 2*sizeof(short), psize*xdim, SLCfile) != psize*xdim) die("error reading SLC file", name);
	for (k=0; k<psize*xdim; k++) { 
		cdata[k].r = (float) dfact * sdata[2*k];
		cdata[k].i = (float) dfact * sdata[2*k +1];
	}
	return(EXIT_SUCCESS);
}
/*---------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* read i2 SLC complex data					*/					
/* SLCfile - SLC file						*/
/* name - name of SLC_file					*/
/* sdata - input i2 data					*/
/* cdata - output complex data (double)				*/
/* xdim - length of row						*/
/*								*/
int read_SLC_short2double(FILE *SLCfile, char *name, short *sdata, dcomplex *cdata, int xdim, int psize, double dfact)
{
long     k;

	if (fread(sdata, 2*sizeof(short), psize*xdim, SLCfile) != psize*xdim) die("error reading SLC file", name);
	for (k=0; k<psize*xdim; k++) {
		cdata[k].r = (double) dfact * sdata[2*k];
		cdata[k].i = (double) dfact * sdata[2*k +1];
	}
	return(EXIT_SUCCESS);
}
/*---------------------------------------------------------------*/
