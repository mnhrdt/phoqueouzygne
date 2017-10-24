/*-------------------------------------------------------*/
#include <math.h>
#include "gmtsar.h"
#include "siocomplex.h"
#include "xcorr.h"
/* Cain, S. C. and M. M. Hayat, Projection-Based Image Registration in the
 * Presence of Fixed-Pattern Noise, 
 * IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 10, NO. 12, DECEMBER 2001
 * p. 1860
*/
int xcorr_1d(struct FCOMPLEX *, struct FCOMPLEX *, float *, int, int, int *);
/*-------------------------------------------------------------------------------*/
void do_freq_proj_corr(struct xcorr xc, int iloc)
{
	int	nxc, nyc, maxlag,ifac;
	int	xoff, yoff;
	int 	i, j, ii, ix, iy;
	float	row[256], col[256], fix, fiy;
	float 	max_corr_x, max_corr_y, max_corr;
	struct FCOMPLEX row1[256], col1[256];
	struct FCOMPLEX row2[256], col2[256];
	struct FCOMPLEX ccin[256], ccout[4*256];

	for (i=0; i<256; i++) {
		row1[i].r = col1[i].r = row2[i].r = col2[i].r = 0.0;
		row1[i].i = col1[i].i = row2[i].i = col2[i].i = 0.0;
		}

	/* d1, c1 is the master					*/
	/* d2, c2 is the slave					*/
	if (debug) print_complex(xc.c1, xc.npy, xc.npx, 1);
	if (debug) print_complex(xc.c2, xc.npy, xc.npx, 1);

	/* sum over rows and columns */
	for (i=0; i<xc.npy; i++){
		for (j=0; j<xc.npx; j++){
			row1[j].r += xc.c1[i*xc.npx+j].r;
			row2[j].r += xc.c2[i*xc.npx+j].r;
			col1[j].r += xc.c1[i+j*xc.npx].r;
			col2[j].r += xc.c2[i+j*xc.npx].r;
			}
		}

	maxlag = xc.npx / 2;
	maxlag = 32;
	/* do 1D xcorr along rows and columns */
	xcorr_1d(row1, row2, row, xc.npx, maxlag, &xoff);
	xcorr_1d(col1, col2, col, xc.npy, maxlag, &yoff);


	max_corr = calc_time_corr(xc, -1 * xoff, -1 * yoff);
	if (verbose) fprintf(stderr," (freq proj) jpeak %f xoffset %d corr %4.2lf %4.2f\n", (float) xoff, xc.x_offset, row[xoff+xc.npx/2], max_corr);
	if (verbose) fprintf(stderr," (freq proj) ipeak %f yoffset %d corr %4.2lf %4.2f\n", (float) yoff, xc.y_offset, col[yoff+xc.npy/2], max_corr);

	xc.loc[iloc].xoff = -1 * xoff;
	xc.loc[iloc].yoff = -1 * yoff;
	xc.loc[iloc].corr = max_corr;

/*
	for (i=-8;i<8;i++) fprintf(stderr," %d %4.3f %4.3f \n", i,row[i+xoff+xc.npx/2], col[i+yoff+xc.npy/2]);
*/
	ifac = 4;
	max_corr_x = max_corr_y = -99999.99;

	for (i=0; i<xc.npx; i++) {
		ccin[i].r = row[i]; 
		ccin[i].i = 0.0; 
		}

	fft_interpolate_1d(ccin, xc.npx, ccout, ifac);

	for (i=0; i<ifac*xc.npx; i++) {
		if (ccout[i].r > max_corr_x) {
			max_corr_x = ccout[i].r;
			ix = (ifac*xc.npx)/2 - i;
			}
		}

	for (i=0; i<xc.npy; i++) {
		ccin[i].r = col[i]; 
		ccin[i].i = 0.0; 
		}

	fft_interpolate_1d(ccin, xc.npy, ccout, ifac);

/*
	for (i=-8;i<8;i++) fprintf(stderr," %d %4.3f  \n", i,ccout[i-ifac*xoff+ifac*xc.npy/2].r);
*/
	for (i=0; i<ifac*xc.npy; i++) {
		if (ccout[i].r > max_corr_y) {
			max_corr_y = ccout[i].r;
			iy = (ifac*xc.npy)/2 - i;
			}
		}

	fix = ix / (float) ifac;
	fiy = iy / (float) ifac;

}
/*-------------------------------------------------------------------------------*/
int xcorr_1d(struct FCOMPLEX  *x, struct FCOMPLEX *y, float *c, int n, int maxlag, int *xoff)
{
int	i,j, k, lag;
float	denom, nom, n1, n2, max;
denom = n1 = n2 = 0.0;
max = -9999.0;

	/* calculate denom */
	for (i=0; i<n; i++){
		n1 += x[i].r*x[i].r;
		n2 += y[i].r*y[i].r;
		}
	denom = sqrtf(n1*n2);

	for (lag=-1*maxlag; lag<=maxlag; lag++){
		nom = 0.0;
		for (i=0; i<n; i++){
			k = i + lag;
			if (k < 0) k = k + n;
			if (k > n-1) k = k - n;
			nom += x[i].r*y[k].r;
			}	

		c[lag+n/2] = nom / denom;
		if (c[lag+n/2] > max){
			max = c[lag+n/2];
			*xoff = lag;
			}
		}

/*
	for (i=-1*maxlag; i<maxlag; i++) {
		fprintf(stdout," row %d %f ", i, c[i+n/2]);
		if (i == *xoff) fprintf(stdout,"**\n");
		if (i != *xoff) fprintf(stdout,"\n");
		}
*/

}
/*-------------------------------------------------------------------------------*/
