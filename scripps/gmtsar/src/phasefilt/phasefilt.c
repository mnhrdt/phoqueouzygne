#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmtsar.h"
#include "gmt.h"
#include "netcdf.h"
/*-------------------------------------------------------------------------------------	*/
/*
 *	implementation of adaptive non-linear phase filter based on:				
 * 	Goldstein, R. M. and C. L. Werner, 1998, Radar interferogram filtering 
 *	for geophysical applications, (25),21, 4035-4038.
 *
 *	Baran, I., M. P. Stewart, B. M. Kampes, Z. Perski, and P. Lilly, 2003
 *	A modification of the Golstein radar interferogram filter,
 *	IEE trans. geoscience and remote sensing, (41), 9, 2114-2118.
 *
 *	rjm, SDSU, may 2010
 *	comments: 
 *	- handled edges in a kludgey way
 *	- added memory allocation checks
 *
 *	future:
 *	- added option to calculate power spectra with a taper
 *	- (synthetics suggested a slight improvement from reduced sidelobes)
 *	
 *      Dec 2010 - bug in weighting fixed [apply to complex]
 *     	default psize changed to 32 
 *      '-complex_out' option added to write out filtered real and imag
 *
 *-------------------------------------------------------------------------------------	*/
int verbose;
int debug;

int calc_corr(char *, struct GRD_HEADER *, char *, struct GRD_HEADER *, int, int, float *, float *);
int apply_pspec(int, int, float, struct FCOMPLEX *, struct FCOMPLEX *);

char *USAGE  = "\n USAGE:\nphasefilt -imag imag.grd -real real.grd [-alpha alpha][-psize size][-amp1 amp1.grd -amp2 amp2.grd][-diff][-v]\n"
" applies Goldstein adaptive filter to phase [output: filtphase.grd]\n"
" or applies modified Goldstein adaptive filter to phase [output: filtphase.grd, corrfilt.grd]\n"
"-imag [required] GMT format file of imaginary component\n"
"-real [required] GMT format file of real component\n"
"-alpha 	exponent for filter - usually between 0.0 and 1.5 (0.0 should not filter).\n"
"		default: 0.5	[Goldstein filter] (anything above 1.0 may be excessive)\n"
"		alpha < 0 will set alpha = (1 - coherence) [modified Goldstein]\n"
"-psize 	patch size for filtering. Must be power of two.\n"
"		default: 32\n"
"-amp1 		GMT format file of amplitude image of image 1. Needed (and applies) modified filter.\n"
"-amp2 		GMT format file of amplitude image of image 2. Needed (and applies) modified filter.\n"
"-diff 		Calculate difference between input phase and output phase.\n"
"-complex_out	Write out filtered real and imaginary (filtphase_real.grd and filtphase_imag.grd)\n"
"-v 		Verbose.\n"
"-debug 	Debug.\n"
"\n"
"example:\n"
"phasefilt -imag imag.grd -real real.grd -alpha 0.5\n";

int main(int argc, char **argv){
int	i, j, ii, jj, k1, k2;
int	xdim, ydim;
int	nxp, nyp; 
int	fdir, bdir, psize, corrflag, dflag, comflag;
float	alpha, pcorr, swgt;
float	*outphase, *wgt, *amp, *corr, *diff, *ftmp;
double 	x_inc, y_inc;
struct 	FCOMPLEX *data, *fdata, *patch0, *patch1;
char	sre[256], sim[256];
char	amp1[256], amp2[256];

/* GMT stuff */
int 	argc2 = 1;
char 	*argv2[2] = {"dummy",0};
struct 	GRD_HEADER grd_imag, grd_real;
struct 	GRD_HEADER grd_amp1, grd_amp2;

	argc2 = GMT_begin (argc2, argv2);

	debug = 0;
	verbose = 0;
	corrflag = 0;
	comflag = 0;
	if (argc < 3) die(USAGE,"");

	/* defaults */
	psize = 32;		/* size of patch  # changed from 64 */
	alpha = 0.5;		/* exponent */
	dflag = 0;		/* write out difference */
	parse_command_line(argv, argc, USAGE, sre, sim, &alpha, &psize, amp1, amp2, &dflag, &comflag);

	/* fft direction */
	fdir = -1;
	bdir = 1;

	/* patch size 		*/
	/* currently square 	*/
	nxp = psize;
	nyp = psize;

	/* read dimensions 	*/
	read_file_hdr(sim, &grd_imag, sre, &grd_real, &xdim, &ydim);

	if ((data = (struct FCOMPLEX *) malloc(xdim * ydim * sizeof(struct FCOMPLEX))) == NULL) die("error allocating memory","");
	if ((fdata = (struct FCOMPLEX *) malloc(xdim * ydim * sizeof(struct FCOMPLEX))) == NULL) die("error allocating memory","");
	if ((patch0 = (struct FCOMPLEX *) malloc(nxp * nyp * sizeof(struct FCOMPLEX))) == NULL) die("error allocating memory","");
	if ((patch1 = (struct FCOMPLEX *) malloc(nxp * nyp * sizeof(struct FCOMPLEX))) == NULL) die("error allocating memory","");

	if ((wgt = (float *) malloc(nxp * nyp * sizeof(float))) == NULL) die("error allocating memory","");
	if ((amp = (float *) malloc(xdim * ydim * sizeof(float))) == NULL) die("error allocating memory","");
	if ((outphase = (float *) malloc(xdim * (ydim+nyp) * sizeof(float))) == NULL) die("error allocating memory","");

	/* array for difference */
	if (dflag) {if ((diff = (float *) malloc(xdim * ydim * sizeof(float))) == NULL) die("error allocating memory","");}

	/* array for writing filtered real and imag  */
	if (comflag) {if ((ftmp = (float *) malloc(xdim * ydim * sizeof(float))) == NULL) die("error allocating memory","");}

	for (i=0; i<xdim*ydim; i++) outphase[i] = fdata[i].r = fdata[i].i = 0.0;

	/* read data 	*/
	read_data(sim, &grd_imag, sre, &grd_real, data, amp);

	/* calculate correlation 	*/
	if (alpha < 0) { 
		corr = (float *) malloc(xdim * ydim * sizeof(float));
		calc_corr(amp1, &grd_amp1, amp2, &grd_amp2, xdim, ydim, amp, corr); 
		corrflag = 1;
		}

	if (verbose) {
		if (corrflag == 1) fprintf(stderr,"phasefile: using coherence-dependent alpha (1 - gamma)\n");
		if (corrflag == 0) fprintf(stderr,"phasefilt: constant alpha (%6.2f)\n", alpha);
		}

	/* create weights for each patch 				*/
	/* each patch overlaps each other by half and summed 		*/
	/* ideally, total wgt for each pixl = 1				*/
	/* except at edges...						*/
	make_wgt(wgt, nxp, nyp);

	for (ii=0; ii<=(ydim-nyp); ii+=nyp/2){
		for (jj=0; jj<=(xdim-nxp); jj+=nxp/2){

			pcorr = 0.0;
			swgt = 0.0;
			for (i=0; i<nyp; i++){
				for (j=0; j<nxp; j++) {
					k1 = (ii+i)*xdim + jj + j;
					k2 = i*nxp+j;
					patch0[k2].r = data[k1].r;
					patch0[k2].i = data[k1].i;
					if (corrflag) {
						pcorr += wgt[k2]*corr[k1];
						swgt  += wgt[k2];
						}
					}
				}

			/* set alpha to 1.0 - coherence 	*/
			/* Baran et al., 2003			*/
			if (corrflag) alpha = 1.0 - pcorr/swgt;

			cfft2d(&nxp, &nyp, patch0, &fdir);

			apply_pspec(nxp, nyp, alpha, patch0, patch1);

			cfft2d(&nxp, &nyp, patch1, &bdir);

			for (i=0; i<nyp; i++){
				for (j=0; j<nxp; j++) {
					k1 = (ii+i)*xdim +jj+j;
					k2 = i*nxp+j;
					fdata[k1].r = fdata[k1].r + wgt[k2]*patch1[k2].r;
					fdata[k1].i = fdata[k1].i + wgt[k2]*patch1[k2].i;
					}
				}
			}
		}

	for (i=0; i<xdim*ydim; i++) outphase[i] = atan2f(fdata[i].i,fdata[i].r);
	
	write_grdfile("filtphase.grd", "phasefilt", "phase", outphase, x_inc, y_inc, xdim, ydim);

	if (corrflag) write_grdfile("filtcorr.grd", "phasefilt", "corr", corr, x_inc, y_inc, xdim, ydim);

	if (comflag) {
		for (i=0; i<xdim*ydim; i++)  ftmp[i] = fdata[i].r; 
		write_grdfile("filtphase_real.grd", "phasefilt", "real", ftmp, x_inc, y_inc, xdim, ydim);
		for (i=0; i<xdim*ydim; i++)  ftmp[i] = fdata[i].i; 
		write_grdfile("filtphase_imag.grd", "phasefilt", "imag", ftmp, x_inc, y_inc, xdim, ydim);
		}

	if (dflag) {
		for (i=0; i<xdim*ydim; i++) diff[i] = atan2f(data[i].i,data[i].r) - outphase[i];
		write_grdfile("filtdiff.grd", "phasefilt", "diff", diff, x_inc, y_inc, xdim, ydim, verbose);
		}

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int write_grdfile(char *fname, char *prog, char *type, float *data, double x_inc, double y_inc, int xdim, int ydim, int verbose)
{
FILE *fout;

	if (verbose) fprintf(stderr," writing %s \n", fname);

	x_inc = y_inc = 1.0;

	if ((fout = fopen(fname,"w")) == NULL) die("cannot open ",fname);
	create_GMT_binary_hdr(xdim, ydim, x_inc, y_inc, prog, type, fout);
	fwrite(data, sizeof(float), xdim*ydim, fout);
	fclose(fout);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int read_file_hdr(char *sim, struct GRD_HEADER *grd_imag, char *sre, struct GRD_HEADER *grd_real, int *xdim, int *ydim)
{
	if (GMT_read_grd_info(sim, grd_imag)) die("error reading file",sim);
	if (GMT_read_grd_info(sre, grd_real)) die("error reading file",sre);

	if ((grd_imag->nx != grd_real->nx) || (grd_imag->ny != grd_real->ny)) die("dimensions not equal!","");

	*xdim = grd_imag->nx;
	*ydim = grd_imag->ny;

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int read_data(char *imname, struct GRD_HEADER* imhdr, char *rename, struct GRD_HEADER* rehdr, struct FCOMPLEX *cdata, float *amp)
{
long	n, i;
float	*im, *re;

	n = imhdr->nx * imhdr->ny;

	im = (float *) malloc(n * sizeof(float));
	re = (float *) malloc(n * sizeof(float));

	if (GMT_read_grd(imname, imhdr, im, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) die("error reading data",imname);
	if (GMT_read_grd(rename, rehdr, re, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) die("error reading data",rename);

	for (i=0; i<n; i++){
		cdata[i].r = re[i];
		cdata[i].i = im[i];
		amp[i] = sqrtf(re[i]*re[i] + im[i]*im[i]);
		}

	free((char *) re);
	free((char *) im);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int make_wgt(float *wgt, int nxp, int nyp)
{
int	i, j;
float	wx,wy;

	for (i=0; i<nyp/2; i++){
		wy = 1.0 - fabsf((float) (i) - (nyp/2.0 - 1)) / (nyp/2.0 - 1);
		for (j=0; j<nxp/2; j++){
			wx = 1.0 - fabsf((float) (j) - (nxp/2.0 - 1)) / (nxp/2.0 - 1);
			wgt[i*nxp+j] = wgt[(i+1)*nxp-j-1] = wx*wy;
			wgt[(nyp-i-1)*nxp+j] = wgt[(nyp-i)*nxp-j-1] = wx*wy;
			}
		}

	if (debug) print_patch(nxp, nyp, wgt);

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
/* the classic Goldstein filter							*/
int apply_pspec(int m, int n, float alpha, struct FCOMPLEX *in, struct FCOMPLEX *out)
{
int 	i;
float   wgt;

	if (alpha < 0) die("alpha < 0; something is rotten in Denmark","");
	/* pow(x,a/2) == pow(sqrt(x),a) */  
	for (i=0; i<m*n; i++) {
		wgt = powf((in[i].r*in[i].r + in[i].i*in[i].i),alpha/2.0);
		out[i].r = wgt*in[i].r;	
		out[i].i = wgt*in[i].i;	
		}

	return(EXIT_SUCCESS);
}
/*-------------------------------------------------------------------------------------	*/
int calc_corr(char *amp1, struct GRD_HEADER *grd_amp1, char *amp2, struct GRD_HEADER *grd_amp2, int xdim, int ydim, float *amp, float *corr)
{
	int i,n;
	int xdim2, ydim2;
	float *a1, *a2;
	float a;

	n = xdim * ydim;

	read_file_hdr(amp1, grd_amp1, amp2, grd_amp2, &xdim2, &ydim2);

	if ((xdim != xdim2) || (ydim != ydim2)) die("amp files are different size than real and imag files","");

	a1 = (float *) malloc(xdim * ydim * sizeof(float));
	a2 = (float *) malloc(xdim * ydim * sizeof(float));

	if (GMT_read_grd(amp1, grd_amp1, a1, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) die("error reading data",amp1);
	if (GMT_read_grd(amp2, grd_amp2, a2, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) die("error reading data",amp2);

	for (i=0; i<n; i++) {
		a = a1[i] * a2[i];
		if (a > 0.0) {
			corr[i] = amp[i] / sqrtf(a);
		} else { 
			corr[i] = 0.0;
		}
		if (corr[i] < 0.0) corr[i] = 0.0;
		if (corr[i] > 1.0) corr[i] = 1.0;
		}

return(EXIT_SUCCESS);
}
/*--------------------------------------------------------------------*/
int parse_command_line(char **a, int na, char *USAGE, char *sre, char *sim, float *alp, int *ps, char *amp1, char *amp2, int *dflag, int *comflag)
{
int	n;
int	flag[4];

for (n=0; n<4; n++) flag[n] = 0;

for(n = 1; n < na; n++) {
	if (!strcmp(a[n], "-v")) {
		verbose = 1;
		/* only args after verbose (if set) will be echoed */
	} else if (!strcmp(a[n], "-real")) {
		n++;
		if (n == na) die (" no option after -real!\n","");
		strcpy(sre, a[n]);
		flag[0] = 1;
		if (verbose) fprintf(stderr, "real %s \n", sre);
	} else if (!strcmp(a[n], "-imag")) {
		n++;
		if (n == na) die (" no option after -imag!\n","");
		strcpy(sim, a[n]);
		flag[1] = 1;
		if (verbose) fprintf(stderr, "imag %s \n", sim);
	} else if (!strcmp(a[n], "-alpha")) {
		n++;
		if (n == na) die (" no option after -alpha!\n","");
		*alp = strtod(a[n],NULL);
		if (verbose) fprintf(stderr, "alpha %f \n", *alp);
	} else if (!strcmp(a[n], "-amp1")) {
		n++;
		if (n == na) die (" no option after -amp1!\n","");
		strcpy(amp1, a[n]);
		flag[2] = 1;
		if (verbose) fprintf(stderr, "amp1 %s \n", amp1);
	} else if (!strcmp(a[n], "-amp2")) {
		n++;
		if (n == na) die (" no option after -amp2!\n","");
		strcpy(amp2, a[n]);
		flag[3] = 1;
		if (verbose) fprintf(stderr, "amp2 %s \n", amp2);
	} else if (!strcmp(a[n], "-psize")) {
		n++;
		if (n == na) die (" no option after -psize!\n","");
		*ps = atoi(a[n]);
		if (verbose) fprintf(stderr, "patch size %d \n", *ps);
	} else if (!strcmp(a[n], "-diff")) {
		*dflag = 1;
		if (verbose) fprintf(stderr, "calculate difference \n");
	} else if (!strcmp(a[n], "-complex_out")) {
		*comflag = 1;
	} else if (!strcmp(a[n], "-debug")) {
		verbose = 1;
		debug = 1;
	} else {
		fprintf(stderr," %s *** option not recognized ***\n\n",a[n]);
		fprintf(stderr,"%s",USAGE);
		exit(1);
		}
	}

	/*
	flag[0]	real 
	flag[1]	imag 
	flag[2]	amp1 	only needed if modified filter
	flag[3]	amp2	only needed if modified filter 
	*/

	/* check for real and imag files are set */
	if (flag[0] == 0) die("real file not set","");
	if (flag[1] == 0) die("imag file not set","");

	/* if amp files are defined set alpha = -1 to denote coherence-based alpha */
	if ((flag[2] == 1) && (flag[3] == 1)) *alp = -1.0; 

	if ((flag[2] == 1) && (flag[3] == 0)) die("amp1 set but not amp2 (needed for modified filter)","");
	if ((flag[2] == 0) && (flag[3] == 1)) die("amp2 set but not amp1 (needed for modified filter)","");

	/* if alpha < 0 check that both amp files are set */
	if ((*alp < 0) && (flag[2] == 0)) die("amp1 file not set (needed for modified filter)","");
	if ((*alp < 0) && (flag[3] == 0)) die("amp2 file not set (needed for modified filter)","");

return(EXIT_SUCCESS);
}
/*--------------------------------------------------------------------*/
