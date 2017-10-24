/***************************************************************************
 * phasediff reads two complex SAR_SLC images and computes the phase       * 
 * difference of the images removing the effects of the curved earth       * 
 * and optionally the topography. Orbit information is provided in the     * 
 * second PRM-file.                                                        *
 * Model phase can also be removed as an option.                           *
 * The images should be matched to a level that will                       *
 * produce interference fringes and must be the same size.                 *
***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell  and Evelyn J. Price                        *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/7/95                                                        *
 * Rewritten: Rob Mellors (San Diego State University)                     *
 * Date   :  11/7/09                                                       *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 *                                                                         *
 * 11/18/96     Changed sign of earth flattening phase to fit with signs   *
 *              of baseline estimates.                                     *
 * 11/18/96     Changed to read GIPS headers to gather parameters.         *
 * 11/21/96     Changed to update baseline estimate.                       *
 * 06/19/97     Changed to NOT update baseline parameters                  *
 * 06/25/97     Changed to change baseline along the frame                 *
 * 03/10/99     Changed to remove reference phase due to topography        *
 * 04/06/01     Changed to remove model or other known phase               *
 * 06/21/02     Changed to scale the interferogram by 1/sqrt(amp)          *
 * 08/16/06     Changed to use the topo_ra to shift the range              *
 *              coordinate of the repeat image for long baselines.         *
 * 07/02/08     provide ability to shift topo_ra                           *
 * 07/11/08     Added spacecraft height start and end to make phase        *
 *              continuous across patches                                  *
 * 09/18/08     Added TEC_start and TEC_end to account for the ionosphere  *
 * 08/18/09     Added small correction to the phase for the elevation-     *
 *              dependent range shift, important for long baselines        *
 * 11/07/09     Code rewritten by RJM to read PRM-files and write grd-files*
 * 04/24/10     changed calc_phase subroutine based on pdiff.c             * 
 *              changed the phase calculation part to match pdiff.c        *
 *              this change is to account for slight range change due to   * 
 *              topography, which is accounted for in pdiff.c              *
 * 10/23/10     changed calc_phase to calc_drho and completely rewrote     *
 *              the topographic phase correction to use all the nonlinear  *
 *              terms.
 ***************************************************************************/

#include "gmtsar.h"
#include "lib_functions.h"

char    *USAGE = "\nUsage: "
"phasediff ref.PRM rep.PRM [-topo topo_ra.grd] [-model modelphase.grd]\n (topo_ra and model in GMT grd format)\n";

void calc_drho(int, double *, double *, double, double, double, double, double, double *);
void get_prm(struct PRM *, char *);
void print_prm_params(struct PRM, struct PRM);
void fix_prm_params(struct PRM *, char *);
void read_optional_args(int, char **, struct PRM *, int *, struct PRM *, int *);
void calc_average_topo(double *, int, int, float *);
int read_SLC_short2float(FILE *, char *, short *, fcomplex *, int, int, double);
FILE *create_GMT_binary_float(struct PRM, char *, char *);

int main (int argc, char **argv)
{
int	j, k, istart;
int	topoflag, modelflag;
int	xdim, ydim;		/* size of SLC file */
int	ydim_start;		/* start of SLC filesize */
int 	xdimt, ydimt, xt,yt;	/* size of topo file, increment */
int 	xdimm, ydimm, xm, ym;	/* size of model file, increment */
short	*d1, *d2;		/* pointers to input data files */
double  *real, *imag, *xs, *shft, *ss, *as, *topo2;
double  *range, *drho, *drho0;
float 	*optrr, *optri, *topo, *model;
double	drange, dt, tspan, time, time2;
double  ht0, htc, htf, dht, ddht, height;
double 	alpha, cnst, pha, avet;
double	B, Bh, Bv, dBh, dBv, ddBh, ddBv;
double	Bh0, Bhc, Bhf, Bv0, Bvc, Bvf;
double  ys,test;
double  xdect, ydect, xdecm, ydecm;
double  ymint, xmint, yinct, xinct, rdumt; 
double  yminm, xminm, yincm, xincm, rdumm; 
char    title[128];
FILE	*SLCfile1, *SLCfile2;
FILE	*realfile, *imagfile;
fcomplex *intfp, *iptr1, *iptr2, pshif;

struct PRM p1, p2, tp, mp;	

	verbose = 0;
	topoflag = modelflag = 0;
	xdect = ydect = 1.0;
	xdecm = ydecm = 1.0;
	avet = 0.0;
	ydim_start = 0;

	if (argc < 3) die (USAGE,"");

	/* read prm file into two pointers */
	get_prm(&p1, argv[1]);
	get_prm(&p2, argv[2]);

	if(verbose) fprintf(stderr,"near range: %lf %lf \n",p1.near_range, p2.near_range); 

	if (argc > 3) read_optional_args(argc, argv, &tp, &topoflag, &mp, &modelflag);

	if (debug) print_prm_params(p1, p2);
	if (p2.baseline_start < -9000) die("baseline < -9000 not set ?","");

	/* near_range, SC_clock_start, and SC_clock_stop need to be changed */
	fix_prm_params(&p1, argv[1]);
	fix_prm_params(&p2, argv[2]);

	if(verbose) fprintf(stderr,"near range: %lf %lf \n",p1.near_range, p2.near_range);
		
	/* open SLC files */
	if ((SLCfile1 = fopen(p1.SLC_file,"r")) == NULL) die("Can't open SLCfile",p1.SLC_file); 
	if ((SLCfile2 = fopen(p2.SLC_file,"r")) == NULL) die("Can't open SLCfile",p2.SLC_file); 

	/* set width and length */
	/* check dimensions of the two SLC files */ 
	if (p1.num_rng_bins == p2.num_rng_bins) { 
		xdim = p1.num_rng_bins;
	}
	else {
		die("The dimensions of range do not match", "");
	}
	
	if (p1.num_patches * p1.num_valid_az == p2.num_patches * p2.num_valid_az) {
                ydim = p1.num_patches * p1.num_valid_az;
        }
	else {
                die("The dimensions of azimuth do not match", "");
	}
	fprintf(stderr, " xdim %d, ydim %d \n", xdim, ydim);
	
	/* set heights */
	htc = p1.ht; 
	ht0 = p1.ht_start;
	htf = p1.ht_end;
			
	/* open output files */
	realfile = create_GMT_binary_float(p1, "real.grd", "real");
	imagfile = create_GMT_binary_float(p1, "imag.grd", "imag");

	/* allocate memory */
	drho = (double *) malloc(xdim * sizeof(double));
	drho0 = (double *) malloc(xdim * sizeof(double));
	range = (double *) malloc(xdim * sizeof(double));

	optrr = (float *) malloc(xdim * sizeof(float));
	optri = (float *) malloc(xdim * sizeof(float));

	intfp = (fcomplex *) malloc(xdim * sizeof(fcomplex));
	iptr1 = (fcomplex *) malloc(xdim * sizeof(fcomplex));
	iptr2 = (fcomplex *) malloc(xdim * sizeof(fcomplex));

	d1 = (short *) malloc(2 * xdim * sizeof(short));
	d2 = (short *) malloc(2 * xdim * sizeof(short));

	shft = (double *) malloc(xdim * sizeof(double));
	xs   = (double *) malloc(xdim * sizeof(double));
	topo2 = (double *) malloc(xdim * sizeof(double));
	imag = (double *) malloc(xdim * sizeof(double));
	real = (double *) malloc(xdim * sizeof(double));
	ss = (double *) malloc(xdim * sizeof(double));
	as = (double *) malloc(xdim * sizeof(double));

	/* open and read topo file and allocate memory */
	if (topoflag) {
		readgrdsize_(&xdimt, &ydimt, tp.input_file);
		if (xdim % xdimt != 0 || ydim % ydimt != 0)
			die("The dimension SLC must be multiplication factor of the topo_ra", tp.input_file); 
		topo = (float *) malloc(xdimt*ydimt*sizeof(float));
		readgrd_(topo, &xdimt, &ydimt, &ymint, &xmint, &yinct, &xinct, &rdumt, title, tp.input_file);				
		if(verbose) fprintf(stderr,"\n%s %s %d %d\n", title, tp.input_file, xdimt, ydimt);
		if(verbose) fprintf(stderr,"\n%f %f %f %f %f\n", xmint, ymint, xinct, yinct, rdumt); 

		/* xinct or yinct may be negative */
		xdect = abs(xinct);
		ydect = abs(yinct);
		
		/* calculate the average  and remove the average from the topography */

		calc_average_topo(&avet, xdimt, ydimt, topo);
		if(verbose) fprintf(stderr," read topo file: average %f \n", avet);
	} else {topo = NULL;}

	/* open and read the model file and allocate the memory */
	
	if (modelflag) {
		readgrdsize_(&xdimm, &ydimm, mp.input_file);
		if (xdim % xdimm != 0 || ydim % ydimm != 0)
			die("The dimension SLC must be multiplication factor of the modelphase", mp.input_file); 
		model = (float *) malloc(xdimm*ydimm*sizeof(float));
		readgrd_(model, &xdimm, &ydimm, &yminm, &xminm, &yincm, &xincm, &rdumm, title, mp.input_file);
		xdecm = abs(xincm);
		ydecm = abs(yincm);

	} else {model = NULL;}
	
        /*   compute the time span and the time spacing */

	tspan = 86400.*fabs(p2.SC_clock_stop - p2.SC_clock_start);
        dt = tspan/(ydim-1);
	if (tspan < 0.01 || p2.prf < 0.01) die("check sc_clock_start, _end, or prf","");

	/* setup the default parameters */

	drange = SOL/(2.0*p2.fs);
	alpha = p2.alpha_start*PI/180.0;
	cnst = -4.0*PI/p2.lambda;
	for (k=0;k<xdim;k++){
		range[k]=p1.near_range+k*drange;
		topo2[k]=0.;
		xs[k] = k;
	}

	/* calculate initial baselines */
	Bh0 = p2.baseline_start*cos(p2.alpha_start*PI/180.0);
	Bv0 = p2.baseline_start*sin(p2.alpha_start*PI/180.0);
	Bhf = p2.baseline_end*cos(p2.alpha_end*PI/180.0);
	Bvf = p2.baseline_end*sin(p2.alpha_end*PI/180.0);

        /* first case is quadratic baseline model, second case is default linear model */

        if(p2.baseline_center != NULL_DOUBLE ||  p2.alpha_center != NULL_DOUBLE){
		
		Bhc = p2.baseline_center*cos(p2.alpha_center*PI/180.0);
		Bvc = p2.baseline_center*sin(p2.alpha_center*PI/180.0);
		dBh = (-3.*Bh0 + 4*Bhc -Bhf)/tspan;
		dBv = (-3.*Bv0 + 4*Bvc -Bvf)/tspan;
		ddBh = (2.*Bh0 - 4*Bhc + 2*Bhf)/(tspan*tspan);
		ddBv = (2.*Bv0 - 4*Bvc + 2*Bvf)/(tspan*tspan);

        }
        else {
		dBh = (Bhf - Bh0)/tspan;
		dBv = (Bvf - Bv0)/tspan;
		ddBh = ddBv = 0.0;
        }

	/* calculate height increment */
	dht = (-3.*ht0 + 4*htc -htf)/tspan;
	ddht = (2.*ht0 - 4*htc + 2*htf)/(tspan*tspan);

	/* revise params in accordance with first_line 	*/
	if ((p1.first_line > 0) && (p2.first_line > 0 )) ydim_start = p1.first_line - 1; 

	/* now go through all the rows 		*/
	for(j=ydim_start;j<(ydim+ydim_start);j++){

		/* read data from complex i2 SLC 	*/
	 	read_SLC_short2float(SLCfile1, p1.SLC_file, d1, &iptr1[0], xdim, 1, DFACT);
	 	read_SLC_short2float(SLCfile2, p2.SLC_file, d2, &iptr2[0], xdim, 1, DFACT);

		yt = j/ydect;   /* for topo_ra */
		ym = j/ydecm;  /* for modelphase */

                /* calculate the change in baseline and height along the frame */
                time = j * dt;
                time2 = time * time;
		Bh = Bh0 + dBh*time + ddBh*time2;
		Bv = Bv0 + dBv*time + ddBv*time2;
		B  = sqrt(Bh*Bh + Bv*Bv);
		alpha = atan2(Bv,Bh);
		height = ht0 + dht*time + ddht * time2;

		for (k=0;k<xdim;k++) { 
			shft[k]=0.;
			if(topoflag) {
				xt=k/xdect;
				topo2[k]= topo[xt+xdimt*yt];
			}
		}

		/* calculate the combined earth curvature and topography correction */
		calc_drho(xdim,range,topo2,avet,p1.RE,height,B,alpha,drho);

		/* loop over range to make topographic and model phase corrections */
		for (k=0; k<xdim; k++){
			intfp[k] = iptr1[k];
			pha=cnst*drho[k];
			if (modelflag) { 
				xm = k/xdecm; /* xm is increment for model phase. note they are all integers */
				pha = pha - model[xm+xdimm*ym]; 
			}
			pshif = Cexp(pha);
			intfp[k] = Cmul(intfp[k],pshif);
		}
		
		/* shift the range of the repeat image to improve image matching for very long baselines > 1000 m */
		if ((topoflag > 0)  && (p2.baseline_start > 1000.0)) { 

		/* compute the range change with no topography so the range shift can be determined for the spline */
			calc_drho(xdim,range,shft,avet,p1.RE,height,B,alpha,drho0);

	        	for(k=0;k<xdim;k++){
				shft[k]=(drho0[k]-drho[k])/drange;
                		real[k]=intfp[k].r;
               	 		imag[k]=intfp[k].i;
        		}
	 	/* shift the real part */
			spline_(&istart,&xdim, xs, real, ss, as);
		        for(k=0; k<xdim; k++) {
                		ys = xs[k] + shft[k];
                		evals_(&istart,&ys, &xdim, xs, real, ss, &test);
                		intfp[k].r = test;
			}
		/* shift imaginary part  */
        		spline_(&istart,&xdim, xs, imag, ss, as);
        		for(k=0; k<xdim; k++) {
                		ys = xs[k] + shft[k];
                		evals_(&istart,&ys, &xdim, xs, real, ss, &test);
                		intfp[k].i = test;
			}
		}

		/* make interferogram */
		for(k=0;k<xdim;k++){
			iptr2[k] = Conjg(iptr2[k]);
      			intfp[k] = Cmul(intfp[k],iptr2[k]);
			optrr[k] =  intfp[k].r;
			optri[k] =  intfp[k].i; 
		}

		/* write out data (row j) */
		fwrite(optrr, sizeof(float), xdim, realfile);
		fwrite(optri, sizeof(float), xdim, imagfile);
	}

	fclose(realfile);
	fclose(imagfile);

	return(EXIT_SUCCESS);
}


/*--------------------------------------------------------------*/
void calc_drho(int xdim, double *range, double *topo, double avet, double re,
		double height, double B, double alpha, double *drho)
{
int k;
double rho,sint,cost,cosa,sina;
double term1,term2,c,c2,ret,ret2;

	sina = sin(alpha);
	cosa = cos(alpha);
	c = re + height;
	c2 = c*c;
	for(k=0; k<xdim; k++){

/* compute the look angle using equation (C26) in Appendix C */
		rho = range[k];
		ret = re+avet+topo[k];
		ret2 = ret*ret;
		cost = ((rho*rho + c2 - ret2)/(2.*rho*c));
		if( cost >= 1.) die("calc_drho","cost >= 0"); 
		sint = sqrt(1. - cost*cost);

/* compute the range change using equation (c23) in Appendic C */
 		term1 = -B*(sint*cosa-cost*sina);
		term2 = B*B*(cost*cosa+sint*sina)/(2.*rho);
		drho[k] = term1 + term2; 
  
/* New (Eric Lindsey, April 2015): compute the range change using the full nonlinear equation */
/*  		term1 = rho*rho + B*B - 2*rho*B*(sint*cosa-cost*sina);
		drho[k] = -rho + sqrt(term1);*/
	}
}

/*--------------------------------------------------------------*/
void calc_average_topo(double *avet, int xdimt, int ydimt, float *topo)
{
double 	sumt;
int    	k, nsum;

	sumt = 0.0;
	nsum = 0;

/* compute the average topography, save the value, and remove it from the topography */
	for (k=0; k < xdimt*ydimt;k++) {
		sumt += topo[k];
		nsum++;
	}
	*avet = sumt/nsum;
	for (k=0; k < xdimt*ydimt;k++) {
		topo[k] = topo[k] - *avet;
	}
	if (verbose) fprintf(stderr," mean topo %lf\n",*avet);
}
