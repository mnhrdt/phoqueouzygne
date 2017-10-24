/***************************************************************************/
/* offset_topo reads  an amplitude image and a topo_ra grid as well as     */
/* an initial guess on how to shift the topo_ra to match the master        */
/* The program uses cross correlation to estimate the refined shift to     */
/* make the topo_ra match the amplitude image more exactly.                */
/* There is an option to output a new shifted topo_ra.                     */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell and Xiaopeng Tong                           *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/22/08                                                       *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE                                                                    *
 ***************************************************************************/

# include "gmt.h"
# include "netcdf.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <string.h>
# define max(a,b) a > b ? a : b
# define min(a,b) ((a) < (b) ? (a) : (b))
/*## include "image_sio.h"
## include "siocomplex.h"*/

int main (int argc, char **argv)
{
	int i,j,k,i1,j1,k1,ibufsizea,ibufsizet;
	int is,js,ns;
	int ni,nj,ntot;
	int xshft, yshft, ib = 200;  /* ib is the width of the edge of the images not used for corr. must be > 2 */
	double ra,rt,avea,avet,suma,sumt,sumc,corr,denom;
	double maxcorr=-9999.;
	int imax,jmax;
	char *file1, *file2, *file3;
	int argc2 = 1;
	char *argv2[2] = {"dummy",0};
	float *ramp, *rtopo, *rtopo_shft;
	struct GRD_HEADER grd_amp, grd_topo;

/* execute GMT_begin */
   	argc2 = GMT_begin (argc2, argv2);

/* get the information from the command line */
	if(argc < 6){
		printf("\nUsage: offset_topo amp_master.grd topo_ra.grd rshift ashift ns [topo_shift.grd] \n \n");
		printf("   amp_master.grd - amplitude image of master \n");
		printf("   topo_ra.grd    - topo in range/azimuth coordinates of master \n");
		printf("   rshift         - guess at integer range shift \n");
		printf("   ashift         - guess at integer azimuth shift \n");
		printf("   ns             - integer search radius \n");
		printf("   topo_shift.grd - shifted topo_ra - optional, will be shifted by rshift, ashift \n \n");
		exit(-1);
	}

	file1 = argv[1];
	file2 = argv[2];
	xshft = atoi(argv[3]);
	yshft = atoi(argv[4]);
	ns = atoi(argv[5]);
	if(argc >= 7) file3 = argv[6];

/* read the header of the amplitude and topo images */
	if (GMT_read_grd_info (file1, &grd_amp)) {
		fprintf (stderr, "Error opening file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if (GMT_read_grd_info (file2, &grd_topo)) {
		fprintf (stderr, "Error opening file %s\n", file2);
		exit (EXIT_FAILURE);
	}

/* make sure the dimensions match */
        if(grd_amp.nx != grd_topo.nx) {
		fprintf (stderr, "file dimensions do not match \n");
		exit (EXIT_FAILURE);
	}

/* allocate the memory for the two files */
        ibufsizet = grd_amp.nx * grd_topo.ny;
        ibufsizea = grd_amp.nx * grd_amp.ny;
	if((ramp = (float *) malloc(ibufsizea*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for amplitude image.\n");
		exit(-1);
        }
	if((rtopo = (float *) malloc(ibufsizet*sizeof(float))) == NULL){
		fprintf(stderr,"Sorry, couldn't allocate memory for topo_ra image.\n");
		exit(-1);
        }
	if(argc >= 7) {
		if((rtopo_shft = (float *) malloc(ibufsizet*sizeof(float))) == NULL){
			fprintf(stderr,"Sorry, couldn't allocate memory for topo_shift image.\n");
			exit(-1);
        }
	}

/* read the two grids */
	if(GMT_read_grd(file1, &grd_amp, ramp, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) {
		fprintf (stderr, "Error reading file %s\n", file1);
		exit (EXIT_FAILURE);
	}
	if(GMT_read_grd(file2, &grd_topo, rtopo, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE )) {
		fprintf (stderr, "Error reading file %s\n", file2);
		exit (EXIT_FAILURE);
	}
        ni = min(grd_amp.ny,grd_topo.ny);
	fprintf(stderr," %d %d %d \n",ni,grd_amp.ny,grd_topo.ny);
        nj=grd_topo.nx;

/* compute average */
	ntot=0;
	suma=0;
	sumt=0.;
	for (i=0;i<ni;i++){
		for(j=0;j<nj;j++){
			k=i*nj+j;
			ntot++;
			suma=suma+ramp[k];
			sumt=sumt+rtopo[k];
		}
	}
	avea=suma/ntot;
	avet=sumt/ntot;

/*  compute the normalized cross correlation function */
	for(is=-ns+yshft; is<ns+1+yshft; is++){
	for(js=-ns+xshft; js<ns+1+xshft; js++){
		ntot=0;
		sumc=0.;
		suma=0.;
		sumt=0.;
		for (i=0+ib;i<ni-ib;i++){
			i1=i-is;
			for(j=0+ib;j<nj-ib;j++){
				j1=j-js;
				k=i*nj+j;
				k1=i1*nj+j1;
				if(i1 >=0 && i1 < ni && j1 >=0 && j1 < nj) {
				ntot++;
				ra=ramp[k]-avea;
                                rt=rtopo[k1+1]-rtopo[k1-1];
				sumc=sumc+ra*rt;
				suma=suma+ra*ra;
				sumt=sumt+rt*rt;
				}
			}
		}
		corr=0;
		denom=suma*sumt;
		if(denom > 0.) corr=sumc/sqrt(denom);
		/*printf(" rshift = %d  ashift = %d  correlation = %f\n",js,is,corr);*/
                if(corr > maxcorr) {
			maxcorr=corr;
			imax=is;
			jmax=js;
		}
	}
	}
	printf(" optimal: rshift = %d  ashift = %d  max_correlation = %f\n",jmax,imax,maxcorr);

/* write the shifted topo phase file */
	if(argc >= 7){
		for (i=0;i<ni;i++){
			i1=i-imax; 
			for(j=0;j<nj;j++){
				j1=j-jmax;  
				k=i*nj+j;
				k1=i1*nj+j1;
				rtopo_shft[k]=0.;
				if(i1 >=0 && i1 < ni && j1 >=0 && j1 < nj) rtopo_shft[k]=rtopo[k1];
			}
		}

/*   write the shifted grd-file */
	GMT_write_grd(file3, &grd_topo, rtopo_shft, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE);		
	free((char *) rtopo_shft);
	}
	free((char *) ramp);
	free((char *) rtopo);
   return(0);
}
