/***************************************************************************
 * bperp reads the PRM file of a repeat pass file and extracts the         *
 * baseline information.  Then it computes a complete array of             *
 * perpendicular baseline that depends on both range and azimuth.          *
 * This baseline array is used for taking sums and differences of          *
 * interferograms                                                          *
 ***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  7/7/95                                                        *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *     Modified from ihBperp to read and write grd files.                  *
 *                              Xiaopeng Nov 24 2010                       *
 * DATE                                                                    *
 ***************************************************************************/

#include"gmtsar.h"
#include"lib_functions.h"
#include<math.h>


char *USAGE = "Usage: bperp rep.PRM phase.grd bperp.grd \n\n"
"   rep.PRM    -  input repeat PRM file that provides the interferometric baseline \n"   
"   phase.grd  -  input phase grd file that provides the dimensions of the grd file \n"   
"   bperp.grd  -  output grd file of the perpendicular baseline\n\n"; 

void get_prm(struct PRM *, char *);
void fix_prm_params(struct PRM *, char *);

int 
main(int argc, char **argv)
{
struct      PRM prm;
float       *phase, *bp, *theta;
double      drange; 
double      ymin, xmin, yinc, xinc, rland, rdum; 
double      dBh, dBv, Bh0, Bhf, Bv0, Bvf, B, alpha, Bh, Bv, dht, height;
char        grdfile[128], outfile[128], title[128]; 
int         xdim, ydim;
int         j, k;  

	if(argc < 4) die("\n", USAGE); 

	/* open PRM file */ 
	get_prm(&prm, argv[1]);

	/* near_range, SC_clock_start, and SC_clock_stop need to be changed */
        fix_prm_params(&prm, argv[1]);

	/* open phase grd file to get the grd file size */ 
	strcpy(grdfile,argv[2]); 
	readgrdsize_(&xdim, &ydim, grdfile);

	/* prepare the output bperp grd file */ 
	strcpy(outfile,argv[3]);

	/* allocate memory */ 
	phase = (float *) malloc(ydim * xdim * sizeof(float));
	bp = (float *) malloc(ydim * xdim * sizeof(float));
	theta = (float *) malloc(xdim * sizeof(float));

	readgrd_(phase, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rdum, title, grdfile);

/* calculate drange */

        drange = abs(xinc)*SOL/(2.0*prm.fs);

/* compute starting and ending baselines in rectangular co-ordinates */

        dBh = 0.;
        dBv = 0.;
        Bh0 = prm.baseline_start*cos(prm.alpha_start*PI/180.);
        Bhf = prm.baseline_end*cos(prm.alpha_end*PI/180.);
        Bv0 = prm.baseline_start*sin(prm.alpha_start*PI/180.);
        Bvf = prm.baseline_end*sin(prm.alpha_end*PI/180.);
        dBh = (Bhf - Bh0)/ydim;  
        dBv = (Bvf - Bv0)/ydim;    
        dBh = (Bhf - Bh0)/ydim;  
        dBv = (Bvf - Bv0)/ydim;    

/* calculate height increment if available */
        dht = 0.0;        

        if (prm.ht_start > 0.0 && prm.ht_end > 0.0) {
	  dht = (prm.ht_end - prm.ht_start)/ydim;
	}

        for(j=0;j<ydim;j++){

/* change the baseline and alpha along the satellite track */

        Bh = Bh0 + dBh*j;
        Bv = Bv0 + dBv*j;
        B  = sqrt(Bh*Bh + Bv*Bv);
        alpha = atan2(Bv,Bh);
	height = prm.ht_start + dht*j;

        calc_theta(xdim,prm.near_range,drange,prm.RE,height,alpha,theta);

/* compute the perpendicular baseline */
          for(k=0;k<xdim;k++){
                bp[j*xdim+k] = B*cos((double)theta[k]-alpha);
          }
	}
	

/* write the grd file */ 
	strcpy(title, "perpendicular baseline");
	writegrd_(bp, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rland, &rdum, title, outfile);

return(EXIT_SUCCESS);
} 

calc_theta(xdim,range0,drange,re,height,alpha,theta)
int xdim;
double range0,drange,re,height,alpha;
float *theta;
{
        int k;
        double rho,eta;
        double c,c2,re2;

        c = re + height;
        c2 = c*c;
        re2 = re*re;

        for(k=0;k<xdim;k++){
                rho = range0 + k*drange;
                eta = ((rho*rho + c2 - re2)/(2.*rho*c));
		if ( eta >=1. ) die("eta >= 0","");
                theta[k]=acos(eta);
        }
}
	
