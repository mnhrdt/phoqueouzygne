/***************************************************************************/
/* phase2topo reads residual phase and computes residual topography.       */
/***************************************************************************/

/***************************************************************************
 * Creator:  David T. Sandwell                                             *
 *           (Scripps Institution of Oceanography)                         *
 * Date   :  10 July, 1999                                                 *
 ***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 * DATE: Nov 9 2010 - modified to read and write grd files                 *
 *                  - use a linear relationship between the stacked        *
 *                    residual phase and DEM corrections                   *
 *                  - this linear relationship is valid because the        *
 *                    deviations from the SRTM is small (a few meters)     *
 *                                                         by Xiaopeng     *
 ***************************************************************************/

#include"gmtsar.h"
#include"lib_functions.h"
#include<math.h> 

char *USAGE = "Usage: phase2topo master.PRM topo_in.grd res_phase.grd topo_out.grd \n \n"
"    master.PRM    - master PRM files used for mapping \n"
"    topo_in.grd   - name of input topography in the radar co-ordinates of the master. \n"
"    res_phase.grd - name of input phase per unit baseline\n"
"    topo_out.grd  - name of output corrected topography in the radar co-ordinates of the master. \n\n"
"    Note the residual phase should be scaled by the perpendicular baseline (see bperp). \n";

void calc_phase2topo(int, double, double, double, double, float *, float *);

int
main(int argc, char **argv)
{
int     k,j;
int	xdim, ydim, xdimp, ydimp, ixdec; 
float   *topo, *scale, *res, *phase; 
double  drange, cnst; 
double  ymin, xmin, yinc, xinc, rland, rdum; 
struct  PRM prm; 
char    topofilename[128], phasefilename[128], title[128];
char    new_topo[128];

	rland = 999;
	debug = 0;

	if (argc < 5) die("\n", USAGE);

	/* open PRM file */ 
	get_prm(&prm, argv[1]);

	/* open grd rat file */ 
	strcpy(topofilename,argv[2]);
	readgrdsize_(&xdim, &ydim, topofilename);

	/* prepare the residual phase filename */ 
	strcpy(phasefilename,argv[3]);
	readgrdsize_(&xdimp, &ydimp, phasefilename);

	/* make sure the topo and residual phase have the same dimensions */
	if (xdimp != xdim || ydimp != ydim) { 
		die("\n", "topo and residue phase should have the same dimensions. \n"); 
        }

	/* prepare the output filename */
	strcpy(new_topo, argv[4]); 


	/* allocate the memory for the arrays */
	topo = (float *) malloc(ydim * xdim * sizeof(float));
	phase = (float *) malloc(ydim * xdim * sizeof(float));
	scale = (float *) malloc(xdim * sizeof(float));
	res = (float *) malloc(xdim * sizeof(float));


	ixdec = (int) prm.num_rng_bins/xdim; 
	drange = ixdec*SOL/(2.0*prm.fs);
	cnst = -prm.lambda/4.0/M_PI;

	readgrd_(topo, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rdum, title, topofilename);
	readgrd_(phase, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rdum, title, phasefilename);

	for (j=0; j<ydim; j++){
		for (k=0; k< xdim; k++) { 
			scale[k] = topo[j*xdim + k];
			res[k] = phase[j*xdim + k];
		}
		calc_phase2topo(xdim, prm.near_range, drange, prm.RE, prm.ht, scale, res); 
		for (k=0; k< xdim; k++) { 
			if ( phase[j*xdim + k] != rdum )  
			topo[j*xdim + k] = cnst * scale[k] + topo[j*xdim + k];
		}
	}

	/* write the grd file */ 
	rdum = rland * 1;
	writegrd_(topo, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rland, &rdum, title, new_topo);

return(EXIT_SUCCESS);
}

/*-------------------------------------------------------------*/
void calc_phase2topo(int xdim, double range0, double drange, double re, double height, float *scale, float *res)
{
int     k;
double rho,etat;
double sint,c,c2,ret,ret2;

        c = re + height;
	c2 = c * c;

        for (k=0; k<xdim; k++){
                ret = re+scale[k];
		ret2 = ret * ret;
                rho = range0 + k*drange;
                etat = ((rho*rho + c2 - ret2)/(2.*rho*c));
                if ( etat >=1. ) die("eta >= 0","");

                sint = sqrt(1. - etat*etat);
                scale[k] = res[k] * rho * c * sint / ret;

/*              printf("ret = %f res[k] = %f rho = %f c = %f sint = %f scale = %f \n", ret, res[k], rho, c, sint, scale[k]); */
                }
}
/*-------------------------------------------------------------*/

