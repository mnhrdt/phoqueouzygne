/*-------------------------------------------------------*/
#include "gmtsar.h"
#include "xcorr.h"
/*-------------------------------------------------------------------------------*/
void do_highres_corr(struct xcorr xc, int iloc)
{
int	i, j, ic, jc;
int	nx, ny, nx2, ny2, ifc;
float	max_corr, ipeak, jpeak, sub_xoff, sub_yoff;

	ifc = xc.interp_factor;

	/* size of complex version of correlation 	*/
	/* must be power of two 			*/
	nx = xc.n2x;
	ny = xc.n2y;

	/* size of interpolated matrix			*/
	nx2 = ifc*nx;
	ny2 = ifc*ny;

	/* remove calculated offset from 1 pixel resolution 	*/
	/* factor of 2 on xoff for range interpolation		*/
	jc = (xc.nxc / 2) - nx / 2 - (int) xc.loc[iloc].xoff;
	ic = (xc.nyc / 2) - ny / 2 - (int) xc.loc[iloc].yoff;

	/* copy values from correlation to complex 	*/
	/* use values centered around highest value	*/
	for (i=0; i<ny; i++) {
		for (j=0; j<nx; j++) {
			xc.md[i*nx+j].r = powf(xc.corr[(ic+i)*xc.nxc + (jc+j)], 0.25);
			xc.md[i*nx+j].i = 0.0;
			}
		}

	if (debug) print_complex(xc.md, nx, ny, 1);

	fft_interpolate_2d(xc.md, ny, nx, xc.cd_exp, ny2, nx2, ifc);

	if (ifc <=4) print_complex(xc.cd_exp, nx2, ny2, 1);

	/* find maximum in interpolated matrix		*/
	max_corr = -1.0;
	ipeak = jpeak = -1;
	for (i=0; i<ny2; i++) {
		for (j=0; j<nx2; j++) {
			if (xc.cd_exp[i*nx2+j].r > max_corr) {
				max_corr = xc.cd_exp[i*nx2+j].r;
				ipeak = i - (ny2 / 2);	
				jpeak = j - (ny2 / 2);	
				}
			}
		}

	/* fft interpolation */
	sub_xoff = jpeak / (float) (ifc);
	sub_yoff = ipeak / (float) (ifc);

	if (debug) {
		fprintf(stderr," highres [ri %d ifc %d](%4.1f %4.1f) (nx %d ny %d) jpeak %f ipeak %f : %f %f : %4.2f\n", 
		xc.ri, ifc, xc.loc[iloc].xoff, xc.loc[iloc].yoff, nx2, ny2, jpeak, ipeak, sub_xoff, sub_yoff,  max_corr);
		}

	xc.loc[iloc].xfrac = sub_xoff;
	xc.loc[iloc].yfrac = sub_yoff;
}
/*-------------------------------------------------------------------------------*/
