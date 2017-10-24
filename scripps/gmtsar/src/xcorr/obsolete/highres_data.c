/*-------------------------------------------------------*/
#include "gmtsar.h"
#include "xcorr.h"
/*-------------------------------------------------------------------------------*/
void do_highres_data(struct xcorr xc, int iloc)
{
int	i, j, nx, ny, nx_exp, ny_exp;
int	ic, jc, ii, jj, ifc;
float	xoff, yoff;
double	max_corr, cave;

	ii = jj = -99;

	ifc = xc.interp_factor;
	nx = xc.n2x;
	ny = xc.n2y;
	nx_exp = nx * ifc;
	ny_exp = ny * ifc;

	for (i=0; i <nx*ny; i++) xc.md[i].r = xc.md[i].i = xc.sd[i].r = xc.sd[i].i = 0.0;
	for (i=0; i <nx_exp*ny_exp; i++) xc.md_exp[i].r = xc.md_exp[i].i = xc.sd_exp[i].r = xc.sd_exp[i].i = 0.0;

	jc = (int) (xc.loc[iloc].xoff - xc.x_offset);
	ic = (int) (xc.loc[iloc].yoff - xc.y_offset);

	/* oversample data */
	for (i=0; i<ny; i++){
		for (j=0; j<nx; j++){
			xc.md[i*nx + j].r = (float) xc.d1[(3*xc.ysearch/2 - ny/2 + i - ic)*xc.npx + (3*xc.xsearch/2 - nx/2 + j - jc)];
			xc.sd[i*nx + j].r = (float) xc.d2[(3*xc.ysearch/2 - ny/2 + i)*xc.npy + (3*xc.xsearch/2 - nx/2 + j)];
			}
		}

	if (debug) print_int(xc.d1, xc.npy, xc.npx, 0);
	if (debug) print_int(xc.d2, xc.npy, xc.npx, 0);

	fft_interpolate_2d(xc.md, ny, nx, xc.md_exp, ny_exp, nx_exp, xc.interp_factor);
	fft_interpolate_2d(xc.sd, ny, nx, xc.sd_exp, ny_exp, nx_exp, xc.interp_factor);

	if (debug) print_complex(xc.md_exp, nx_exp, ny_exp, 1);
	if (debug) print_complex(xc.sd_exp, nx_exp, ny_exp, 1);

	fft_multiply(ny_exp, nx_exp, xc.md_exp, xc.sd_exp, xc.cd_exp);

	if (debug) print_complex(xc.cd_exp, nx_exp, ny_exp, 1);

	cave = 0.0;
	max_corr = 0.0;
	for (i=0; i<ny_exp; i++){
		for (j=0; j<nx_exp; j++){
/*
	for (i=(ny_exp/2 - ny_exp/2); i<(ny_exp/2 + ny_exp/2); i++){
		for (j=(nx_exp/2 - nx_exp/2); j<(nx_exp/2 + nx_exp/2); j++){
*/
			xc.interp_corr[i*nx_exp+j] = (Cabs(xc.cd_exp[i*nx_exp+j])*Cabs(xc.cd_exp[i*nx_exp+j])) / ((double) ny_exp * nx_exp);
			cave += xc.interp_corr[i*nx_exp+j];
			if (xc.interp_corr[i*nx_exp+j] > max_corr) {
				max_corr = xc.interp_corr[i*nx_exp+j];
				ii = i;
				jj = j;
				}
			}
		}

	cave /= (1.0 * nx_exp * ny_exp);

	yoff = (ii - ny_exp/2) / ((float) (1 * xc.interp_factor));
	xoff = (jj - nx_exp/2) / ((float) (1 * xc.interp_factor));

/*
	fprintf(stderr," high resolution: off (%d %d) diff (%d %d) ii, jj (%d %d) yoff xoff %6.3f %6.3f \n", xc.y_offset, xc.x_offset, ic, jc, ii, jj, yoff, xoff);
*/

	/* add to results */
	xc.loc[iloc].xoff -= xoff;
	xc.loc[iloc].yoff -= yoff;

	for (i=0; i<nx_exp*ny_exp;i++) xc.interp_corr[i] = 100.0 *(xc.interp_corr[i] / max_corr);

	if (debug) print_patch(xc.interp_corr, ny_exp, nx_exp);

	for (i=0; i<ny_exp; i++){
		for (j=0; j<nx_exp; j++) fprintf(stderr," %4.2f ",xc.interp_corr[i*ny_exp+j]);
		fprintf(stderr,"\n");
		}

	fprintf(stderr,"\n");

}
/*-------------------------------------------------------------------------------*/
