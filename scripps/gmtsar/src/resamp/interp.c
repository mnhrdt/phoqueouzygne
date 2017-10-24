/************************************************************************
* nearest, bilinear, and bicubic interpolations                         *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>
#include "gmtsar.h"

void nearest(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	int i, j, k;

	/* compute the indices of the upper left corner */

	j = (int)(ras[0] + 0.5);
	i = (int)(ras[1] + 0.5);
	k = 2*xdims*i + 2*j;

	/* use the nearest point if it is within the bounds of the slave array */

	if(i < 0 || i >= ydims || j < 0 || j >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {
	  sout[0] = s_in[k];
	  sout[1] = s_in[k+1];
	}
}

void bilinear(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da, real, imag;
	int k00, k01, k10, k11;
	int i0, j0;
	int nclip;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

	/* compute the indices of the 4 corners */
       
	k00 = 2*xdims*i0     + 2*j0;
	k01 = 2*xdims*i0     + 2*(j0+1);
	k10 = 2*xdims*(i0+1) + 2*j0;
	k11 = 2*xdims*(i0+1) + 2*(j0+1);

	/* do the interpolation if all 4 corners are within the bounds of the slave array */

	if(i0 < 0 || i0 >= (ydims-1) || j0 < 0 || j0 >= (xdims-1)) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {
        real = s_in[k00] * (1.0 - da) * (1.0 - dr)
             + s_in[k10] * (da)       * (1.0 - dr)
             + s_in[k01] * (1.0 - da) * (dr) 
             + s_in[k11] * (da)       * (dr);
        if((int)fabs(real) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(real + 0.5);
        imag = s_in[k00+1] * (1.0 - da) * (1.0 - dr)
             + s_in[k10+1] * (da)       * (1.0 - dr)
             + s_in[k01+1] * (1.0 - da) * (dr) 
             + s_in[k11+1] * (da)       * (dr);
        if((int)fabs(imag) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(imag + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}
void bicubic(double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da;
	double rdata[16], idata[16], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

	/* make sure all 4 corners are within the bounds of the slave array */

	if((i0-1) < 0 || (i0+2) >= ydims || (j0-1) < 0 || (j0+2) >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {

	/* safe to do the interpolation */

	for (i=0; i<4; i++){
		for (j=0; j<4; j++){
		k = i*4 +j;
		kk = 2*xdims*(i0-1+i)  + 2*(j0-1+j);
		rdata[k] = s_in[kk];
		idata[k] = s_in[kk+1];
		}
	}

	/* interpolate the real and imaginary data */

	bicubic_one(rdata, idata, dr, da, cz);

        if((int)fabs(cz[0]) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(cz[0] + 0.5);
        if((int)fabs(cz[1]) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}

void bisinc (double *ras, short *s_in, int ydims, int xdims, short *sout)
{
	double dr, da, ns2=NS/2-1;
	double rdata[NS*NS], idata[NS*NS], cz[2];
	int i, j, k, kk;
	int i0, j0;
	int nclip;

	/* compute the residual offsets */
	nclip = 0;
	j0 = (int)floor(ras[0]);
	i0 = (int)floor(ras[1]);
        dr = ras[0] - (double)j0;
        da = ras[1] - (double)i0;
	if(dr < 0. || dr > 1. || da < 0. || da > 1) fprintf(stderr," dr or da out of bounds %f %f \n",dr,da);

	/* make sure all 4 corners are within the bounds of the slave array */

	if((i0-ns2) < 0 || (i0+ns2+1) >= ydims || (j0-ns2) < 0 || (j0+ns2+1) >= xdims) {
	  sout[0] = 0;
	  sout[1] = 0;
	}
	else {

	/* safe to do the interpolation */

	for (i=0; i<NS; i++){
		for (j=0; j<NS; j++){
		k = i*NS +j; 
		kk = 2*xdims*(i0-ns2+i)  + 2*(j0-ns2+j);
		rdata[k] = s_in[kk];
		idata[k] = s_in[kk+1];
		}
	}

	/* interpolate the real and imaginary data */

	sinc_one(rdata, idata, dr, da, cz);

        if((int)fabs(cz[0]) > I2MAX) nclip = nclip + 1;
	sout[0] = (short)clipi2(cz[0] + 0.5);
        if((int)fabs(cz[1]) > I2MAX) nclip = nclip + 1;
	sout[1] = (short)clipi2(cz[1] + 0.5);
	}
	/*if(nclip > 0) fprintf(stderr," %d integers were clipped \n",nclip);*/
}
