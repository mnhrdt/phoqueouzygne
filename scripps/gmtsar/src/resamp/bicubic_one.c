/************************************************************************
* bi-cubic interpolation algorithm modified from GMT                    *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>

double cubic_kernel( double, double);

void bicubic_one(double *rdata, double *idata, double x, double y, double *cz)
{
	int i, j, ij;
	double wx[4], wy[4];
        double arg, w, wsum, rsum, isum;
	double a = -0.3;

	/* These weights are based on the cubic convolution kernel, see for example
	   http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
	   These weights include a free parameter (a).
	*/

	for (i = 0; i < 4; i++) {
		arg = fabs(x + 1 - i);
		wx[i] = cubic_kernel(arg,a);
		arg = fabs(y + 1 - i);
		wy[i] = cubic_kernel(arg,a);
	}

	rsum = isum = wsum = 0.0;
	ij = 0;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
		w = wx[i] * wy[j];
		rsum += rdata[ij+i] * w;
		isum += idata[ij+i] * w;
		wsum += w;
		}
		ij += 4;
	}
	if(wsum <= 0.0) printf(" error wsum is zero \n"); 
	cz[0] = rsum/wsum;
	cz[1] = isum/wsum;
}
