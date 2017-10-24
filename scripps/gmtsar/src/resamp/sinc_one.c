/************************************************************************
* sinc function interpolation                                           *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/28/13                                                     *
************************************************************************/
#include <math.h>
#include <stdio.h>
#include "gmtsar.h"

double sinc_kernel( double );

void sinc_one(double *rdata, double *idata, double x, double y, double *cz)
{
	int i, j, ij, ns2 = NS/2-1 ;
	double wx[NS], wy[NS];
        double arg, w, wsum, rsum, isum;

	for(i = 0; i < NS ; i++){
		arg = fabs(x + ns2 - i);
		wx[i] = sinc_kernel(arg);
		arg = fabs(y + ns2 - i);
		wy[i] = sinc_kernel(arg);
	}

	rsum = isum = wsum = 0.0;
	ij = 0;
	for (j = 0; j < NS; j++) {
		for (i = 0; i < NS; i++) {
		w = wx[i] * wy[j];
		rsum += rdata[ij+i] * w;
		isum += idata[ij+i] * w;
		wsum += w;
		}
		ij += NS;
	}
	if(wsum <= 0.0) printf(" error wsum is zero \n"); 
	cz[0] = rsum/wsum;
	cz[1] = isum/wsum;
}
