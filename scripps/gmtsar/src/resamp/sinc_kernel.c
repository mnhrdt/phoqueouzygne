/************************************************************************
  computes sinc function kernel for interpolation
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/28/13                                                     *
************************************************************************/
#include <math.h>
#define PI 3.1415926535897932

double sinc_kernel (double x)
{
double arg, f;

	arg = fabs(PI*x);
        if(arg > 0.){
		f=sin(arg)/arg;
	}
	else {
		f=1.;
	}
        return (f);
}
