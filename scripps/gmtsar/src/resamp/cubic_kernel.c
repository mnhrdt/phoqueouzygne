/************************************************************************
  kernel computes a bi-cubic spline kernel using the formula given at 
  the following web page
  http://undergraduate.csse.uwa.edu.au/units/CITS4241/Handouts/Lecture04.html
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
#include <stdio.h>

double cubic_kernel(double arg, double a)

/* note arg must be positive and a must be between -3 and 0. */
{

double arg2, arg3, f;

	arg2=arg*arg;
	arg3=arg2*arg;
	if(arg <= 1.){
		f=(a+2)*arg3-(a+3)*arg2+1.;
	}
	else if(arg <=2.){
		f=a*arg3-5*a*arg2+8*a*arg-4*a;
	}
	else {
		f=0.;
	}
        return (f);
}

