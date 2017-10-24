/************************************************************************
    program to test the sinc interpolation routine
************************************************************************/
#include <stdio.h>
#include <time.h>

main(argc,argv)
int argc;
char *argv[];
{
	double rdata[64],idata[64],cz[2];
        double x, y;
	int i, j;

/* put the value 1 in various array elements and check output */

	for(i=0;i<64;i++){
		rdata[i]=0.0;
		idata[i]=0.0;
	}
	rdata[27]=1.;
/*	rdata[28]=1.;
	rdata[35]=1.;
	rdata[36]=1.;*/

/* loop over a very dense array of x and y */

	for(j=0;j<101;j++){
		x=(double)j/100.;
	for(i=0;i<101;i++){
		y=(double)i/100.;
	sinc_one(rdata,idata,x,y,cz);
        printf(" %f %f %f \n",x,y,cz[0]);
	}
	}
}

