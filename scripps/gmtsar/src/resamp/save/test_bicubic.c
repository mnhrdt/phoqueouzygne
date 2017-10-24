/************************************************************************
    program to test the bicubic interpolation routines
************************************************************************/
#include <stdio.h>
#include <time.h>

main(argc,argv)
int argc;
char *argv[];
{
	double rdata[16],idata[16],cz[2];
        double x, y;
	int i, j;

/* put the value 1 in various array elements and check output */

	for(i=0;i<16;i++){
		rdata[i]=0.0;
		idata[i]=0.0;
	}
	rdata[5]=1.;
	rdata[6]=1.;
	rdata[9]=1.;
	rdata[10]=1.;

/* loop over a very dense array of x and y */

	for(j=0;j<101;j++){
		x=(double)j/100.;
	for(i=0;i<101;i++){
		y=(double)i/100.;
	bicubic_one(rdata,idata,x,y,cz);
        printf(" %f %f %f \n",x,y,cz[0]);
	}
	}
}

