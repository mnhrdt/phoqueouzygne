/************************************************************************
    program to test speed and accuracy of FFT routines

	Sun Performance Library

************************************************************************/
/************************************************************************
* David T. Sandwell		(Scripps Institution of Oceanography)   *
* Date   : 12/23/96						        *
************************************************************************/
/************************************************************************
* Modification History:                                                	*
*								        *
* Date                                                                  *
*                                                                       *
************************************************************************/ 
#include <stdio.h>
#include <time.h>
#include "../../../include/soi.h"

main(argc,argv)
int argc;
char *argv[];
{
        int ndata=32768,fsgn,ier,dir;
	int n,i,k,ntimes=2048;
	double delr,dmin,dmax;
	fcomplex *fdata,*fft_vec,RCmul(),Cmul(),Conjg();

        printf("compute speed and accuracy of FFT\n");


/* allocate memory */

	if((fft_vec = (fcomplex *) malloc(ndata*sizeof(fcomplex))) == NULL){
	  fprintf(stderr,"Sorry, couldn't allocate memory\n");
	  exit(-1);
	}

/*  fill the array */
	for(i=0;i<ndata;i++){
          fft_vec[i].r = 0.;
          fft_vec[i].i = 0.;
	}
        fft_vec[0].r =  1.;
        fft_vec[0].i = -1.;
        fft_vec[1].r = -1.;
        fft_vec[1].i =  1.;
        fft_vec[2].r =  1.;
        fft_vec[2].i = -1.;

/*  print the array */

	printf("\n first 10 elements \n");
	for(i=0;i<10;i++){
	  printf(" %e %e \n",fft_vec[i].r,fft_vec[i].i);
	}

/* do the FFTs */

        dir=0;
          cfft1d_(&ndata,fft_vec,&dir);

	printf(" doing the fft ndata, ntimes, %d %d \n ", ndata, ntimes);
 	clock();   
    	for(k=0;k<ntimes;k++){
          dir=-1;
          cfft1d_(&ndata,fft_vec,&dir);

          dir=1;
            cfft1d_(&ndata,fft_vec,&dir);
	}

/*  print the array */

	printf("\n first 10 elements \n");
	for(i=0;i<10;i++){
	  printf(" %e %e \n",fft_vec[i].r,fft_vec[i].i);
	}
	print_time((float) clock()/CLOCKS_PER_SEC); 

	free((char *) fft_vec);
}

print_time(timer)
float timer;
{
	int min;
	float sec;

	min = (int) timer/60;
	sec = timer - ((float) min)*60.0;

	fprintf(stderr,"\nProcessing Elapsed time: %d min %.2f sec\n\n",min,sec);
}
