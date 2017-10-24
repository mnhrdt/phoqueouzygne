/************************************************************************
* conv2d is a subroutine to perform a 2-D convolution of a data array   *
* and a filter.  Both arrays are floats.                                *
************************************************************************/
/************************************************************************
* Creator: David T. Sandwell    (Scripps Institution of Oceanography    *
* Date   : 04/20/98                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*                                                                       *
* Date                                                                  *
************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
# define min(x, y) (((x) < (y)) ? (x) : (y))
# define max(x, y) (((x) > (y)) ? (x) : (y))

int conv2d(rdat,ni,nj,filt,nif,njf,fdat,ic,jc,rnorm)

int   *ni,*nj;  	/* size of rdat */
int   *nif,*njf;	/* size of filt, both numbers should be odd */
int   *ic,*jc;		/* location of filter */
float *rdat;		/* array of unfiltered data */
float *filt;		/* filter array */
float *fdat;		/* single filtered point */
float *rnorm;		/* sum of filter coefficients and used in convolution */

{
   int 	nif2,njf2,iflt,jflt;
   int 	i,j,i0,i1,j0,j1;
   float	filter;

/*  check the filter lengths */

   nif2 = (int)(*nif/2.);
   njf2 = (int)(*njf/2.);
   if(2*nif2 == *nif || 2*njf2 == *njf) {
     fprintf(stderr," nif njf %d %d should be odd \n",*nif,*njf);
     exit(-1);
   }

/*  compute the bounds of the convolution */

   i0 = max(0,*ic-nif2);
   i1 = min(*ni,*ic+nif2);
   j0 = max(0,*jc-njf2);
   j1 = min(*nj,*jc+njf2);

   *fdat = 0.;
   *rnorm = 0.;

/*  do it - but don't use any NaN's */

   for(i=i0;i<i1+1;i++) {
     iflt = i-*ic+nif2;
     for(j=j0;j<j1+1;j++) {
       jflt = j-*jc+njf2;
       filter = filt[jflt+*njf*iflt];
       *fdat  = *fdat  + filter*rdat[j+*nj*i];
       *rnorm = *rnorm + filter;
     }
   }
   return(0);
}
