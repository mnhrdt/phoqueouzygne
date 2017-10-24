/************************************************************************
* trans_col calls a 1-D fft on columns of data.				*
************************************************************************/
/************************************************************************
* Creator: Evelyn J. Price	(Scripps Institution of Oceanography)	*
* Date   : 11/18/96							*
************************************************************************/
/************************************************************************
* Modification History							*
* 									*
* Date									*
*									*
************************************************************************/

#include "../../include/soi.h"
#include <stdlib.h>
#include <math.h>

trans_col(xnum,ynum,data)
int xnum,ynum;
fcomplex **data;
{
	int i,k,dir;
	fcomplex *fft_vec;

	if((fft_vec = (fcomplex *) malloc(ynum*sizeof(fcomplex))) == NULL){
	  fprintf(stderr,"trans_col: Can't allocate memory for fft_vec.\n");
	  return(-1);
	}

	for(i=0;i<xnum;i++){
	  for(k=0;k<ynum;k++){
	      fft_vec[k] = data[k][i];
	  }
          dir = -1;
	  cfft1d_(&ynum,fft_vec,&dir);
	  for(k=0; k<ynum;k++){
	    data[k][i] = fft_vec[k];
	  }
	}
	free((char *) fft_vec);
}
