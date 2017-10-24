/************************************************************************
* strech performs a range shift depending on range using interpolation  *
************************************************************************/
/************************************************************************
* Creator: Meng Wei, Scripps Institution of Oceanography)     		*
* Date   : 09/29/06                                                     *
************************************************************************/
/************************************************************************
* Modification History                                                  *
*                                                                       *
* Date                                                                  *
************************************************************************/

#include "../../../include/soi.h"
#include "../../../include/siocomplex.h"
#include <math.h>

aastretch(fdata,ipatch,nrows,num_valid_az,num_rng_bins,coef)
int ipatch,nrows,num_valid_az,num_rng_bins;
fcomplex **fdata;
float coef;
{	
	fcomplex *vec_ipt,*vec_out;
	int i,j,n,an,dx,low_ind;
	float *x;
/* if abs(coef) is too big, error */
	if(coef*num_valid_az>=10||coef*num_valid_az<=-10){
	  fprintf(stderr,"The a_stretch_a is too large.\n");
          exit(-1);
        }
/* allocate memory for vectors */
	low_ind = (nrows-num_valid_az)/2;
	if((vec_ipt = (fcomplex *) malloc((nrows-low_ind)*sizeof(fcomplex))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for vec_input.\n");
          exit(-1);
        }
	if((vec_out = (fcomplex *) malloc((nrows-low_ind)*sizeof(fcomplex))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for vec_output.\n");
          exit(-1);
        }
	if((x = (float *) malloc((nrows-low_ind)*sizeof(float))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for shift.\n");
          exit(-1);
        }
/* calculate shift for per colomn */
	for(j=0;j<num_rng_bins;j++){
		n = 0;
		for(i=0;i<nrows-low_ind;i++){	
			vec_ipt[i].r = fdata[i+low_ind][j].r;
			vec_ipt[i].i = fdata[i+low_ind][j].i;
		}
		for(i=0;i<nrows-low_ind;i++){
		  if(n<nrows-low_ind-1){
                	x[i] = ((ipatch-1)*num_valid_az + i)*(1-coef);
			an = (int)((ipatch-1)*num_valid_az*(1-coef))+ n;
		/*	if(j==1)fprintf(stderr,"i=%d,xi=%f,an=%d \n",i,x[i],an);  */
                	if(i==0){
			   vec_out[n].r = vec_ipt[i].r;
                           vec_out[n].i = vec_ipt[i].i;
                           n = n + 1;
                        }
                	else{
			  dx = (int)x[i]-(int)x[i-1];
                          if(dx==1){
                            vec_out[n].r = (vec_ipt[i-1].r*(x[i]-an)+vec_ipt[i].r*(an-x[i-1]))/(x[i]-x[i-1]);
                            vec_out[n].i = (vec_ipt[i-1].i*(x[i]-an)+vec_ipt[i].i*(an-x[i-1]))/(x[i]-x[i-1]);
                            n = n + 1;
                          }
                          else if(dx==2){
                             vec_out[n].r = (vec_ipt[i-1].r*(x[i]-an)+vec_ipt[i].r*(an-x[i-1]))/(x[i]-x[i-1]);
                             vec_out[n].i = (vec_ipt[i-1].i*(x[i]-an)+vec_ipt[i].i*(an-x[i-1]))/(x[i]-x[i-1]);
                             vec_out[n+1].r = (vec_ipt[i-1].r*(x[i]-an-1)+vec_ipt[i].r*(an+1-x[i-1]))/(x[i]-x[i-1]);
                             vec_out[n+1].i = (vec_ipt[i-1].i*(x[i]-an-1)+vec_ipt[i].i*(an+1-x[i-1]))/(x[i]-x[i-1]);
                             n = n + 2;
                          }
			  else{}
			}
		   }
        	}
		for(i=0;i<nrows-low_ind;i++){
                   fdata[i+low_ind][j].r = vec_out[i].r;
                   fdata[i+low_ind][j].i = vec_out[i].i;
                 }
	}
/*	fprintf(stderr,"aastretch once \n");  */
	free((fcomplex *) vec_ipt);
	free((fcomplex *) vec_out);
	free((float *) x);
}
