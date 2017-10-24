#include "../../../include/soi.h"
#include "../../../include/siocomplex.h"
#include <math.h>

main(argc,argv)
{
	int n,i,j,nrows,num_valid_az,num_rng_bins,num_patches,tnrows,ipatch,low_ind,hi_ind,inum,nnrows,nn_valid_az;
	int ashft,nend,aashft,aaashft;
	fcomplex **fdata,**tdata;
	float coef;
	FILE *fopen(),*fpo,*fponew;
	coef = -0.13;
	tnrows = 74;
	nrows = 24;
	num_valid_az = 10;
	num_rng_bins = 10;
	num_patches=6;
	nnrows = num_patches*num_valid_az+nrows-num_valid_az;
	nn_valid_az = num_patches*num_valid_az;
	nlooks = 1;
/* open output file*/
	if((fpo = fopen("file_old","a")) == NULL){
          fprintf(stderr,"Bad output file pointer.\n");
          exit(-1);
        }
	if((fponew = fopen("file_new","w")) == NULL){
          fprintf(stderr,"Bad file_orig pointer.\n");
          exit(-1);
        }
/* allocate memory */
	if((tdata=(fcomplex **) malloc(tnrows*sizeof(fcomplex))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for input data.\n");
          exit(-1);
        }
        for(i=0;i<tnrows;i++){
          if((tdata[i]=(fcomplex *) malloc(num_rng_bins*sizeof(fcomplex)))==NULL){
            fprintf(stderr,"sorry, couldn't allocate memory for input data.\n");
            exit(-1);
          }
        }
	if((fdata=(fcomplex **) malloc(nrows*sizeof(fcomplex))) == NULL){
          fprintf(stderr,"Sorry, couldn't allocate memory for input data 2.\n");
          exit(-1);
        }
        for(i=0;i<nrows;i++){
          if((fdata[i]=(fcomplex *) malloc(num_rng_bins*sizeof(fcomplex)))==NULL){
            fprintf(stderr,"sorry, couldn't allocate memory for input data 2.\n");
            exit(-1);
          }
        }
	if(nlooks!=1){
          fprintf(stderr,"Bad nlooks. \n");
          exit(-1);
        }
/* make fdata */
	for(j=0;j<num_rng_bins;j++){
                for(i=0;i<tnrows;i++){
                        tdata[i][j].r = i+j;
                        tdata[i][j].i = sqrt(i+j);
                }
        }
	low_ind = (nrows-num_valid_az)/2;
        hi_ind = (nrows+num_valid_az)/2;

	for(ipatch=1;ipatch<=num_patches;ipatch++){
	    for(j=0;j<num_rng_bins;j++){
		for(i=0;i<nrows;i++){
			fdata[i][j].r = tdata[(ipatch-1)*num_valid_az+i][j].r;
			fdata[i][j].i = tdata[(ipatch-1)*num_valid_az+i][j].i;
	/*		if(j==1)fprintf(stderr,"%d %f\n",i,fdata[i][j].r); */
		}
	    }		
	    aastretch(fdata,ipatch,nrows,num_valid_az,num_rng_bins,coef);
	/*    nend = hi_ind - (num_valid_az-(int)(num_valid_az*(1-coef))); */
	    aashft = ipatch*num_valid_az - (int)(ipatch*num_valid_az*(1-coef));
	    aaashft = num_valid_az - (int)(num_valid_az*(1-coef));
	    if(ipatch==1) ashft  = num_valid_az - (int)(num_valid_az*(1-coef));
	    else if(ipatch==num_patches){
		ashft = num_valid_az - (int)(num_valid_az*(1-coef));;
	        for(i=2;i<num_patches;i++){
	   	    ashft = ashft + num_valid_az-(int)(i*num_valid_az*(1-coef))+(int)((i-1)*num_valid_az*(1-coef));
		}
		ashft = -ashft;
	    }	 
	    else{
		ashft = num_valid_az-(int)(ipatch*num_valid_az*(1-coef))+(int)((ipatch-1)*num_valid_az*(1-coef));
	    }
	    nend = hi_ind - ashft;
	    fprintf(stderr,"iashft=%d shft=%d ashft=%d nend = %d \n",aashft,aaashft,ashft,nend); 
	    pfdata(fdata,low_ind,nend,num_rng_bins,fponew);
	}
	
	aastretch(tdata,1,nnrows,nn_valid_az,num_rng_bins,coef);
	nend = nnrows-low_ind;
	pfdata(tdata,low_ind,nend,num_rng_bins,fpo);
	fclose(fpo);
	fclose(fponew);
}

pfdata(fdata,nstart,nend,num_rng_bins,fpo)
fcomplex **fdata;
int nstart,nend,num_rng_bins;
FILE *fpo;
{
	int i,j;
	for(i=nstart;i<nend;i++){
		for(j=0;j<num_rng_bins;j++){
			fprintf(fpo,"%f %f ",fdata[i][j].r,fdata[i][j].i);
			if(j==num_rng_bins-1)fprintf(fpo,"\n");
/*			if(j==num_rng_bins-1)fprintf(stderr,"%d \n",i);	*/
		}
	}
}
