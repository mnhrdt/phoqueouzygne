/*****************************************************************************************
 *  Program to do InSAR time-series analysis.                                            *
 *  Use Small Baseline Subset (SBAS) algorithm.                                          *
 ***************************************************************************************** 
 * Creator: Xiaopeng Tong and David Sandwell 						 *
 *          (Scripps Institution of Oceanography)					 *
 * Date: 12/23/2012  									 *
 ****************************************************************************************/ 
/*****************************************************************************************
 *  Modification history: 								 *
 *  08/31/2013 debug the program							 *
 *  03/20/2014 add DEM error, mean velocity 						 *
 *  03/22/2014 add correlation, use weighted least-squares  				 *
 *  04/01/2014 add seasonal term			  				 *
 *  08/19/2014 allocate memory with 1D array instead of multiple malloc			 *
 *  08/19/2014 do not require the velocity curve go through origin 			 *
 *  08/19/2014 remove seasonal term	 			 			 *
 *  08/19/2014 fix temporal smoothing	 			 			 *
 ****************************************************************************************/

/* Reference: 
P. Berardino, G. Fornaro, R. Lanari, and E. Sansosti, “A new algorithm
for surface deformation monitoring based on small baseline differential
SAR interferograms,” IEEE Trans. Geosci. Remote Sensing, vol. 40, pp.
2375–2383, Nov. 2002.

Schmidt, D. A., and R. Bürgmann(2003),
Time-dependent land uplift and subsidence in the Santa Clara valley, 
California, from a large interferometric synthetic aperture radar data set, 
J. Geophys. Res., 108, 2416, doi:10.1029/2002JB002267, B9.
*/

/* Use DGELSY to solve the equations */
/* Calling DGELSY using column-major order */

# include <stdio.h>
# include <string.h>
# include <math.h>
# include <stdlib.h>
# include "gmtsar.h"
# define Malloc(type,n) ((type*)malloc((n)*sizeof(type)))
# define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
# ifdef DEBUG
# define checkpoint() printf("Checkpoint at line %d in file %s\n", __LINE__, __FILE__) 
# else 
# define checkpoint()
# endif

void dgelsy_(const int *m, const int *n,const int *nrhs,double *G, const int *lda,double* b, const int *ldb,int *jpvt,const double *rcond, const int *rank,double *work,const int *lwork,const int *info); 
void die (char *s1 , char *s2);

char *USAGE = " \n\nUSAGE: sbas intf.tab scene.tab N S xdim ydim [-smooth sf] [-wavelength wl] [-incidence inc] [-range -rng] [-rms] [-dem]\n\n"
" input: \n"
"intf.tab		--  list of unwrapped (filtered) interferograms:\n"
"   format:   unwrap.grd  corr.grd  ref_id  rep_id  B_perp \n"
"scene.tab		--  list of the SAR scenes in chronological order\n"
"   format:   scene_id   number_of_days \n"
"   note:     the number_of_days is relative to a reference date \n"
"N             		--  number of the interferograms\n"
"S             		--  number of the SAR scenes \n"
"xdim and ydim 		--  dimension of the interferograms\n"
"-smooth sf		--  smoothing factors, default=0 \n"
"-wavelength wl		--  wavelength of the radar wave (m) default=0.236 \n"
"-incidence theta 	--  incidence angle of the radar wave (degree) default=37 \n" 
"-range rng 		--  range distance from the radar to the center of the interferogram (m) default=866000 \n" 
"-rms 			--  output RMS of the data misfit grids (mm): rms.grd\n" 
"-dem 			--  output DEM error (m): dem.grd \n\n" 
" output: \n"
"disp_##.grd   		--  cumulative displacement time series (mm) grids\n"
"vel.grd 		--  mean velocity (mm/yr) grids \n\n"
" example:\n"
"   sbas intf.tab scene.tab 88 28 700 1000 \n\n";

int main(int argc, char **argv)
{
	/* define variables */

	char ** gfile, ** cfile;
	char uplo[1]="U";
        char title[64],outfile[36],tmp1[36],tmp2[36],tmp3[36];
	int i,j,k,p,m,n,nrhs=1,info,rank,xdim,lwork,lworkM,ydim,xin,yin;
	int nM=4,ldaM,ldbM;
        int N,S;
        int ldb,lda,*flag,*jpvt,*jpvtM,*H,*L;
	int flag_rms=0,flag_dem=0;
	float *phi,*grdin,*corin,sf,*disp,*res,*dem,*bperp,*vel,*amp,*shift;
	float *var,old,new;
	double *G,*A,*Gs,*d,*ds,*work,*workM,*time,*gM,*dM;
	double rcond=1e-3,pred,rng,wl,theta,scale;
	double xmin,ymin,xinc,yinc,rdum,rland;
	double sumxx,sumxy,sumx,sumy;
	FILE *infile, *datefile;
	
	if (argc < 7) die("\n",USAGE);

	/* default parameters are for ALOS-1 */
	/* use range at the center of the image*/
	rng=866000;
	/* wavelength of the radar wave*/
	wl=0.236;
	/* incidence angle at the center of the image*/
	theta=37;  
	/* smoothing factor */
	sf=0;

	/* parse the commands */
	if ((infile = fopen(argv[1],"r")) == NULL) die("Can't open file",argv[1]); 
        if ((datefile = fopen(argv[2],"r")) == NULL) die("Can't open file",argv[2]);
	N = atoi(argv[3]);
        S = atoi(argv[4]);
        xdim = atoi(argv[5]);
        ydim = atoi(argv[6]);

	fprintf(stderr,"\n");
	
	for (i=7;i<argc;i++) {
		if (!strcmp(argv[i],"-smooth")) { 
			i++;
			if (i==argc) die("no option after -smooth!\n","");
			sf=atof(argv[i]);
		} else if (!strcmp(argv[i],"-wavelength")) {
			i++;
			if (i==argc) die("no option after -wavelength!\n","");
			wl=atof(argv[i]);
		} else if (!strcmp(argv[i],"-incidence")) {
			i++;
			if (i==argc) die("no option after -incidence! \n","");
			theta=atof(argv[i]);
		} else if (!strcmp(argv[i],"-range")) {
			i++;
			if (i==argc) die("no option after -range \n","");
			rng=atof(argv[i]);
		} else if (!strcmp(argv[i],"-rms")) {
			flag_rms=1;
			fprintf(stderr,"compute RMS misfit\n");
		} else if (!strcmp(argv[i],"-dem")) {
			flag_dem=1;
			fprintf(stderr,"compute DEM error\n");
		} else {
                fprintf(stderr," %s *** option not recognized ***\n\n",argv[n]);
                fprintf(stderr,"%s",USAGE);
                exit(1);
		}
	}

        fprintf(stderr,"\nsetting smoothing to %7.3lf \n",sf);
        fprintf(stderr,"setting radar wavelength to %7.3lf m \n",wl);
        fprintf(stderr,"setting radar incidence angle to %7.3lf degree \n", theta);
        fprintf(stderr,"setting range to %9.3lf m \n", rng);

        scale=4.0*M_PI/wl/rng/sin(theta/180.0*M_PI);

	m = N+S-2;
	n = S; 

	lwork=max(1,m*n+max(m*n,nrhs)*16);
	lda=max(1,m);
	ldb=max(1,max(m,n));
	ldaM=max(1,S);
	ldbM=max(1,max(S,nM));
	lworkM=max(1,S*nM+max(S*nM,1)*16);

/* memory allocation */

	if ((jpvt = Malloc(int,n)) == NULL) die("memory allocation!","jpvt");
	if ((jpvtM = Malloc(int,nM)) == NULL) die("memory allocation!","jpvtM");
        if ((work = Malloc(double, lwork)) == NULL) die("memory allocation!","work");
        if ((workM = Malloc(double, lworkM)) == NULL) die("memory allocation!","workM");
        if ((d = Malloc(double,ldb)) == NULL) die("memory allocation!","d");
        if ((ds = Malloc(double,N)) == NULL) die("memory allocation!","ds");
        if ((dM = Malloc(double,ldbM)) == NULL) die("memory allocation!","dM");
        if ((bperp = Malloc(float,N)) == NULL) die("memory allocation!","bperp");
        if ((gfile = Malloc(char *,N)) == NULL) die("memory allocation!","gfile");
        if ((cfile = Malloc(char *,N)) == NULL) die("memory allocation!","cfile");
	for (i=0;i<N;i++) {
		if ((gfile[i] = Malloc(char,36)) == NULL) die("memory allocation!","gfile[i]");
		if ((cfile[i] = Malloc(char,36)) == NULL) die("memory allocation!","cfile[i]");
	}
        if ((L = Malloc(int,S)) == NULL) die("memory allocation!","L");
        if ((time = Malloc(double, S)) == NULL) die("memory allocation!","time");
	if ((H = Malloc(int,N*2)) == NULL) die("memory allocation!","H");
        if ((G = Malloc(double, m*n)) == NULL) die("memory allocation!","G");
        if ((A = Malloc(double, m*n)) == NULL) die("memory allocation!","A");
        if ((Gs = Malloc(double, N*n)) == NULL) die("memory allocation!","Gs");
        if ((gM = Malloc(double, S*nM)) == NULL) die("memory allocation!","gM");
        if ((grdin = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","grdin");
        if ((corin = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","corin");
        if ((flag = Malloc(int, xdim * ydim)) == NULL) die("memory allocation!","flag");
        if ((dem = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","dem");
        if ((res = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","res");
        if ((vel = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","vel");
        if ((amp = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","amp");
        if ((shift = Malloc(float, xdim * ydim)) == NULL) die("memory allocation!","shift");
        if ((phi = Malloc(float, N*xdim*ydim)) == NULL) die("memory allocation!","phi");
        if ((var = Malloc(float, N*xdim*ydim)) == NULL) die("memory allocation!","var");
        if ((disp = Malloc(float, S*xdim*ydim)) == NULL) die("memory allocation!","disp");

/* initialization */
        for (i=0;i<m*n;i++) G[i]=0;
        for (i=0;i<n*N;i++) Gs[i]=0;

        for (k=0;k<ydim;k++) {
                for (j=0;j<xdim;j++) {
                        flag[j+xdim*k]=0;
                        res[j+xdim*k]=0;
			dem[j+xdim*k]=0;
			for (p=0;p<S;p++) { 
				disp[p*xdim*ydim+j*ydim+k]=1e-6; /* to differ from rland */
			}
                }
        }

	printf("read table file ...\n");
        /* read in scene.tab */
        i=0;
        while(fscanf(datefile, "%s %s", &tmp1[0],&tmp2[0]) != EOF) {
                L[i]=atoi(tmp1);
		time[i]=atof(tmp2);
                i=i+1;
        }
        if (i != S) die("S and number of the SAR scenes don't match!","");
        fprintf(stderr,"number of SAR scenes is %d \n", S);
	for (i=S-1;i>=0;i--) time[i]=time[i]-time[0];

	/* read in intf.tab */
	i=0;
        while(fscanf(infile, "%s %s %s %s %s", gfile[i],cfile[i],&tmp1[0],&tmp2[0],&tmp3[0]) != EOF) {
		H[i*2+0]=atoi(tmp1);
		H[i*2+1]=atoi(tmp2);
		bperp[i]=atof(tmp3);
		i=i+1;
	}
        if (i != N) die("N and number of interferograms don't match!","");
        fprintf(stderr,"number of interferograms is %d \n", N);


        /* read in N 2-dimensional grd file into 3D array */
	printf("read phase and correlation grids ...\n");
        for (i=0;i<N;i++) {
                readgrdsize_(&xin,&yin,cfile[i]);
                if (xin != xdim || yin != ydim) die("dimension don't match!",cfile[i]);
		readgrdsize_(&xin,&yin,gfile[i]);
                if (xin != xdim || yin != ydim) die("dimension don't match!",gfile[i]);
                readgrd_(corin,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rdum,title,cfile[i]);
                readgrd_(grdin,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rdum,title,gfile[i]);
                for (k=0;k<ydim;k++) {
                        for (j=0;j<xdim;j++) {
                                phi[i*xdim*ydim+ydim*j+k] = grdin[j+k*xdim];
                                if (fabs(grdin[j+k*xdim] - rdum) < 1e-6) {
					flag[j+xdim*k] = 1;
				}
				if (corin[j+k*xdim]>=1e-2 && corin[j+k*xdim]<=0.99) {
				/* Rosen et al., 2000 IEEE */
					var[i*xdim*ydim+ydim*j+k] = sqrt((1.0-corin[j+k*xdim]*corin[j+k*xdim])/(corin[j+k*xdim]*corin[j+k*xdim])); 
				}
				else if (corin[j+k*xdim]<1e-2) {
					var[i*xdim*ydim+ydim*j+k] = 99.99;
				}
				else {
					var[i*xdim*ydim+ydim*j+k] = 0.1;
				}
                        }
                }
	}

	/* fill the G matrix */
        printf("fill the G matrix ...\n");
        for (i=0;i<N;i++) {
                for (j=0;j<S-1;j++) {
                        if (L[j]>=H[i*2+0] && L[j]<H[i*2+1]) {
                                G[i+m*j]=1;
                                Gs[i+N*j]=1;
                        }
                }
                G[i+m*(n-1)]=bperp[i]*scale;
                Gs[i+N*(n-1)]=bperp[i]*scale;
        }

	/* scale the smoothing by time */
        for (i=0;i<S-2;i++) {
                G[i+N+i*m]=sf/(time[i+1]-time[i]);
                G[i+N+(i+1)*m]=-sf/(time[i+2]-time[i+1]);
        }


        /* save G matrix to A as it will get destroyed*/
        for (i=0;i<m*n;i++) A[i]=G[i];


        /* loop over xdim by ydim pixel */
	printf("run least-squares problem over %d by %d pixel ...\n",xdim,ydim);
        for (k=0;k<ydim;k++) {
                for (j=0;j<xdim;j++) {
                        /* check the dummy value of the grd file */
                        if (flag[j+xdim*k] != 1) {
				for (i=0;i<m;i++) { 
					d[i]=0;
                                	if (i<N) { 
						d[i] = (double)phi[i*xdim*ydim+ydim*j+k]/var[i*xdim*ydim+ydim*j+k];
						ds[i] = d[i];
					}
				}

			        for (i=0;i<m;i++) { 
					for (p=0;p<n;p++) { 
						if (i<N) {
							G[i+m*p]=A[i+m*p]/var[i*xdim*ydim+ydim*j+k]; 
						}
						else {
							G[i+m*p]=A[i+m*p]; 
						}
					}
				}
  
                                dgelsy_(&m,&n,&nrhs,G,&lda,d,&ldb,jpvt,&rcond,&rank,work,&lwork,&info);

/*
				dgelss_(&m,&n,&nrhs,G,&lda,d,&ldb,sv,&rcond,&rank,work,&lwork,&info);
*/
				if (info != 0) printf("warning! input has an illegal value\n");
				if (j==0 && k==0) { 
					if (rank==n) { 
						printf("matrix is full rank: %d\n\n", rank);
					}
					else if (rank<n) { 
						printf("matrix is rank-deficient: %d\n\n", rank);
					}
				}
                                for (i=0;i<S;i++) { 
					for (p=0;p<i;p++) { 
						disp[i*xdim*ydim+j*ydim+k]=disp[i*xdim*ydim+j*ydim+k]+d[p];
					}
					disp[i*xdim*ydim+j*ydim+k]=-79.58*wl*disp[i*xdim*ydim+j*ydim+k];
				}

				if (flag_dem == 1) dem[j+xdim*k]=d[n-1];

				if (flag_rms == 1) {
					new=0;
					old=0;
				/* check the WRMS reduction */
                        		for (i=0;i<N;i++) {
                                		pred=0;
                                		for (p=0;p<n;p++) {
                                        		pred=pred+d[p]*Gs[i+p*N]/var[i*xdim*ydim+ydim*j+k];
                                		}
                                		new=new+(ds[i]-pred)*(ds[i]-pred);
                                		old=old+ds[i]*ds[i];
                        		}	
                        		res[j+xdim*k]=(sqrt(old)-sqrt(new))/sqrt(old);
				}

				/* fitting a straight line by least squares */
                                sumxy=0;
                                sumxx=0;
				sumx=0;
				sumy=0;
                                for (i=0;i<S;i++) {
                                       	sumxy=sumxy+time[i]*disp[i*xdim*ydim+j*ydim+k];
                                       	sumxx=sumxx+time[i]*time[i];
					sumy=sumy+disp[i*xdim*ydim+j*ydim+k];
					sumx=sumx+time[i];
                               	}
                                vel[j+xdim*k]=(S*sumxy-sumx*sumy)/(S*sumxx-sumx*sumx)*365.0;
                        }
                        else {
                                for (i=0;i<S;i++) disp[i*xdim*ydim+j*ydim+k] = rdum;
				vel[j+xdim*k] = rdum;
				if (flag_rms == 1) res[j+xdim*k] = rdum;
				if (flag_dem == 1) dem[j+xdim*k] = rdum;
                        }
                }
        }

	/* write output */
	printf("write output ...\n");
        strcpy(title,"displacement time series (mm)");
        for (i=0;i<S;i++) {
                strcpy(outfile,"disp_");
                for (k=0;k<ydim;k++) {
                        for (j=0;j<xdim;j++) {
                                grdin[j+k*xdim]=disp[i*xdim*ydim+j*ydim+k]; 
                      }
                }
                sprintf(tmp1,"%02d",i+1);
                strcat(outfile,tmp1);
                strcat(outfile,".grd");
                writegrd_(grdin,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rland,&rdum,title,outfile);
        }

	if (flag_rms == 1) {
	        strcpy(title,"WRMS reduction %");
		sprintf(outfile,"rms.grd");
		writegrd_(res,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rland,&rdum,title,outfile);
	}

	if (flag_dem == 1) {
		strcpy(title,"DEM error (m)");
        	sprintf(outfile,"dem.grd");
        	writegrd_(dem,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rland,&rdum,title,outfile);
	}

	strcpy(title,"mean velocity (mm/yr)");
	sprintf(outfile,"vel.grd");
	writegrd_(vel,&xin,&yin,&ymin,&xmin,&yinc,&xinc,&rland,&rdum,title,outfile);


        /* free memory */

        for (i=0;i<N;i++) {
		free(gfile[i]);
		free(cfile[i]);
        }
        free(phi);
        free(var);
	free(gfile);
	free(cfile);
        free(disp);
	free(G);
	free(A);
	free(Gs);
        free(H);
	free(d);
	free(ds);
	free(L);
	free(res);
	free(vel);
	free(amp);
	free(shift);
	free(time);
        free(flag);
        free(bperp);
        free(dem);
        free(grdin);
        free(corin);
        free(work);
        free(workM);
	free(jpvt);
	free(jpvtM);
	free(gM);
	free(dM);

	fclose(infile);
	fclose(datefile);

        return(0);
}

