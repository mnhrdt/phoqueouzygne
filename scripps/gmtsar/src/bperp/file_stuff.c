#include "gmtsar.h"
#include "lib_functions.h"

/*--------------------------------------------------------------*/
FILE *create_GMT_binary_float(struct PRM p1, char *filename, char *string)
{
int	nx, ny;
double	x_inc, y_inc;
FILE	*datafile;

	if (verbose) fprintf(stderr,"create_GMT_binary_float\n");
	/* attempt to open data and hdr files */
	if ((datafile = fopen(filename,"w")) == NULL) die("cannot open ",filename);

	/* write out information to header file */
	nx = p1.num_rng_bins;
	ny  = p1.num_patches * p1.num_valid_az;

	x_inc = 1.0;
	y_inc = 1.0;

	create_GMT_binary_hdr(nx, ny, x_inc, y_inc, "phasediff", string, datafile);

	return(datafile);
}

/*--------------------------------------------------------------*/
void print_prm_params(struct PRM p1, struct PRM p2)
{
	fprintf(stderr," SLC 1: num_rng_bins %d num_lines %d \n",p1.num_rng_bins, p1.num_lines);
	fprintf(stderr," SLC 2: num_rng_bins %d num_lines %d \n",p2.num_rng_bins, p2.num_lines);
	fprintf(stderr," lambda %f \n",p2.lambda);
	fprintf(stderr," baseline_start %f \n",p2.baseline_start);
	fprintf(stderr," B_end %f \n",p2.baseline_end);
	fprintf(stderr," alpha_start %f \n",p2.alpha_start);
	fprintf(stderr," alpha_end %f \n",p2.alpha_end);
	fprintf(stderr," near_range %f \n",p2.near_range);
	fprintf(stderr," rng_samp_rate %f \n",p2.fs);
	fprintf(stderr," sc_clock_start %f \n",p2.SC_clock_start);
	fprintf(stderr," sc_clock_stop %f \n",p2.SC_clock_stop);
	fprintf(stderr," prf %f \n",p2.prf);
}

/*--------------------------------------------------------------*/
void fix_prm_params(struct PRM *p, char *s)
{
double delr;

	delr = SOL/p->fs/2.0;

	/* these are from prm2gips */
	p->near_range = p->near_range + (p->st_rng_bin - p->chirp_ext + p->rshift-1)*delr;
	p->SC_clock_start = p->SC_clock_start + p->ashift/(p->prf*86400.0) + (p->nrows-p->num_valid_az)/(2.0*p->prf*86400);
	p->SC_clock_stop  = p->SC_clock_start + (p->num_valid_az*p->num_patches)/(p->prf*86400.0);

}
/*--------------------------------------------------------------*/
/* xdim, ydim are size of SLC files 				*/
/* xdimd, ydimd are size of topo (or model) file 		*/
/* xdimd and ydimd must be integer factor of xdim, ydim		*/
/* xdim/xdimd = xdec						*/
/* ydim/ydimd = ydec						*/
/* xdim 	SLC file x dimension				*/
/* ydim		SLC file y dimension				*/
float *read_aux_binary_grd(char *input_file_name, int xdim, int ydim, int *xdimd, int *ydimd, double *xdec, double *ydec)
{
float	*ptr;
long 	nitems;
FILE	*fin;

	fin = fopen(input_file_name,"r");
	if (fin == NULL) die("error opening ",input_file_name);

	read_GMT_binary_dimensions(fin, xdimd, ydimd);
	nitems = (*xdimd) * (*ydimd);

	if (verbose) fprintf(stderr," done reading file %s (%d by %d)\n", input_file_name, *xdimd, *ydimd); 

	if ((ptr = (float *) malloc((*xdimd) * (*ydimd) * sizeof(float))) == NULL) 
					die("problem allocating memory",":topo");

	/* check dimensions */
	*ydec = floor((double) (ydim)/(*ydimd));
	*xdec = floor((double) (xdim)/(*xdimd));

	if (verbose) fprintf(stderr," file %s ( %d by %d => decimation %6.2lf by %6.2lf)\n"
		,input_file_name, *xdimd, *ydimd, *xdec, *ydec);

	if(((*xdec)*(*xdimd) != xdim) | ((*ydec)*(*ydimd) != ydim)){ 
		fprintf(stderr," SLC file %d by %d ; topo file %d by %d \n", xdim, ydim, *xdimd, *ydimd);
		die("dimensions must be integral to the SLC dimensions",input_file_name); 
		}

	fread(ptr, sizeof(float), nitems, fin);

	fclose(fin);

return(ptr);
}
/*--------------------------------------------------------------*/
