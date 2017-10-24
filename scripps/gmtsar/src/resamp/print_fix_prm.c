#include "gmtsar.h"
#include "lib_functions.h"

/*--------------------------------------------------------------*/
void print_prm_params(struct PRM p1, struct PRM p2)
{
	fprintf(stderr," SLC 1: num_rng_bins %d num_lines %d \n",p1.num_rng_bins, p1.num_lines);
	fprintf(stderr," SLC 2: num_rng_bins %d num_lines %d \n",p2.num_rng_bins, p2.num_lines);
	fprintf(stderr," lambda %f \n",p2.lambda);
	fprintf(stderr," near_range %f \n",p2.near_range);
	fprintf(stderr," rng_samp_rate %f \n",p2.fs);
	fprintf(stderr," sc_clock_start %f \n",p2.SC_clock_start);
	fprintf(stderr," sc_clock_stop %f \n",p2.SC_clock_stop);
	fprintf(stderr," prfm %f \n",p1.prf);
	fprintf(stderr," prfs %f \n",p2.prf);
	fprintf(stderr," rshift %f \n",p2.rshift+p2.sub_int_r);
	fprintf(stderr," ashift %f \n",p2.ashift+p2.sub_int_a);
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
