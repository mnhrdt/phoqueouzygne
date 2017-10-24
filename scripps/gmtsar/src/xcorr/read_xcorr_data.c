#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmtsar.h"
#include "xcorr.h"

void read_complex_short(FILE *, struct FCOMPLEX *, int, int, int, short *);
void read_real_float(FILE *, struct FCOMPLEX *, int, int, int, float *);

/*-------------------------------------------------------*/
void read_xcorr_data(struct xcorr xc, int iloc)
{
int	iy, ishft;
short	*tmp_m,*tmp_s;
float	*tmp2_m,*tmp2_s;

	tmp_m = (short *) malloc(2*xc.m_nx*sizeof(short));		/* whole line */
	tmp2_m = (float *) malloc(2*xc.m_nx*sizeof(float));		/* whole line */
	tmp_s = (short *) malloc(2*xc.s_nx*sizeof(short));		/* whole line */
	tmp2_s = (float *) malloc(2*xc.s_nx*sizeof(float));		/* whole line */

	/* set locations and read data for master 	*/
	/* read whole line at correct y offset		*/
	iy = xc.loc[iloc].y - xc.npy/2;
	
	if (debug) fprintf(stderr," reading data from master at y = %d and %d items\n", iy, xc.m_nx);

	if (xc.format == 0) read_complex_short(xc.data1, xc.d1, iy, xc.npy, xc.m_nx, tmp_m);
	if (xc.format == 1) read_real_float(xc.data1, xc.d1, iy, xc.npy, xc.m_nx, tmp2_m);


	/* set locations and read data for slave */
	ishft = (int) xc.loc[iloc].y * xc.astretcha;
	iy = xc.loc[iloc].y + xc.y_offset + ishft - xc.npy/2;

	if (debug) fprintf(stderr," reading data from slave at y = %d and %d items\n", iy, xc.s_nx);
	if (xc.format == 0) read_complex_short(xc.data2, xc.d2, iy, xc.npy, xc.s_nx, tmp_s);
	if (xc.format == 1) read_real_float(xc.data2, xc.d2, iy, xc.npy, xc.s_nx, tmp2_s);

	free ((char *) tmp_m);
	free ((char *) tmp2_m);
	free ((char *) tmp_s);
	free ((char *) tmp2_s);
}
/*-------------------------------------------------------*/
void read_complex_short(FILE *f, struct FCOMPLEX *d, int iy, int npy, int nx, short *tmp)
{
	long num_to_seek;
	int	i, j;
	
	num_to_seek = 2*iy*nx*sizeof(short); 
	fseek(f, num_to_seek, SEEK_SET);	/* from beginning */

	/* need to read two parts of complex numbers */
	for (i=0; i<npy; i++){
		fread(&tmp[0],2*sizeof(short), nx, f);	/* read whole line */

		/* read into complex float */
		for (j=0; j<nx; j++) {
			d[i*nx+j].r = (float) tmp[2*j];
			d[i*nx+j].i = (float) tmp[2*j + 1];
			}
		}
}
/*-------------------------------------------------------*/
void read_real_float(FILE *f, struct FCOMPLEX *d, int iy, int npy, int nx, float *tmp)
{
        long num_to_seek;
        int     i, j;

        num_to_seek = iy*nx*sizeof(float);
        fseek(f, num_to_seek, SEEK_SET);        /* from beginning */

        /* need to read two parts of complex numbers */
        for (i=0; i<npy; i++){
                fread(&tmp[0],sizeof(float), nx, f);  /* read whole line */

                /* read into complex float */
                for (j=0; j<nx; j++) {
                        d[i*nx+j].r = (float) tmp[j];
                        d[i*nx+j].i = 0.0;
                        }
                }
}
