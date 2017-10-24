#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gmtsar.h"
#include "xcorr.h"

/*-------------------------------------------------------*/
/* determine pixel location of each point for testing 	*/
/* use master						*/
/* pixel locations for slave includes x, y offset 	*/
void get_locations(struct xcorr *xc)
{
	int	i, j, n;

/* need to space more than offset from side */
/* also add two more range locations then below we don't use the first and last range locations */
/* these first and last locations are unreliable because the data have been extended */
	xc->x_inc = (xc->m_nx - 2*(xc->xsearch + xc->nx_corr))/(xc->nxl + 3);
	xc->y_inc = (xc->m_ny - 2*(xc->ysearch + xc->ny_corr))/(xc->nyl + 1);

/* add an extra column on memory just in case */
	xc->loc = malloc(xc->nyl*(xc->nxl+1)*sizeof(struct locs));

	n = 0;

	for (j=1;j<=xc->nyl; j++){
/* remove the first and last of each range line to avoid regions of low SAR amplitude and unreliable correlation*/
		for (i=2;i<=xc->nxl+1; i++){

			xc->loc[n].x = xc->npx + i*xc->x_inc;
			xc->loc[n].y = xc->npy + j*xc->y_inc;

			n++;
			}
		}

	xc->nlocs = n;

	fprintf(stderr," locations  n %d nx %d nyl %d nxl %d x_inc %d y_inc %d\n", n, xc->m_nx, xc->nyl, xc->nxl, xc->x_inc, xc->y_inc);
}
/*-------------------------------------------------------*/
