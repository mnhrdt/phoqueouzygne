/************************************************************************
* ram2ras maps range and azimuth from a master image location into the  *
* corresponding range and azimuth location of the slave image.          *
************************************************************************/
/************************************************************************
* Creator: David Sandwell       (Scripps Institution of Oceanography)   *
* Date   : 03/22/13                                                     *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*                                                                       *
* DATE                                                                  *
************************************************************************/
#include "gmtsar.h"
#include <stdio.h>
#include <stdlib.h>

void ram2ras(struct PRM ps, double *ram, double *ras)

{

double dr;

	dr = SOL/ps.fs/2.0;

	/* this is the range coordinate */
	ras[0] = ram[0] + ((ps.rshift+ps.sub_int_r)+ram[0]*ps.stretch_r+ram[1]*ps.a_stretch_r);

	/* this is the azimuth coordinate */
	ras[1] = ram[1] + ((ps.ashift+ps.sub_int_a)+ram[0]*ps.stretch_a+ram[1]*ps.a_stretch_a);
}

