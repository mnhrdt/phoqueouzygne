/************************************************************************
* rng_cmp performs the range compression operation on raw radar echos.  *
*	using a precomputed range reference function.			*
************************************************************************/
/************************************************************************
* Creator: Evelyn J. Price	(Scripps Institution of Oceanography)	*
* Date   : 11/18/96							*
************************************************************************/
/************************************************************************
* Modification History							*
* 									*
* Date									*
************************************************************************/ 

#include "../../include/soi.h"
#include "../../include/siocomplex.h"
#include <math.h>

rng_cmp(ranfft,data,ref)
int ranfft;
fcomplex *data, *ref;
{
	int i, dir;
        dir = -1;
	cfft1d_(&ranfft,data,&dir);
	for(i=0;i<ranfft;i++){
	  data[i] = Cmul(ref[i],data[i]);
	}
        dir = 1;
	cfft1d_(&ranfft,data,&dir);
}
