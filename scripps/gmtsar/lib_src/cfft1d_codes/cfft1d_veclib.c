/************************************************************************
* cfft1d is a subroutine used to call and initialize veclib FFT         *
* in a Mac OS X computer.                                               *
************************************************************************/
/************************************************************************
* Creator: Robert Kern          (Scripps Institution of Oceanography    *
* Date   : 12/2005                                                      *
************************************************************************/
/************************************************************************
* Modification history:                                                 *
*                                                                       *
* DATE                                                                  *
************************************************************************/
#include <Accelerate/Accelerate.h>
/*#include <vecLib/vecLib.h>*/
#include <stdio.h>
#include <stdlib.h>

void cfft1d_(int* np, DSPComplex* c, int* dir);
void cfft1d_cleanup_();

static int n = 0; 
static int log2n;
static FFTSetup setup;
static DSPSplitComplex d;
static float scale;
static int inited = 0;

void cfft1d_(int* np, DSPComplex* c, int* dir)
{
    if (*dir == 0) return;
    if (n != *np) {
        cfft1d_cleanup_();
        n = *np;
        for (log2n=1; (1<<log2n)<*np; log2n++) {}

        d.realp = (float*)malloc(n*sizeof(float));
        d.imagp = (float*)malloc(n*sizeof(float));

        setup = vDSP_create_fftsetup(log2n, 0);
        scale = 1.0/n;
        inited = 1;
    }
    
    vDSP_ctoz(c, 2, &d, 1, n);

    vDSP_fft_zip(setup, &d, 1, log2n, (*dir==-1 ? FFT_FORWARD : FFT_INVERSE));

    vDSP_ztoc(&d, 1, c, 2, n);

    if (*dir == 1) {
        vDSP_vsmul((float*)c, 1, &scale, (float*)c, 1, 2*n);
    }
}

void cfft1d_cleanup_()
{
    if (inited) {
        free(d.realp);
        free(d.imagp);
        vDSP_destroy_fftsetup(setup);
        inited = 0;
    }
}
