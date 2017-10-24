/*  program to test the fft shift routine */
#include "../../include/soi.h"
#include "../../include/siocomplex.h"
#include <math.h>
void main() {

int k,nd=1024;
fcomplex *datai, *datao;
double arg;

datai  = (fcomplex *) malloc(nd*sizeof(fcomplex));
datao  = (fcomplex *) malloc(nd*sizeof(fcomplex));

/*  fill the data array with a gaussian function */

  for (k=0; k<nd; k++){
    datai[k].r=0.0;
    datai[k].i=0.0;
    arg=(k-nd/2)*(k-nd/2)/4.;
    datai[k].r=(k-nd/2)*exp(-arg);
    datai[k].i=-(k-nd/2)*exp(-arg);
    datao[k].r=datai[k].r;
    datao[k].i=datai[k].i;
  }
  shift(nd,datao,50.5);
  for (k=0; k<nd; k++){
    printf(" %d %lf %lf %lf %lf \n",k,datai[k].r,datai[k].i,datao[k].r,datao[k].i);
  }
}
