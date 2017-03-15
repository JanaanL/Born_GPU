#ifndef SINC_BOB_H
#define SINC_BOB_H 1
#include "stdio.h"

class sinc_bob{

  public:
    sinc_bob(){ table=0;}
    sinc_bob(int nt, int nsinc);
    void mksinc (float *sinc,int lsinc,float d,float *space);
    void toep (int m,float *r,float *f,float *g,float *a);
    ~sinc_bob(){ delete_sinc();}
    void delete_sinc();
    float *return_tab(float v){

      int itab=(int)(v/dtab+.5);
      if(itab==ntab) itab=ntab-1;
       return table[itab];
    }
    float **table;
    int ntab,nsinc;
    float dtab;
    
};
#endif
