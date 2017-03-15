#ifndef CGSTEP_H
#define CGSTEP_H 1
#include<my_vector.h>
#include<step.h>
#include<stdio.h>
class cgstep: public step{
  public:
    cgstep(){ s=0; ss=0; };
    virtual bool step_it(int iter, my_vector *x, my_vector *g, my_vector *rr, my_vector *gg, double *sc);
    virtual void alloc_step(my_vector *mod,my_vector *dat);
    void forget_it(){ 
      forget=true;
      s->zero(); ss->zero();  
      small=1e-20;
    }
 
    ~cgstep(){
        if(s!=0) delete s;
        if(ss!=0) delete ss;
        s=0;
        ss=0;
     }

  private:
    my_vector *s, *ss;
    bool forget;
    double small;
  };


#endif
