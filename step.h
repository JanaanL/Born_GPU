#ifndef STEP_H
#define STEP_H 1
#include<my_vector.h>

class step{
  public:
    step(){};
    virtual void alloc_step(my_vector *mod, my_vector *dat){}
    virtual bool step_it(int iter, my_vector *x, my_vector *g, my_vector *rr, my_vector *gg,double *sc){return false;}
};

#endif

