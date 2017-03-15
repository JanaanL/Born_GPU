#ifndef DERIV_3D_H
#define DERIV_3D_H 1
#include "my_operator.h"
#include "source_func_3d.h"
class deriv: public my_operator{
  public:
    deriv(hypercube_float *mod, hypercube_float *dat){
       set_domain(mod); set_range(dat);
       };

    ~deriv(){};
    virtual bool forward(bool add, my_vector *model, my_vector *data,int iter=0);
    virtual bool adjoint(bool add, my_vector *model, my_vector *data,int iter=0);

 

};
#endif
