#ifndef tmute_H
#define tmute_H 1
#include "my_operator.h"
#include "oc_float.h"
#include "source_func_3d.h"
class tmute: public my_operator{
  public:
    tmute(float t, float v, oc_float *mod, oc_float *dat){
       set_domain(mod); set_range(dat);
       t_0=t; vmute=v;
       };

    ~tmute(){};
    virtual bool forward(bool add, my_vector *model, my_vector *data);
    virtual bool adjoint(bool add, my_vector *model, my_vector *data);

   float t_0,vmute;
  
};
#endif
