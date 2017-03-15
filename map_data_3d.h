#ifndef MAP_DATA_3D_H
#define MAP_DATA_3D_H 1
#include "my_operator.h"
#include "source_func_3d.h"
#include "float_3d.h"
class map_data: public my_operator{
  public:
    map_data(){};
    map_data(int npt,float *s_x,float *s_y, float *s_z, int *locs,hypercube_float *model,
      hypercube_float *dom, hypercube_float *ran,int ntbig);

    virtual bool forward(bool add, my_vector *model, my_vector *data,int iter=0);
    virtual bool adjoint(bool add, my_vector *model, my_vector *data,int iter=0);
    ~map_data(){delete [] scale;}


    float *scale;
    int *map;
    int ntb;

};
#endif
