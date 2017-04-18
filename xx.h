#ifndef RTM_OP_3D_H
#define RTM_OP_3D_H 1
#include "float_1d.h"
#include "my_operator.h"
#include "data_rtm_3d.h"
#include "image_rtm_3d.h"
#include "fd_prop_3d.h"
#include "source_func_3d.h"
class rtm_zero_op: public my_operator, public fd_prop{
  public:
    rtm_zero_op(){basic_init_op();}
    rtm_zero_op(param_func *par,vel_fd_3d *vel, source_func *source_func, data_rtm_3d *data_insert, image_rtm_3d *image, float aper, bool verb, bool encode, int *rand_vec, bool do_src=true,bool redo_src=false);
      void create_source_fields();
      void migrate(image_rtm_3d *m,data_rtm_3d *d);

    virtual bool forward(bool add, my_vector *model, my_vector *data, int iter);
    virtual bool adjoint(bool add, my_vector *model, my_vector *data, int iter);
    ~rtm_zero_op(){delete_rtm_op();}
    void delete_rtm_op(){};
    void create_transfer_sinc_data( int nsinc);
    void create_random_trace(int nshots, int *rand_vec);
    bool encode;

  private:
    data_rtm_3d *data;
    image_rtm_3d *image;
    int *rand_vec;
   
    int jtd;
    float dtd;
    bool redo;
    int base;
    float *slice_p0;


};
#endif
