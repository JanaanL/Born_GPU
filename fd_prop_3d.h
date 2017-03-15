#ifndef FDPROP_3D_H
#define FDPROP_3D_H 1
#include "vel_fd_3d.h"
#include "source_func_3d.h"
#include "param_func.h"
class fd_prop{
  public:
   fd_prop(){};

    bool  calc_stability(std::vector<float> *ds,int n);
    ~fd_prop(){delete_fd_op();}
    void delete_fd_op(){};
    void set_vel(vel_fd_3d *v){ vel=v;}

    void set_verb(bool v){ verb=v;}
    void set_bounds(int n1, int n2, int n3){ nboundt=n1; nbound=n2, nbound_y=n3;}
    void set_fd_basics(param_func *par, source_func *source_func, float ap,bool v);
    void create_transfer_sinc_source(int nsinc);

    vel_fd_3d *vel;
    float dt;
    float vmax;
    int nt;
    int nz,nx,ny;
    int nboundt,nbound,nbound_y;
    source_func *source;
    float aper;
    int jts;
    float dtw;
    bool verb;
 


};
#endif
