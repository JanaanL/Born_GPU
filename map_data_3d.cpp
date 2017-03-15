#include "map_data_3d.h"
#include "math.h"

map_data::map_data(int npt,float *s_x,float *s_y,float *s_z, int *locs, hypercube_float *model,hypercube_float *dom, hypercube_float *ran,int ntbig){

  scale=new float[npt];
  map=locs;
  int h=0;
  ntb=ntbig;
  int ia=0;  
  int n1=model->get_axis(1).n;
  int n2=model->get_axis(2).n;
  

  for(int i4=0; i4 < npt; i4++){

    int ix=(int)s_x[i4];
    int iz=(int)s_z[i4];
    int iy=(int)s_y[i4];

    float sm=0;
    float yn=s_y[i4]-iy;
    float xn=s_x[i4]-ix;

    int ia=i4;
    float zn=s_z[i4]-iz;
    zn=0.;
    scale[ia]=expf(-zn*zn-xn*xn-yn*yn);
    map[ia]=ix+iy*n1+iz*n1*n2;

  }
  fprintf(stderr,"CHECK %d %f %f %f\n",map[npt-1],s_z[npt-1],s_x[npt-1],s_y[npt-1]);
  set_domain(dom);
  set_range(ran);
}


bool map_data::adjoint(bool add, my_vector *model, my_vector *data,int iter){

  float_3d *d=(float_3d*) data;
  float_3d *m=(float_3d*) model;
  if(!add) m->zero();

  int nt=d->get_axis(1).n;
  int n2=d->get_axis(2).n;
  int n3=d->get_axis(3).n;

  for(int i3=0; i3 < n3; i3++){
    for(int i2=0; i2 < n2; i2++){
      for(int it=0; it < nt; it++){
        m->vals[ntb*i2+it+4+i3*ntb*n2]+=d->vals[it+i2*nt+i3*nt*n2];
      }
    }
  }

}

bool map_data::forward(bool add, my_vector *model, my_vector *data,int iter){

  float_3d *d=(float_3d*) data;
  float_3d *m=(float_3d*) model;
  if(!add) d->zero();

  int nt=d->get_axis(1).n;
  int n2=d->get_axis(2).n;
  int n3=d->get_axis(3).n;

  for(int i3=0; i3 < n3; i3++){
    for(int i2=0; i2 < n2; i2++){
      for(int it=0; it < nt; it++){
        d->vals[it+i2*nt+i3*nt*n2]+=m->vals[ntb*i2+it+4+i3*ntb*n2];
      }
    }
  }

}

