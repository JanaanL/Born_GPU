#ifndef WAVELET_RTM_3D_H
#define WAVELET_RTM_3D_H 1
#include<sregf.h>
#include<hypercube_float.h>
#include "wavefield_insert_3d.h"
class source_func: public wavefield_insert_3d{
  public:
    source_func(){};
    source_func(char *tag);
    


    void set_source_file(char *tg){
      tag=tg;
    }
    ~source_func(){ }
    hypercube_float *create_domain(int ishot);
    virtual int get_points(bool e){if(e); return 0;}
    void set_compute_size(hypercube *dom, float aper,int nbt,int nb,int nby, int fat, int blocksize);
    int y_points(){ return ay.n;}
    int x_points(){ return ax.n;}
    int z_points(){ return az.n;}
    virtual void get_source_func(hypercube_float *domain, int ishot,int nts, int *ilocs, float *vals){
      if(domain==0); if(ishot==0); if(nts==0); if(ilocs==0); if(vals==0);
    }
    float get_dt(){return dt;}
    float dt;
    std::vector<float> sx,sz,sy;
    int nx,nz,ny,jt,jts,jtd;
    float aper;
    axis az,ax,ay;
    int nboundt, nbound, nbound_y;
    std::string tag;
};

class wavelet_source_func: public source_func{
  public:
    wavelet_source_func(){};
    wavelet_source_func(char *tag);
    void set_sz(float s_z){ sz.push_back(s_z);}
    virtual int get_points(bool e);
    virtual void get_source_func(hypercube_float *domain, int ishot,int nts, int *ilocs, float *vals);
    ~wavelet_source_func(){delete wavelet;}
    void set_sources_axes(float s_z, axis src_axis1, axis src_axis2);

    hypercube_float *wavelet;

    
};
class wavefield_source_func: public source_func{
  public:
    wavefield_source_func(){};
    wavefield_source_func(char *tag);
    virtual int get_points(bool e);
    void set_sources_depth(float s_z){
    sz.push_back(s_z);
    fprintf(stderr,"JUST SET %f  \n",s_z);
    }
    virtual void get_source_func(hypercube_float *domain, int ishot,int nts, int *ilocs, float *vals);
 ~wavefield_source_func(){delete wavefield;}
    hypercube_float *wavefield;
    
};

#endif
