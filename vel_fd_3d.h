#include  "sregf.h"
#include "hypercube_float.h"
#ifndef VEL_RTM_3D_H 
#define VEL_RTM_3D_H 1
class vel_fd_3d: public sregf{
  public:
    vel_fd_3d(){fld=0;dens=0; fprintf(stderr,"default constructor \n");}
    vel_fd_3d(char *tag);
    vel_fd_3d(char *tag,char *tg1);
    hypercube_float *rand_subspace(int irand, int nrandt, int nrandblr, int randblr2, float dt, hypercube_float *space);

    float get_min_samp(){
       float dx=fld->get_axis(1).d;
       float dy=fld->get_axis(2).d;
       float dz=fld->get_axis(3).d;
       if(dz < dx && dz < dy ) return dz;
       if(dy < dx ) return dz;
       return dx;
    }
        float get_max_samp(){
       float dx=fld->get_axis(1).d;
       float dy=fld->get_axis(2).d;
       float dz=fld->get_axis(3).d;
       if(dz > dx && dz > dy ) return dz;
       if(dy > dx) return dz;
       return dx;
    }
    void set_old(){ old=true;}
    void set_mid(){mid=true;}
    bool has_dens(){
      if(dens==0) return false;
      return true;
    }
    float max_vel(){
      return (float) fld->my_max();
    }
    float min_vel(){
      return (float) fld->my_min();
    }
    axis get_axis(int iax){
      return fld->get_axis(iax);
    }
    void set_zero(){
       zero=true;
    
    }
   ~vel_fd_3d(){
   fprintf(stderr,"in delete 1 \n");
       if(fld!=0) delete fld;
       fld=0;

   }
   
   hypercube_float *fld;
   hypercube_float *dens;
   sregf *dfile;
   bool old;
   bool mid,zero;




};
#endif
