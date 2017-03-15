#ifndef IMAGE_RTM_3D_H
#define IMAGE_RTM_3D_H 1
#include "vel_fd_3d.h"   
#include "float_3d.h"
#include "oc_float.h"

class image_rtm_3d: public oc_float{
  public:
   image_rtm_3d();
   image_rtm_3d(char *tg,vel_fd_3d *v);
   
   ~image_rtm_3d(){
     if(image!=0) delete image;
   }
   void write_volume();
   void write_final_volume();
   void add_image(hypercube_float *img);
   float_3d *extract_sub(axis a1, axis a2, axis a3);

   void set_source_file(oc_float *ptr);
   hypercube_float *image, *mute;
   oc_float *myf;

};
#endif
