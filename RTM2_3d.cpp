#include<sep_params.h>
#include<sregf.h>
#include <seplib.h>
#include <vel_fd_3d.h>
#include <data_rtm_3d.h>
#include <image_rtm_3d.h>
#include "source_func_3d.h"
#include "rtm_zero_op_3d.h"
#include "gpu_funcs_3d.h"

 main(int argc, char **argv){
 
   sep_params pars=sep_params(argc,argv);
   data_rtm_3d data=data_rtm_3d("data",&pars);
   source_func wavelet=source_func("wavelet"); 

   float src_depth=pars.get_float("src_depth",0.);
   vel_fd_3d vel=vel_fd_3d("velocity");
   image_rtm_3d image=image_rtm_3d("image",&vel);

   axis asx=data.get_axis(4);
   axis asy=data.get_axis(5);
   wavelet.set_sources_3d(src_depth,asx,asy);

   image.zero();
   image.set_source_file(&image);

   oc_float *d2=new oc_float("dtest.H",&data);
   data.set_source_file(d2);

   int n_gpus=pars.get_int("n_gpus",1);
   setup_cuda(n_gpus,argc,argv);

   float aper=pars.get_float("aper",8.);

   int *rand_vec=0;
   bool encode=false;

   rtm_zero_op op=rtm_zero_op(&pars,&vel,&wavelet,&data,&image,aper,true,encode,rand_vec,true);


   op.adjoint(false,&image,&data,1);

   return 0;
 }
       






