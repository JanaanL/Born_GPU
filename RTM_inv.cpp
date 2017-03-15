#include<sep_params.h>
#include<sregf.h>
#include <seplib.h>
#include <vel_fd_3d.h>
#include <data_rtm_3d.h>
#include <image_rtm_3d.h>
#include "source_func_3d.h"
#include "rtm_zero_op_3d.h"
#include "gpu_funcs_3d.h"
#include "cgstep.h"
#include "lin_solver.h"
#include "tmute.h"

 main(int argc, char **argv){
 
   sep_params pars=sep_params(argc,argv);
   data_rtm_3d data=data_rtm_3d("data",&pars);
   source_func wavelet=source_func("wavelet"); 
   axis asx=data.get_axis(4);
   axis asy=data.get_axis(5); //added
   int niter=pars.get_int("niter",5);
   float src_depth=pars.get_float("src_depth",0.);
   wavelet.set_sources_3d(src_depth,asx,asy);
   vel_fd_3d vel=vel_fd_3d("velocity");
   image_rtm_3d image=image_rtm_3d("image",&vel);
  // image.read_all("itest.H",image.image);
   //image.write_all(image.tagit.c_str(),image.image);
   image.zero();

   image.set_source_file(&image);
   //image.write_volume();
 ////  image.write_description("image");
   oc_float *d2=new oc_float("dtest.H",&data);
   data.set_source_file(d2);

   int n_gpus=pars.get_int("n_gpus",0);
   setup_cuda(n_gpus,argc,argv);

   float aper=pars.get_float("aper",8.);
//   int *rand_vec=0;
   int *rand_vec=0;//new int[data.get_axis(4).n*data.get_axis(5).n];

   bool encode=false;
   rtm_zero_op op=rtm_zero_op(&pars,&vel,&wavelet,&data,&image,aper,true,encode,rand_vec,true,false);
   //fprintf(stderr,"CREATING TEMP data \n");
   //op.forward(false,&image,d2);
   //oc_float *i2=new oc_float("iout.H",image.image);
   //i2->zero();
   //fprintf(stderr,"CREATING TEMP model\n");
  // op.adjoint(false,i2,d2);
  // seperr("");
  // op.dot_test(true);
fprintf(stderr,"Iterations: %d \n",niter);

   cgstep *st=new cgstep();

   simple_solver *solv=new simple_solver(st,&data,&op);

   solv->set_verbose(3);
   solv->solve(niter);

   oc_float *mine=(oc_float*) solv->return_model();
   image.tagit="image";
   image.add(mine);

   return 0;

 }






