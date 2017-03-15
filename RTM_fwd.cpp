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
   int wavefield=pars.get_int("wavefield_source",0);
   fprintf(stderr,"WHAT IS WAVEFIELD %d \n",wavefield);
   
   source_func *wavelet;
      float src_depth=pars.get_float("src_depth",0.);

   if(wavefield==1){
     fprintf(stderr,"initiailizing a wavefield\n");
     wavefield_source_func  *ws=new wavefield_source_func("wavefield"); 

     ws->set_sources_depth(src_depth);

     wavelet=ws;

   }
   else{
     wavelet_source_func  *ws=new wavelet_source_func("wavelet"); 
          fprintf(stderr,"initiailizing a wavelet\n");

     axis asx=data.get_axis(4);
     axis asy=data.get_axis(5); //added
     ws->set_sources_axes(src_depth,asx,asy);
     wavelet=ws;
   }




   vel_fd_3d vel=vel_fd_3d("velocity");

   image_rtm_3d image=image_rtm_3d("image",&vel);

//   image.read_all("image",image.image);
   //image.write_all(image.tagit.c_str(),image.image);
//   image.zero();
   image.set_source_file(&image);
   //image.write_volume();
//   image.write_description("image");

   oc_float *d2=new oc_float("dataout",&data);
   data.set_source_file(d2);

   int n_gpus=pars.get_int("n_gpus",1);
   setup_cuda(n_gpus,argc,argv);
   float aper=pars.get_float("aper",8.);
   int *rand_vec=0;
   bool encode=false;

   rtm_zero_op op=rtm_zero_op(&pars,&vel,wavelet,&data,&image,aper,true,encode,rand_vec,false);
//   rtm_zero_op op=rtm_zero_op(&pars,&vel,&wavelet,&data,&image,aper,true/*,true,false*/);
   //fprintf(stderr,"CREATING TEMP data \n");
   op.forward(false,&image,d2,1);
   //oc_float *i2=new oc_float("iout.H",image.image);
   //i2->zero();
   //fprintf(stderr,"CREATING TEMP model\n");
  // op.adjoint(false,i2,d2);
  // seperr("");
//   op.dot_test(true);
//   seperr("");
//   op.adjoint(false,&image,&data);

//  cgstep *st=new cgstep();
//  oc_float *wt_space=data.clone();
//  tmute tm=tmute(pars.get_float("tmute",.3),pars.get_float("vmute",1500.),wt_space,&data);
//  simple_solver *solv=new simple_solver(st,&data,&op);

//  solv->add_equation(dat,mat);
  
//  solv->set_verbose(3);
//  solv->solve(pars.get_int("niter"));
  
//  oc_float *mine=(oc_float*) solv->return_model();
//  image.tagit="image";
//  image.add(mine);
   delete wavelet;

 //  image.write_volume();
   return 0;
 }
       






