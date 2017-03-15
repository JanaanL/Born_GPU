#include "rtm_zero_op_3d.h"
#include "sinc_bob.h"
#include "gpu_funcs_3d.h"
#include <math.h>
#include "deriv_3d.h"
#include "laplac_3d.h"
#include "map_data_3d.h"
rtm_zero_op::rtm_zero_op(param_func *par,vel_fd_3d *vel_3d, source_func *source_func, data_rtm_3d *dat, image_rtm_3d *img,float ap,bool v,bool enc,int *r_vec, bool do_src,bool redo_src){
   std::vector<float> ds; 
   data=dat;
   ds.push_back(dat->dt);

   ds.push_back(source_func->get_dt());

   set_vel(vel_3d);

   calc_stability(&ds,data->nt);
   set_fd_basics(par,source_func,ap,v);
   set_ntblock(50);
   image=img;

   set_domain(image);
   set_range(data);
   dtd=data->dt;
   jtd=(int)(dtd/dt);
   nt=(data->nt-1)*jtd+1;
   base=0;
   redo=redo_src;
   create_transfer_sinc_data(8);

   encode=enc;
   rand_vec=r_vec;

   if(do_src && !redo) create_source_fields();
}

void rtm_zero_op::create_transfer_sinc_data(int nsinc){
  jtd=(int)((data->dt)/dt);
  sinc_bob myr(jtd,nsinc);
  transfer_sinc_table_d(nsinc,jtd,myr.table);
}

void rtm_zero_op::create_source_fields(){

  int nshots=data->nshots();
  auxinout("source_fields");

  int npts,nt_big;
  nt_big=9+data->nt; 
  npts=source->get_points(encode);
  int *locs=new int[npts];
  float *vals=new float[npts*(nt_big)];

  int seed=0;
  //nshots=1;

  for(int ishot=0; ishot< nshots; ishot++){

     if(verb) fprintf(stderr,"Forward propagating shot %d of %d \n",ishot,nshots);
     hypercube_float *src_p0=source->create_domain(ishot);
     hypercube_float *src_p1=src_p0->clone();

     source->get_source_func(src_p0,ishot,nt_big,locs,vals);

     transfer_source_func(npts,nt_big,locs,vals);
     seed=locs[3] - ((int) locs[3]/10000)*10000;

     hypercube_float *vloc=vel->rand_subspace(seed,nbound,nbound,nbound,dt,src_p0);

     transfer_vel_func1(nx,ny,nz,vloc->vals);

     time_t startwtimes, endwtimes;
     startwtimes = time(&startwtimes);

     source_prop(nx,ny,nz,false,true,src_p0->vals,src_p1->vals,jts,npts,nt);

     endwtimes = time(&endwtimes);

     fprintf(stderr,"Finished propagating wavefield, in %f \n",difftime(endwtimes, startwtimes));
     srite("source_fields",src_p0->vals,nz*nx*ny*4);
     srite("source_fields",src_p1->vals,nz*nx*ny*4);


     delete vloc;
     delete src_p0;
     delete src_p1;
   }
   delete [] locs;
   delete [] vals;
}

bool rtm_zero_op::forward(bool add, my_vector *mvec, my_vector *dvec, int iter){


  data_rtm_3d *d=(data_rtm_3d*) dvec;
  image_rtm_3d *m=(image_rtm_3d*) mvec;

  if(!add) dvec->zero();

  int npts_s=source->get_points(encode);

  data->set_source_file(d);
  image->set_source_file(m);

  int nt_big=8+data->nt;

  int npts=data->get_points();
  int rec_nx=data->get_points1();
  int rec_ny=data->get_points2();

  int *locs=new int[npts];
  int *locs_s=new int[npts_s];
  float *vals=new float[npts_s*(nt_big)];

  float *s_z=new float[npts],*s_x=new float[npts],*s_y=new float[npts];

  axis ax_pt1(rec_nx);
  axis ax_pt2(rec_ny);

  float_3d *shot=new float_3d(data->get_axis(1),data->get_axis(2),data->get_axis(3));
  float_3d *shot_deriv=new float_3d(data->get_axis(1),data->get_axis(2),data->get_axis(3));

  axis ax_t=axis(nt_big);
  float_3d *rec_func=new float_3d(ax_t,ax_pt1,ax_pt2);
  int nshots=data->nshots();
  //nshots=1;
  for(int ishot=0; ishot<nshots; ishot++){

    if(verb) fprintf(stderr,"Forward shot %d of %d\n",ishot,data->nshots());
    hypercube_float *src_p0=source->create_domain(ishot);

    hypercube_float *src_p1=src_p0->clone();
 
   axis ad1=src_p0->get_axis(1);
    axis ad2=src_p0->get_axis(2);
    axis ad3=src_p0->get_axis(3);

    source->get_source_func(src_p0,ishot,nt_big,locs_s,vals);

    transfer_source_func(npts_s,nt_big,locs_s,vals);

    int seed=locs_s[3] - ((int) locs_s[3]/10000)*10000;

    hypercube_float *vrand=vel->rand_subspace(seed,nbound,nbound,nbound,dt,src_p0);
    hypercube_float *vnone=vel->rand_subspace(seed,0,0,0,dt,src_p0);

    hypercube_float *img=image->extract_sub(ad1,ad2,ad3);

    transfer_vel_func1(nx,ny,nz,vrand->vals);
    transfer_vel_func2(nx,ny,nz,vnone->vals);

    data->get_source_func(src_p0,ishot,s_x,s_y,s_z,8,nt_big,shot);

    map_data mapit(npts,s_x,s_y,s_z,locs,src_p0,rec_func,shot_deriv,nt_big);

    rec_func->zero();

    transfer_receiver_func(rec_nx,rec_ny,nt_big,locs,rec_func->vals/*,s_z,s_x,s_y*/);

    hypercube_float *img_lap=img->clone();

    laplac lap=laplac(img_lap,img);

    lap.forward(false,img_lap,img);

    //lap.dot_test(true);

    rtm_forward(ad1.n,ad2.n,ad3.n,jtd,img->vals,rec_func->vals,npts_s,nt,nt_big,rec_nx,rec_ny);
    //rtm_forward(ad1.n,ad2.n,ad3.n,jtd,img->vals,rec_func->vals,npts_s,nt,nt_big,rec_nx,rec_ny);

    deriv der=deriv(shot_deriv,shot);

    mapit.forward(false,rec_func,shot_deriv);

    der.forward(false,shot_deriv,shot);


   srite("buf1.H",rec_func->vals,rec_func->get_n123()*4);
      srite("buf2.H",shot_deriv->vals,shot_deriv->get_n123()*4);
   srite("buf3.H",shot->vals,shot->get_n123()*4);

    data->add_data(ishot,shot);
fprintf(stderr,"SDS 1 \n");
    delete src_p0; delete src_p1; delete img; delete img_lap;
    delete vrand; delete vnone;
fprintf(stderr,"SDS 21 \n");

  }
  delete [] s_z;
  delete [] s_x;
  delete [] s_y;
fprintf(stderr,"SDS 31 \n");

  delete [] locs_s;
  delete [] vals;
fprintf(stderr,"SDSa 41 \n");

  delete [] locs;
  fprintf(stderr,"SDS b41 \n");

  delete rec_func;
fprintf(stderr,"SDS 41c \n");

  delete shot;
fprintf(stderr,"SDS 41 \n");

  delete shot_deriv;
fprintf(stderr,"SDS 15 \n");

  return true;

}

bool rtm_zero_op::adjoint(bool add, my_vector *mvec, my_vector *dvec, int iter){

  if(redo){
     create_source_fields();
  }
  if(!add) {
    mvec->zero();
  }

  data_rtm_3d *d=(data_rtm_3d*) dvec;
  image_rtm_3d *m=(image_rtm_3d*) mvec;
  int npts_s=source->get_points(encode);
  data->set_source_file(d);

  image->set_source_file(m);
  int nt_big=8+data->nt;
  int npts=data->get_points();
  int rec_nx=data->get_points1();
  int rec_ny=data->get_points2();

  int *locs=new int[npts];
  int *locs_s=new int[npts_s];
  float *vals=new float[npts_s*(nt_big)];

  float *s_z=new float[npts],*s_x=new float[npts],*s_y=new float[npts];
  axis ax_pt1(rec_nx);
  axis ax_pt2(rec_ny);

  float_3d *shot=new float_3d(data->get_axis(1),data->get_axis(2),data->get_axis(3));
  float_3d *shot_deriv=new float_3d(data->get_axis(1),data->get_axis(2),data->get_axis(3));
  axis ax_t=axis(nt_big);

  float_3d *rec_func=new float_3d(ax_t,ax_pt1,ax_pt2);
  int seed=0;

  sseek("source_fields",0,0);
  int nts=(int)(nt/jts)+1;
  int nshots=data->nshots();

  //nshots=1;

  for(int ishot=0; ishot<nshots; ishot++){
    if(verb) fprintf(stderr,"Adjoint shot %d of %d \n",ishot,data->nshots());
    hypercube_float *src_p0=source->create_domain(ishot);
    hypercube_float *src_p1=src_p0->clone();
    axis ad1=src_p0->get_axis(1);
    axis ad2=src_p0->get_axis(2);
    axis ad3=src_p0->get_axis(3);

    source->get_source_func(src_p0,ishot,nt_big,locs_s,vals);
    data->get_source_func(src_p0,ishot,s_x,s_y,s_z,8,nt_big,shot);

    transfer_source_func(npts_s,nt_big,locs_s,vals);

    seed=locs_s[3] - ((int) locs_s[3]/10000)*10000;

    hypercube_float *vrand=vel->rand_subspace(seed,nbound,nbound,nbound,dt,src_p0);
    //hypercube_float *vnone=vel->rand_subspace(seed,0,0,0,dt,src_p0);
    hypercube_float *img=src_p0->clone();

    transfer_vel_func1(nx,ny,nz,vrand->vals);
    //transfer_vel_func2(nx,ny,nz,vnone->vals);

    deriv der=deriv(shot_deriv,shot);
    der.adjoint(false,shot_deriv,shot);

    map_data mapit(npts,s_x,s_y,s_z,locs,src_p0,rec_func,shot_deriv,nt_big);
    mapit.adjoint(false,rec_func,shot_deriv);

    sreed("source_fields",src_p1->vals,4*ad1.n*ad2.n*ad3.n);
    sreed("source_fields",src_p0->vals,4*ad1.n*ad2.n*ad3.n);

    transfer_receiver_func(rec_nx,rec_ny,nt_big,locs,rec_func->vals/*,s_z,s_x,s_y*/);

    time_t startwtimem, endwtimem;
    startwtimem = time(&startwtimem);

    rtm_adjoint(ad1.n,ad2.n,ad3.n,jtd,src_p0->vals,src_p1->vals,img->vals,npts_s,nt/*,src,recx*/);

    endwtimem = time(&endwtimem);

    fprintf(stderr,"Adjoint, wall clock time = %f \n", difftime(endwtimem, startwtimem));

    hypercube_float *img_lap=img->clone();
    laplac lap=laplac(img_lap,img);
    lap.adjoint(false,img_lap,img);
    //lap.dot_test(true);

    //image->add_image(img);
    image->add_image(img_lap);

    delete src_p0; delete src_p1; delete img; delete img_lap;
    delete vrand;
  }

  image->write_volume();

  delete [] s_z;
  delete [] s_x;
  delete [] s_y;
  delete [] locs_s;
  delete [] vals;
  delete [] locs;
  delete rec_func;
  delete shot;
  delete shot_deriv;

  return true;
}
