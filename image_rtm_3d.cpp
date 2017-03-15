#include "image_rtm_3d.h"
#include "deriv_3d.h"

image_rtm_3d::image_rtm_3d(char *tg,vel_fd_3d *v){
  tagit=tg;
  image=new hypercube_float(v->fld);
  std::vector<axis> axes; axes.push_back(v->get_axis(1));
  axes.push_back(v->get_axis(2));

  axes.push_back(v->get_axis(3));
  init_nd(axes);

  auxinout(tg);

  write_description(tg);

  image->zero();
//  hypercube_float *mute=(hypercube_float*)image->clone_zero();
//  sreed("mute",mute->vals,4*get_axis(1).n*get_axis(2).n*get_axis(3).n);
}

void image_rtm_3d::set_source_file(oc_float *ptr){
  myf=ptr;
  sseek(myf->tagit.c_str(),0,0);
  read_all(myf->tagit,image);
}

void image_rtm_3d::add_image(hypercube_float *img){
  axis a1o=get_axis(1);
  axis a2o=get_axis(2);
  axis a3o=get_axis(3);
  axis a1i=img->get_axis(1);
  axis a2i=img->get_axis(2);
  axis a3i=img->get_axis(3);

  int f3=(int)((a3i.o-a3o.o)/a3o.d);
  int f2=(int)((a2i.o-a2o.o)/a2o.d);
  int f1=(int)((a1i.o-a1o.o)/a1o.d);

  for(int i3=0; i3 < a3i.n; i3++){
    if(i3+f3 >=0 && i3+f3 < a3o.n){
      for(int i2=0; i2 < a2i.n; i2++){
        if(i2+f2 >=0 && i2+f2 < a2o.n){
          for(int i1=0; i1 < a1i.n; i1++){
            if(i1+f1>=0 && i1+f1 <a1o.n){
                image->vals[i1+f1+(i2+f2)*a1o.n+(i3+f3)*a1o.n*a2o.n]+=img->vals[i1+i2*a1i.n+i3*a1i.n*a2i.n];
            }
          }
        }
      }
    }
  }

}

float_3d *image_rtm_3d::extract_sub(axis a1i, axis a2i, axis a3i){
  axis a1o=get_axis(1);
  axis a2o=get_axis(2);
  axis a3o=get_axis(3);

  sseek(myf->tagit.c_str(),0,0);
  image->zero();
  read_all(myf->tagit,image);

  float_3d *img=new float_3d(a1i,a2i,a3i);
  img->zero();
/*
  hypercube_float *mute=(hypercube_float*)image->clone_zero();
  sreed("mute",mute->vals,4*a1o.n*a2o.n*a3o.n);
  sseek("mute",0,0);
  for(int h=0; h<image->get_n123(); h++){
    image->vals[h]=image->vals[h]*mute->vals[h];
  }
  */

  int f3=(int)((a3i.o-a3o.o)/a3o.d);
  int f2=(int)((a2i.o-a2o.o)/a2o.d);
  int f1=(int)((a1i.o-a1o.o)/a1o.d);

  for(int i3=0; i3 < a3i.n; i3++){
    if(i3+f3 >=0 && i3+f3 < a3o.n){
      for(int i2=0; i2 < a2i.n; i2++){
        if(i2+f2 >=0 && i2+f2 < a2o.n){
          for(int i1=0; i1 < a1i.n; i1++){
            if(i1+f1>=0 && i1+f1 <a1o.n){
              img->vals[i1+(i2)*a1i.n+(i3)*a1i.n*a2i.n]=image->vals[i1+f1+(i2+f2)*a1o.n+(i3+f3)*a1o.n*a2o.n];
            }
          }
        } 
      }
    }
  }

  return img;
}

void image_rtm_3d::write_volume(){
 // hypercube_float *mute=(hypercube_float*)image->clone_zero();
  //sreed("mute",mute->vals,4*image->get_axis(1).n*image->get_axis(2).n*image->get_axis(3).n);
  //sseek("mute",0,0);
//  for(int h=0; h<image->get_n123(); h++){
  //  image->vals[h]=image->vals[h]*mute->vals[h];
 // }
  sseek(myf->tagit.c_str(),0,0);
  write_all(myf->tagit,image);
srite("im2_out",image->vals,4*image->get_axis(1).n*image->get_axis(2).n*image->get_axis(3).n);

}

void image_rtm_3d::write_final_volume(){

  sseek("image",0,0);
  write_all("image",image);

}

