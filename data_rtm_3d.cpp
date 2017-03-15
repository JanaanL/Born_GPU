#include "data_rtm_3d.h"
#include <math.h>

data_rtm_3d::data_rtm_3d(char *tg,param_func *par){

  // Here we set up the data vectors, relevant axes and set to zero

  tag_init(tg);

  axis a=get_axis(1);
  max_time=a.o+a.d*(a.n-1);
  max_time=max_time;
  dt=a.d;
  nt=a.n;
  rec_locs=0; data=0;
  dsamp=dt;

  sz=par->get_float("rec_depth",0.);
 fprintf(stderr,"sz is %f \n",sz);

}

int data_rtm_3d::get_points(){
  // Simply get the surface grid size
  return get_axis(2).n*get_axis(3).n;
}

int data_rtm_3d::get_points1(){
  // Return x axis length
  return get_axis(2).n;
}

int data_rtm_3d::get_points2(){
  // Return y axis length
  return get_axis(3).n;
}

void data_rtm_3d::add_data(int ishot, hypercube_float *dat){

  int n1=get_axis(1).n;
 int  n2=get_axis(2).n;
  int n3=get_axis(3).n;
  sseek_block(myf->tagit.c_str(),get_axis(2).n*ishot*get_axis(3).n,get_axis(1).n*4,0);
  hypercube_float *tmp=(hypercube_float*)dat->clone_zero();
  for(long long i3=0; i3 < get_axis(3).n; i3++)
    sreed(myf->tagit.c_str(),tmp->vals+i3*n1*n2,n1*n2*4);

  for(long long i=0; i <  tmp->get_n123(); i++) tmp->vals[i]+=dat->vals[i];

  sseek_block(myf->tagit.c_str(),get_axis(2).n*ishot*get_axis(3).n,get_axis(1).n*4,0);
  for(long long i3=0; i3 < get_axis(3).n; i3++)
    srite(myf->tagit.c_str(),tmp->vals+i3*n1*n2,n1*n2*4);

  delete tmp;

}

void data_rtm_3d::get_source_func(hypercube_float *domain, int ishot, float *s_x, float *s_y, float *s_z, int nsinc,  int nts,hypercube_float *time){

  // Read the data and set up the source geometry

  axis a1=domain->get_axis(1);
  axis a2=domain->get_axis(2);
  axis a3=domain->get_axis(3);
  axis at=get_axis(1);

  // Seek to shot position and read
  sseek_block(myf->tagit.c_str(),get_axis(3).n*get_axis(2).n*ishot,get_axis(1).n*4,0);
  sreed(myf->tagit.c_str(),time->vals,get_axis(1).n*get_axis(2).n*get_axis(3).n*4);

  int i=0;

  // Set up the regular source geometry
  for(int i3=0; i3< (get_axis(3).n); i3++){
    for(int i2=0; i2< get_axis(2).n; i2++, i++){
      s_z[i]=(sz-a3.o)/a3.d;
      s_x[i]=(get_axis(2).o+get_axis(2).d*i2-a1.o)/a1.d;
      s_y[i]=(get_axis(3).o+get_axis(3).d*i3-a2.o)/a2.d;
    }
  }
  fprintf(stderr,"first receiver depth is %f \n",s_z[0]);
  sseek(myf->tagit.c_str(),0,0);
}

/*void data_rtm_3d::get_source_func_encode(hypercube_float *domain, int ishot, bool encode, int *rvec, float *s_z, float *s_x, float *s_y, int nsinc,  int nts,hypercube_float *time){

  // Read the data and set up the source geometry

  axis a1=domain->get_axis(1);
  axis a2=domain->get_axis(2);
  axis a3=domain->get_axis(3);
  axis at=get_axis(1);

  float sx=get_axis(4).o+get_axis(4).d*ishot;
  float sy=get_axis(5).o+get_axis(5).d*ishot;
  int ns=get_axis(4).n*get_axis(5).n;

  fprintf(stderr,"  RTM Data Details, encode=%d, nshots=%d \n",encode,ns);

  int i=0;

  // Encode and sum all shots into one supershot
  hypercube_float *tmp=(hypercube_float*)time->clone_zero();
  for(int is=0; is<ns; is++){
    sseek_block(myf->tagit.c_str(),get_axis(3).n*get_axis(2).n,get_axis(1).n*4*is,0);
    sreed(myf->tagit.c_str(),tmp->vals,get_axis(1).n*get_axis(2).n*get_axis(3).n*4);

    for(int k=0; k < (int) tmp->get_n123(); k++) time->vals[k] += rvec[is]*tmp->vals[k]/ns;
      fprintf(stderr,"    Data summing, %d %d \n",is,rvec[is]);


  }
fprintf(stderr,"Shots summed \n");
    // Set up the regular source geometry RIGHT NOW THIS IS THE SAME FOR ALL SHOTS
    for(int i3=0; i3< (get_axis(3).n); i3++){
      for(int i2=0; i2< get_axis(2).n; i2++, i++){
        s_z[i]=(sz-a1.o)/a1.d;
        s_x[i]=(get_axis(2).o+get_axis(2).d*i2-a2.o)/a2.d;
        s_y[i]=(get_axis(3).o+get_axis(3).d*i3-a3.o)/a3.d;
      }
    }

  sseek(myf->tagit.c_str(),0,0);

}*/


