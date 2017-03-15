#include<oc_float.h>
#include<math.h>
#include <cstdlib>
#include <string.h>
#define BUFSZ 1000000
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )

oc_float::oc_float(char *tag,std::vector<axis> axes,bool alloc){
  this->init_nd(axes);
  auxinout(tag);
  if(alloc) {
    zero_file();
    }
 
   name="oc_float";
   tagit=tag;
   write_description(tag);
   temp=false;
}
oc_float::oc_float(std::string tag,hypercube *hyper){

  int n=hyper->get_ndim();
  std::vector<axis> axes=hyper->return_axes(n);
  this->init_nd(axes);
   temp=false;

   tagit=tag;
   write_description(tag);
   auxclose(tag.c_str());
   auxinout(tag.c_str());
   name="oc_float";
}
oc_float::oc_float(char *tag){
  tag_init(tag);
   name="oc_float";
      temp=false;

}
void oc_float::scale(double r){
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    for(int i=0; i < blk; i++) buf[i]=buf[i]*r;
    sseek(tag,-blk*4,1);
    srite(tag,buf,blk*4);
    done+=block;
  }
  delete [] buf;
    sseek(tag,0,0);

}
void oc_float::allocate(){
  zero_file();
}
void oc_float::deallocate(){
  auxclose(tagit.c_str());
  auxout(tagit.c_str());
  float a;srite(tagit.c_str(),&a,4);
  auxclose(tagit.c_str());
  auxinout(tagit.c_str());
  
}
oc_float *oc_float::clone(bool alloc, std::string tag){
  int ndims=this->get_ndim();

  std::vector<axis> axes;

  
  std::string val;
  
  if(tag=="NONE") val=make_temp();
  else val=tag;
  oc_float *tmp=new oc_float(val,this);

  if(alloc) {
    tmp->allocate();
    tmp->add(this);
  }
 return tmp;
}
void oc_float::set_val(double val){
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];

  for(int i=0; i < BUFSZ; i++) buf[i]=val;
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    srite(tag,buf,blk*4);
    done+=block;
  }
  
  sseek(tag,0,0);
  delete [] buf;
}
void oc_float::normalize(float val){
  double mx=my_max();
  double sc=val/mx;
  scale(sc);
}
void oc_float::info(char *str,int level){
 double sm=0.,mymin,mymax;
 int imin=0,imax=0;
 
  char tag[1024]; strcpy(tag,tagit.c_str());
  auxclose(tag);
  auxinout(tag);
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ],*hold=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    if(done==0) memcpy(hold,buf,BUFSZ*4);
    if(done==0) {mymin=mymax=buf[0];}
    for(int i=0; i < blk;i++){
      if(mymin > buf[i]) {mymin=buf[i]; imin=i;}
      if(mymax < buf[i]) {mymax=buf[i]; imax=i;}
      sm+=buf[i]*buf[i];
    }
    //for(int i=0; i < blk; i++) buf[i]=buf[i]*r;
    done+=block;
    srite(str,buf,block*4);
  }
  delete [] buf; 
  if(1==3){
 fprintf(stderr,"N123 %d \n",(int)get_n123());
 fprintf(stderr,"    NAME=%s TYPE=%s tagit=%s\n",str,name.c_str(),tagit.c_str());

  for(int i=0; i < get_ndim(); i++) fprintf(stderr,"    n%d=%d",i+1,get_axis(i+1).n);
  fprintf(stderr,"\n       N=%d min(%d)=%g max(%d)=%g RMS=%g   \n",
      (int)get_n123(), imin, mymin,imax, mymax, sqrt(sm)/(1.0*get_n123()));
  long long print;
  if(level!=0){
  if(level<1)  print=get_n123();
  else print=MIN(level,BUFSZ);
  int ic=0;
   for(long long i=0; i <MIN(BUFSZ,sz); i++){
     if(fabs(hold[i])>.0000001 && ic <print){
      // fprintf(stderr,"val %d %f \n",(int)i,hold[i]);
       ic+=1;
       }
   } 
  }
  }
  delete [] hold;
}
my_vector *oc_float::clone_vec(){
  //my_vector *m=(my_vector*) this->clone();
  oc_float *m=this->clone();
  return (my_vector*)m;
}
my_vector *oc_float::clone_space(){
  oc_float *m=this->clone(false);
  return (my_vector*)m;
}
void oc_float::random(){
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
//srand ( time(NULL) );
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    for(int i=0; i < blk; i++) buf[i]=(float)rand()/(float)RAND_MAX-.5;
    srite(tag,buf,block*4);
    done+=block;
  }
  delete [] buf;
  sseek(tag,0,0);

}
double oc_float::dot(my_vector *other){
  double ret=0;
  check_same(other);
  oc_float *o=(oc_float*) other;
  char tag2[1024]; strcpy(tag2,o->tagit.c_str());
  sseek(tag2,0,0);
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  float *buf2=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    if(strcmp(tag,tag2)==0) memcpy(buf2,buf,sizeof(float)*blk);
    else sreed(tag2,buf2,block*4);
    for(int i=0; i < blk; i++){
      ret+=(double)buf[i]*(double)buf2[i];
      //fprintf(stderr,"check %d %f %f %g \n",i,buf[i],buf2[i],ret);
    }
    done+=block;
  }
  fprintf(stderr,"DOT %s %s %g \n",tag,tag2,ret);
  delete [] buf; delete [] buf2;
  sseek(tag,0,0);

  return ret;
}
void  oc_float::mult(my_vector *other){
 
  double ret=0;
  check_same(other);
  oc_float *o=(oc_float*) other;
  char tag2[1024]; strcpy(tag2,o->tagit.c_str());
  sseek(tag2,0,0);
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  float *buf2=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    sreed(tag2,buf2,block*4);
    for(int i=0; i < blk; i++)buf[i]=buf[i]*buf2[i];
    sseek(tag,-blk*4,1);
    srite(tag,buf,blk*4);
    done+=block;
  }
  delete [] buf; delete [] buf2;
  sseek(tag,0,0);

 
}
void  oc_float::scale_add(double sc1, my_vector *v1, double sc2, my_vector *v2){
 
  
  double ret=0;
  check_same(v1); check_same(v2);
  oc_float *ver1=(oc_float*) v1,*ver2=(oc_float*)v2;
  char tag2[1024]; strcpy(tag2,ver1->tagit.c_str());
  sseek(tag2,0,0);
  char tag[1024]; strcpy(tag,ver2->tagit.c_str());
  sseek(tag,0,0);
  char tago[1024]; strcpy(tago,tagit.c_str());
  sseek(tago,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  float *buf2=new float[BUFSZ];
  float *out=new float[BUFSZ];
  
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    sreed(tag2,buf2,block*4);
    for(int i=0; i < blk; i++)out[i]=buf[i]*sc1+sc2*buf2[i];
    srite(tago,out,blk*4);
    done+=block;
  }
  delete [] buf; delete [] buf2; delete [] out;

   sseek(tag,0,0);

}
double  oc_float::sum(){
  double my_s=0;
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    for(int i=0; i < blk; i++) my_s+=buf[i];
    done+=block;
  }
    sseek(tag,0,0);

  delete [] buf;
  return my_s;
}
double oc_float::my_min(){
  double my_min;
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    if(done==0) my_min=buf[0];
    for(int i=0; i < blk; i++) buf[i]=MIN(buf[i],my_min);
    done+=block;
  }
    sseek(tag,0,0);

  delete [] buf;
}
double oc_float::my_max(){
  double my_max;
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,block*4);
    if(done==0) my_max=buf[0];
    for(int i=0; i < blk; i++) buf[i]=MAX(buf[i],my_max);
    done+=block;
  }
    sseek(tag,0,0);

  delete [] buf;
}
void oc_float::scale_add(const double scale_me, my_vector *other, const double scale_other){

  check_same(other);
  oc_float *o=(oc_float*) other;
  char tag2[1024]; strcpy(tag2,o->tagit.c_str());

  sseek(tag2,0,0);
//fprintf(stderr,"scale_add \n");
  char tag[1024]; strcpy(tag,tagit.c_str());
//  fprintf(stderr,"READING FROM %s %s to %s \n",tag,tag2,tag);
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  float *buf2=new float[BUFSZ];
//fprintf(stderr,"scale_add %d %lld\n",BUFSZ,(long long int) sz);

  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    sreed(tag,buf,blk*4);
    sreed(tag2,buf2,blk*4);
    for(int i=0; i < blk; i++)buf[i]=buf[i]*scale_me+scale_other*buf2[i];
    sseek(tag,-blk*4,1);
    srite(tag,buf,blk*4);
    done+=block;
  }
    sseek(tag,0,0);

  delete [] buf; delete [] buf2;

}
void oc_float::add(my_vector *other){
  scale_add(1.,other,1.);
}
bool oc_float::check_match(const my_vector *v2){
   if(-1==v2->name.find("oc_float")){
     fprintf(stderr,"vector not oc_float");
     return false;
   }
   oc_float *h2=(oc_float*) v2;

   if(get_ndim_g1()!=h2->get_ndim_g1()){
     fprintf(stderr,"vectors not the same number of dimensions\n");
     return false;
   }

   for(int i=0; i < get_ndim_g1(); i++){
     axis a1=get_axis(i+1),a2=h2->get_axis(i+1);
     if(a1.n!=a2.n){
       fprintf(stderr,"vectors axis=%d not the same number of samples %d,%d \n",
         i+1,a1.n,a2.n);
       return false;
     }

  if(fabs((a1.o-a2.o)/a1.d) > .01){
       fprintf(stderr,"vectors axis=%d not the same origin %f,%f \n",
         i+1,a1.o,a2.o);
       return false;
     }  if(a1.n!=a2.n){
       fprintf(stderr,"vectors axis=%d not the same sampling %f,%f \n",
         i+1,a1.d,a2.d);
       return false;
     }     
   }

   return true;
}
void oc_float::zero_file(){
  char tag[1024]; strcpy(tag,tagit.c_str());
  sseek(tag,0,0);
  long sz=get_n123(), done=0,block;
  int blk;
  float *buf=new float[BUFSZ];
  for(int i=0; i < BUFSZ;i++) buf[i]=0;
  while(done < sz){
    block=MIN(sz-done,BUFSZ); blk=(int)block;
    srite(tag,buf,block*4);
    done+=block;
  }
  delete [] buf;
  sseek(tag,0,0);

}
std::string oc_float::make_temp(){
  char temp_file[4096];
  strcpy(temp_file,"TEMP_XXXXXX"); mkstemp(temp_file);
  std::string tmp=temp_file;
  return tmp;

}
