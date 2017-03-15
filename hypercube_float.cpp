#include<hypercube_float.h>
#include<math.h>
#include <cstdlib>
#include <string.h>
#define MIN(a,b) ( ((a)<(b)) ? (a):(b) )

hypercube_float::hypercube_float(std::vector<axis> axes,bool alloc){
  this->init_nd(axes);
  if(alloc) {
    this->vals=new float[this->get_n123()];
    }
  else this->vals=0;
   name="hypercube_float";

}
hypercube_float::hypercube_float(hypercube *hyper){
  int n=hyper->get_ndim();
  std::vector<axis> axes=hyper->return_axes(n);
  this->init_nd(axes);
  this->vals=new float [this->get_n123()];
   name="hypercube_float";
}
hypercube_float::hypercube_float(std::vector<axis> axes, float *vals){

  this->init_nd(axes);
  this->vals=new float[this->get_n123()];
  this->set(vals);
   name="hypercube_float";
}
hypercube_float *hypercube_float::clone(bool alloc){
 int ndims=this->get_ndim();
std::vector<axis> axes;
 for(int i=0; i< ndims; i++) axes.push_back(get_axis(i+1));

 hypercube_float *tmp=new hypercube_float(axes,alloc);
 
 if(alloc && this->vals!=0) {
   for(int i=0; i< get_n123(); i++) tmp->vals[i]=this->vals[i];
 }
 else if(alloc){
    for(int i=0; i< get_n123(); i++) tmp->vals[i]=0;
 }
 return tmp;
}
void hypercube_float::add(hypercube_float *vec){

  for(int i=0; i < get_n123(); i++) vals[i]+=vec->vals[i];
}

void hypercube_float::set(float *array){

       memcpy((void*)this->vals,(const void*) array, 
         sizeof(float)*this->get_n123());

}
void hypercube_float::set_val(double val){
  for(int i=0; i < this->get_n123(); i++) this->vals[i]=val;
}
void hypercube_float::normalize(float val){
 float mymax=vals[0];

 for(int i=1; i< get_n123(); i++){
 
   if(mymax < vals[i]) mymax=vals[i];
 }
 float sc;
 sc=val/mymax;

 for(int i=0; i< get_n123(); i++){
   vals[i]=vals[i]*sc;
   
  }

}
void hypercube_float::info(char *str,int level){
 double sm=0.,mymin,mymax;
 mymax=mymin=vals[0];
 int imin=0,imax=0;
fprintf(stderr,"N123 %d \n",(int)get_n123());
 fprintf(stderr,"    NAME=%s TYPE=%s \n",str,name.c_str());
 for(int i=0; i < get_n123(); i++){
   if(mymin > vals[i]) {mymin=vals[i]; imin=i;}
   if(mymax < vals[i]) {mymax=vals[i]; imax=i;}
   sm+=vals[i]*vals[i];
  }

  for(int i=0; i < get_ndim(); i++) fprintf(stderr,"    n%d=%d",i+1,get_axis(i+1).n);
  for(int i=0; i < get_ndim(); i++) fprintf(stderr,"    o%d=%f",i+1,get_axis(i+1).o);
  for(int i=0; i < get_ndim(); i++) fprintf(stderr,"    d%d=%f",i+1,get_axis(i+1).d);
  fprintf(stderr,"\n       N=%d min(%d)=%g max(%d)=%g RMS=%g   \n",
      (int)get_n123(), imin,
mymin,
imax,
mymax,
sqrt(sm)/(1.0*get_n123()));
  long long print;
  if(level!=0){
  if(level<1)  print=get_n123();
  else print=MIN(level,(int)get_n123());
   for(long long i=0; i < print; i++){
     //fprintf(stderr,"val %d %f \n",(int)i,vals[i]);
  } 
  }

 
}
my_vector *hypercube_float::clone_vec(){
  //my_vector *m=(my_vector*) this->clone();
  hypercube_float *m=this->clone();
  return (my_vector*)m;
}
my_vector *hypercube_float::clone_space(){
  //my_vector *m=(my_vector*) this->clone();
  hypercube_float *m=this->clone(false);

  return (my_vector*)m;
}
void hypercube_float::random(){
  for(int i=0; i < get_n123(); i++){
    vals[i]=(float)rand()/(float)RAND_MAX-.5;
  }
}
double hypercube_float::dot(my_vector *other){
  double ret=0;
  check_same(other);
  hypercube_float *o=(hypercube_float*) other;
  for(long long i=0; i < get_n123(); i++){
    ret+=(double)this->vals[i]*(double)o->vals[i];
  }
//  fprintf(stderr,"IN DOT %g %d \n",ret,get_n123());
  return ret;
}
void  hypercube_float::mult(my_vector *other){
 
  check_same(other);
  hypercube_float *o=(hypercube_float*) other;
  for(long long i=0; i < get_n123(); i++){
   vals[i]=vals[i]*o->vals[i];
  }
 
}

void  hypercube_float::scale_add(double sc1, my_vector *v1, double sc2, my_vector *v2){
 
  check_same(v1);
    check_same(v2);
  hypercube_float *v_1=(hypercube_float*) v1;
    hypercube_float *v_2=(hypercube_float*) v2;
  for(long long i=0; i < get_n123(); i++){
   vals[i]=v_1->vals[i]*sc1+v_2->vals[i]*sc2;
  }
 
}
double  hypercube_float::sum(){
 
  double my_s=0;
  for(long long i=0; i < get_n123(); i++){
   my_s+=vals[i];
  }
 return my_s;
}

void hypercube_float::take_min(my_vector*other,my_vector *change){
 check_same(other);
       
 hypercube_float *c;
  hypercube_float *o=(hypercube_float*) other;
   if(change!=0){
       check_same(change);
       c=(hypercube_float*) change;
     }
  for(long long i=0; i < get_n123(); i++){
    if(o->vals[i] < vals[i]){
      vals[i]=o->vals[i];
       if(change!=0) c->vals[i]=1;
    }
  }
 
 }
double hypercube_float::my_min(){
  double my_min=vals[0];
  for(long long i=1; i < get_n123(); i++){
    if(vals[i]< my_min) my_min=vals[i];
  }
  return my_min;
}
 double hypercube_float::my_max(){
  double my_max=vals[0];
  for(long long i=1; i < get_n123(); i++){
    if(vals[i]>my_max) my_max=vals[i];
  }
  return my_max;
}


 void hypercube_float::take_max(my_vector *other,my_vector *change){
     check_same(other);
hypercube_float *c;
     if(change!=0){
       check_same(change);
       c=(hypercube_float*) change;
     }
  hypercube_float *o=(hypercube_float*) other;
  for(long long i=0; i < get_n123(); i++){
    if(o->vals[i] > vals[i]) {
      vals[i]=o->vals[i];
      if(change!=0) c->vals[i]=1;
    }
  }
    
    }
 void hypercube_float::scale_add(const double scale_me, my_vector *other, const double scale_other){
  check_same(other);
  hypercube_float *o=(hypercube_float*) other;
  for(int i=0; i< get_n123();i++){ this->vals[i]=this->vals[i]*scale_me+
    scale_other*o->vals[i];
  }
}
void hypercube_float::add(my_vector *other){
  check_same(other);
  hypercube_float *o=(hypercube_float*) other;
  for(int i=0; i < get_n123(); i++) this->vals[i]+=o->vals[i];
}
bool hypercube_float::check_match(const my_vector *v2){
   if(-1==v2->name.find("hypercube_float")){
     fprintf(stderr,"vector not hypercube_float");
     return false;
   }
   hypercube_float *h2=(hypercube_float*) v2;

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
