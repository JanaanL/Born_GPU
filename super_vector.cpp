#include "super_vector.h"
#include "param_func.h"


super_vector::super_vector(my_vector *v1, my_vector *v2,bool sp){
  super_init("super_vector");
  just_space=sp;
fprintf(stderr,"s_v 1 \n");

    if(just_space){

  vecs.push_back(v1->clone_space());
  vecs.push_back(v2->clone_space());
  }
  else{
     vecs.push_back(v1->clone_vec());
  vecs.push_back(v2->clone_vec());

  }
}
super_vector::super_vector(my_vector *v1, my_vector *v2, my_vector *v3,bool sp){
  super_init("super_vector");
 just_space=sp;
fprintf(stderr,"s_v 2 \n");

   if(just_space){

  vecs.push_back(v1->clone_space());
  vecs.push_back(v2->clone_space());
   vecs.push_back(v3->clone_space());
   }
   else{
  vecs.push_back(v1->clone_vec());
  vecs.push_back(v2->clone_vec());
   vecs.push_back(v3->clone_vec());
   }
}
super_vector::super_vector(std::vector<my_vector*> vs,bool sp){
  super_init("super_vector");
 just_space=sp;
fprintf(stderr,"s_v 3 \n");

   if(just_space){

  for(int i=0; i < (int)vs.size(); i++) vecs.push_back(vs[i]->clone_space());
  }
  else{
    for(int i=0; i < (int)vs.size(); i++) vecs.push_back(vs[i]->clone_vec());

  }
}
void super_vector::add_vector(my_vector *v1){
fprintf(stderr,"s_v 4 \n");

  if(just_space){
  vecs.push_back(v1->clone_space());  
  }
  else{
    vecs.push_back(v1->clone_vec());  

  }
}
my_vector *super_vector::return_vector(int ivec){
  //should give an error condition
fprintf(stderr,"s_v 5 \n");

  return vecs[ivec];
}
void super_vector::scale(const double num){
fprintf(stderr,"s_v 6 \n");

  for(int i=0; i < (int)vecs.size(); i++) vecs[i]->scale(num);
}
void super_vector::scale_add(const double mes,  my_vector *vec, const double other) {
fprintf(stderr,"s_v 7 \n");

  if(!check_match(vec)) 
    fprintf(stderr,"vectors don't match scale_add\n");
  super_vector *v=(super_vector*) vec;
  for(int i=0; i <(int) vecs.size(); i++) vecs[i]->scale_add(mes,
  v->vecs[i],other);
}
void super_vector::set_val(double  val){
fprintf(stderr,"s_v 8 \n");

  for(int i=0; i < (int) vecs.size(); i++) vecs[i]->set_val(val);
}
void super_vector::add(my_vector *vec){
fprintf(stderr,"s_v 9 \n");

 if(!check_match(vec)) 
    fprintf(stderr,"vectors don't match add\n");
    super_vector *v=(super_vector*) vec;
 for(int i=0; i <(int) vecs.size(); i++) vecs[i]->add(v->vecs[i]);
}
double super_vector::dot(my_vector *vec){
fprintf(stderr,"s_v 10 \n");

if(!check_match(vec)) 
    fprintf(stderr,"vectors don't match dot\n");
    super_vector *v=(super_vector*) vec;
 double dd=0.,xx;
   for(int i=0; i <(int) vecs.size(); i++) {
       xx=vecs[i]->dot(v->vecs[i]);
       
       dd+=xx;
     }
  return dd;

}
void super_vector::info(char *str,int level){
  fprintf(stderr,"super vector %s \n",str);
fprintf(stderr,"s_v 11 \n");

 for(int i=0; i < (int)vecs.size();i++){
   char temp_ch[128]; sprintf(temp_ch,"  portion: %d",i);
   vecs[i]->info(temp_ch,level);
 }
}
void super_vector::random(){
  for(int i=0; i < (int) vecs.size(); i++){
    vecs[i]->random();
  }

}
void super_vector::inverse_hessian(my_vector *vec){

if(!check_match(vec)) 
    fprintf(stderr,"vectors don't match scale_add\n");
    super_vector *v=(super_vector*) vec;
for(int i=0; i < (int) vecs.size();i++) vecs[i]->inverse_hessian(v->vecs[i]);

}
bool super_vector::check_match(const my_vector *v2){

  if(-1==v2->name.find("super_vector")) {
    fprintf(stderr,"not a super vector \n");
    return false;
  }
  super_vector *v=(super_vector*) v2;
  if(v->vecs.size()!=vecs.size()){
    fprintf(stderr,"number of vectors (%d %d) not the same \n",(int)v->vecs.size(),(int)vecs.size());
    return false;
  }

  for(int i=0; i <(int) vecs.size(); i++){
     if(!vecs[i]->check_match(v->vecs[i])){
       fprintf(stderr,"check matched failed for vector %d \n",i);
       return false;
     }
   }

  return true;
}
my_vector *super_vector::clone_it(bool space){
fprintf(stderr,"s_v 12 \n");

  super_vector *vec= new  super_vector(vecs,space);
  return (my_vector*) vec;
}
