 #include<sreg.h>
#include <math.h>
#include <string.h>
using std::string;
 
 sreg::sreg(string tag){
   this->tag_init(tag);
 }
 sreg::sreg(sreg *sfile){
    this->file_init(sfile);

 }
 sreg::sreg(std::vector<axis> axes){
   this->set_axes(axes);

 }

 void sreg::tag_init(string tag){
   std::vector<axis> mya;
   int n=1;
   int nt;
  float dt,ot;
  char junk[124];
   while(0==sep_get_data_axis_par(tag.c_str(),&n,&nt,&ot,&dt,junk)){
     n++;
     if(nt==1 && fabsf(ot)<.0001 && fabsf(dt-1.)< .00001) break;
    }
   n--;

   for(int i=1; i <= n; i++) {
     mya.push_back(read_axis(tag,i));
   }
   tagit=tag;
   this->init_nd(mya);
   title=get_string("title",tag);
 }
 void sreg::file_init(sreg *sfile){

   std::vector<axis> mya;
   for(int i=0; i < sfile->get_ndim(); i++) mya.push_back(sfile->get_axis(i+1));
   this->init_nd(mya);


 }
 void sreg::write_description(string tag){
   for(int i=1; i <=this->get_ndim();i++) this->write_axis(tag,i,this->get_axis(i));
 }
 
 axis sreg::read_axis(string tag, int iax){
    axis a;
    char l[128];
    if(0!=sep_get_data_axis_par(tag.c_str(),&iax,&a.n,&a.o,&a.d,l)) 
     seperr("trouble getting axis %d from tag %s \n",iax,l);
    a.label=l;
    return a;
  }
  
 void sreg::write_axis(string tag, int iax, axis a){
    char l[128];
    strcpy(l,a.label.c_str());
    if(0!=sep_put_data_axis_par(tag.c_str(),&iax,&a.n,&a.o,&a.d,l)) 
      seperr("trouble writing axis %d to %s\n",iax,tag.c_str());
 }
 
  int sreg::get_int(const char *arg) const{
    int ret;

    if(0==auxpar((char*)arg,"d",&ret,(char*)tagit.c_str())) seperr("Couldn't grab %s\n",arg);
    return ret;
  }
  int sreg::get_int(const char *arg, int def) const
   {
    int ret;

    int ierr=auxpar((char*)arg,"d",&ret,(char*)tagit.c_str());
     if(ierr==0) ret=def;
    

    return ret;
  }
  float sreg::get_float(const char *arg) const
  {
    float ret;
    if(0==auxpar((char*)arg,"f",&ret,(char*)tagit.c_str())) seperr("Couldn't grab %s\n",arg);
    return ret;
  }
  float sreg::get_float(const char *arg, float def) const
  {
    float ret;
    if(0==auxpar((char*)arg,"f",&ret,(char*)tagit.c_str())) ret=def;
    return ret;
  }
   std::string sreg::get_string(const std::string arg) const
   {
    char x[1024];
    std::string ret;
    if(0==auxpar((char*)arg.c_str(),"s",&x,(char*)tagit.c_str())) seperr("Couldn't grab %s\n",arg.c_str());
    ret=x;
    return ret;
  }
     std::string sreg::get_string(std::string arg, std::string def) const
     {
    char x[1024];
    std::string ret;
    int iret=auxpar((char*)arg.c_str(),"s",x,(char*)tagit.c_str());
    if(iret==0) ret=def;
    else ret=x;
    return ret;
  };
