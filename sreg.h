#ifndef SREG_H
#define SREG_H 1
#include<axis.h>
#include<hypercube.h>
#include<string>
extern "C" {
#include "sep3d.h"
}
class sreg: public hypercube
{
  public:  
    sreg(){} //Default constructor
    sreg(std::string tag);
    sreg(sreg *sep_file);
    sreg(std::vector<axis> axes);
    void file_init(sreg *sep_file);
    void beg_file(std::string tag){ sseek(tag.c_str(),0,0);}
    void tag_init(std::string tag);
    void write_description(std::string tag);
     int get_int(const char *arg) const;
    int get_int(const char *arg, int def) const;
    float get_float(const char *arg, float def) const;
    float get_float(const char *arg) const;
 //   std::string get_string(const char *arg) const;
//    std::string get_string(const char *arg, char *def) const;
    
    int get_int(const std::string arg) const {return this->get_int(arg.c_str());}
    int get_int(const std::string arg, int def) const {
  
       return this->get_int(arg.c_str(),def);}
    float get_float(const std::string arg, float def) const {return this->get_float(arg.c_str(),def);}
    float get_float(const std::string arg) const {return this->get_float(arg.c_str());}
    std::string get_string(const std::string arg) const ;
    std::string get_string(const std::string arg, std::string def) const;
    void set_title(std::string tit){title=tit;}
    std::string get_title(){ return title;}
    

   
     ~sreg(){}
  
    std::string tagit,title;
  private:

    axis read_axis(std::string tag, int iax);
    void write_axis(std::string tag, int iax, axis mya);

};
#endif


