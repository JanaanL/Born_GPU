#ifndef PARAM_FUNC_H
#define PARAM_FUNC_H 1
#include<string>
#include<stdbool.h>


class param_func {
  public:
  
    param_func();
    virtual int get_int(const char *arg) const { if(arg==0); return 0;}
    virtual int get_int(const char *arg, int def) const {if(arg==0 && def==0); return 0;}
    virtual float get_float(const char *arg, float def)const { if(arg==0 && def==0.);return 0.;}
    virtual float get_float(const char *arg) const{ if(arg==0); return 0.;}
    virtual std::string get_string(const char *arg)const{if(arg==0); return "A";}
    virtual std::string get_string(const char *arg, const char *def)const{if(arg==0 && def==0);
      return "A";}
    virtual std::string get_string(const std::string arg)const{
       if(arg==""); return "A";}
    virtual std::string get_string(const std::string arg, std::string def)const{
      if(arg=="" && def=="") ; return "A";}
    virtual int get_int(const std::string arg)const{if(arg=="");return 0;}
    virtual int get_int(const std::string arg, int def)const{ if(def==0 && arg=="");return 0;}
    virtual int *get_ints(const std::string arg, int nmax)const{
       if(arg=="" && nmax==0) ;
      int *n=new int[1]; return n; }
    virtual int *get_ints(const std::string  arg, int nmax, int *def)const{
     if(arg=="" && nmax==0 && def==0) ;
      int *n=new int[1]; return n; }
    virtual float get_float(const std::string arg, float def)const{
       if(arg=="" && def==0) ;return 0.;}
    virtual float get_float(const std::string arg)const{
    if(arg=="");
      return 0.;}
    virtual int *get_ints(const char *arg, int nmax) const {
      if(arg=="" && nmax==0);
      int *n=new int[1]; return n; }
    virtual int *get_ints(const char *arg, int nmax, int *def) const {
      if(arg=="" && nmax==0 && def==0) ;
     int *n=new int[1]; return n; }
    virtual bool is_input_file(const std::string tag)const{if(tag=="");return false;}
    virtual bool is_inout_file(const std::string tag)const{if(tag=="");return false;}
    virtual void error(const std::string aa) const {if(aa=="");}
    virtual ~param_func(){};
 };
#endif

