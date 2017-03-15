#ifndef oc_float_H
#define oc_float_H 1
//#include<cg.h>
#include "sregf.h"
#include<axis.h>
#include<hypercube.h>
#include<my_vector.h>

class oc_float: public sregf, public my_vector{

  public:
    oc_float(){};   
    virtual my_vector *clone_vec();
    virtual my_vector *clone_space();
    void check_same(my_vector *other){ if(other==0);}
    double dot(my_vector *other);
    virtual void scale_add(const double mes,  my_vector *vec, const double other);
    void add(my_vector *other);
    oc_float(char *tag,std::vector<axis> axes,bool alloc=true);
    oc_float(std::string tag,hypercube *hyper);
    oc_float(char *tag);
    virtual void random();
    virtual void scale(double r);
    void zero_file();
    virtual void mult(my_vector *vec);
    virtual double sum();
    virtual void scale_add(const double sc1, my_vector *v1,double sc2, my_vector *v2);
    virtual double my_min();
    virtual double my_max();
    oc_float *clone(bool alloc=true, std::string tag="NONE");
    void set_val(double val);
    void normalize(float val);
    void allocate();
    void init_ndf(std::vector<axis> ax){ init_nd(ax); allocate();}
    void deallocate();
    virtual ~oc_float(){
     //if(temp)  this->deallocate();
    } 
    std::string make_temp();
    virtual bool check_match(const my_vector *v2);
    void info(char *str,int level=0);
    private:
      bool temp;
 };
 
 #endif
 

