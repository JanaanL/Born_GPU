#ifndef SUPER_VECTOR_H
#define SUPER_VECTOR_H 1
#include "my_vector.h"
#include <vector>
class super_vector: public my_vector{
  public:
    super_vector(){};
    super_vector(my_vector *v1, my_vector *v2,bool sp=false);
    super_vector(my_vector *v1, my_vector *v2, my_vector *v3,bool sp=false);
    super_vector(std::vector<my_vector*> vs,bool sp=false);
    void add_vector(my_vector *v1);
    my_vector *return_vector(int ivec);
    virtual void scale(const double num);
    virtual void scale_add(const double mes,  my_vector *vec, const double other) ;
    virtual void set_val(double  val);
    virtual void add(my_vector *vec);
    virtual double dot(my_vector *vec);
    virtual void info(char *str,int level=0);
    virtual void random();
    virtual void inverse_hessian(my_vector *vec);
    my_vector *clone_vec(){clone_it(false);}
    my_vector *clone_space(){ clone_it(true);}
    my_vector *clone_it(bool alloc);
  
    virtual bool check_match(const my_vector *v2);
    void super_init(std::string n){ name=n;}
    
    ~super_vector(){
      for(int i=0; i < (int)vecs.size(); i++) delete vecs[i];
      vecs.clear();
    }
    
    
    bool just_space;
    std::vector<my_vector*> vecs;



};

#endif
