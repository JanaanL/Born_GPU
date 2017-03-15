#ifndef BASIC_SOLVER_H
#include<my_vector.h>
#include<my_operator.h>
#include<cgstep.h>
#include <vector>


class lin_solver{
  public:
    lin_solver(){}
    
    
    
    ~lin_solver(){clean_up();}
    void set_verbose(int iv);
    bool solve(int niter);
    void init_solve(step *s);
    
  

    void create_wt_op(my_operator *op, my_operator *wt);
    void create_solver_vecs();
    my_vector *init_rr(my_vector *data, my_operator *op, my_operator *wt, my_vector *m0);
    virtual void update_model(my_vector *mod);
    virtual my_vector *return_model();
    
  
    my_vector *g,*m,*gg,*rr;
    my_operator *oper,*wt_op,*iop;
    int verb;
    bool have_op,have_wt;
    step *st;
    
  private:
    void clean_up();
    double scale;
    
    
    

};
class simple_solver:public lin_solver{
  public:
    simple_solver(){}
    simple_solver(step *s,my_vector *data, my_operator *op,
      my_operator *wt=0, my_vector *m0=0);
    ~simple_solver(){};
  

  private:
      void clean_up();

};
class reg_solver:public lin_solver{
  public:
    reg_solver(){};
    reg_solver(step *s,my_vector *data, my_operator *op,my_operator *reg, float eps,
       my_operator *wt=0, my_vector *m=00);
 
    ~reg_solver(){};
    
  private:
    void clean_up();
};
class prec_solver:public lin_solver{
  public:
    prec_solver(){};
    prec_solver(step *s,my_vector *data, my_operator *op,my_operator *reg, float eps,
     my_operator *wt=0, my_vector *m0=0);

    ~prec_solver(){clean_up();}
      virtual void update_model(my_vector *mod);
      virtual my_vector *return_model();

  private:
    void clean_up();
    my_operator *pop,*prec_chain_op;
};
#endif
