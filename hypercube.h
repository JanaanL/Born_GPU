#ifndef HYPERCUBE_H
#define HYPERCUBE_H 1
#include<axis.h>
#include <vector>

class hypercube{

  public:
    hypercube(){}
       
    void set_axes(std::vector<axis> axes);
    void set_axis(int idim, axis myaxis);
    axis get_axis(int idim);
    long long get_n123(){return n123;}
    std::vector<int> return_ns();
     void deallocate(){
         axes.clear();
       }
     ~hypercube(){
       this->deallocate();
      }
    void init_nd(std::vector<axis> axes);
    std::vector<axis> return_axes(int nmax) const;
    int get_ndim(){return axes.size();}
    int get_ndim_g1();
    std::vector<axis> get_axes();
    bool same_size(hypercube *other);
    
  protected:
    long long n123;
    std::vector<axis> axes;
 
};
#endif

