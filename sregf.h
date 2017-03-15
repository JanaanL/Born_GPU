#ifndef SREGF_H
#define SREGF_H 1
#include<sreg.h>
#include<hypercube_float.h>
#include<vector>


class sregf:public sreg{

  public: 
    sregf(){} //Default constructor
    sregf(std::string tag){this->tag_init(tag);}
    sregf(sreg *sep_file){
      this->file_init(sep_file);}
    sregf(std::vector<axis> ax){
   
    this->init_nd(ax);
    
    }
  
    void read_all(std::string tag,hypercube_float *cube);
    void read_next(std::string tag, hypercube_float *cube);
    void write_all(std::string tag, hypercube_float *cube);
    void write_next(std::string tag, hypercube_float *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, hypercube_float *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, std::vector<int> fwind, hypercube_float *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_float *cube);
    void write_slice(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_float *cube);
    void write_slice(std::string tag,std::vector<int> nwind,  hypercube_float *cube);
    void write_slice(std::string tag, std::vector<int> nwind, std::vector<int> fwind, hypercube_float *cube);
  
     void read_subset(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_float *cube);
     void write_subset(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_float *cube);
  
    private:
 };
 
 #endif


