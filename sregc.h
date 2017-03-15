#ifndef sregc_H
#define sregc_H 1
#include<sreg.h>
#include<hypercube_complex.h>
#include<vector>


class sregc:public sreg{

  public: 
    sregc(){} //Default constructor
    sregc(std::string tag){this->tag_init(tag);}
    sregc(sreg *sep_file){
      this->file_init(sep_file);}
    sregc(std::vector<axis> ax){
   
    this->init_nd(ax);
    
    }
  
    void read_all(std::string tag,hypercube_complex *cube);
    void read_next(std::string tag, hypercube_complex *cube);
    void write_all(std::string tag, hypercube_complex *cube);
    void write_next(std::string tag, hypercube_complex *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, hypercube_complex *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, std::vector<int> fwind, hypercube_complex *cube);
    void read_slice(std::string tag,  std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_complex *cube);
    void write_slice(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_complex *cube);
    void write_slice(std::string tag,std::vector<int> nwind,  hypercube_complex *cube);
    void write_slice(std::string tag, std::vector<int> nwind, std::vector<int> fwind, hypercube_complex *cube);
  
     void read_subset(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_complex *cube);
     void write_subset(std::string tag, std::vector<int> nwind, std::vector<int> fwind, std::vector<int> jwind, hypercube_complex *cube);
  
    private:
 };
 
 #endif


