#include<axis.h>
#include<math.h>
extern "C" {
#include "seplib.h"
}

using std::string;

axis::axis(int n, float o, float d, string label){
  string unit="Undefined";
  this->basic_set(n,o,d,label,unit);
}
axis::axis(int n, float o, float d, string label, string unit){
  this->basic_set(n,o,d,label,unit);
  }
axis::axis(int n, float o, float d){  
  string unit="Undefined",label="Undefined";
  this->basic_set(n,o,d,label,unit);
}
bool axis::same_axis(axis ax){
   bool match=true;  
  if(this->d < 1e-7) seperr("d of axis is to small");
  if(this->n !=ax.n) match=false;
  if(fabs((this->o-ax.o)/this->d) > .01) match=false;
  if(fabs(this->d/ax.d-1.) > .01) match=false;
  return match;
}
axis::axis(int n){
     string unit="Undefined",label="Undefined";
  this->basic_set(n,0.,1,label,unit);
 }
  void axis::basic_set(int n, float o, float d, string label, string unit)
  {
  this->n=n;
  this->o=o;
  this->d=d;
  this->label=label;
  this->unit=unit;
  }
 

 
