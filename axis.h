#ifndef AXIS_H
#define AXIS_H 1
#include<iostream>
#include<string>
#include<stdio.h>

class axis{
  public:
    axis(){}  //Default constructor does nothing
    axis(int n, float o, float d, std::string label, std::string unit);
    axis(int n, float o, float d, std::string label);
    axis(int n, float o, float d);
    axis(int n);
    void print_info(int i){
      fprintf(stderr, "AXIS %d %d %f %f \n",i,n,o,d);
     }

    bool same_axis(axis ax);
    int n;
    float o,d;
    std::string label,unit;
    
  ~axis(){}
  private:
    void basic_set(int n, float o, float d, std::string label, std::string unit);
    
    };
#endif



