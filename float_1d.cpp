#include<float_1d.h>

float_1d::float_1d(axis a1)
{
  base_init_1df(a1);
}
float_1d::float_1d(axis a1, float *array)
{
 this->base_init_1df(a1);
 this->set(array);
 }
 void float_1d::base_init_1df(axis a1)
 {
   std::vector<axis> a; a.push_back(a1);
   this->init_ndf(a);
   }
  float_1d::float_1d(int n1)
  {
    axis a1=axis(n1);

    this->base_init_1df(a1);
  }
  float_1d::float_1d(int n1, float *array)
  {
    
    axis a1=axis(n1);

    this->base_init_1df(a1);
    this->set(array);
  }
   


