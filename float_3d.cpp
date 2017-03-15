#include<float_3d.h>

float_3d::float_3d(axis a1, axis a2, axis a3)
{
  base_init_3df(a1,a2,a3);
}
float_3d::float_3d(axis a1, axis a2, axis a3, float *array)
{
 this->base_init_3df(a1,a2,a3);
 this->set(array);
 }
 void float_3d::base_init_3df(axis a1, axis a2, axis a3)
 {
   std::vector<axis> axes;
   axes.push_back(a1);
   axes.push_back(a2);
   axes.push_back(a3);
   this->init_nd(axes);
   this->allocate();
   }
  float_3d::float_3d(int n1, int n2, int n3)
  {
    axis a1=axis(n1);
    axis a2=axis(n2); 
    axis a3=axis(n3);
    this->base_init_3df(a1,a2,a3);
  }
  float_3d::float_3d(int n1, int n2, int n3,float *array)
  {
    
    axis a1=axis(n1);
    axis a2=axis(n2); 
    axis a3=axis(n3);
    this->base_init_3df(a1,a2,a3);
    this->set(array);
  }
   


