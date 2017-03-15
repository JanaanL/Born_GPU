#include "tmute.h"
#include "oc_float.h"

bool tmute::adjoint(bool add, my_vector *model, my_vector *data){
  oc_float *dat=(oc_float*) data;
  oc_float *mod=(oc_float*) model;
  if(!add) mod->zero();
  axis a1=dat->get_axis(1);
  axis a2=dat->get_axis(2);
  axis a3=dat->get_axis(3);
  axis a4=dat->get_axis(4);
  axis a5=dat->get_axis(5);
  std::vector<axis> as;
  as.push_back(a1); as.push_back(a2); as.push_back(a3);
  hypercube_float *d=new hypercube_float(as);
  hypercube_float *m=new hypercube_float(as);
  sseek(dat->tagit.c_str(),0,0);
  sseek(mod->tagit.c_str(),0,0); 
  for(int i5=0; i5 < a5.n ;i5++){
   for(int i4=0; i4 < a4.n ;i4++){

    sreed(dat->tagit.c_str(),d->vals,a1.n*a2.n*a3.n*4);
    sreed(mod->tagit.c_str(),m->vals,a1.n*a2.n*a3.n*4);

    for(int i3=0; i3 < a3.n; i3++){

     for(int i2=0; i2 < a2.n; i2++){
        float t0=t_0+(a2.o+a2.d*i2)/vmute;
        for(int i1=(int)((t0-a1.o)/a1.d); i1 < a1.n; i1++){
          m->vals [i3*a1.n*a2.n+i2*a1.n+i1]+=d->vals[i3*a1.n*a2.n+i2*a1.n+i1];
        }
     }
    }
    sseek(mod->tagit.c_str(),-a1.n*a2.n*a3.n*4,1);

    srite(mod->tagit.c_str(),m->vals,a1.n*a2.n*a3.n*4);
   }
  }
  delete d; delete m;

}
bool tmute::forward(bool add, my_vector *model, my_vector *data){
  oc_float *dat=(oc_float*) data;
  oc_float *mod=(oc_float*) model;
  axis a1=dat->get_axis(1);
  axis a2=dat->get_axis(2);
  axis a3=dat->get_axis(3);
  axis a4=dat->get_axis(4);
  axis a5=dat->get_axis(5);
  std::vector<axis> as;as.push_back(a1); as.push_back(a2); as.push_back(a3);
  hypercube_float *d=new hypercube_float(as);
  hypercube_float *m=new hypercube_float(as);
  
  if(!add) dat->zero();
  sseek(dat->tagit.c_str(),0,0);
  sseek(mod->tagit.c_str(),0,0);
  for(int i5=0; i5 < a5.n ;i5++){
   for(int i4=0; i4 < a4.n ;i4++){
    for(int i3=0; i3 < a3.n; i3++){
      sreed(dat->tagit.c_str(),d->vals,a1.n*a2.n*a3.n*4);
      sreed(mod->tagit.c_str(),m->vals,a1.n*a2.n*a3.n*4);
      for(int i2=0; i2 < a2.n; i2++){
         float t0=t_0+(a2.o+a2.d*i2)/vmute;
         for(int i1=(int)((t0-a1.o)/a1.d); i1 < a1.n; i1++){
           d->vals [i3*a1.n*a2.n+i2*a1.n+i1]+=m->vals[i3*a1.n*a2.n+i2*a1.n+i1];
         }
      }
    }
    sseek(dat->tagit.c_str(),-a1.n*a2.n*a3.n*4,1);
    srite(dat->tagit.c_str(),d->vals,a1.n*a2.n*a3.n*4);
   }
  }

  delete d; delete m;
}

