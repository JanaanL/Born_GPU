#include "laplac_3d.h"
#include "float_3d.h"

bool laplac::adjoint(bool add, my_vector *model, my_vector *data,int iter){
  float_3d *d=(float_3d*) data;
  float_3d *m=(float_3d*) model;
  if(!add) m->zero();

  axis a1=d->get_axis(1);
  axis a2=d->get_axis(2);
  axis a3=d->get_axis(3);
  int i=0;

  for(int i3=0; i3 < a3.n-1; i3++){
    for(int i2=0; i2 < a2.n; i2++){
      for(int i1=0; i1 < a1.n; i1++){
	i=i3*a1.n*a2.n+i2*a1.n+i1;
	//m->vals[i]=-6*d->vals[i]+d->vals[i+1]+d->vals[i-1]+d->vals[i+a1.n]+d->vals[i-a1.n]+d->vals[i+a1.n*a2.n]+d->vals[i-a1.n*a2.n];
	m->vals[i]+=d->vals[i+a1.n*a2.n]-d->vals[i];
      }
    }
  }

  /*for(int i3=0; i3 < a3.n; i3++){
    m->vals [i3*a1.n*a2.n]+= d->vals[i3*a1.n*a2.n+1] - 2*d->vals[i3*a1.n*a2.n] + d->vals[i3*a1.n*a2.n-1];
    for(int i2=0; i2 < a2.n; i2++){
      m->vals [i2*a1.n+i3*a1.n*a2.n]+=d->vals[i2*a1.n+i3*a1.n*a2.n+1] - 2*d->vals[i2*a1.n+i3*a1.n*a2.n ] +d->vals[i2*a1.n+i3*a1.n*a2.n-1];
      for(int i1=0; i1 < a1.n; i1++){
         m->vals [i3*a1.n*a2.n+i2*a1.n+i1]+=d->vals[i1+i2*a1.n+i3*a1.n*a2.n+1] - 2*d->vals[i1+i2*a1.n+i3*a1.n*a2.n]+ d->vals[i1+i2*a1.n+i3*a1.n*a2.n-1];
      }
    }
  }*/

  /*for(int i3=0; i3 < a3.n; i3++){
    m->vals [i3*a1.n*a2.n]+=(d->vals[i3*a1.n*a2.n+1] - d->vals[i3*a1.n*a2.n])/a1.d;
    for(int i2=0; i2 < a2.n; i2++){
      m->vals [i2*a1.n+i3*a1.n*a2.n]+=(d->vals[i2*a1.n+i3*a1.n*a2.n+1] - d->vals[i2*a1.n+i3*a1.n*a2.n ])/a1.d;
      for(int i1=1; i1 < a1.n; i1++){
         m->vals [i3*a1.n*a2.n+i2*a1.n+i1]+=(d->vals[i1+i2*a1.n+i3*a1.n*a2.n+1] - d->vals[i1+i2*a1.n+i3*a1.n*a2.n])/a1.d;
	}
    }
  }*/

}

bool laplac::forward(bool add, my_vector *model, my_vector *data,int iter){
  hypercube_float *d=(float_3d*) data;
  hypercube_float *m=(float_3d*) model;

  axis a1=d->get_axis(1);
  axis a2=d->get_axis(2);
  axis a3=d->get_axis(3);
  if(!add) d->zero();
  int i=0;

  for(int i3=0; i3 < a3.n-1; i3++){
    for(int i2=0; i2 < a2.n; i2++){
      for(int i1=0; i1 < a1.n; i1++){
	i=i3*a1.n*a2.n+i2*a1.n+i1;
	//d->vals[i]=-6*m->vals[i]+m->vals[i+1]+m->vals[i-1]+m->vals[i+a1.n]+m->vals[i-a1.n]+m->vals[i+a1.n*a2.n]+m->vals[i-a1.n*a2.n];
	d->vals[i+a1.n*a2.n]+=m->vals[i];
	d->vals[i]-=m->vals[i];
      }
    }
  }

  /*for(int i3=0; i3 < a3.n; i3++){
    d->vals [1+i3*a1.n*a2.n]+= m->vals [i3*a1.n*a2.n];
    d->vals [i3*a1.n*a2.n]-= 2*m->vals [i3*a1.n*a2.n];
    d->vals [-1+i3*a1.n*a2.n]+= m->vals [i3*a1.n*a2.n];
    for(int i2=0; i2 < a2.n; i2++){
      d->vals [1+i2*a1.n+i3*a1.n*a2.n]+= m->vals [i2*a1.n+i3*a1.n*a2.n];
      d->vals [i2*a1.n+i3*a1.n*a2.n]-= 2*m->vals [i2*a1.n+i3*a1.n*a2.n];
      d->vals [-1+i2*a1.n+i3*a1.n*a2.n]+= m->vals [i2*a1.n+i3*a1.n*a2.n];
      for(int i1=0; i1 < a1.n; i1++){
         d->vals [i1+i2*a1.n+i3*a1.n*a2.n]-= 2*m->vals [i3*a1.n*a2.n+i2*a1.n+i1];
         d->vals [i1+i2*a1.n+i3*a1.n*a2.n+1]+= m->vals [i3*a1.n*a2.n+i2*a1.n+i1];
         d->vals [i1+i2*a1.n+i3*a1.n*a2.n-1]+= m->vals [i3*a1.n*a2.n+i2*a1.n+i1];
       }
    }
  }*/

  /*for(int i3=0; i3 < a3.n; i3++){
    d->vals [1+i3*a1.n*a2.n]+= m->vals [i3*a1.n*a2.n]/a1.d;
    d->vals [i3*a1.n*a2.n]-= m->vals [i3*a1.n*a2.n]/a1.d;
    for(int i2=0; i2 < a2.n; i2++){
      d->vals [1+i2*a1.n+i3*a1.n*a2.n]+= m->vals [i2*a1.n+i3*a1.n*a2.n]/a1.d;
      d->vals [i2*a1.n+i3*a1.n*a2.n]-= m->vals [i2*a1.n+i3*a1.n*a2.n]/a1.d;
      for(int i1=1; i1 < a1.n; i1++){
         d->vals [i1+i2*a1.n+i3*a1.n*a2.n]-= m->vals [i3*a1.n*a2.n+i2*a1.n+i1]/a1.d;
         d->vals [i1+i2*a1.n+i3*a1.n*a2.n+1]+= m->vals [i3*a1.n*a2.n+i2*a1.n+i1]/a1.d;
       }
    }
  }*/

}

