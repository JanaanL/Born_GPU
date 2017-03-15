#ifndef I_OP_H
#define I_OP_H 1
#include "my_operator.h"
class i_op: public my_operator{
  public:
  i_op(my_vector *mod, my_vector *dat);
  virtual bool forward(bool add, my_vector *model, my_vector *data);
  virtual bool adjoint(bool add, my_vector *model, my_vector *data);






};
#endif
