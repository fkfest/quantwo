#include "expression.h"


const SlotType* Expression::add(const SlotType& slottype)
{
  SlotTypes::iterator it = _slottypes.begin();
  it = _slottypes.insert(it,slottype);
  return &(*it);
}

const Tensor* Expression::add(const Tensor& tensor)
{
  TensorsSet::iterator it = _tensors.begin();
  it = _tensors.insert(it,tensor);
  return &(*it);
}


std::ostream & operator << (std::ostream& o, const Expression& exp) {
  o << "---- decl" << std::endl;
  // index spaces
  const SlotTypes& sts = exp.slottypes();
  SlotTypes::const_iterator ists;
  _foreach(ists,sts){
    o << *ists << std::endl;
  }
  // tensors....
  const TensorsSet& ts = exp.tensors();
  TensorsSet::const_iterator its;
  _foreach(its,ts){
    o << "tensor: " << *its << std::endl;
  }
  
  // contractions...
  
  return o;
}