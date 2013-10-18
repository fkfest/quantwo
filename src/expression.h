#ifndef Expression_H
#define Expression_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <stdint.h>
#include "globals.h"
#include "utilities.h"
#include "types.h"
#include "tensor.h"
#include "action.h"


class Expression {
public:
  Expression(){};
  const SlotTypes& slottypes() const { return _slottypes; };
  const TensorsSet& tensors() const { return _tensors; };
  const SlotType * add(const SlotType& slottype); 
  const Tensor * add(const Tensor& tensor); 
  
//private:
  SlotTypes _slottypes;
  TensorsSet _tensors;
  std::list<Contraction> _contractions;
  std::list<Summation> _summations;
  
};

//! output operator for expression
std::ostream & operator << (std::ostream& o, const Expression& exp);

#endif