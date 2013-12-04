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

class Expression;

// term in terms of DiagramTensors 
// TODO: move to a new file
class Diagram {
public:
  Diagram(){};
  // search for the best contraction order
//  std::vector<Action*> 
  void binarize(Expression& expr) const;
  // all slot types in this diagram
  SlotTs _slottypes;
  // all tensors in diagram, including the "vacuum tensor" (_tensor[0])
  std::vector<DiagramTensor> _tensors;
};

//! output operator for diagrams
std::ostream & operator << (std::ostream& o, const Diagram& d);

class Expression {
public:
  Expression();
  const SlotTypes& slottypes() const { return _slottypes; };
  const TensorsSet& tensors() const { return _tensors; };
  const SlotType * add( const SlotType& slottype ); 
  const Tensor * add( const Tensor& tensor ); 
  //! new name for a tensor. TODO: Reuse some intermediate names!
  std::string newname( const Symmetries& syms, const Cuts& cuts );
  
//private:
  SlotTypes _slottypes;
  TensorsSet _tensors;
  std::list<Contraction> _contractions;
  std::list<Summation> _summations;
  std::string _lastname;
  // diagrammatic representation  
  std::list<Diagram> _diagrams;
};

//! output operator for expression
std::ostream & operator << (std::ostream& o, const Expression& exp);

#endif
