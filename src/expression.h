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
// max number of tensors in a contraction(needed in binarize)
const uint MAXNTENS = 12;
// term in terms of DiagramTensors 
// TODO: move to a new file
class Diagram {
public:
  Diagram() : _fac(1) {};
  // create an intermediate tensor from a contraction of ten1 and ten2.
  DiagramTensor newTensor( const DiagramTensor& ten1, const DiagramTensor& ten2, std::string name = "" ) const;
  // contraction cost of ten1 and ten2 to res (res has to be created before!)
  Cost contractionCost( const DiagramTensor& ten1, const DiagramTensor& ten2, const DiagramTensor& res ) const;
  // search for the best contraction order
//  std::vector<Action*> 
  void binarize(Expression& expr) const;
  // generates an expression-tensor from a diagram-tensor
  Tensor exprTensor( const DiagramTensor& ten ) const;
  // generates an expression-contraction from a diagrammatic contraction R=AB
  Contraction exprContraction( const DiagramTensor& tenA, const DiagramTensor& tenB, const DiagramTensor& tenR, 
                               const Tensor * pA, const Tensor * pB ) const;
  // generates an expression summation from a diagrammatic "summation" R = a*A
  Summand exprSummation( const DiagramTensor& tenA, Factor fac, const DiagramTensor& tenR, const Tensor * pA ) const;
  // transforms to tensors and intermediates using bitmasks from binarize-function
  // if accumulate = true: adds the tensor as a new summand 
  const Tensor * transform2Expr( Expression& expr, const Array<DiagramTensor>& inters, const Array<std::bitset<MAXNTENS> >& order,
                       std::bitset<MAXNTENS> bt, bool accumulate = false ) const;
  // add tensor
  const DiagramTensor * add( DiagramTensor dten, const Tensor * pTen = 0, bool pushfront = false );
  // all slot types in this diagram
  SlotTs _slottypes;
  // all tensors in diagram, including the "vacuum tensor", i.e., the "result" (_tensor[0])
  Array<DiagramTensor> _tensors;
  // all cuts in diagram
  Cuts _cuts;
  Factor _fac;
};

//! output operator for diagrams
std::ostream & operator << (std::ostream& o, const Diagram& d);

class Expression {
public:
  Expression();
  const SlotTypes& slottypes() const { return _slottypes; };
  const TensorsSet& tensors() const { return _tensors; };
  const std::set< const Tensor *>& residualtensors() const { return _residuals; };
  const SlotType * add( const SlotType& slottype ); 
  // if accumulate = true: add as a new summand
  const Tensor * add( const Tensor& tensor, bool accumulate = false );
  const Action * add( const Action * pAction );
  // add residual tensor
  void addresidual( const Tensor * pRes ) {_residuals.insert(pRes);};
  // finds residual that corresponds to res (or creates it new) and adds a summand
  const Tensor * add2residual( const Tensor& res, const Summand& sumd );
  // searches for the same tensor, if considerprops==false does not consider symmetry and cuts
  const Tensor * find( const Tensor& tensor, bool considerprops = true ) const;
  //! new name for a tensor. TODO: Reuse some intermediate names!
  std::string newname( const Symmetries& syms, const Cuts& cuts );
  
//private:
  SlotTypes _slottypes;
  TensorsSet _tensors;
  std::set< const Tensor * > _residuals;
  std::list<Contraction> _contractions;
  std::list<Summation> _summations;
  std::string _lastname;
  // diagrammatic representation  
  std::list<Diagram> _diagrams;
};

// prints contractions recursively
void print_code(std::ostream& o, const Tensor& ten);
//! output operator for expression
std::ostream & operator << (std::ostream& o, const Expression& exp);

#endif
