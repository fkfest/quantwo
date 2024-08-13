#ifndef Diagram_H
#define Diagram_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <stack>
#include <algorithm>
#include <assert.h>
#include <stdint.h>
#include "globals.h"
#include "utilities.h"
#include "types.h"
#include "tensor.h"
#include "action.h"
#include "expression.h"

typedef std::map<std::string,std::string> PerMap;
class Expression;

// max number of tensors in a contraction(needed in binarize)
const uint MAXNTENS = 12;
class Diagram {
public:
  Diagram() : _fac(1) {};
  // create an intermediate tensor from a contraction of ten1 and ten2.
  DiagramTensor newTensor( const DiagramTensor& ten1, const DiagramTensor& ten2, std::string name = "" ) const;
  // contraction cost of ten1 and ten2 to res (res has to be created before!)
  Cost contractionCost( const DiagramTensor& ten1, const DiagramTensor& ten2, const DiagramTensor& res ) const;
  // search for the best contraction order
  void binarize(Expression& expr) const;
  // generates an expression-tensor from a diagram-tensor
  Tensor exprTensor( const DiagramTensor& ten ) const;
  // generates an expression-contraction from a diagrammatic contraction R=AB
  Contraction exprContraction( const DiagramTensor& tenA, const DiagramTensor& tenB, const DiagramTensor& tenR,
                               const Tensor * pA, const Tensor * pB ) const;
  // generates an expression summation from a diagrammatic "summation" R = a*A
  Summation exprSummation( const DiagramTensor& tenA, const DiagramTensor& tenR, const Tensor * pA ) const;
  // transforms to tensors and intermediates using bitmasks from binarize-function
  const Tensor * transform2Expr( Expression& expr, const Array<DiagramTensor>& inters, const Array<std::bitset<MAXNTENS> >& order,
                       std::bitset<MAXNTENS> bt ) const;
  // add tensor
  const DiagramTensor * add( DiagramTensor dten, const Tensor * pTen = 0, bool pushfront = false );
  bool isresidual(const DiagramTensor& dten) const;
  bool equal(const Diagram& diag) const;
  void createPermMap(const Array<std::string>& aslots, const Array<std::string>& bslots);
  void permute(Array<std::string>& slots);
  // all slot types in this diagram
  SlotTs _slottypes;
  // all tensors in diagram, including the "vacuum tensor", i.e., the "result" (_tensor[0])
  Array<DiagramTensor> _tensors;
  // all cuts in diagram
  Cuts _cuts;
  Factor _fac;
  PerMap _permmap;
};

//! output operator for diagrams
std::ostream & operator << (std::ostream& o, const Diagram& d);

#endif