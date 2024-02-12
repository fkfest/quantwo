#ifndef Expression_H
#define Expression_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
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
  void binarize(Expression& expr) const;
  // generates an expression-tensor from a diagram-tensor
  Tensor exprTensor( const DiagramTensor& ten ) const;
  // generates an expression-contraction from a diagrammatic contraction R=AB
  Contraction exprContraction( const DiagramTensor& tenA, const DiagramTensor& tenB, const DiagramTensor& tenR,
                               const Tensor * pA, const Tensor * pB ) const;
  // generates an expression summation from a diagrammatic "summation" R = a*A
  Summation exprSummation( const DiagramTensor& tenA, const DiagramTensor& tenR, const Tensor * pA ) const;
  // transforms to tensors and intermediates using bitmasks from binarize-function
  // if accumulate = true: adds the tensor as a new summand
  const Tensor * transform2Expr( Expression& expr, const Array<DiagramTensor>& inters, const Array<std::bitset<MAXNTENS> >& order,
                       std::bitset<MAXNTENS> bt, bool accumulate = false ) const;
  // add tensor
  const DiagramTensor * add( DiagramTensor dten, const Tensor * pTen = 0, bool pushfront = false );
  bool isresidual(const DiagramTensor& dten) const;
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
  // finds residual that corresponds to res (or creates it new) and adds an action
  const Tensor * add2residual( const Tensor& res, const Action * pAct );
  // searches for the same tensor, if considerprops==false does not consider symmetry and cuts
  const Tensor * find( const Tensor& tensor, bool considerprops = true ) const;
  //! new name for a tensor. TODO: Reuse some intermediate names!
  std::string newname( const Symmetries& syms, const Cuts& cuts );
  // print Julia TensorOperations code
  void printjulia(std::ofstream& out) const;
  std::string elemconame(const std::string& name, const SlotTs& slottypes) const;
  // sort diagram list with elemcocompare_diags compare function
  void elemcosort_diags();

//private:
  SlotTypes _slottypes;
  TensorsSet _tensors;
  std::set< const Tensor * > _residuals;
  std::list<Contraction> _contractions;
  std::list<Summation> _summations;
  std::string _lastname;
  // diagrammatic representation
  std::list<Diagram> _diagrams;
  std::map<std::string,std::string> _internames;
};

// compare function to sort diagrams by containing integral names according to order provided in names vector
static bool elemcocompare_diags(const Diagram& diagA, const Diagram& diagB){
  std::set<std::string> sorted;
  std::vector<std::string> names = {"oovv","OOVV","oOvV",
                                    "d_voov","d_VOOV","d_vOoV",
                                    "d_oooo","d_OOOO","d_oOoO",
                                    "d_vovo","d_VOVO","d_vOvO",
                                    "d_vvoo","d_VVOO","d_vVoO",
                                    "d_vvvv","d_VVVV","d_vVvV",
                                    "d_ovvo","d_OVVO","d_oVvO",
                                    "d_oVoV"};
  std::vector<std::string>::iterator posend = names.begin();
  for( std::vector<std::string>::iterator name = names.begin(); name != names.end(); name++ ){
    if ( diagA._tensors[1].name() == *name ){
      posend += std::distance(names.begin(),name);
      for( std::vector<std::string>::iterator it = names.begin(); it != posend; ++it ){
        if ( *it == diagB._tensors[1].name() ) return false;
      }
      return true;
    }
  }
  return false;
}

// prints contractions recursively
void print_code(std::ostream& o, const Tensor& ten);
//! output operator for expression
std::ostream & operator << (std::ostream& o, const Expression& exp);

#endif
