#ifndef Expression_H
#define Expression_H

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
#include "diagram.h"

class Diagram;

class Expression {
public:
  Expression();
  struct comp_pTen{
    bool operator()(const Tensor* lhs, const Tensor* rhs) const{
      if (*lhs < *rhs) return true;
      else if (*rhs < *lhs) return false;
      else return false;
    }
  };
  const SlotTypes& slottypes() const { return _slottypes; };
  const TensorsList& tensors() const { return _tensors; };
  const std::set< const Tensor *, comp_pTen >& residualtensors() const { return _residuals; };
  const SlotType * add( const SlotType& slottype );
  const Tensor * add( const Tensor& tensor);
  const Action * add( const Action * pAction );
  // add residual tensor
  void addresidual( const Tensor * pRes ) {_residuals.insert(pRes);};
  // finds residual that corresponds to res (or creates it new) and adds an action
  const Tensor * add2residual( const Tensor& res, const Action * pAct );
  // searches for the same tensor, if considerprops==false does not consider symmetry and cuts
  const Tensor * find( const Tensor& tensor, bool considerprops = true ) const;
  //! new name for a tensor. TODO: Reuse some intermediate names!
  std::string newname( const Symmetries& syms, const Cuts& cuts );
  void equalDiagrams();
  // print Julia TensorOperations code
  void printjulia(std::ofstream& out) const;
  // print Julia tensor load and drop
  void printjulia(std::ofstream& out, const std::string& tensorname, std::stack<std::string>& LIFO) const;
  void printfac(std::ofstream& out, const Factor& fac) const;
  //! penalize virtuals in tensoropt calls
  std::string juliacost(const std::vector<SlotTs>& slottypes, const Array<std::string>& resslots, 
                        const Array<std::string>& aslots, const Array<std::string>& bslots,
                        const Array<std::string>& cslots) const;
  std::string elemconame(const std::string& name, const SlotTs& slottypes) const;
  //! sort diagram list with elemcocompare_diags compare function
  void elemcosort_diags();
  //! compare function to sort diagrams by integral names according to order provided in the function
  static bool elemcocompare_diags(const Diagram& diagA, const Diagram& diagB);
  //! calculate a number specific to the "external orbital structure" (number and commutative position) in aslots 
  uint extorb(const Array<std::string>& resslots, const Array<std::string>& aslots) const;

//private:
  SlotTypes _slottypes;
  TensorsList _tensors;
  std::set< const Tensor *,comp_pTen > _residuals;
  std::list<Contraction> _contractions;
  std::list<Summation> _summations;
  std::string _lastname;
  // diagrammatic representation
  std::list<Diagram> _diagrams;
  std::map<std::string,std::string> _internames;
};


// prints contractions recursively
void print_code(std::ostream& o, const Tensor& ten);
//! output operator for expression
std::ostream & operator << (std::ostream& o, const Expression& exp);

#endif
