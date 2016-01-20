#ifndef Equation_H
#define Equation_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>
#include "utilities.h"
#include "product.h"
#include "term.h"
#include "matrices.h"
#include "globals.h"
#include "lexic.h"

/*!
    Equation
 */
class Equation : public Sum<Term, TFactor> {
public:
  
};

/*!
    lexic Equation
 */
class LEquation {
public:
  // add string
  LEquation & operator += (const Lelem & lel) { _eqn.add(lel); return *this; };
  // clear _eqn
  void reseteq() { _eqn = LelString();};
  // extract expression (remove parentheses and sums)
  bool extractit();
  // transform LelString/ to Sum<Term>, if excopsonly: do only pure excitation operators (and bra/ket)
  bool do_sumterms(bool excopsonly=false);
  // add new operator
  void addnewop(const std::string& name, const LelString& oper){ _newops[name] = oper; };
  LelString eqn() const { return _eqn;};
  // get sum of terms
  Sum<Term,TFactor> sumterms() const { return _sumterms;};
private:
  // add connections to term, and term to _sumterms
  void addterm(Term& term, bool plus, lui beg, lui end, 
               Product< long int > const & indxoperterm, bool excopsonly);
  // handle bra/ket
  Oper handle_braket(Lelem const & lel, Term & term, bool excopsonly=false);
  // handle explicit excitation index (like ^{ab}_{ij})
  Oper handle_explexcitation(Term& term, std::string const & name, bool dg, bool excopsonly=false, bool phi=true);
  // handle excitation index
  Oper handle_excitation( Term& term, std::string const & name, bool dg, int lmel = 0, bool excopsonly=false );
  // handle factor
  TFactor handle_factor(Lelem const & lel) const;
  // handle operator
  Oper handle_operator(Lelem const & lel, Term & term, bool excopsonly=false);
  // handle name, ups and downs of operators 
  // (name, excitation class, name additions, dagger, non-conserving character, orbital types...)
  // returns true if found up or down
  bool handle_namupdown(std::string& name, short& excl, std::string& nameup, std::string& namedown, bool& dg, int& lmel, 
                        std::vector< Product<Orbital::Type> >& orbtypes, const std::string& lelnam ) const;
  // returns true if explicit given orbtypes found
  bool handle_orbtypes(std::vector< Product<Orbital::Type> >& orbtypes, const std::string& string) const;
  // handle sum
  void handle_sum(Lelem const & lel, Term & term) const;
  // handle parameter
  void handle_parameters(Term& term, bool excopsonly = false);
  // handle permutation
  Permut handle_permutation(Lelem const & lel) const;
  // reset term (set to Term() and reset lastorbs)
  void reset_term(Term& term) const;
  // correct used orbitals
  void correct_orbs(Term& term, const Product<Orbital>& orbs);
  
  LelString _eqn;
  Equation _sumterms;
  // save names and corresponding info of pure excitation and deexciation operators
  typedef std::map<std::string, LExcitationInfo> LExcitationMap;
  LExcitationMap _excops;
  // save parameters and sums in term
  LelString _paramterm, _sumsterm;
  // connections "map" in _eqn (starts from 1)
  // positive: connected; negative: disconnected
  // the earlier the connection comes the more important it is:
  // e.g. ((-2,-3),(2,3,4)) --> elements 2 and 3 are disconnected, and element 4 is connected to 2 and/or 3
  ConnectionsMap _connections;
  //custom operators from \newop
  NewOpMap _newops;
};

std::ostream & operator << (std::ostream & o, LEquation const & inp);



#endif
