#ifndef Lexic_H
#define Lexic_H
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
#include "equation.h"

/*!
    Lexic elements
*/
class Lelem {
  public:
  // enumerate lexic
  enum Lex {Bra, Ket, LPar, RPar, Oper, Param, Num, Frac, Plus, Minus, Times, Div, Sum, Perm };
  // enumerate types expressions in parantheses (Normal, Connected, Disconnected, ...)
  enum Conn {Normal, Connect, Disconnect }; 
  // constructor from name and Lex
  Lelem (std::string name, Lex lex, Conn conn=Normal);
  std::string name() const;
  Lex lex() const;
  Conn conn() const;
  // return if the bra was already expanded
  bool expandedbra() const;
  // return same element with a status "bra is expanded"
  Lelem braexpanded() const;
  // check equality
  bool operator == (Lelem const & lel) const;
  private:
  std::string _name;
  Lex _lex;
  Conn _conn;
  bool _expandedbra;
};

std::ostream & operator << (std::ostream & o, Lelem const & lel);

typedef std::map< Orbital::Type,Orbital > TOrb4Type;
/*
 *  excitation info (for \mu_i etc) 
 */
class LExcitationInfo {
public:
  LExcitationInfo() {};
  LExcitationInfo( const TOrb4Type& orbs, short exclass, Matrices::Spinsym spin, int termpos = -1 ) :
    _orbs4excops(orbs), _exccls(exclass), _spinsymexcs(spin), _posexcopsterm(termpos){};
    
  void reset_term_info() { _posexcopsterm = -1; };
  // save indices that the excitation (and deexcitation) operators have got.
  TOrb4Type _orbs4excops;
  // save excclass
  short _exccls;
  // save spinsymmetry
  Matrices::Spinsym _spinsymexcs;
  // save position of a Matrix of the _excops in term, if -1: the _excops is not present in this term
  int _posexcopsterm;
};

/*!
    lexic Equation
 */
class LEquation {
public:
  // add string
  LEquation & operator += (const Lelem & lel) { _eqn *= lel; return *this; };
  // clear _eqn
  void reseteq() { _eqn = Product<Lelem>();};
  // extract expression (remove parentheses and sums)
  bool extractit();
  // transform Product<Lelem> to Sum<Term>, if excopsonly: do only pure excitation operators (and bra/ket)
  bool do_sumterms(bool excopsonly=false);
  // add new operator
  void addnewop(const std::string& name, const Product<Lelem>& oper){ _newops[name] = oper; };
  Product<Lelem> eqn() const { return _eqn;};
  // get sum of terms
  Sum<Term,TFactor> sumterms() const { return _sumterms;};
private:
  // get position of the closing bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  lui closbrack(Product<Lelem> const & eqn, lui ipos) const;
  // get position of the opening bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  lui openbrack(Product<Lelem> const & eqn, lui ipos) const;
  // add connections: beg and end are positions of opening and closing parantheses in aterm
  Product<long int> addconnections(const Product< Lelem >& aterm, lui beg, lui end) const;
  // expand all newops
  Product<Lelem> expandnewops(Product<Lelem> const & eqn) const;
  // find the end position of current element in aterm, if bk==true: bra and ket are treated as brackets
  lui elem(Product<Lelem> const & aterm, lui beg, bool bk=false) const;
  // find the end position of current term
  lui term(Product<Lelem> const & eqn, lui beg) const { return elem(eqn,beg,true); };
  // expand equation
  Product<Lelem> expandeqn(Product<Lelem> const & eqn, std::vector< Product<long int> > & connections) const;
  // expand a term
  Product<Lelem> expandterm(Product<Lelem> const & aterm, std::vector< Product<long int> > & connections) const;
  // expand parantheses pair
  Product<Lelem>  expandpar(Product<Lelem> const & aterm, lui beg, std::vector< Product<long int> > & connections) const;
  // test if eqn is completely expanded
  bool expanded(Product<Lelem> const & eqn) const;
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
  
  Product<Lelem> _eqn;
  Equation _sumterms;
  // save names and corresponding info of pure excitation and deexciation operators
  typedef std::map<std::string, LExcitationInfo> LExcitationMap;
  LExcitationMap _excops;
  // save parameters and sums in term
  Product<Lelem> _paramterm, _sumsterm;
  // connections "map" in _eqn (starts from 1)
  // positive: connected; negative: disconnected
  // the earlier the connection comes the more important it is:
  // e.g. ((-2,-3),(2,3,4)) --> elements 2 and 3 are disconnected, and element 4 is connected to 2 and/or 3
  std::vector< Product<long int> > _connections;
  //custom operators from \newop
  typedef std::map< std::string, Product<Lelem> > NewOpMap;
  NewOpMap _newops;
};

std::ostream & operator << (std::ostream & o, LEquation const & inp);

#endif
