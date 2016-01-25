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

typedef std::map< Orbital::Type,Orbital > TOrb4Type;
/*
 *  excitation info (for \mu_i etc) 
 */
class LExcitationInfo {
public:
  LExcitationInfo() {};
  LExcitationInfo( const TOrb4Type& orbs, std::vector<OrbitalTypes> orbts, 
                   short exclass, bool dgr, Matrices::Spinsym spin, int termpos = -1 ) :
    _orbs4excops(orbs),_orbtypes(orbts), _exccls(exclass), _dg(dgr), 
    _spinsymexcs(spin), _posexcopsterm(termpos){};
    
  void reset_term_info() { _posexcopsterm = -1; };
  TOrb4Type & orbs4excops() { return _orbs4excops; };
  const std::vector<OrbitalTypes> & orbtypes() const { return _orbtypes; };
  short exccls() const { return _exccls; };
  bool dg() const { return _dg;};
  int lmel() const { assert(_orbtypes.size() == 2); return _orbtypes[1].size()-_orbtypes[0].size(); };
  Matrices::Spinsym spinsymexcs() const { return _spinsymexcs;};
  int posexcopsterm() const { return _posexcopsterm;};
  void set_posexcopsterm(int pos){ _posexcopsterm = pos;};
private:
  // indices that the excitation (and deexcitation) operators have got.
  TOrb4Type _orbs4excops;
  // orbital types (first set for occ and second set for virt indices )
  std::vector<OrbitalTypes> _orbtypes;
  // excclass
  short _exccls;
  // dagger
  bool _dg;
  // spinsymmetry
  Matrices::Spinsym _spinsymexcs;
  // position of a Matrix of the _excops in term, if -1: the _excops is not present in this term
  int _posexcopsterm;
};
/*
 *  excitation map
 */
class LExcitationMap : public std::map<std::string, LExcitationInfo> {
public:
  // get excitation info corresponding to excitation index in name.
  // if not there yet, handle excitation index and add to the map. 
  // Return iterator to this excitation
  LExcitationMap::iterator get_add( std::string const & name, int lmel = 0 );
  
  // set lastorbs in globalterm (if smaller)
  void set_lastorbs(const Product< Orbital >& orbs, Spin::Type spintype);
  // correct used orbitals
  void correct_orbs(const Product<Orbital>& orbs);
  // term to get lastorbs
  const Term& orbsterm() const { return _globalterm;};
private:  
  // global "term" to generate orbitals for excitations 
  Term _globalterm;
};
/* 
 * name parsing
   handle name, ups and downs of operators 
   (name, excitation class, name additions, dagger, non-conserving character, orbital types...)
 */
struct LParsedName {
  enum Type {
    Name = 0x000,
    Lmel = 0x001,
    Nameadd = 0x002,
    Dg = 0x004,
    Orbs = 0x008,
    Excl = 0x01,
    Excitation = 0x02,
    Orbtypes = 0x04
  };
  int lmel;
  std::string name, nameadd, excitation;
  bool foundsscipt, dg;
  //occ (superscript) and virt (subscript) orbitals
  Product<Orbital> occ,virt;
  short int excl;
  std::vector<OrbitalTypes> orbtypes;
  Product<Orbital> orbs;
  LParsedName() : lmel(0),dg(false),excl(0){};
  // parse namein for name, subscript and superscript
  // if try2set=Name - stop after setting the name
  LParsedName( const std::string& namein, uint try2set );
private:
  void parse_superscript( const std::string& up,  uint try2set );
  void parse_subscript( const std::string& down,  uint try2set );
  bool gen_orbtypes(const std::string& string);
};

/*!
    lexic Equation
 */
class LEquation {
public:
  // add string
  LEquation & operator += (const Lelem & lel) { _eqn.add(lel); return *this; };
  // clear _eqn
  void reseteq() { _eqn = LelString(); _connections = ConnectionsMap(); };
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
  // handle sum
  void handle_sum(Lelem const & lel, Term & term) const;
  // handle parameter
  void handle_parameters(Term& term, bool excopsonly = false);
  // handle permutation
  Permut handle_permutation(Lelem const & lel) const;
  // reset term (set to Term() and reset lastorbs)
  void reset_term(Term& term) const;
  
  LelString _eqn;
  Equation _sumterms;
  // save names and corresponding info of pure excitation and deexciation operators
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
