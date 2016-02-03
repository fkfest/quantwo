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
#include "operators.h"
#include "globals.h"
#include "lexic.h"

/*!
    Equation
 */
class Equation : public Sum<Term, TFactor> {
public:
  
};

typedef std::map< Orbital::Type,Orbital > TOrb4Type;
class LExcitationMap;
/*
 *  excitation info (for \mu_i etc) 
 */
class LExcitationInfo {
  friend class LExcitationMap;
public:
  LExcitationInfo() {};
  LExcitationInfo( const Product<Orbital>& orbs, int lmelec, Matrices::Spinsym spin ) :
    _orbs(orbs), _lmel(lmelec), _spinsymexcs(spin){};
    
  // excitation class
  short exccls(bool dg = false) const{ return (_orbs.size()-lmel(dg))/2; };
  // electron non-conserving number
  int lmel(bool dg = false) const { return dg ? -_lmel : _lmel; };
  // spin symmetry of the excitation
  Matrices::Spinsym spinsymexcs() const { return _spinsymexcs;};
  // return orbitals of the excitation
  Product<Orbital> orbitals(bool dg = false) const;
private:
  // orbitals for the excitation
  // FIXME will replace most of the other entries here!
  Product<Orbital> _orbs;
  // electron non-conserving number
  int _lmel;
  // spinsymmetry
  Matrices::Spinsym _spinsymexcs;
};
/*
 *  excitation map FIXME \mu_1 vs \mu_1^dg = should have the same orbitals?
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
    Lmel = 0x001,       // +-1
    Nameadd = 0x002,    // \snam{blabla} or just blabla
    Dg = 0x004,         // \dg
    Orbs = 0x008,       // ij_1
    Excl = 0x01,        // 1, 2, 3...
    Excitation = 0x02,  // \mu_1, \nu_3
    Orbtypes = 0x04     // ^{ii}_{ta}
  };
  int lmel;
  std::string name, nameadd, excitation;
  bool foundsscipt, dg;
  //occ (superscript) and virt (subscript) orbitals
  Product<Orbital> occ,virt;
  short int excl;
  // orbital types (first set for occ and second set for virt indices )
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
  Oper handle_excitation( std::string const & name, bool dg, int lmel = 0, bool excopsonly=false );
  // handle factor
  TFactor handle_factor(Lelem const & lel) const;
  // handle operator
  Oper handle_operator(Lelem const & lel, Term & term, bool excopsonly=false);
  // handle sum
  Product<Orbital> handle_sum(const Lelem& lel);
  // handle parameter
  Matrices handle_parameter(Lelem const & lel);
  // handle permutation
  Permut handle_permutation(Lelem const & lel) const;
  // reset term (set to Term() and reset lastorbs)
  void reset_term(Term& term) const;
  
  LelString _eqn;
  Equation _sumterms;
  // save names and corresponding info of pure excitation and deexciation operators
  LExcitationMap _excops;
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
