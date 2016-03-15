/// \brief Creation and comparison of "unique" graphs (diagrams) 
 
#ifndef UNIGRAPH_H
#define UNIGRAPH_H

#include <vector>
#include <list>
#include "types.h"
#include "globals.h"
#include "term.h"
#include "evertices.h"
#include "matrix.h"

class Permutation : public std::map<uint,uint> {
  typedef std::map<uint,uint> Base;
public: 
  using Base::Base;
  bool operator < (const Permutation& perm) const { 
    if ( this->size() < perm.size() ) return true;
    if ( perm.size() < this->size() ) return false;
    return static_cast<const Base&>(*this) < static_cast<const Base&>(perm);
  };
};

/// The graphs have a canonical order of operators (neglecting orbital names),
/// sets of equivalent vertices for permutational symmetry handling, 
/// and a connection vector (each element is a pointer to the next vertex connected to the
/// current vertex by an incomming line (i.e., a creator for the current vertex))
class UniGraph {
public:
  UniGraph(const Term& term);
  // gen term
  Term gen_term();
  // ordered matrices
  Product<Matrix> ordmats() const;
  const Order& connections() const { return _vertconn;};
  const Equivalents& equivals() const { return _equivs;};
  const PermVertices& eqperms() const { return _eqperms;};
  const PermVertices& eqperm_from() const { return _eqperm_from;};
  // check equality of two graphs
  bool is_equal(const UniGraph& ug) const;
  // search for the minimal _vertconn using _equivs and _eqperms
  void minimize();
  // generate permutation using orbitals from ug (have to be generated before!)
  std::pair<Permut,TFactor> permutation( const UniGraph& ug ) const; 
private:
  // apply permutations to minimize the connections-vector further
  // calls gen_perms to set _perms
  void apply_eqperms(Order& connections, const Order& vertorder,
                     PermVertices& ep, PermVertices& epf, PermVertices& epfo);
  // set _perms
  void gen_perms(const PermVertices& from, const PermVertices& to);
  // canonical order of matrices
  Order _matsord;
  // vertices connections v[0]=1 --> a line from vertex 0 to 1
  Order _vertconn;
  // additional permutations to get to vertconn
  Permutation _perms;
  // store orbital types according to vertconn
  OrbitalTypes _orbtypes;
  // store external orbitals
  Product<Orbital> _extorbs_crea, _extorbs_anni;
  // store new orbitals (for generation of permutations)
  Array<Orbital> _neworbs;
  // external vertices
  JointVertices _extvertices;
  // vectors of equivalent (joint) vertices for permutational symmetry
  // T_2 (0,1) T_2 (2,3) --> [ (0,1)(2,3) ][ (0)(1) ][ (2)(3) ] 
  Equivalents _equivs;
  // vertices for allowed permutations (e.g. external lines etc) and vertices connected to those 
  // (needed because of non-conserving stuff)
  PermVertices _eqperms, _eqperm_from;
  const Term * pTerm;
};

std::ostream & operator << (std::ostream & o, const UniGraph& ug);

#endif