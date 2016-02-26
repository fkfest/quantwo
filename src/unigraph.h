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

/// The graphs have a canonical order of operators (neglecting orbital names),
/// sets of equivalent vertices for permutational symmetry handling, 
/// and a connection vector (each element is a pointer to the next vertex connected to the
/// current vertex by an incomming line (i.e., a creator for the current vertex))
class UniGraph {
public:
  typedef std::vector<JointVertices> PermVertices;
  UniGraph(const Term& term);
  // ordered matrices
  Product<Matrix> ordmats() const;
  const Order& connections() const { return _vertconn;};
  const Equivalents& equivals() const { return _equivs;};
  const PermVertices& eqperms() const { return _eqperms;};
  const PermVertices& eqperm_from() const { return _eqperm_from;};
  // search for the minimal _vertconn using _equivs and _eqperms
  void minimize();
private:
  typedef std::map<uint,uint> Permutation;
  // canonical order of matrices
  Order _matsord;
  // vertices connections v[0]=1 --> a line from vertex 0 to 1
  Order _vertconn;
  // additional permutations to get to vertconn
  Permutation _perms;
  // store orbital types according to vertconn for a safety check
  OrbitalTypes _orbtypes;
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