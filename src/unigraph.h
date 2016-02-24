/// \brief Creation and comparison of "unique" graphs (diagrams) 
 
#ifndef UNIGRAPH_H
#define UNIGRAPH_H

#include <vector>
#include <list>
#include "types.h"
#include "globals.h"
#include "term.h"
#include "evertices.h"
#include "matrices.h"

/// The graphs have a canonical order of operators (neglecting orbital names),
/// sets of equivalent vertices for permutational symmetry handling, 
/// and a connection vector (each element is a pointer to the next vertex connected to the
/// current vertex by an incomming line (i.e., a creator for the current vertex))
class UniGraph {
public:
  UniGraph(const Term& term);
  // ordered matrices
  Product<Matrices> ordmats() const;
  const Order& connections() const { return _vertconn;};
  const Equivalents& equivals() const { return _equivs;};
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
  typedef std::vector<Permutation> Permutations;
  // vectors of equivalent (joint) vertices for permutational symmetry
  // T_2 (0,1) T_2 (2,3) --> [ (0,1)(2,3) ][ (0)(1) ][ (2)(3) ] 
  Equivalents _equivs;
  // allowed permutations (e.g. external lines etc)
  Permutations _eqperms;
  const Term * pTerm;
};

std::ostream & operator << (std::ostream & o, const UniGraph& ug);

#endif