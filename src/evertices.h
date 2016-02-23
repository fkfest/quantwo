/// Equivalence of vertices for UniGraph
#ifndef EVERTICES_H
#define EVERTICES_H

#include <vector>
#include <list>
#include <assert.h>
#include "types.h"
#include "globals.h"
#include "utilities.h"

typedef std::vector<uint> JointVertices;
struct EquiVertices : public std::vector<JointVertices> {
  EquiVertices() : std::vector<JointVertices>(){};
  void add( const JointVertices& jv ) { 
    assert(size()==0 || begin()->size() == jv.size()); 
    push_back(jv); };
  void add( uint vert ) {
    assert(size()==0 || begin()->size() == 1 );
    push_back(JointVertices(1,vert));
  };
};

/// T_2 (0,1) T_2 (2,3) --> [ (0,1)(2,3) ][ (0)(1) ][ (2)(3) ] 
typedef std::list<EquiVertices> Equivalents;

std::ostream & operator << (std::ostream & o, const EquiVertices& ev);
std::ostream & operator << (std::ostream & o, const Equivalents& eqv);
#endif