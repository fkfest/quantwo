/// Equivalence of vertices for UniGraph
#ifndef EVERTICES_H
#define EVERTICES_H

#include <vector>
#include <list>
#include <assert.h>
#include "types.h"
#include "arrays.h"
#include "globals.h"
#include "utilities.h"

typedef Array<uint> JointVertices;
struct EquiVertices : public std::vector<JointVertices> {
  EquiVertices() : std::vector<JointVertices>(){};
  void add( const JointVertices& jv ) { 
    assert(size()==0 || begin()->size() == jv.size()); 
    push_back(jv); };
  void add( uint vert ) {
    assert(size()==0 || begin()->size() == 1 );
    push_back(JointVertices(1,vert));
  };
  // next permutation in vord using JointVertices
  // relies on the fact that vord[j->front()] < vord[i->front()] is a valid check
  bool next_permutation( Order& vord )
  {
    assert( size() > 1 );
    const_iterator i = begin(), j;
    uint vordi = vord[i->front()];
    for ( ++i ; i != end(); ++i ) {
      uint vordj = vordi;
      vordi = vord[i->front()];
      if ( vordj < vordi ) break;
    }
    if ( i != end() ) { 
      for ( j = begin(); vordi <= vord[j->front()]; ++j ) {}
      swap(vord,i,j);
      reverse(vord,begin(),i);
      return true;
    }
    reverse(vord,begin(),end());
    return false;
  };
  // swap vertices in vord from JointVertices sets it and jt 
  void swap( Order& vord, const_iterator it, const_iterator jt )
  {
    for ( uint i = 0; i < it->size(); ++i )
      std::swap(vord[(*it)[i]],vord[(*jt)[i]]);
  };
  // reverse vertices in vord from JointVertices first to last
  void reverse( Order& vord, const_iterator first, const_iterator last )
  {
    while ((first!=last)&&(first!=--last)) {
      swap (vord,first,last);
      ++first;
    }
  };

};

/// T_2 (0,1) T_2 (2,3) --> [ (0,1)(2,3) ][ (0)(1) ][ (2)(3) ] 
typedef std::vector<EquiVertices> Equivalents;

struct PermVertices : public std::vector<JointVertices> {
  PermVertices(): std::vector<JointVertices>(){};
  // this = ord(pvert)
  PermVertices(const PermVertices& pvert, const Order& ord): std::vector<JointVertices>() {
    this->resize(pvert.size());
    PermVertices::iterator itpv = this->begin();
    for( const JointVertices& pvs: pvert){
      itpv->reserve(pvs.size());
      for ( const uint& pv: pvs)
        itpv->push_back(ord[pv]);
      ++itpv;
    }
  };
  // this = ord(pvert) (this should have proper sizes already)
  void set_with_order(const PermVertices& pvert, const Order& ord) {
    assert( this->size() == pvert.size() );
    PermVertices::iterator itpv = this->begin();
    for( const JointVertices& pvs: pvert){
      assert( itpv->size() == pvs.size() );
      JointVertices::iterator it = itpv->begin();
      for ( const uint& pv: pvs) {
        *it = ord[pv];
        ++it;
      }
      ++itpv;
    }
  }
  void sort() {
    for( JointVertices& pvs: *this){
      std::sort(pvs.begin(),pvs.end());
    }
  }
  void minpermute(Order& connections, const PermVertices& orig) {
    assert( this->size() == orig.size() );
    PermVertices::const_iterator itpvso = orig.begin();
    Order temp(connections.size());
    for( JointVertices& pvs: *this){
      // resort
      InsertionSort(&connections.front(),&pvs.front(),pvs.size());
      assert(pvs.size() == itpvso->size());
      for ( uint i = 0; i < pvs.size(); ++i ) {
        temp[i] = connections[pvs[i]];
      }
      for ( uint i = 0; i < itpvso->size(); ++i ) {
        connections[(*itpvso)[i]] = temp[i];
      }
      ++itpvso;
    }
  };
};

std::ostream & operator << (std::ostream & o, const EquiVertices& ev);
std::ostream & operator << (std::ostream & o, const Equivalents& eqv);
#endif