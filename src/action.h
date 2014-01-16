#ifndef Action_H
#define Action_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <stdint.h>
#include "globals.h"
#include "utilities.h"
#include "types.h"
#include "tensor.h"

// base class for actions on the tensors
class Action {
public:
  Action(){};
  virtual ~Action(){};
  virtual Cost cost( Cost mincost ) = 0;
};


// R = fac * A B
class Contraction : public Action {
public:
  Contraction() : p_A(0), p_B(0), //p_R(0), 
                  _fac(1), _cost(-1){};
  Contraction( const Tensor& a, const Tensor& b, //const Tensor& r, 
               const Slots& AinB, const Slots& BinA,
               const Slots& AinR, const Slots& RinA,
               const Slots& BinR, const Slots& RinB,
               const Factor& fac = 1 );
  ~Contraction(){};
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
  void print( std::ostream& o, const Tensor& res ) const;
//private:
  const Tensor *p_A, *p_B; //, *p_R;
  Factor _fac;
  // lists same slots in tensors, e.g., _AinB: slots in tensor A that will be contracted with tensor B 
  Slots 
    _AinB, _BinA,
    _AinR, _RinA,
    _BinR, _RinB;
  // save cost
  Cost _cost;
};

typedef std::list<Tensor> TensorsSet;
typedef Array<Slots> Slots4Ten;
typedef Array<Actions> ActionsSets;

// R = \sum_i fac_i "Contraction"
class Summation : public Action {
public:
  Summation() {};
  Summation( const ActionsSets& pActs );
  ~Summation(){};
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
//private:
  ActionsSets _pActs;
  // save cost
  Cost _cost;
};




#endif