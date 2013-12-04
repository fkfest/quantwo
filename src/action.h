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
  Contraction() : p_A(0), p_B(0), p_R(0), _fac(1), _cost(-1){};
  Contraction( const Tensor& a, const Tensor& b, const Tensor& r, 
               const Slots& AinB, const Slots& BinA,
               const Slots& AinR, const Slots& RinA,
               const Slots& BinR, const Slots& RinB,
               const Factor& fac = 1 );
  ~Contraction(){};
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
//private:
  const Tensor *p_A, *p_B, *p_R;
  Factor _fac;
  // lists same slots in tensors
  Slots 
    _AinB, _BinA,
    _AinR, _RinA,
    _BinR, _RinB;
  // save cost
  Cost _cost;
};

typedef std::set<Tensor> TensorsSet;
typedef Array<const Tensor*> TensorPointers;

// R = \sum_i fac_i A_i
class Summation : public Action {
public:
  typedef Array<Factor> Factor4Ten;
  typedef Array<Slots> Slots4Ten;
  Summation() : p_R(0) {};
  Summation( const TensorPointers& pA, const Tensor& r,
             const Slots4Ten& AinR, const Slots4Ten& RinA,
             const Factor4Ten& facs ) :
             p_A(pA), p_R(&r), _facs(facs),
             _AinR(AinR), _RinA(RinA) {};
  ~Summation(){};
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
//private:
  TensorPointers p_A;
  const Tensor * p_R;
  Factor4Ten _facs;
  Slots4Ten
    _AinR, _RinA;
  // save cost
  Cost _cost;
};




#endif