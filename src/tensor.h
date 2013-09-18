#ifndef Tensor_H
#define Tensor_H

#include <string>
#include <set>
#include <iostream>
#include <assert.h>
#include "globals.h"
#include "types.h"
#include "product.h"
#include "orbital.h"
#include "inpline.h" // for name-handling
#include "sum.h"

typedef uint64_t Length;

class SlotType {
public:
  SlotType(Length n) : nIndices(n) {};
private:
  Length nIndices;
};

typedef std::list<SlotType> SlotTypes;
typedef std::vector<uint> Slots;

class Cut {
  enum Type {
    NoType = 0,
    // cut according to a list of indices (pairlist, tripleslist)
    List = 1,
    // cut according to a set of other indices (pairdomains etc)
    Domain = 2,
    // forced triangular mode (neglect all non-triangular blocks)
    Triang = 10
  };
  // cut strength
  enum Strength{
    DefSth = 0,
    // distant (e.g. strong+close+weak+distant pairs)
    Dist = 1,
    // weak (e.g. strong+close+weak pairs)
    Weak = 2,
    // close (e.g. strong+close pairs)
    Close = 3,
    // strong (e.g. strong pairs only)
    Strong = 4,
    // no strength, needed for search
    NoSth = 100
  };
  Cut::Type cutType;
  Cut::Strength 
      cutStrength,
      // strength of the parent defining cut
      PDCStrength;
  // type of the cut ( 1 = single; 2 = pair; 3 = triple ... )
  // eq. or larger (for united domains/lists) than number of cut (defining) slots in cutSlots and defSlots.
  uint nSlotCutType;
  Slots cutSlots;
  Slots defSlots;
};

//slot types
typedef std::vector<const SlotType*> SlotTs;
typedef std::vector<Cut> Cuts;

class Tensor {
public:
  Tensor( const SlotTs& slots_ ) : slots(slots_) {};
  Tensor( const SlotTs& slots_, const Cuts& cuts_ ) : slots(slots_), cuts(cuts_) {};
  bool operator == ( const Tensor& Ten ) const;
private:
  SlotTs slots;
  Cuts cuts;
};

typedef double Factor;
typedef double Cost;

// R = fac * A B
class Contraction {
public:
  Cost cost() const;
  // will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost ) const;
private:
  Factor fac;
  // lists same slots in tensors
  Slots 
    AinB, BinA,
    AinR, RinA,
    BinR, RinB;
  Tensor *pA, *pB, *pR;
  
};

// R = \sum_i fac_i A_i
class Summation {
  
typedef std::vector<Factor> Factor4Ten;
typedef std::vector<Slots> Slots4Ten;
typedef std::vector<Tensor*> TensorPointers;

public:
private:
  Factor4Ten fac;
  Slots4Ten
    AinR, RinA;
  Tensor *pR;
  TensorPointers pA;
};

class Expression {
  
};

#endif