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
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
  void print( std::ostream& o, const Tensor& res ) const;
  std::string fingerprint(const Tensor& ten) const;
//private:
  inline static std::set<std::string> _printed;
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

// R = \sum_i fac_i A_i
class Summation : public Action {
public:
  Summation() : p_A(0), _fac(1), _cost(-1){};
  Summation( const Tensor& a, //const Tensor& r,
               const Slots& AinR, const Slots& RinA, const Factor& fac =1) : p_A(&a), _fac(fac), _AinR(AinR), _RinA(RinA) {};
  // for mincost > 0: will return either actual cost if it's smaller than mincost, or (mincost + 1)
  Cost cost( Cost mincost = -1 );
  void print( std::ostream& o, const Tensor& res ) const;
  const Tensor * p_A;
  Factor _fac;
  Slots _AinR, _RinA;
//private:
  // save cost
  Cost _cost;
};

void slotNames4Refs(Array<std::string>& xslotnames, Array<std::string>& yslotnames,
                           std::map<SlotType::Type,std::string>& oldnames, const Slots& sXinY, const Slots& sYinX,
                           const SlotTs& xslottypes, const SlotTs& yslottypes );

void fillFreeSlotNames(Array<std::string>& xslotnames, std::map<SlotType::Type,std::string>& oldnames,
                              const Tensor& xten);


#endif

