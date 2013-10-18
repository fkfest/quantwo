#include "action.h"

Contraction::Contraction(const Tensor& a, const Tensor& b, const Tensor& r, 
                         const Slots& AinB, const Slots& BinA, 
                         const Slots& AinR, const Slots& RinA, const Slots& BinR, const Slots& RinB, const Factor& fac) :
  p_A(&a), p_B(&b), p_R(&r), _fac(fac), _AinB(AinB), _BinA(BinA), _AinR(AinR), _RinA(RinA), _BinR(BinR), _RinB(RinB), _cost(-1)
{
  assert( _AinB.size() == _BinA.size() );
  assert( _AinR.size() == _RinA.size() );
  assert( _BinR.size() == _RinB.size() );
  assert( _AinB.size()+_AinR.size() >= p_A->slots().size() );
  assert( _BinA.size()+_BinR.size() >= p_B->slots().size() );
  
}


Cost Contraction::cost(Cost mincost)
{
  if ( _cost >= 0 ) return _cost;
  _cost = 0;
  if ( p_B == 0){
    if ( p_A ){
      // scaling by a factor: cost = number of elements in A
    }
  } else {
    assert( p_A );
    // multiplication: cost = nA*nL*nB (take care of locality!)
  }
  
  return _cost;
}

Cost Summation::cost(Cost mincost)
{
  Cost cst = 0;
  return cst;
}
