#include "action.h"

Contraction::Contraction(const Tensor& a, const Tensor& b, //const Tensor& r,
                         const Slots& AinB, const Slots& BinA,
                         const Slots& AinR, const Slots& RinA, const Slots& BinR, const Slots& RinB, const Factor& fac) :
  p_A(&a), p_B(&b), //p_R(&r),
  _fac(fac), _AinB(AinB), _BinA(BinA), _AinR(AinR), _RinA(RinA), _BinR(BinR), _RinB(RinB), _cost(-1)
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
      if ( std::abs(_fac - 1) > Numbers::verysmall ) {
        // scaling by a factor: cost = number of elements in A
      }
    }
  } else {
    assert( p_A );
    // multiplication: cost = nA*nL*nB (take care of locality!)
    if ( _cost > mincost ) {
      // do less work, we don't need to calculate the whole cost
      // if at some point we know that the cost is larger than the previous minimal cost...
    }
  }
  return _cost;
}

void slotNames4Refs(Array<std::string>& xslotnames, Array<std::string>& yslotnames,
                           std::map<SlotType::Type,std::string>& oldnames, const Slots& sXinY, const Slots& sYinX,
                           const SlotTs& xslottypes, const SlotTs& yslottypes )
{
  assert( sXinY.size() == sYinX.size() );
  assert( xslotnames.size() == xslottypes.size() );
  assert( yslotnames.size() == yslottypes.size() );
  for ( uint iSt = 0; iSt < sXinY.size(); ++iSt ){
    uint
      iSlotX = sXinY[iSt],
      iSlotY = sYinX[iSt];
    assert( iSlotX < xslotnames.size() );
    assert( iSlotY < yslotnames.size() );
    assert( xslottypes[iSlotX] == yslottypes[iSlotY] );
    if ( xslotnames[iSlotX] == "" ){
      if ( yslotnames[iSlotY] == "" ){
        const SlotType * pslott = xslottypes[iSlotX];
        xslotnames[iSlotX] = yslotnames[iSlotY] = oldnames[pslott->type()] = pslott->name(oldnames[pslott->type()]);
      } else {
        xslotnames[iSlotX] = yslotnames[iSlotY];
      }
    } else if ( yslotnames[iSlotY] == "" ){
      yslotnames[iSlotY] = xslotnames[iSlotX];
    }
    assert( xslotnames[iSlotX] == yslotnames[iSlotY] );
  }
}
void fillFreeSlotNames(Array<std::string>& xslotnames, std::map<SlotType::Type,std::string>& oldnames,
                              const Tensor& xten)
{
  const SlotTs& xslottypes = xten.slots();
  assert( xslotnames.size() == xslottypes.size() );
  for ( uint iSlotX = 0; iSlotX < xslotnames.size(); ++iSlotX ){
    if ( xslotnames[iSlotX] == "" ){
      // has to be a phantom slot!
      assert( xten.phantomSlot(iSlotX) );
      const SlotType * pslott = xslottypes[iSlotX];
      xslotnames[iSlotX] = oldnames[pslott->type()] = pslott->name(oldnames[pslott->type()]);
    }
  }
}

std::string Contraction::fingerprint(const Tensor& ten) const
{
  std::string fingerprint;
  for(auto elem : _RinA){
    fingerprint += std::to_string(elem);}
  for(auto elem : _AinR){
    fingerprint += std::to_string(elem);}
  for(auto elem : _RinB){
    fingerprint += std::to_string(elem);}
  for(auto elem : _BinR){
    fingerprint += std::to_string(elem);}
  for(auto elem : _AinB){
    fingerprint += std::to_string(elem);}
  for(auto elem : _BinA){
    fingerprint += std::to_string(elem);}
  for(auto elem : ten.slots()){
    fingerprint += std::to_string(elem->type());}
  for(auto elem : p_A->slots()){
    fingerprint += std::to_string(elem->type());}
  return fingerprint;
}

void Contraction::print(std::ostream& o, const Tensor& res) const
{
  std::map<SlotType::Type,std::string> slotnames;
  Array<std::string>
      resslots(res.slots().size()),
      aslots(p_A->slots().size()),
      bslots(p_B->slots().size());

    slotNames4Refs(resslots,aslots,slotnames,_RinA,_AinR,res.slots(),p_A->slots());
    slotNames4Refs(resslots,bslots,slotnames,_RinB,_BinR,res.slots(),p_B->slots());
    slotNames4Refs(aslots,bslots,slotnames,_AinB,_BinA,p_A->slots(),p_B->slots());

    fillFreeSlotNames(resslots,slotnames,res);
    fillFreeSlotNames(aslots,slotnames,*p_A);
    fillFreeSlotNames(bslots,slotnames,*p_B);

    if( ! (_printed.find(res.name()) != _printed.end())){
      o << "." << res.name() << "[" << resslots << "] ";

      if ( _fac < 0 )
        o << "-= ";
      else
        o << "+= ";
      if ( std::abs(std::abs(_fac) - 1) > Numbers::verysmall ) o << std::abs(_fac) << "*";
      o << p_A->name() << "[" << aslots << "] ";
      o << p_B->name() << "[" << bslots << "]";
      o << std::endl;
    }
    if ( res.type() == "A" ) _printed.insert(res.name());
}

Cost Summation::cost(Cost mincost)
{
  Cost cst = 0;
  return cst;
}

void Summation::print(std::ostream& o, const Tensor& res) const
{
  std::map<SlotType::Type,std::string> slotnames;
  Array<std::string> 
      resslots(res.slots().size()),
      aslots(p_A->slots().size());

  slotNames4Refs(resslots,aslots,slotnames,_RinA,_AinR,res.slots(),p_A->slots());
  fillFreeSlotNames(resslots,slotnames,res);

  std::map<SlotType::Type,std::string> slotnamesA(slotnames);
  fillFreeSlotNames(aslots,slotnamesA,*(p_A));

  o << "." << res.name() << "[" << resslots << "] ";
  if ( _fac < 0 )
    o << "-= ";
  else
    o << "+= ";
  if ( std::abs(std::abs(_fac) - 1) > Numbers::verysmall ) o << std::abs(_fac) << " ";
  o << p_A->name() << "[" << aslots << "] " << std::endl;
}