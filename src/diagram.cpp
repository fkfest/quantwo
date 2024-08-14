#include "diagram.h"

const DiagramTensor* Diagram::add(DiagramTensor dten, const Tensor* pTen, bool pushfront)
{
  if (pTen){
    if (pTen->_cuts.size() > 0) error( "add cuts handling in Diagram", "Diagram::add");
    dten._name = pTen->_name;
    dten._syms = pTen->_syms;
    dten._phantomSlots = pTen->_phantomSlots;
  }
  if ( pushfront){
    //residual tensor
    if( pTen->slots().size() > 4 ){
      dten._connect.slotref[0] = 5;
      dten._connect.slotref[1] = 4;
      dten._connect.slotref[2] = 3;
      dten._connect.slotref[3] = 2;
      dten._connect.slotref[4] = 1;
      dten._connect.slotref[5] = 0;
    }
    else if( pTen->slots().size() > 2 ){
      dten._connect.slotref[0] = 3;
      dten._connect.slotref[1] = 2;
      dten._connect.slotref[2] = 1;
      dten._connect.slotref[3] = 0;
    }
    else{
      dten._connect.slotref[0] = 1;
      dten._connect.slotref[1] = 0;
    }
    _tensors.push_front(dten);
    return &(_tensors.front());
  }
  else{
    _tensors.push_back(dten);
    return &(_tensors.back());
  }
}

DiagramTensor Diagram::newTensor( const DiagramTensor& ten1, const DiagramTensor& ten2, std::string name ) const
{
  DiagramTensor dT(name);
//  std::bitset<MAXNINDICES>
//    connectionmask = ten1._connect.bitmask&ten2._connect.bitmask;
  dT._connect.bitmask = ten1._connect.bitmask^ten2._connect.bitmask;

// set slotref
// the slottypes in _slottypes are in canonical order, therefore we
// simply have a reverse order of slots here.
  uint nSlots = dT._connect.bitmask.count();
  dT._connect.slotref.resize(nSlots);
  for ( uint ist = 0; ist < nSlots; ++ist ){
    dT._connect.slotref[nSlots-ist-1] = ist;
  }
// TODO: set symmetry. may be also cuts?

  return dT;
}

Tensor Diagram::exprTensor(const DiagramTensor& ten) const
{
  SlotTs slots;
  std::bitset<MAXNINDICES> bitm = ten._connect.bitmask;
  slots.resize(bitm.count());
  for ( uint ipos = 0, ist = 0; ist < slots.size() && ipos < _slottypes.size(); ++ipos, bitm >>= 1 ){
    if ( bitm[0] ) {
      uint iSlot = ten._connect.slotref[bitm.count()-1];
      slots[iSlot] = _slottypes[ipos];
      ++ist;
    }
  }
  Cuts cuts;
  if ( _cuts.size() > 0 )
    error("generate cuts for expr Tensor");
  std::string name(ten._name);
  if(ten._connect.bitmask == _tensors[0]._connect.bitmask){
    if( ten.type() == "f" ) name ="f";
    else if( ten.type() != "I" ) name ="R";
  }
  return Tensor( slots, ten._syms, cuts, name );
}

bool Diagram::isresidual(const DiagramTensor& dten) const
{
  assert(dten._connect.bitmask.count() < 7);
  if (dten._connect.bitmask.count() == 2){
  assert(dten._connect.bitmask.count() < 5);
    if(_tensors[0]._connect.bitmask.count() > 2){//double and single residual in _tensors
      if(dten._connect.bitmask == _tensors[1]._connect.bitmask ) return true;
      else return false;}
    else{//only single residual in _tensors (and therefore at first place)
      if(dten._connect.bitmask == _tensors[0]._connect.bitmask ) return true;
      else return false;
    }
  }
  else if(dten._connect.bitmask.count() == 4){
    if(dten._connect.bitmask == _tensors[0]._connect.bitmask ) return true;
    else return false;
  }
  else if(dten._connect.bitmask.count() == 6){
    if(dten._connect.bitmask == _tensors[0]._connect.bitmask ) return true;
    else return false;
  }
  else{
    error("Something is wrong in Diagram::isresidual()");
    return 0;
  }
}

const Tensor * Diagram::transform2Expr(Expression& expr, const Array< DiagramTensor >& inters, const Array< std::bitset< MAXNTENS > >& order,
                             std::bitset< MAXNTENS > bt ) const
{
  unsigned long ibt = bt.to_ulong();
  const Action * pAct = 0;
  if ( bt.count() > 1 ) {
    // get parents
    std::bitset<MAXNTENS>
        bt1 = order[ibt],
        bt2 = bt^bt1;
    const Tensor
      * pTen1 = transform2Expr(expr,inters,order,bt1),
      * pTen2 = transform2Expr(expr,inters,order,bt2);
    unsigned long
        ibt1 = bt1.to_ulong(),
        ibt2 = bt2.to_ulong();
    // action...
    Contraction contr = exprContraction(inters[ibt1],inters[ibt2],inters[ibt],pTen1,pTen2);
    pAct = expr.add(&contr);
  }
  assert( bt.count() > 0 );
  // add the intermediate
  Tensor ten(exprTensor(inters[ibt]));
  if(inters[ibt].name() == "" && isresidual(inters[ibt])){
    return expr.add2residual(ten,pAct);
  }
  else{
    ten.add(pAct);
    return expr.add(ten);
  }
}

void Diagram::binarize(Expression& expr) const
{
  uint nmats = _tensors.size();
  // first tensor is the residual tensor
  if ( nmats > 0 ) --nmats;
  assert( nmats <= MAXNTENS );
  // Dynamical programming: First we calculate the cost of simplest intermediates, then of
  // more complex etc. and save the cheapest way and cost of their evaluation
  Array<uint> mat_idx;
  for ( uint i = 0; i < nmats; ++i ) mat_idx.push_back(i);
  Array<Cost> cost;
  Array< std::bitset<MAXNTENS> > order;
  Array<DiagramTensor> inters;
  cost.resize(1<<nmats);
  order.resize(1<<nmats);
  inters.resize(1<<nmats);

  std::bitset<MAXNTENS> bt, bt1, bt2;
  for ( uint i = 0; i < mat_idx.size(); ++i ){
    bt.reset();
    bt[mat_idx[i]] = true;
    cost[bt.to_ulong()] = 0;
    order[bt.to_ulong()] = 0;
    assert( mat_idx[i] < _tensors.size() );
    // copy tensors, first tensor in _tensors is the residual tensor
    inters[bt.to_ulong()] = _tensors[mat_idx[i]+1];
  }
  uint nsteps=0;
  // number of simple tensors in intermediate tensors
  for ( uint L = 2; L <= mat_idx.size(); ++L ){
    Array<uint>::iterator itmat = mat_idx.begin() + L;
    // consider all possible intermediate tensors of size L
    do {
      Array<uint>::iterator
        itmbeg = mat_idx.begin(),
        itmend = mat_idx.begin(),
        itm;
      ++itmend;
      bt.reset();
      for ( uint i = 0; i < L; ++i ) bt[mat_idx[i]] = true;
      for ( itm = mat_idx.begin(), bt1.reset(), bt2 = bt; itm != itmend; ++itm ) {
        bt1[*itm] = true;
        bt2[*itm] = false;
      }
      unsigned long
        ibt = bt.to_ulong(),
        ibt1 = bt1.to_ulong(),
        ibt2 = bt2.to_ulong();
      cost[ibt] = MAXCOST;
      inters[ibt] = newTensor(inters[ibt1],inters[ibt2]);
      // number of simple tensors in the first parent of the current intermediate tensor
      for ( uint nn = 1; nn <= L/2; ++nn, ++itmend ){
        // removes the rest redundancy for even L
        if ( nn*2 == L ) ++itmbeg;
        // consider all possible "first" parents
        do {
//          xout << Array<uint>(mat_idx.begin(),itmend) << " x " << Array<uint>(itmend,itmat) << std::endl;
          for ( itm = mat_idx.begin(), bt1.reset(), bt2 = bt; itm != itmend; ++itm ) {
            bt1[*itm] = true;
            bt2[*itm] = false;
          }
          ibt1 = bt1.to_ulong(),
          ibt2 = bt2.to_ulong();
          Cost qq = cost[ibt1] + cost[ibt2];
          if ( qq < cost[ibt] ) {
            qq += contractionCost(inters[ibt1],inters[ibt2],inters[ibt]);
            if ( qq < cost[ibt] ) {
              cost[ibt] = qq;
              order[ibt] = bt1;
            }
          }
          ++nsteps;
        } while(next_combination(itmbeg,itmend,itmat));
      }
    } while(next_combination(mat_idx.begin(),itmat,mat_idx.end()));
  }
  assert( _tensors[0]._connect.bitmask == inters[bt.to_ulong()]._connect.bitmask );
  // residual tensor
  Tensor res = exprTensor(_tensors[0]);
  if(inters.size() == 2){//we have R=a*A
    Tensor ten(exprTensor(inters[1]));
    const Tensor *pTen = expr.add(ten);
    Summation sum = exprSummation(inters[1],_tensors[0],pTen);
    const Action * pAct = expr.add(&sum);
    expr.add2residual(res,pAct);
  }
  else{//we have R=a*A*B*..
    //recursive calls of transform2Expr inside depending on "relations" of residual tensor bt
    transform2Expr(expr,inters,order,bt);
  }
  for ( TensorsList::iterator it = expr._tensors.begin(); it != expr._tensors.end(); ++it) {
    if ( it->equal(res) ) {
      const Tensor * pRes = &(*it);
      expr.addresidual(pRes);
    }
  }
//  LOOP STRUCTURE USING bitset next_combination
//  nsteps = 0;
//  std::bitset<MAXNTENS> xx, xxt;
//  Array<uint> tt1,tt2;
//  for ( uint L = 2; L <= nmats; ++L ){
//    for ( uint i = nmats-L; i < nmats; ++i ) xx[i] = true;
//    do {
//      tt.clear();
//      for ( uint i = 0; i < nmats; ++i )
//        if ( xx[i] ) tt.push_back(i);
//      for ( uint nn = 1; nn <= L/2; ++nn ){
//        uint end = L;
//        if ( nn*2 == L ) --end;
//        for ( uint i = L-nn; i < end; ++i ) xxt[i] = true;
//        do {
//          tt1.clear();
//          tt2.clear();
//          for ( uint i = 0; i < end; ++i )
//            if ( xxt[i] )
//              tt1.push_back(tt[i]);
//            else
//              tt2.push_back(tt[i]);
//          if ( end != L ) tt1.push_back(tt[end]);
//          ++nsteps;
//        } while(next_combination(xxt));
//      }
//    } while(next_combination(xx));
//  }

}

bool DiagramTensor::equalestimate( const DiagramTensor& ten ) const{
  uint diff = 0;
  for( uint i = 0; i < this->_connect.slotref.size(); ++i ){
    if( this->_connect.slotref[i] != ten._connect.slotref[i] ) diff++;
  }
  if( diff > 2 ) return false;
  if((this->_connect.bitmask ^ ten._connect.bitmask).count() > 4) return false;
  return true;
}

bool Diagram::equal(const Diagram& diag) const{
  if( this->_tensors.size() != diag._tensors.size() ) return false;
  if( this->_slottypes.size() != diag._slottypes.size() ) return false;
  //assuming _tensors[1] is electron integral
  if( this->_tensors[1]._connect.bitmask != diag._tensors[1]._connect.bitmask ) return false;
  if( this->_tensors[1]._connect.slotref != diag._tensors[1]._connect.slotref ) return false;
  for( uint i = 0; i < this->_tensors.size(); ++i ){
    if(_tensors[i].slotTypeLetters(this->_slottypes) != diag._tensors[i].slotTypeLetters(diag._slottypes)) return false;
  }
  for ( uint i = 0; i < this->_tensors.size(); ++i ){
    if( !this->_tensors[i].equalestimate(diag._tensors[i]) ) return false;
  }
  return true;
}

void Diagram::createPermMap(const Array<std::string>& aslots, const Array<std::string>& bslots){
  assert(aslots.size() == bslots.size());                                           
  for(uint i = 0; i < aslots.size(); i++){                                        
    if( aslots[i] != bslots[i] ){                                                   
      _permmap[aslots[i]] = bslots[i];                                                 
    }                                                                           
  }  
}

void Diagram::permute(Array<std::string>& slots){
  for( uint i = 0; i < slots.size(); i++ ){                                  
    PerMap::iterator it = _permmap.find(slots[i]);                             
    if ( it != _permmap.end()){                                                    
      slots[i] = _permmap[slots[i]];                                         
    }                                                                           
  }             
}

std::ostream & operator << (std::ostream& o, const Diagram& d) {
  o << "Diagram: {";
  o << d._fac;
  for ( uint i = 0; i < d._tensors.size(); ++i ){
    o << d._tensors[i].name() << "["<< d._tensors[i].slotTypeLetters(d._slottypes) << "]";
    o << "(" << d._tensors[i]._connect.bitmask << ")";
    o << "(" << d._tensors[i]._connect.slotref << ")";
  }
  o << "}";
  return o;
}
