#include "expression.h"

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

Cost Diagram::contractionCost( const DiagramTensor& ten1, const DiagramTensor& ten2, const DiagramTensor& res ) const
{
  // find contracted slots
  std::bitset<MAXNINDICES>
    connectionmask = ten1._connect.bitmask&ten2._connect.bitmask;
  assert( (connectionmask&res._connect.bitmask) == 0 );
  assert( res._connect.bitmask == (ten1._connect.bitmask^ten2._connect.bitmask));
  xout << "connectionmask: " << connectionmask << std::endl;

  return 1;
}

//std::vector<Action*> 
void Diagram::binarize(Expression& expr) const
{
  const uint MAXNTENS = 12;
  
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
  xout << "in " << nsteps << " steps" << std::endl;
  bt.reset();
  for ( uint i = 0; i < mat_idx.size(); ++i ) bt[i] = true;
  xout << "recursive cost: " << cost[bt.to_ulong()] << std::endl;
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

// returns a character or a string in curly brackets
static std::string ReadAndAdvance(std::size_t& ipos, const std::string& s)
{
  assert( ipos < s.size() );
  std::size_t first = ipos;
  ++ipos;
  if ( s[first] == '{' ){
    // closing bracket
    ipos = s.find_first_of( '}', first + 1 );
    assert( ipos != std::string::npos );
    ++ipos;
  }
  return s.substr(first,ipos-first);
}

Expression::Expression()
{
  // add default tensors
  const TsPar& deftensors = Input::sPars["tensor"];
  TsPar::const_iterator idt;
  _foreach(idt,deftensors) {
    SlotTs sts;
    std::size_t 
      ibra = idt->first.find_first_of('['),
      iket = idt->first.find_first_of(']');
    // generate tensor name
    std::string name = idt->first.substr(0,ibra);
    if ( ibra == std::string::npos || iket == std::string::npos ){
      if ( ibra != std::string::npos || iket != std::string::npos )
        error("Missmatch in [ ] in the default tensor definition");
    } else {
      // slottypes
      for ( std::size_t ii = ibra+1; ii < iket; ){
        SlotTypes::iterator ist = _slottypes.begin();
        ist = _slottypes.insert( ist, SlotType(ReadAndAdvance(ii,idt->first)) );
        sts.push_back(&(*ist));
      }
    }
    Tensor tens(sts,name);
    tens.CreateCutFromDesc(idt->second);
    _tensors.insert(tens);
  }
}

const SlotType* Expression::add(const SlotType& slottype)
{
  SlotTypes::iterator it = _slottypes.begin();
  it = _slottypes.insert(it,slottype);
  return &(*it);
}

const Tensor* Expression::add(const Tensor& tensor)
{
  TensorsSet::iterator it = _tensors.begin();
  it = _tensors.insert(it,tensor);
  return &(*it);
}

std::string Expression::newname(const Symmetries& syms, const Cuts& cuts)
{
  const std::string& names = Input::sPars["syntax"]["intermeds"];
  if ( _lastname.size() == 0 ) {
    _lastname = names[0];
    return _lastname;
  }
  for ( int i = _lastname.size()-1; i >= 0; --i) {
    uint indx = names.find(_lastname[i]);
    if(indx==std::string::npos)
      error("Something wrong with tensor names","Expression::newname");
    else if ( indx < names.size()-1 ) { // not the last possible name
      _lastname[i]=names[indx+1];
      return _lastname;
    }
    _lastname[i]=names[0];//set to first letter
  }
  //add another letter
  _lastname += names[0];
  return _lastname;

}

std::ostream & operator << (std::ostream& o, const Diagram& d) {
  o << "Diagram: {";
  for ( uint i = 0; i < d._tensors.size(); ++i ){
    o << d._tensors[i].name() << "["<< d._tensors[i].slotTypeLetters(d._slottypes) << "]";
    o << "(" << d._tensors[i]._connect.bitmask << ")";
  }
  o << "}";
  return o;
}

std::ostream & operator << (std::ostream& o, const Expression& exp) {
  o << "---- decl" << std::endl;
  // index spaces
  const SlotTypes& sts = exp.slottypes();
  SlotTypes::const_iterator ists;
  _foreach(ists,sts){
    o << *ists << std::endl;
  }
  // tensors....
  const TensorsSet& ts = exp.tensors();
  TensorsSet::const_iterator its;
  _foreach(its,ts){
    o << "tensor: " << *its << std::endl;
  }
  
  // contractions...
  
  return o;
}
