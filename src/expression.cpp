#include "expression.h"
#include "sum.h"

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

static void setContractionSlots( Slots& sXinY, Slots& sYinX, const DiagramTensor& tenX, const DiagramTensor& tenY ) {
  std::bitset<MAXNINDICES> overlap = tenX._connect.bitmask&tenY._connect.bitmask;

  sXinY.resize(overlap.count());
  sYinX.resize(sXinY.size());
  for ( uint ipos = 1, ist = 0; ist < sXinY.size(); ++ipos, overlap >>= 1 ){
    if ( overlap[0] ) {
      uint iStRef = (tenX._connect.bitmask >> ipos).count();
      sXinY[ist] = tenX._connect.slotref[iStRef];
      iStRef = (tenY._connect.bitmask >> ipos).count();
      sYinX[ist] = tenY._connect.slotref[iStRef];
      ++ist;
    }
  }
}
Contraction Diagram::exprContraction(const DiagramTensor& tenA, const DiagramTensor& tenB, const DiagramTensor& tenR,
                                     const Tensor * pA, const Tensor * pB ) const
{
  Slots
    sAinB, sBinA,
    sAinR, sRinA,
    sBinR, sRinB;
  setContractionSlots(sAinB,sBinA,tenA,tenB);
  setContractionSlots(sAinR,sRinA,tenA,tenR);
  setContractionSlots(sBinR,sRinB,tenB,tenR);
  if(pA->name() == "T"){
    if(_tensors[0]._connect.bitmask == tenR._connect.bitmask)
      return Contraction(*pB,*pA,sBinA,sAinB,sBinR,sRinB,sAinR,sRinA,_fac);
    else
      return Contraction(*pB,*pA,sBinA,sAinB,sBinR,sRinB,sAinR,sRinA);
  }
  else{
    if(_tensors[0]._connect.bitmask == tenR._connect.bitmask)
      return Contraction(*pA,*pB,sAinB,sBinA,sAinR,sRinA,sBinR,sRinB,_fac);
    else
      return Contraction(*pA,*pB,sAinB,sBinA,sAinR,sRinA,sBinR,sRinB);
  }
}

Summation Diagram::exprSummation(const DiagramTensor& tenA, const DiagramTensor& tenR, const Tensor* pA) const
{
  Slots
    sAinR, sRinA;
  setContractionSlots(sAinR,sRinA,tenA,tenR);
  return Summation(*pA,sAinR,sRinA,_fac);
}

Cost Diagram::contractionCost( const DiagramTensor& ten1, const DiagramTensor& ten2, const DiagramTensor& res ) const
{
  Cost cost = 1;
  // find contracted slots
  std::bitset<MAXNINDICES>
    connectionmask = ten1._connect.bitmask&ten2._connect.bitmask,
    allindices = ten1._connect.bitmask|ten2._connect.bitmask;
  assert( (connectionmask&res._connect.bitmask) == 0 );
  assert( res._connect.bitmask == (ten1._connect.bitmask^ten2._connect.bitmask));
  if ( _cuts.size() ) {
    error("implement cost function with cuts","Diagram::contractionCost");
  } else {
    std::bitset<MAXNINDICES> allin(allindices);
    for ( uint ipos = 0; allin.any(); ++ipos, allin>>=1 ) {
      if ( allin[0] ) {
        cost *= _slottypes[ipos]->length();
      }
    }
  }
  return cost;
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
  bool explspin = Input::iPars["prog"]["explspin"];
  std::string tensor;
  // add default tensors
  if(explspin)
    tensor = "xtensor";
  else
    tensor = "tensor";
  const TsPar& deftensors = Input::sPars[tensor];
  for (const auto& dt: deftensors) {
    SlotTs sts;
    std::size_t
      ibra = dt.first.find_first_of('['),
      iket = dt.first.find_first_of(']');
    // generate tensor name
    std::string name = dt.first.substr(0,ibra);
    if ( ibra == std::string::npos || iket == std::string::npos ){
      if ( ibra != std::string::npos || iket != std::string::npos )
        error("Missmatch in [ ] in the default tensor definition");
    } else {
      // slottypes
      for ( std::size_t ii = ibra+1; ii < iket; ){
        SlotTypes::iterator ist = _slottypes.begin();
        ist = _slottypes.insert( ist, SlotType(ReadAndAdvance(ii,dt.first)) );
        sts.push_back(&(*ist));
      }
    }
    Tensor tens(sts,name);
    tens.CreateCutFromDesc(dt.second);
    _tensors.push_back(tens);
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
  for ( TensorsList::iterator it = _tensors.begin(); it != _tensors.end(); ++it) {
    if ( it->equal(tensor) ) {
      return &(*it);
    }
  }
  _tensors.push_back(tensor);
  if ( tensor.name() == "" ){
    const Contraction * pContr = dynamic_cast< const Contraction * >(_tensors.back().parents().back());
    if( pContr ){
      if ( _internames.count((*pContr).fingerprint(tensor)) > 0 ){
          _tensors.back()._name = _internames[(*pContr).fingerprint(tensor)];}
      else{
          _tensors.back()._name = newname(tensor.syms(),tensor.cuts());
          _internames[(*pContr).fingerprint(tensor)] = _tensors.back()._name;
      }
    }
    else
      _tensors.back()._name = newname(tensor.syms(),tensor.cuts());
  }
  return &(_tensors.back());
}

const Action* Expression::add(const Action* pAction)
{
  const Contraction * pContr = dynamic_cast< const Contraction * >( pAction );
  if ( pContr ) {
    _contractions.push_back(*pContr);
    return &(_contractions.back());
  }
  const Summation * pSum = dynamic_cast< const Summation * >( pAction );
  assert( pSum );
  _summations.push_back(*pSum);
  return &(_summations.back());
}

const Tensor* Expression::add2residual(const Tensor& res, const Action * pAct)
{
  for ( TensorsList::iterator it = _tensors.begin(); it != _tensors.end(); ++it) {
    if ( it->equal(res) ) {
      it->add(pAct);
      return &(*it);
    }
  }
  error("Something is wrong: residual not found!","Expression::add2residual");
  return 0;
}

const Tensor* Expression::find(const Tensor& tensor, bool considerprops) const
{
  for ( TensorsList::const_iterator it = _tensors.begin(); it != _tensors.end(); ++it)
    if ( it->equal(tensor,considerprops) ) return &(*it);
  return 0;
}


std::string Expression::newname(const Symmetries& syms, const Cuts& cuts)
{
  const std::string& names = Input::sPars["syntax"]["intermeds"];
  if ( _lastname.size() == 0 ) {
    _lastname = names[0];
    return _lastname;
  }
  for ( int i = _lastname.size()-1; i >= 0; --i) {
    std::size_t indx = names.find(_lastname[i]);
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
  o << d._fac;
  for ( uint i = 0; i < d._tensors.size(); ++i ){
    o << d._tensors[i].name() << "["<< d._tensors[i].slotTypeLetters(d._slottypes) << "]";
    o << "(" << d._tensors[i]._connect.bitmask << ")";
  }
  o << "}";
  return o;
}

void print_action(std::ostream& o, const Action * pAct, const Tensor& ten)
{
  assert( pAct );
  const Contraction * pContr = dynamic_cast< const Contraction * >(pAct);
  if (pContr) {
    print_code(o,*(pContr->p_A));
    print_code(o,*(pContr->p_B));
    pContr->print(o,ten);
  }
  else{
    const Summation * pSum = dynamic_cast< const Summation * >(pAct);
    assert( pSum );
    pSum->print(o,ten);
  }
}

void print_code(std::ostream& o, const Tensor& ten)
{
  if (ten._parents.size() > 0) {
    for(const Action * pAct : ten._parents)
      print_action(o,pAct,ten);
  }
}

bool Expression::elemcocompare_diags(const Diagram& diagA, const Diagram& diagB){
  std::set<std::string> sorted;
  std::vector<std::string> order = {
                                    "d_vvvv","d_VVVV","d_vVvV",
                                    "d_vvvo","d_VVVO","d_vVvO",
                                    "d_vovv","d_VOVV","d_vOvV",
                                    "d_vvov","d_VVOV","d_vVoV",
                                    "d_vvoo","d_VVOO","d_vVoO",
                                    "d_vovo","d_VOVO","d_vOvO",
                                    "d_voov","d_VOOV","d_vOoV",
                                    "d_vooo","d_VOOO","d_vOoO",
                                    "d_oooo","d_OOOO","d_oOoO",
                                    "d_ovoo","d_OVOO","d_oVoO",
                                    "d_oovo","d_OOVO","d_oOvO",
                                    "d_ooov","d_OOOV","d_oOoV",
                                    "d_ovvo","d_OVVO","d_oVvO",
                                    "d_ovov","d_OVOV","d_oVoV",
                                    "oovv","OOVV","oOvV",
                                    "d_ovvv","d_OVVV","d_oVvV"
                                    };
  std::vector<std::string>::iterator posend = order.begin();
  for( std::vector<std::string>::iterator name = order.begin(); name != order.end(); name++ ){
    if ( diagA._tensors[1].name() == *name ){
      posend += std::distance(order.begin(),name);
      for( std::vector<std::string>::iterator it = order.begin(); it != posend; ++it ){
        if ( *it == diagB._tensors[1].name() ) return false;
      }
      return true;
    }
  }
  return false;
}

void Expression::elemcosort_diags()
{
  _diagrams.sort(elemcocompare_diags);
}

void Expression::printjulia(std::ofstream& out) const
{
  std::stack<std::string> LIFO;
  for( std::list<Diagram>::const_iterator diagcit = _diagrams.begin(); diagcit != _diagrams.end(); diagcit++ ){

    std::map<SlotType::Type,std::string> slotnames;
    std::vector<Array<std::string>> slots;
    std::vector<SlotTs> slottypes;
    Diagram diag = *diagcit;

    for( auto dtit = diag._tensors.begin(); dtit != diag._tensors.end(); dtit++ ){
      Tensor ten(diag.exprTensor(*dtit));
      slottypes.push_back(ten._slots);
    }
    
    if ( slottypes.size() == 2 )//R=fac*A
    {
      Slots
        sAinR, sRinA;
    
      Array<std::string>
          resslots(slottypes[0].size()),
          aslots(slottypes[1].size());

      setContractionSlots(sAinR,sRinA,diag._tensors[1],diag._tensors[0]);

      slotNames4Refs(resslots,aslots,slotnames,sRinA,sAinR,slottypes[0],slottypes[1]);

      // print load and drop statements
      printjulia(out, diag._tensors[1].name(), LIFO);

      out << "@tensoropt ";
      out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
      if (std::abs(std::abs(diag._fac) - 1.0) > 1.e-6){
        out << sgnchar(diag._fac) << "= " << std::abs(diag._fac) << " * ";
      }
      else{
        out << sgnchar(diag._fac) << "= ";
      }
      out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "]" << std::endl;
    }
    else if ( slottypes.size() == 3 )//R=fac*A*B
    {
      Slots
        sAinB, sBinA,
        sAinR, sRinA,
        sBinR, sRinB;
    
      Array<std::string>
          resslots(slottypes[0].size()),
          aslots(slottypes[1].size()),
          bslots(slottypes[2].size());

      setContractionSlots(sAinB,sBinA,diag._tensors[1],diag._tensors[2]);
      setContractionSlots(sAinR,sRinA,diag._tensors[1],diag._tensors[0]);
      setContractionSlots(sBinR,sRinB,diag._tensors[2],diag._tensors[0]);

      slotNames4Refs(aslots,bslots,slotnames,sAinB,sBinA,slottypes[1],slottypes[2]);
      slotNames4Refs(resslots,aslots,slotnames,sRinA,sAinR,slottypes[0],slottypes[1]);
      slotNames4Refs(resslots,bslots,slotnames,sRinB,sBinR,slottypes[0],slottypes[2]);

      // print load and drop statements
      printjulia(out, diag._tensors[1].name(), LIFO);

      out << "@tensoropt ";
      out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
      if (std::abs(std::abs(diag._fac) - 1.0) > 1.e-6){
        out << sgnchar(diag._fac) << "= " << std::abs(diag._fac) << " * ";
      }
      else{
        out << sgnchar(diag._fac) << "= ";
      }
      out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
      out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "]" << std::endl;

    }
    else if ( slottypes.size() == 4 )//R=fac*A*B*C
    {
      Slots
        sAinB, sBinA,
        sAinC, sCinA,
        sAinR, sRinA,
        sBinR, sRinB,
        sBinC, sCinB,
        sCinR, sRinC;
    
      Array<std::string>
          resslots(slottypes[0].size()),
          aslots(slottypes[1].size()),
          bslots(slottypes[2].size()),
          cslots(slottypes[3].size());

      setContractionSlots(sAinB,sBinA,diag._tensors[1],diag._tensors[2]);
      setContractionSlots(sAinC,sCinA,diag._tensors[1],diag._tensors[3]);
      setContractionSlots(sAinR,sRinA,diag._tensors[1],diag._tensors[0]);
      setContractionSlots(sBinR,sRinB,diag._tensors[2],diag._tensors[0]);
      setContractionSlots(sBinC,sCinB,diag._tensors[2],diag._tensors[3]);
      setContractionSlots(sCinR,sRinC,diag._tensors[3],diag._tensors[0]);

      slotNames4Refs(resslots,aslots,slotnames,sRinA,sAinR,slottypes[0],slottypes[1]);
      slotNames4Refs(resslots,bslots,slotnames,sRinB,sBinR,slottypes[0],slottypes[2]);
      slotNames4Refs(resslots,cslots,slotnames,sRinC,sCinR,slottypes[0],slottypes[3]);
      slotNames4Refs(aslots,bslots,slotnames,sAinB,sBinA,slottypes[1],slottypes[2]);
      slotNames4Refs(aslots,cslots,slotnames,sAinC,sCinA,slottypes[1],slottypes[3]);
      slotNames4Refs(bslots,cslots,slotnames,sBinC,sCinB,slottypes[2],slottypes[3]);

      // print load and drop statements
      printjulia(out, diag._tensors[1].name(), LIFO);

      out << "@tensoropt ";
      out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
      if (std::abs(std::abs(diag._fac) - 1.0) > 1.e-6){
        out << sgnchar(diag._fac) << "= " << std::abs(diag._fac) << " * ";
      }
      else{
        out << sgnchar(diag._fac) << "= ";
      }
      out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
      out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "] * ";
      out << elemconame(diag._tensors[3].name(),slottypes[3]) << "[" << container2csstring(cslots) << "]" << std::endl;

    }
    else
      error("printjulia not implemented for more than 4 tensors in a diagram");
    if (std::next(diagcit) == _diagrams.end() && ! LIFO.empty())
      out << LIFO.top() << " = " << "nothing" << std::endl;
  }
}

void Expression::printjulia(std::ofstream& out, const std::string& tensorname, std::stack<std::string>& LIFO) const{
  std::string upper = tensorname;
  transform(tensorname.begin(), tensorname.end(), upper.begin(), toupper);
  if ( !LIFO.empty() && LIFO.top() != tensorname ) {
    out << LIFO.top() << " = " << "nothing" << std::endl;  
    LIFO.pop();
  }
  if ( tensorname != "f" ){
    if ( LIFO.empty() ) {
      if (upper == "OOVV" ) //undressed
        out << tensorname << " = " << "ints2(EC,\"" << tensorname << "\")" << std::endl;  
      else //dressed
        out << tensorname << " = " << "load(EC,\"" << tensorname << "\")" << std::endl;  
      LIFO.push(tensorname);
    }
  }
}

std::string Expression::elemconame(const std::string& name, const SlotTs& slottypes) const{
  std::string newname = name;
  bool alpha = false;
  bool beta = false;
  std::vector<std::string> T = {"T1","T2","T3"};
  std::vector<std::string> R = {"R1","R2","R3"};
  assert(slottypes.size() % 2 == 0);
  uint exclevel = (slottypes.size()/2)-1;
  for (const auto it : slottypes ){
    if ( it->type() == SlotType::Type::OccA || it->type() == SlotType::Type::VirtA)
      alpha = true;
    if ( it->type() == SlotType::Type::OccB || it->type() == SlotType::Type::VirtB)
      beta = true;
  }
  if ( alpha && beta ){
    assert(exclevel > 0);
    std::string spinstring;
    if ( exclevel == 1 )
      spinstring = "ab";
    else if (exclevel == 2){
      if (slottypes[1]->type() == SlotType::Type::VirtA)
        spinstring = "aab";
      else if (slottypes[1]->type() == SlotType::Type::VirtB)
        spinstring = "abb";
      else
        error("Unexpected orbital type. Maybe orbital order unexpected?","Expression::elemconame");
    }
    else
      error("Higher than triples not implemented","Expression::elemconame");
    if ( name == "R" ) newname = R[exclevel]+spinstring;
    else if ( name == "T" ) newname = T[exclevel]+spinstring;
  }
  else if ( alpha && !beta ){
    if ( name == "R" ) newname = R[exclevel]+"a";
    else if (name == "T" ) newname = T[exclevel]+"a";
    else if (name == "f" ){
      if ( slottypes[0]->type() == SlotType::OccA && slottypes[1]->type() == SlotType::OccA) newname = "fij";
      else if ( slottypes[0]->type() == SlotType::OccA && slottypes[1]->type() == SlotType::VirtA) newname = "fia";
      else if ( slottypes[0]->type() == SlotType::VirtA && slottypes[1]->type() == SlotType::OccA) newname = "fai";
      else newname = "fab";
    }
  }
  else if (beta && !alpha ){
    if (name == "R" ) newname = R[exclevel]+"b";
    else if (name == "T" ) newname = T[exclevel]+"b";
    else if (name == "f" ){
      if ( slottypes[0]->type() == SlotType::OccB && slottypes[1]->type() == SlotType::OccB) newname = "fIJ";
      else if ( slottypes[0]->type() == SlotType::OccB && slottypes[1]->type() == SlotType::VirtB) newname = "fIA";
      else if ( slottypes[0]->type() == SlotType::VirtB && slottypes[1]->type() == SlotType::OccB) newname = "fAI";
      else newname = "fAB";
    }
  }
  else{
    if (name == "R" ) newname = R[exclevel];
    else if (name == "T" ) newname = T[exclevel];
  }
  return newname;
}

std::ostream & operator << (std::ostream& o, const Expression& exp) {
  o << "---- decl" << std::endl;
  // index spaces
  const SlotTypes& sts = exp.slottypes();
  for (const auto& st: sts){
    o << st << std::endl;
  }
  // tensors....
  const TensorsList& ts = exp.tensors();
  std::set<Tensor> uniquetensortypes;
  for (const auto& t: ts){
    if( (uniquetensortypes.insert(t)).second)
      o << "tensor: " << t << std::endl;
  }
  uniquetensortypes.clear();

  // contractions...
  o << std::endl << "---- code (\"eval_residual\")" << std::endl;
  // print init and save statements for Koeppel's algoopt program
  for (const auto& t: ts){
    if( (t.name() != "T" && t.name()[0] != 'f' && t.type() != "I") && uniquetensortypes.insert(t).second){
      o << "init " << t.name() << "[" << t.slotTypeLetters() << "]" << std::endl;
      if( t.name() == "R" ) o << "save " << t.name() << "[" << t.slotTypeLetters() << "]" << std::endl;
    }
  }
  std::set< const Tensor *, Expression::comp_pTen > residuals = exp.residualtensors();
  if ( residuals.size() == 0 )
    o << "// No residual tensors set!" << std::endl;
  for (const auto& res: residuals) {
    xout << "// Residual: " << *res << std::endl;
    print_code(o,*res);
  }

  return o;
}
