#include "expression.h"
#include "sum.h"

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

std::string Expression::juliacost(const std::vector<SlotTs>& slottypes, const Array<std::string>& resslots, const Array<std::string>& aslots, const Array<std::string>& bslots, const Array<std::string>& cslots) const
{
  std::set<std::string> uniqueindices;
  std::string coststring = "(";
  for(uint i = 0; i < resslots.size(); i++){
    if ( uniqueindices.insert(resslots[i]).second){
      coststring += resslots[i]+"=>";
      if ( slottypes[0][i]->type() == SlotType::Type::Virt  || 
          slottypes[0][i]->type() == SlotType::Type::VirtA || 
          slottypes[0][i]->type() == SlotType::Type::VirtB)
        coststring += "10*x";
      else
        coststring += "x";
      coststring += ",";
    }
  }
  for(uint i = 0; i < aslots.size(); i++){
    if ( uniqueindices.insert(aslots[i]).second){
      coststring += aslots[i]+"=>";
      if ( slottypes[1][i]->type() == SlotType::Type::Virt  || 
          slottypes[1][i]->type() == SlotType::Type::VirtA || 
          slottypes[1][i]->type() == SlotType::Type::VirtB)
        coststring += "10*x";
      else
        coststring += "x";
      coststring += ",";
    }
  }
  for(uint i = 0; i < bslots.size(); i++){
    if ( uniqueindices.insert(bslots[i]).second){
      coststring += bslots[i]+"=>";
      if ( slottypes[2][i]->type() == SlotType::Type::Virt  || 
          slottypes[2][i]->type() == SlotType::Type::VirtA || 
          slottypes[2][i]->type() == SlotType::Type::VirtB)
        coststring += "10*x";
      else
        coststring += "x";
      coststring += ",";
    }
  }
  for(uint i = 0; i < cslots.size(); i++){
    if ( uniqueindices.insert(cslots[i]).second){
      coststring += cslots[i]+"=>";
      if ( slottypes[3][i]->type() == SlotType::Type::Virt  || 
          slottypes[3][i]->type() == SlotType::Type::VirtA || 
          slottypes[3][i]->type() == SlotType::Type::VirtB)
        coststring += "10*x";
      else
        coststring += "x";
      coststring += ",";
    }
  }
  if ( coststring.back() == ',' ) coststring.pop_back();
  coststring += ") ";
  return coststring;
}

uint Expression::extorb(const Array<std::string>& resslots, const Array<std::string>& aslots) const{
  uint extorb = 0;
  Array<uint> locextorb;
  for( Array<std::string>::const_iterator it = resslots.begin(); it != resslots.end(); it++ ){
    Array<std::string>::const_iterator jt = std::find(aslots.begin(), aslots.end(), *it);
    if( jt != aslots.end() ){
      locextorb.push_back(std::distance(aslots.begin(),jt));
      extorb++;
    }
  }
  for(auto&n : locextorb)
    extorb += n;
  return extorb;
}

void Expression::printfac(std::ofstream& out, const Factor& fac) const{
  if (std::abs(std::abs(fac) - 1.0) > 1.e-6)
    out << sgnchar(fac) << "= " << std::abs(fac) << " * ";
  else
    out << sgnchar(fac) << "= ";
}

void Expression::printjulia(std::ofstream& out) const
{
  std::stack<std::string> LIFO;
  bool startedblock = false;
  Array<std::string> primR, primA, primB, primC;
  for( std::list<Diagram>::const_iterator diagcit = _diagrams.begin(); diagcit != _diagrams.end(); diagcit++ ){

    std::map<SlotType::Type,std::string> slotnames;
    std::vector<Array<std::string>> slots;
    std::vector<SlotTs> slottypes,slottypes2;
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
      printfac(out,diag._fac);
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

      if( std::next(diagcit,1) != _diagrams.end()){
        slotnames.clear();
        Diagram nextdiag = *std::next(diagcit,1);
        Slots
          sAinB2, sBinA2,
          sAinC2, sCinA2,
          sAinR2, sRinA2,
          sBinR2, sRinB2;

        for( auto dtit = nextdiag._tensors.begin(); dtit != nextdiag._tensors.end(); dtit++ ){
          Tensor ten(nextdiag.exprTensor(*dtit));
          slottypes2.push_back(ten._slots);
        }
      
        Array<std::string>
            resslots2(slottypes2[0].size()),
            aslots2(slottypes2[1].size()),
            bslots2(slottypes2[2].size());

        setContractionSlots(sAinB2,sBinA2,nextdiag._tensors[1],nextdiag._tensors[2]);
        setContractionSlots(sAinR2,sRinA2,nextdiag._tensors[1],nextdiag._tensors[0]);
        setContractionSlots(sBinR2,sRinB2,nextdiag._tensors[2],nextdiag._tensors[0]);

        slotNames4Refs(resslots2,aslots2,slotnames,sRinA2,sAinR2,slottypes2[0],slottypes2[1]);
        slotNames4Refs(resslots2,bslots2,slotnames,sRinB2,sBinR2,slottypes2[0],slottypes2[2]);
        slotNames4Refs(aslots2,bslots2,slotnames,sAinB2,sBinA2,slottypes2[1],slottypes2[2]);

        // print load and drop statements
        printjulia(out, diag._tensors[1].name(), LIFO);
        if(diag.equalestimate(*(std::next(diagcit,1))) && !startedblock
            && (extorb(resslots,aslots) == extorb(resslots2,aslots2))
            && (extorb(resslots,bslots) == extorb(resslots2,bslots2)))
        {
          startedblock = true;
          primR = resslots;
          primA = aslots;
          primB = bslots;

          out << "@tensoropt ";
          out << "begin" << std::endl;
          out << "X" << "[" << container2csstring(primR) << "] ";
          out << ":= ";
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "]" << std::endl;

          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(primR) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;;
        }
        else if (diag.equalestimate(*(std::next(diagcit,1))) && startedblock
              && (extorb(resslots,aslots) == extorb(resslots2,aslots2))
              && (extorb(resslots,bslots) == extorb(resslots2,bslots2)))
        {
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.permute(resslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;
        }
        else if (!diag.equalestimate(*(std::next(diagcit,1))) && startedblock
            && (extorb(resslots,aslots) == extorb(resslots2,aslots2))
            && (extorb(resslots,bslots) == extorb(resslots2,bslots2)))
        {
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.permute(resslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;
          out << "end" << std::endl;
          startedblock = false;
        }
        else if (startedblock && ( (extorb(resslots,aslots) != extorb(resslots2,aslots2))
              || (extorb(resslots,bslots) != extorb(resslots2,bslots2))
              || !diag.equalestimate(*(std::next(diagcit,1)))))
        {
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.permute(resslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;;
          out << "end" << std::endl;
          startedblock = false;
        }
        else{
          out << "@tensoropt ";
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "]" << std::endl;
        }
      }
      else{
      //print last diagram
        if( startedblock ){
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.permute(resslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;
          out << "end" << std::endl;
          startedblock = false;
        }
        else{
          out << "@tensoropt ";
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "]" << std::endl;
        }
      }
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
      if( std::next(diagcit,1) != _diagrams.end()){
        slotnames.clear();
        Diagram nextdiag = *std::next(diagcit,1);
        Slots
          sAinB2, sBinA2,
          sAinC2, sCinA2,
          sAinR2, sRinA2,
          sBinR2, sRinB2,
          sBinC2, sCinB2,
          sCinR2, sRinC2;
        for( auto dtit = nextdiag._tensors.begin(); dtit != nextdiag._tensors.end(); dtit++ ){
          Tensor ten(nextdiag.exprTensor(*dtit));
          slottypes2.push_back(ten._slots);
        }
      
        Array<std::string>
            resslots2(slottypes2[0].size()),
            aslots2(slottypes2[1].size()),
            bslots2(slottypes2[2].size()),
            cslots2(slottypes2[3].size());

        setContractionSlots(sAinB2,sBinA2,nextdiag._tensors[1],nextdiag._tensors[2]);
        setContractionSlots(sAinC2,sCinA2,nextdiag._tensors[1],nextdiag._tensors[3]);
        setContractionSlots(sAinR2,sRinA2,nextdiag._tensors[1],nextdiag._tensors[0]);
        setContractionSlots(sBinR2,sRinB2,nextdiag._tensors[2],nextdiag._tensors[0]);
        setContractionSlots(sBinC2,sCinB2,nextdiag._tensors[2],nextdiag._tensors[3]);
        setContractionSlots(sCinR2,sRinC2,nextdiag._tensors[3],nextdiag._tensors[0]);

        slotNames4Refs(resslots2,aslots2,slotnames,sRinA2,sAinR2,slottypes2[0],slottypes2[1]);
        slotNames4Refs(resslots2,bslots2,slotnames,sRinB2,sBinR2,slottypes2[0],slottypes2[2]);
        slotNames4Refs(resslots2,cslots2,slotnames,sRinC2,sCinR2,slottypes2[0],slottypes2[3]);
        slotNames4Refs(aslots2,bslots2,slotnames,sAinB2,sBinA2,slottypes2[1],slottypes2[2]);
        slotNames4Refs(aslots2,cslots2,slotnames,sAinC2,sCinA2,slottypes2[1],slottypes2[3]);
        slotNames4Refs(bslots2,cslots2,slotnames,sBinC2,sCinB2,slottypes2[2],slottypes2[3]);

        // print load and drop statements
        printjulia(out, diag._tensors[1].name(), LIFO);
        if(diag.equalestimate(*(std::next(diagcit,1))) && !startedblock 
            && (extorb(resslots,aslots) == extorb(resslots2,aslots2))
            && (extorb(resslots,bslots) == extorb(resslots2,bslots2))
            && (extorb(resslots,cslots) == extorb(resslots2,cslots2)))
        {
          startedblock = true;
          primR = resslots;
          primA = aslots;
          primB = bslots;
          primC = cslots;

          out << "@tensoropt ";
          out << juliacost(slottypes,resslots,aslots,bslots,cslots);
          out << "begin" << std::endl;
          out << "X" << "[" << container2csstring(primR) << "] ";
          out << ":= ";
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "] * ";
          out << elemconame(diag._tensors[3].name(),slottypes[3]) << "[" << container2csstring(cslots) << "]" << std::endl;

          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(primR) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;;
        }
        else if (diag.equalestimate(*(std::next(diagcit,1))) && startedblock
            && (extorb(resslots,aslots) == extorb(resslots2,aslots2))
            && (extorb(resslots,bslots) == extorb(resslots2,bslots2))
            && (extorb(resslots,cslots) == extorb(resslots2,cslots2)))
        {
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.createPermMap(primC,cslots);
          diag.permute(resslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;
        }
        else if (startedblock && ( (extorb(resslots,aslots) != extorb(resslots2,aslots2))
            || (extorb(resslots,bslots) != extorb(resslots2,bslots2))
            || (extorb(resslots,cslots) != extorb(resslots2,cslots2))
            || !diag.equalestimate(*(std::next(diagcit,1)))))
        {
          diag.createPermMap(primA,aslots);
          diag.createPermMap(primB,bslots);
          diag.createPermMap(primC,cslots);
          diag.permute(primR);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;;
          out << "end" << std::endl;
          startedblock = false;
        }
        else{
          out << "@tensoropt ";
          out << juliacost(slottypes,resslots,aslots,bslots,cslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "] * ";
          out << elemconame(diag._tensors[3].name(),slottypes[3]) << "[" << container2csstring(cslots) << "]" << std::endl;
        }
      }
      else{
      //print last diagram
        if( startedblock ){
            diag.createPermMap(primA,aslots);
            diag.createPermMap(primB,bslots);
            diag.createPermMap(primC,cslots);
            diag.permute(primR);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << "X" << "[" << container2csstring(primR) << "]" << std::endl;
          out << "end" << std::endl;
          startedblock = false;
        }
        else{
          out << "@tensoropt ";
          out << juliacost(slottypes,resslots,aslots,bslots,cslots);
          out << elemconame(diag._tensors[0].name(),slottypes[0]) << "[" << container2csstring(resslots) << "] ";
          printfac(out,diag._fac);
          out << elemconame(diag._tensors[1].name(),slottypes[1]) << "[" << container2csstring(aslots) << "] * ";
          out << elemconame(diag._tensors[2].name(),slottypes[2]) << "[" << container2csstring(bslots) << "] * ";
          out << elemconame(diag._tensors[3].name(),slottypes[3]) << "[" << container2csstring(cslots) << "]" << std::endl;
        }
      }
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
        out << tensorname << " = " << "load4idx(EC,\"" << tensorname << "\")" << std::endl;  
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
    else if (name == "f" ){
      if ( slottypes[0]->type() == SlotType::Occ && slottypes[1]->type() == SlotType::Occ) newname = "fij";
      else if ( slottypes[0]->type() == SlotType::Occ && slottypes[1]->type() == SlotType::Virt) newname = "fia";
      else if ( slottypes[0]->type() == SlotType::Virt && slottypes[1]->type() == SlotType::Occ) newname = "fai";
      else newname = "fab";
    }
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
