#include "tensor.h"
#include "action.h"

std::string TensorBase::type() const
{
  if (_name.size() > 4 ) return std::string("I"); //an integral
  if (_name == "T") return std::string("T"); // an amplitude
  if (_name == "R") return std::string("R"); // a residual
  return std::string("A"); //an intermediate
}

// return -1 if not able to read
static int ReadIndexAndAdvance(std::size_t& ipos, const std::string s){
    assert( ipos < s.size() );
    uint
        islot;
    std::size_t
        last = ipos+1,
        iposnew = last;
    if ( s[ipos] == '{' ){
        ++ipos;
        // closing bracket
        last = s.find_first_of( '}', ipos );
        assert( last != std::string::npos );
        iposnew = last+1;
    }
    if ( !str2num<uint>(islot,s.substr(ipos,last-ipos),std::dec) )
        return -1;
    ipos = iposnew;
    return islot;

}

// reads character or a set of characters in curly brackets, i.e. "a" or "{a1}" etc
static std::string ReadNameAndAdvance(std::size_t& ipos, const std::string& s){
    assert( ipos < s.size() );
    std::size_t iposold = ipos;
    if ( s[iposold] == '{' ){
        // closing bracket
        ipos = s.find_first_of( '}', iposold );
        assert( ipos != std::string::npos );
    }
    ++ipos;
    return s.substr(iposold,ipos-iposold);
}

static std::string nextName(const std::string& s, const std::string& oldname){
  std::string name = "";
  std::size_t ipos = 0;
  while( name != oldname && ipos < s.size() ) {
    name = ReadNameAndAdvance(ipos,s);
  }
  if ( ipos >= s.size() )
    error( "add more slot names: "+s,"nextName");
  return ReadNameAndAdvance(ipos,s);
}

SlotType::SlotType(const std::string& lettertype)
{
  TsPar& orbs = Input::sPars["syntax"];
  if ( curlyfind(orbs["occorb"],lettertype) != std::string::npos ) {
    _type = SlotType::Occ;
    _nIndices = Input::iPars["fact"]["nocc"];
    _internalName = "occorb";
  } else if ( curlyfind(orbs["virorb"],lettertype) != std::string::npos ) {
    _type = SlotType::Virt;
    _nIndices = Input::iPars["fact"]["nvir"];
    _internalName = "virorb";
  } else if ( curlyfind(orbs["actorb"],lettertype) != std::string::npos ) {
    _type = SlotType::Act;
    _nIndices = Input::iPars["fact"]["nact"];
    _internalName = "actorb";
  } else if ( curlyfind(orbs["genorb"],lettertype) != std::string::npos ) {
    _type = SlotType::GenT;
    _nIndices = Input::iPars["fact"]["nocc"]+Input::iPars["fact"]["nvir"];
    _internalName = "genorb";
    if ( Input::iPars["prog"]["multiref"] > 0 ) _nIndices += Input::iPars["fact"]["nact"];
  } else {
    error("Unknown letter-type space!","SlotType constructor");
  }
}


bool SlotType::operator<(const SlotType& st) const
{
  if (_nIndices < st._nIndices) return true;
  if (st._nIndices < _nIndices) return false;
  return _type < st._type;
}

std::string SlotType::name(const std::string& oldname) const
{
  TsPar& orbs = Input::sPars["syntax"];
  return nextName(orbs[_internalName],oldname);
}

void Symmetry::canonicalize()
{
  assert( _simSlots.size() == 0 || _simSlots.size() == _symSlots.size() );
  assert( _symSlots.size() > 0 );
  Slots ref;
  ref.identity(_symSlots.size());
//  lui nswaps =
  InsertionSort( &_symSlots[0], &ref[0], ref.size() );
  _symSlots = _symSlots.refarr(ref);
  if ( _simSlots.size() > 0 ) _simSlots = _simSlots.refarr(ref);
//  int ret = 1;
//  if ( nswaps % 2 == 1 ) ret = _sign;
//  return ret;
}

bool Symmetry::operator<(const Symmetry& sym) const
{
  if ( _symSlots.size() < sym._symSlots.size() ) return true;
  if ( sym._symSlots.size() < _symSlots.size() ) return false;
  if ( _simSlots.size() < sym._simSlots.size() ) return true;
  if ( sym._simSlots.size() < _simSlots.size() ) return false;
  for ( uint i = 0; i < _symSlots.size(); ++i ){
    if ( _symSlots[i] < sym._symSlots[i] ) return true;
    if ( sym._symSlots[i] < _symSlots[i] ) return false;
  }
  for ( uint i = 0; i < _simSlots.size(); ++i ){
    if ( _simSlots[i] < sym._simSlots[i] ) return true;
    if ( sym._simSlots[i] < _simSlots[i] ) return false;
  }
  return ( _sign < sym._sign );
}

bool Symmetry::operator==(const Symmetry& sym) const
{
  if ( _symSlots.size() == sym._symSlots.size() &&
       _simSlots.size() == sym._simSlots.size() &&
       _sign == sym._sign ) {
    for ( uint i = 0; i < _symSlots.size(); ++i )
      if ( _symSlots[i] != sym._symSlots[i] ) return false;
    for ( uint i = 0; i < _simSlots.size(); ++i )
      if ( _simSlots[i] != sym._simSlots[i] ) return false;
    return true;
  }
  return false;
}


std::string Cut::strengthLetter() const
{
  switch (_cutStrength){
    case Cut::DefSth:
      return "";
    case Cut::Strong:
      return "s";
    case Cut::Close:
      return "c";
    case Cut::Weak:
      return "w";
    case Cut::Dist:
      return "d";
    case Cut::NoSth:
      return "n";
  }
  return "";
}

bool Cut::operator<(const Cut& cut) const
{
  if ( _cutType < cut._cutType ) return true;
  if ( cut._cutType < _cutType ) return false;
  if ( _nSlotCutType < cut._nSlotCutType ) return true;
  if ( cut._nSlotCutType < _nSlotCutType ) return false;
  if ( _cutSlots.size() < cut._cutSlots.size() ) return true;
  if ( cut._cutSlots.size() < _cutSlots.size() ) return false;
  if ( _defSlots.size() < cut._defSlots.size() ) return true;
  if ( cut._defSlots.size() < _defSlots.size() ) return false;
  for ( uint i = 0; i < _cutSlots.size(); ++i ){
    if ( _cutSlots[i] < cut._cutSlots[i] ) return true;
    if ( cut._cutSlots[i] < _cutSlots[i] ) return false;
  }
  for ( uint i = 0; i < _defSlots.size(); ++i ){
    if ( _defSlots[i] < cut._defSlots[i] ) return true;
    if ( cut._defSlots[i] < _defSlots[i] ) return false;
  }
  return ( _cutStrength < cut._cutStrength );
}

bool Cut::operator==(const Cut& cut) const
{
  if ( _cutType == cut._cutType &&
       _nSlotCutType == cut._nSlotCutType &&
       _cutSlots.size() == cut._cutSlots.size() &&
       _defSlots.size() == cut._defSlots.size() &&
       _cutStrength == cut._cutStrength ) {
    for ( uint i = 0; i < _cutSlots.size(); ++i )
      if ( _cutSlots[i] != cut._cutSlots[i] ) return false;
    for ( uint i = 0; i < _defSlots.size(); ++i )
      if ( _defSlots[i] != cut._defSlots[i] ) return false;
    return true;
  }
  return false;
}

void Cut::canonicalize()
{
  _cutSlots.resort();
  _defSlots.resort();
}


void Canonicalize(SlotTs& sts, Slots& ref)
{
  ref.identity(sts.size());
  InsertionPSortD(&sts[0],&ref[0],sts.size());
  sts = sts.refarr(ref);
}

void Tensor::CreateCutFromDesc(const std::string& desc)
{

  struct loc{
    // return strength or other special modes (as uint(FTensorCut::CutNoSth)+type) (e.g. forced triangular)
    // return nSlots as a count of all slots (including stars)
    static uint SetSlots( Slots& SlotSet, uint& nSlots,
              uint iFirstChar, uint iLastChar, std::string const &s ){
      assert( iLastChar <= s.size() );
      assert( s.size() > iFirstChar );
      uint iret = 0;
      // first char can represent type
      if ( s[iFirstChar] == 'f' ) {
        // list of phantom slots
        iret = uint(Cut::NoSth)+uint(Cut::NoType);
        ++iFirstChar;
      } else if ( s[iFirstChar] == 't' ) {
        // forced triangular mode
        iret = uint(Cut::NoSth)+uint(Cut::Triang);
        ++iFirstChar;
      } else if ( s[iFirstChar] == 's' ) {
        // forced strong list
        iret = Cut::Strong;
        ++iFirstChar;
      } else if ( s[iFirstChar] == 'c' ) {
        // forced strong+close list
        iret = Cut::Close;
        ++iFirstChar;
      } else if ( s[iFirstChar] == 'w' ) {
        // forced strong+close+weak list
        iret = Cut::Weak;
        ++iFirstChar;
      } else if ( s[iFirstChar] == 'd' ) {
        // forced strong+close+weak+distant list
        iret = Cut::Dist;
        ++iFirstChar;
      }
      bool goodstr = true;
      int islot;
      nSlots = 0;
      for ( std::size_t i = iFirstChar; goodstr && i < iLastChar; ++i ){
        if ( s[i] == '*' ) {
          ++nSlots;
          continue;
        }
        islot = ReadIndexAndAdvance(i,s);
        goodstr = ( islot >= 0 );
        if ( goodstr ) {
          ++nSlots;
          SlotSet.push_back(uint(islot));
          --i;
        }
      }
      assert ( goodstr );
      return iret;
    };
  };

  std::size_t
      // find the cut description
      iFirst = desc.find("cut:"),
      iLast;
  if ( iFirst != std::string::npos ) iFirst += 4;
  bool lastcut = false;
  for ( ; iFirst != std::string::npos && !lastcut; ) {
    // find next comma
    iLast = desc.find_first_of( ',', iFirst );
    std::size_t illast = desc.find_first_of( ';', iFirst );
    if ( illast < iLast ) {
      iLast = illast;
      lastcut = true;
    }

    std::string
        SubDesc = desc.substr(iFirst, iLast-iFirst);

    iFirst = (iLast == std::string::npos) ? iLast : iLast + 1;

    if ( SubDesc.empty() ) continue;

    Cut cut;
    cut._cutStrength = Cut::DefSth;
    cut._PDCStrength = Cut::NoSth;
    uint ltype;
    // Find separator of slot sets.
    std::size_t iSep = SubDesc.find_first_of("/");
    if (iSep == std::string::npos){
      // Cut independent of other orbitals ("List"-cut)
      cut._cutType = Cut::List;
      ltype = loc::SetSlots(cut._cutSlots,cut._nSlotCutType,0,SubDesc.size(),SubDesc);
      if ( ltype >= uint(Cut::NoSth)){
        ltype -= uint(Cut::NoSth);
        cut._cutType = (Cut::Type)ltype;
        ltype = 0;
        if ( cut._cutType == Cut::Triang ) {
          // forced triangular
          // TODO
          error("forced triangular is not implemented yet");
          //_cuts.TriOps.push_back(Out.CutOps.size());
        } else {
          assert( cut._cutType == Cut::NoType );
          // add to set of phantom slots
          _phantomSlots.insert(cut._cutSlots.begin(),cut._cutSlots.end());
        }
      } else {
//        Out.ListOps.push_back(Out.CutOps.size());
        cut._cutStrength = (Cut::Strength)ltype;
      }
    } else {
      // Cut dependent on other orbitals ("Domain"-cut)
      cut._cutType = Cut::Domain;
      uint dummy;
      ltype = loc::SetSlots(cut._cutSlots,dummy,0,iSep,SubDesc);
      // no triangular cut or phantom slots here
      assert( ltype < uint(Cut::NoSth) );
      cut._cutStrength = (Cut::Strength)ltype;
      ltype = loc::SetSlots(cut._defSlots,cut._nSlotCutType,iSep+1,SubDesc.size(),SubDesc);
      // no triangular cut or phantom slots here
      assert( ltype < uint(Cut::NoSth) );
      // not explicitly given strength means no strength here
      if ( ltype == 0 ) ltype = uint(Cut::NoSth);
      cut._PDCStrength = (Cut::Strength)ltype;
//      _cuts.DomOps.push_back(Out.CutOps.size());
    }
    if ( cut._cutType != Cut::NoType )
      // all but phantom-slots cuts
      _cuts.push_back( cut );
  }

}

std::string DiagramTensor::slotTypeLetters( const SlotTs& slottypes ) const
{
  assert( _connect.bitmask.size() >= slottypes.size() );
  std::string ret;
  std::vector<std::string> rr;
  rr.resize(_connect.slotref.size());
  for ( uint i = 0; i < slottypes.size(); ++i ){
    if ( _connect.bitmask[i] ) {
      uint ist = _connect.slotref[(_connect.bitmask>>(i+1)).count()];
      assert( ist < rr.size() );
      rr[ist] += slottypes[i]->name();
    }
  }
  for ( uint i = 0; i < rr.size(); ++i ){
    ret += rr[i];
  }
  return ret;
}

void Tensor::add(const Action* pAct)
{
  if (pAct) {
    for (const auto& act: _parents){
      if ( pAct == act ) return;
    }
    const Contraction * pContr = dynamic_cast< const Contraction * >( pAct );
    if(pContr){
      if(find(pContr->p_A))
        insert_action(pAct, pContr->p_A);
      else
        _parents.push_back(pAct);
    }
    else{
      _parents.push_back(pAct);
    }
  }
}

std::string Tensor::slotTypeLetters() const
{
  std::string ret;
  for (const auto& st: _slots){
    ret += st->name();
  }
  return ret;
}

bool Tensor::operator<(const Tensor& ten) const
{
  if ( _name < ten._name ) return true;
  if ( ten._name < _name ) return false;
  if ( _slots.size() < ten._slots.size() ) return true;
  if ( ten._slots.size() < _slots.size() ) return false;
  if ( _syms.size() < ten._syms.size() ) return true;
  if ( ten._syms.size() < _syms.size() ) return false;
  if ( _cuts.size() < ten._cuts.size() ) return true;
  if ( ten._cuts.size() < _cuts.size() ) return false;
  if ( _dummy < ten._dummy ) return true;
  if ( ten._dummy < _dummy ) return false;
  for ( uint i = 0; i < _slots.size(); ++i ){
    if ( _slots[i] < ten._slots[i] ) return true;
    if ( ten._slots[i] < _slots[i] ) return false;
  }
  for ( uint i = 0; i < _syms.size(); ++i ){
    if ( _syms[i] < ten._syms[i] ) return true;
    if ( ten._syms[i] < _syms[i] ) return false;
  }
  for ( uint i = 0; i < _cuts.size(); ++i ){
    if ( _cuts[i] < ten._cuts[i] ) return true;
    if ( ten._cuts[i] < _cuts[i] ) return false;
  }
  return false;
}

bool Tensor::find( const Tensor* p_A){
  for(auto it = _parents.rbegin(); it != _parents.rend(); ++it ){
    const Contraction * pContr = dynamic_cast< const Contraction * >( *(it) );
    if (pContr){
      if(*(pContr->p_A) == *p_A ){
        return true;
      }
    }
  }
  return false;
}

void Tensor::insert_action(const Action* pAct, const Tensor* p_A){
  for(auto it = _parents.rbegin(); it != _parents.rend(); ++it ){
  const Contraction * pContr = dynamic_cast< const Contraction * >( *(it) );
    if (pContr){
      if (*(pContr->p_A) == *p_A ){
        _parents.insert(it.base(),pAct);
        return;
      }
    }
  }
}

bool Tensor::equal(const Tensor& ten, bool considerprops) const
{
  if ( _name == ten._name &&
      _slots.size() == ten._slots.size() &&
      _syms.size() == ten._syms.size()
     ) {
    for ( uint i = 0; i < _slots.size(); ++i )
      if ( _slots[i] != ten._slots[i] ) return false;
    if ( considerprops ) {
      if (_syms.size() != ten._syms.size() ||
          _cuts.size() != ten._cuts.size()) return false;
      for ( uint i = 0; i < _syms.size(); ++i )
        if ( _syms[i] != ten._syms[i] ) return false;
      for ( uint i = 0; i < _cuts.size(); ++i )
        if ( _cuts[i] != ten._cuts[i] ) return false;
    }
    return true;
  }
  return false;
}

std::ostream & operator << (std::ostream& o, const DiagramTensor& t) {
  o << t.name() << "{phantom:" << t.phantom() << " sym:" << t.syms();
  o << "(" << t.connects().bitmask << ")}";
  return o;
}

std::ostream & operator << (std::ostream& o, const SlotType& st)
{
  // introduce a different parameter for the orb-names in ITF?
  TsPar& orbs = Input::sPars["syntax"];
  o << "index-space: ";
  switch (st.type()) {
    case SlotType::Occ:
      o << orbs["occorb"] << ", Closed, c";
      break;
    case SlotType::Virt:
      o << orbs["virorb"] << ", External, e";
      break;
    case SlotType::Act:
      o << orbs["actorb"] << ", Active, a";
      break;
    case SlotType::GenT:
      o << orbs["genorb"] << ", MolOrb, r";
      break;
    case SlotType::AO:
      error("AO space not implemented yet","Printing SlotType");
      o << orbs["aoorb"] << ", BasisAo, A";
      break;
    case SlotType::DF:
      error("DF space not implemented yet","Printing SlotType");
      o << orbs["dforb"] << ", BasisMp2Fit, F";
      break;
    case SlotType::RI:
      error("RI space not implemented yet","Printing SlotType");
      o << orbs["riorb"] << ", BasisRi, R";
      break;
    default:
      error("Unknown space","Printing SlotType");
  }
  o << " // optimization length: " << st.length();
  return o;
}

std::ostream & operator << (std::ostream& o, const Slots& ss) {
  for (const auto& s: ss){
    if (s > 9)
      o << "{" << s << "}";
    else
      o << s;
  }
  return o;
}
std::ostream & operator << (std::ostream& o, const SlotUniqueSet& ss) {
  for (const auto& s: ss){
    if (s > 9)
      o << "{" << s << "}";
    else
      o << s;
  }
  return o;
}
std::ostream & operator << (std::ostream& o, const Symmetry& sym) {
  if ( sym._sign < 0 ) {
    o << "-";
  } else if ( sym._simSlots.empty() ) {
    o << "+";
  }
  o << sym._symSlots;
  if ( !sym._simSlots.empty() ) {
    o << "/" << sym._simSlots;
  }
  return o;
}

std::ostream & operator << (std::ostream& o, const Cut& cut) {
  uint nst = 0;
  switch (cut.type()) {
    case Cut::Domain:
      o << cut.cutSlots() << "/" << cut.strengthLetter();
      o << cut.defSlots();
      nst = cut.defSlots().size();
      break;
    case Cut::List:
      o << cut.strengthLetter() << cut.cutSlots();
      nst = cut.cutSlots().size();
      break;
    default:
      error("Cannot print this cut","Printing cuts");
  }
  for ( uint i = nst; i < cut.nSlotCutType(); ++i ) o << "*";
  return o;
}

std::ostream & operator << (std::ostream& o, const Tensor& t) {
  o << t.name() << "[" << t.slotTypeLetters() << "]";
  o << ", !Create{";
  // type: disk, plain...

  // cut info
  if ( t.cuts().size() > 0 ){
    o << "cut:";
    bool putcomma = false;
    for (const auto& c: t.cuts()){
      if (putcomma) o << ",";
      o << c;
      putcomma = true;
    }
  }
  // symmetry

  o << "}";
  return o;
}

