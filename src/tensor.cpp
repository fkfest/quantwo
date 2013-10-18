#include "tensor.h"

bool SlotType::operator<(const SlotType& st) const
{
  if (_nIndices < st._nIndices) return true;
  if (st._nIndices < _nIndices) return false;
  return _type < st._type;
}

std::string SlotType::typeLetter() const
{
  TsPar& orbs = Input::sPars["syntax"];
  switch (_type) {
    case SlotType::Occ:
      return orbs["occorb"].substr(0,1);
    case SlotType::Virt:
      return orbs["virorb"].substr(0,1);
    case SlotType::Act:
      return orbs["actorb"].substr(0,1);
    case SlotType::GenT:
      return orbs["genorb"].substr(0,1);
    case SlotType::AO:
      error("AO space not implemented yet","Printing SlotType"); 
      return orbs["aoorb"].substr(0,1);
    case SlotType::DF:
      error("DF space not implemented yet","Printing SlotType"); 
      return orbs["dforb"].substr(0,1);
    case SlotType::RI:
      error("RI space not implemented yet","Printing SlotType"); 
      return orbs["riorb"].substr(0,1);
    default:
      error("Unknown space","Printing SlotType");
  }
  return "";
}


std::string Cut::StrengthLetter() const
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


std::string Tensor::SlotTypeLetters() const
{
  std::string ret;
  SlotTs::const_iterator ist;
  _foreach(ist,_slots){
    ret += (*ist)->typeLetter();
  }
  return ret;
}

bool Tensor::operator<(const Tensor& ten) const
{
  if ( _name < ten._name ) return true;
  if ( ten._name < _name ) return false;
  if ( _slots.size() < ten._slots.size() ) return true;
  if ( ten._slots.size() < _slots.size() ) return false;
  if ( _cuts.size() < ten._cuts.size() ) return true;
  if ( ten._cuts.size() < _cuts.size() ) return false;
  for ( uint i = 0; i < _slots.size(); ++i ){
    if ( _slots[i] < ten._slots[i] ) return true;
    if ( ten._slots[i] < _slots[i] ) return false;
  }
  for ( uint i = 0; i < _cuts.size(); ++i ){
    if ( _cuts[i] < ten._cuts[i] ) return true;
    if ( ten._cuts[i] < _cuts[i] ) return false;
  }
  return false;
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
  Slots::const_iterator is;
  _foreach(is,ss){
    if (*is > 9) 
      o << "{" << *is << "}";
    else
      o << *is;
  }
  return o;
}

std::ostream & operator << (std::ostream& o, const Cut& cut) {
  switch (cut.type()) {
    case Cut::Domain:
      o << cut.cutSlots() << "/" << cut.StrengthLetter();
      o << cut.defSlots();
      break;
    case Cut::List:
      o << cut.StrengthLetter() << cut.cutSlots();
      break;
    default:
      error("Cannot print this cut","Printing cuts");
  }
  for ( uint i = cut.defSlots().size(); i < cut.nSlotCutType(); ++i ) o << "*";
  return o;
}

std::ostream & operator << (std::ostream& o, const Tensor& t) {
  o << t.name() << "[" << t.SlotTypeLetters() << "]";
  o << ", !Create{";
  // type: disk, plain...
  
  // cut info
  if ( t.cuts().size() > 0 ){
    o << "cut:";
    Cuts::const_iterator ic;
    bool putcomma = false;
    _foreach(ic,t.cuts()){
      if (putcomma) o << ",";
      o << *ic;
      putcomma = true;
    }
  }
  // symmetry
  
  o << "}";
  return o;
}

