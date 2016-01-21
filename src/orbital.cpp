#include "orbital.h"

Orbital::Orbital() 
{ 
  _type = Orbital::GenT;
  _spin = Spin(Spin::No);
}

Orbital::Orbital(const std::string& name, Electron el)
{
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  if ( isupper((char)name[0]) )
    _spin = Spin(el, spintype);
  else
    _spin = Spin(Spin::No);
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
  
}

Orbital::Orbital(const std::string& name, Orbital::Type type, Electron el)
{
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  if ( isupper((char)name[0]) )
    _spin = Spin(el, spintype);
  else
    _spin = Spin(Spin::No);
  _name=name;
  _name[0] = std::tolower(name[0]);
  _type=type;
}

Orbital::Orbital(const std::string& name, Spin spin)
{
  _spin=spin;
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
}
Orbital::Orbital(const std::string& name, Spin::Type spint, Electron el)
{
  _spin = Spin(el,spint);
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
}
Orbital::Orbital(const std::string& name, Orbital::Type type, Spin spin)
{
  _spin=spin;
  _type=type;
  _name=name;
  _name[0] = std::tolower(name[0]);
}
Orbital::Orbital(const std::string& name, Orbital::Type type, Spin::Type spint, Electron el)
{
  _spin = Spin(el,spint);
  _type=type;
  _name=name;
  _name[0] = std::tolower(name[0]);
}


void Orbital::gentype()
{
  TsPar& orbs = Input::sPars["syntax"];
  if (orbs["virorb"].find(_name[0])!=std::string::npos) {_type=Virt;}
  else if (orbs["occorb"].find(_name[0])!=std::string::npos) {_type=Occ;}
  else if (orbs["genorb"].find(_name[0])!=std::string::npos) {_type=GenT;}
  else if (orbs["actorb"].find(_name[0])!=std::string::npos) {_type=Act;}
  else 
    error("Unknown type of orbital! "+_name,"Orbital::gentype"); 
}

std::string Orbital::name() const
{ return _name;}

Orbital::Type Orbital::type() const
{ return _type;}

Spin Orbital::spin() const
{ return _spin;}
void Orbital::setspin(Spin spin)
{ _spin=spin; }

bool Orbital::operator==(const Orbital& orb) const
{ return _spin==orb._spin && _type==orb._type && _name==orb._name;}
bool Orbital::operator!=(const Orbital& orb) const
{ return !(*this==orb); }

bool Orbital::operator<(const Orbital& orb) const
{
  if (_type<orb._type) return true;
  if (orb._type<_type) return false;
  if (_name<orb._name) return true;
  if (_name>orb._name) return false;
//  assert( _type == orb._type );
  return _spin<orb._spin;
}
std::string Orbital::letname() const
{
  long int iend;
  for ( iend = this->_name.size(); iend > 0 && isdigit(this->_name[iend-1]); --iend ){}
  return _name.substr(0,iend);
}

int Orbital::comp_letname(const Orbital& orb) const
{
  //remove numbers from end
  std::string 
    lname = this->letname(), 
    lnameo = orb.letname();
  if ( lname.size() < lnameo.size() ) return -1;
  else if ( lname.size() > lnameo.size() ) return 1;
  else if ( lname < lnameo ) return -1;
  else if ( lname > lnameo ) return 1;
  return 0;
}
bool Orbital::is_in_set(const TOrbSet& orbset) const
{
  _foreach_cauto(TOrbSet,its,orbset){
    if (this->comp_letname(*its) == 0)
      return true;
  }
  return false;
}
Return Orbital::replace(const Orbital& orb1, const Orbital& orb2, bool smart)
{
  if (*this == orb1) { 
    *this = orb2; 
    if (smart) {
      // restore spin -- will be replaced explicitely 
      this->_spin = orb1._spin;
    }
  }
  if (smart)
    return this->_spin.replace(orb1._spin,orb2._spin);
  else
    return Return::Done;
}
Return Orbital::replace(const Spin& spin1, const Spin& spin2, bool smart)
{
  return this->_spin.replace(spin1,spin2);
}

std::ostream & operator << (std::ostream & o, Orbital const & orb)
{
  if ( orb.spin().type() == Spin::Gen )
    o << char(std::toupper(orb.name().at(0)));
  else
    o << orb.name().at(0);
  if (orb.name().size() > 1)
    o <<"_{" << orb.name().substr(1) << "}";
  if ( orb.spin().type() != Spin::Gen )
    o << orb.spin();
/*  if ( orb.type()==Orbital::Occ ) 
    o << "^occ";
  else if ( orb.type()==Orbital::Virt )
    o << "^virt";
  else
    o << "^gen";
*/
  return o;
}
bool Spin::operator==(const Spin& spin) const
{
  return (_type == spin._type && _el == spin._el);
}

bool Spin::operator<(const Spin& spin) const
{
  if (_el < spin._el) return true;
  if (spin._el < _el) return false;
  return  (_type < spin._type);
}
Return Spin::replace(const Spin& s1, const Spin& s2)
{
  assert( _type != No );
  assert( s2._type != No );
  if ( this->_el == s1._el ) {
    assert( _type == s1._type ); // if not eq --> something wrong with the kroneckers!
    // replace
    if ( _type == s2._type || _type == Gen ){
      *this = s2;
      return Return::Done;
    }
    // delete the term
    if ( (_type == Up && s2._type == Down) ||
         (_type == Down && s2._type == Up) ||
         (_type == GenS && s2._type == GenD) ||
         (_type == GenD && s2._type == GenS) ) 
      return Return::Delete;
    // downgrade the spin-sum
    if ( _type == Up || _type == Down || s2._type == Gen ) {
      return Return::Repeat;
    }
    if (s2._type == Up){
      *this = s2;
      return Return::Done;
    }
    if (s2._type == Down){
      if (_type == GenS){
        *this = s2;
        return Return::Done;
      } else {
        *this = s2;
        return Return::Change_sign;
      }
    }
    error("What am I doing here?","Spin::replace");
  }
  return Return::Done;
}

std::ostream& operator<<(std::ostream& o, const Spin& spin)
{
  switch (spin.type()) {
    case Spin::Up:
      o << "\\alpha";
      break;
    case Spin::Down:
      o << "\\beta";
      break;
    case Spin::GenD:
      o << "\\bar";
    case Spin::GenS:
      o << "\\sigma_{" << spin.el() << "}";
      break;
    default:
      break;
  }
  return o;
}

OrbitalTypes::OrbitalTypes(const std::string& types, lui beg, lui end, bool occ)
{
  if ( beg >= types.size()-1 ) return;
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  lui ipos = beg, ipos1;
  ipos = beg;
  ipos = IL::skip(types,ipos,"{}_^ ");
  while ( (ipos < end) && (ipos1 = IL::nextwordpos(types,ipos,true,false)) != ipos ){//non greedy
    Orbital orb(IL::plainname(types.substr(ipos,ipos1-ipos)),spintype);
    this->push_back(orb.type());
    if ( occ && orb.type() == Orbital::Virt ) {
      xout << "WARNING: Do you really want to have orbital " << orb << " as occupied?" << std::endl;
    } else if ( !occ && orb.type() == Orbital::Occ ) {
      xout << "WARNING: Do you really want to have orbital " << orb << " as virtual?" << std::endl;
    }
    ipos = IL::skip(types,ipos1,"{}_^ ");
  }
}

std::ostream & operator << (std::ostream & o, const Electrons& el)
{
  o << el.name().at(0);
  if (el.name().size()>1)
    o <<"_{" << el.name().substr(1) << "}";
  return o;
}
