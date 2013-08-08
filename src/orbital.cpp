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
Return Orbital::replace(const Orbital& orb1, const Orbital& orb2)
{
  if (*this == orb1) { 
    *this = orb2; 
    // restore spin -- will be replaced explicitely 
    this->_spin = orb1._spin;
  }
  return this->_spin.replace(orb1._spin,orb2._spin);
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


std::ostream & operator << (std::ostream & o, const Electrons& el)
{
  o << el.name().at(0);
  if (el.name().size()>1)
    o <<"_{" << el.name().substr(1) << "}";
  return o;
}
