#include "orbital.h"

Orbital::Orbital() 
{ 
  _type = Orbital::GenT;
  _spin = Spin(Spin::No);
}

Orbital::Orbital(const std::string& name)
{
  bool genspin = isupper((char)name[0]);
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
  if ( genspin )
    _spin = Spin(_name, Spin::GenS);
  else
    _spin = Spin(Spin::No);
  
}

Orbital::Orbital(const std::string& name, Orbital::Type type)
{
  bool genspin = isupper((char)name[0]);
  _name=name;
  _name[0] = std::tolower(name[0]);
  _type=type;
  if ( genspin )
    _spin = Spin(_name, Spin::GenS);
  else
    _spin = Spin(Spin::No);
}

Orbital::Orbital(const std::string& name, Spin spin)
{
  _spin=spin;
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
}
Orbital::Orbital(const std::string& name, Spin::Type spint)
{
  _name=name;
  _name[0] = std::tolower(name[0]);
  gentype();
  _spin = Spin(name,spint);
}
Orbital::Orbital(const std::string& name, Orbital::Type type, Spin spin)
{
  _spin=spin;
  _type=type;
  _name=name;
  _name[0] = std::tolower(name[0]);
}
Orbital::Orbital(const std::string& name, Orbital::Type type, Spin::Type spint)
{
  _type=type;
  _name=name;
  _name[0] = std::tolower(name[0]);
  _spin = Spin(name,spint);
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

std::ostream & operator << (std::ostream & o, Orbital const & orb)
{
  if ( orb.spin().type() == Spin::GenS )
    o << char(std::toupper(orb.name().at(0)));
  else
    o << orb.name().at(0);
  if (orb.name().size() > 1)
    o <<"_{" << orb.name().substr(1) << "}";
  if ( orb.spin().type() != Spin::GenS )
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
  return (_type == spin._type && _name == spin._name);
}

bool Spin::operator<(const Spin& spin) const
{
  if (_name < spin._name) return true;
  if (spin._name < _name) return false;
  return  (_type < spin._type);
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
      o << "\\sigma_{" << spin.name() << "}";
      break;
    default:
      break;
  }
  return o;
}


std::ostream & operator << (std::ostream & o, const Electron& el)
{
  o << el.name().at(0);
  if (el.name().size()>1)
    o <<"_{" << el.name().substr(1) << "}";
  return o;
}
