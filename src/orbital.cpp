#include "orbital.h"

const std::string Orbital::virt_def="abcdefgh";
const std::string Orbital::occ_def="ijklmno";
const std::string Orbital::gen_def="pqrstuvwxyz";

Orbital::Orbital() {}

Orbital::Orbital(std::string name)
{
  _spin=(isupper((char)name[0]) ? GenS : No);
  name[0] = std::tolower(name[0]);
  if (virt_def.find(name[0])!=std::string::npos) {_type=Virt;}
  else if (occ_def.find(name[0])!=std::string::npos) {_type=Occ;}
  else if (gen_def.find(name[0])!=std::string::npos) {_type=GenT;}
  else 
    error("Unknown type of orbital! "+name,"Orbital::Orbital"); 
  _name=name;
}

Orbital::Orbital(std::string name, Orbital::Type type)
{
  _spin=(isupper((char)name[0]) ? GenS : No);
  name[0] = std::tolower(name[0]);
  _type=type;
  _name=name;
}

Orbital::Orbital(std::string name, Orbital::Spin spin)
{
  _spin=spin;
  name[0] = std::tolower(name[0]);
  if (virt_def.find(name[0])!=std::string::npos) {_type=Virt;}
  else if (occ_def.find(name[0])!=std::string::npos) {_type=Occ;}
  else if (gen_def.find(name[0])!=std::string::npos) {_type=GenT;}
  else 
    error("Unknown type of orbital! "+name,"Orbital::Orbital"); 
  _name=name;
}

Orbital::Orbital(std::string name, Orbital::Type type, Orbital::Spin spin)
{
  _spin=spin;
  _type=type;
  name[0] = std::tolower(name[0]);
  _name=name;
}

std::string Orbital::name() const
{ return _name;}

Orbital::Type Orbital::type() const
{ return _type;}

Orbital::Spin Orbital::spin() const
{ return _spin;}
void Orbital::setspin(Orbital::Spin spin)
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
  return _spin<orb._spin;
}

std::ostream & operator << (std::ostream & o, Orbital const & orb)
{
  if ( orb.spin()==Orbital::GenS )
    o << char(std::toupper(orb.name().at(0)));
  else
    o << orb.name().at(0);
  if (orb.name().size()>1)
    o <<"_{" << orb.name().substr(1) << "}";
  if ( orb.spin()==Orbital::Up )
    o << "\\alpha";
  else if ( orb.spin()==Orbital::Down )
    o << "\\beta";
/*  if ( orb.type()==Orbital::Occ ) 
    o << "^occ";
  else if ( orb.type()==Orbital::Virt )
    o << "^virt";
  else
    o << "^gen";
*/
  return o;
}

