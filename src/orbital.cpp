#include "orbital.h"

Orbital::Orbital()
{
  _type = Orbital::NoType;
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
Orbital::Orbital(Orbital::Type type)
{
  TsPar& orbs = Input::sPars["syntax"];
  _spin = Spin(Spin::No);
  _type = type;
  switch( type ){
    case Occ:
      _name = orbs["occorb"][0];
      break;
    case Virt:
      _name = orbs["virorb"][0];
      break;
    case GenT:
      _name = orbs["genorb"][0];
      break;
    case Act:
      _name = orbs["actorb"][0];
      break;
    default:
      Error("Unknown type of orbital!");
  }
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

std::string Orbital::letname() const
{
  lui iend;
  for ( iend = this->_name.size(); iend > 0 && isdigit(this->_name[iend-1]); --iend ){}
  return _name.substr(0,iend);
}
void Orbital::replace_letname(const std::string& newname)
{
  lui iend;
  for ( iend = this->_name.size(); iend > 0 && isdigit(this->_name[iend-1]); --iend ){}
  _name.replace(0,iend,newname);
}
void Orbital::add_prime()
{
  std::string orbname, up, down;
  IL::nameupdown(orbname,up,down,_name);
  up += "\\prime";
  _name = orbname;
  if (!down.empty()) _name += "_{"+down+"}";
  _name += "^{"+up+"}";
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
  for (const auto& orb: orbset){
    if (this->comp_letname(orb) == 0)
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
  (void)smart; //unused variable
  return this->_spin.replace(spin1,spin2);
}

std::ostream & operator << (std::ostream & o, TOrbSet const & orbset)
{
  bool explspin = Input::iPars["prog"]["explspin"];
  for( auto it : orbset ){
    if (explspin) o << "{";
    o << it;
    if (explspin) o << "}";
  }
  return o;
}

std::ostream & operator << (std::ostream & o, Orbital const & orb)
{
  std::string orbname, up, down;
  IL::nameupdown(orbname,up,down,orb.name());
  if ( orb.spin().type() == Spin::Gen )
    o << char(std::toupper(orbname.at(0)));
  else
    o << orbname.at(0);

  if (orbname.size() > 1 || !down.empty())
    o <<"_{" << orbname.substr(1) << down << "}";
  if ( !up.empty() )
    o << "^{" << up << "}";
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

Spin::Type Spin::totype(const std::string& spinpar)
{
  if ((spinpar == "\\alpha") || (spinpar == "\\Up")){return Spin::Up;}
  if ((spinpar == "\\beta") || (spinpar == "\\Down")){return Spin::Down;}
  error("Could not identify spintype.");
  return Spin::Up;
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
      // fall through
    case Spin::GenS:
      o << "\\sigma_{" << spin.el() << "}";
      break;
    default:
      break;
  }
  return o;
}

OrbitalTypes::OrbitalTypes(const std::string& types, bool occ)
{
  if ( types.empty() ) return;
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  lui ipos, ipos1;
  ipos = IL::skip(types,0,"{}_^ ");
  while ( (ipos < types.size()) && (ipos1 = IL::nextwordpos(types,ipos,true,false)) != ipos ){//non greedy
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

std::ostream & operator << (std::ostream & o, const OrbitalTypes& orbt)
{
  for (const auto& ot: orbt ){
    o << Orbital(ot);
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
