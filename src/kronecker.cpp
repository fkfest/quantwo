#include "kronecker.h"

Kronecker::Kronecker(Orbital orb1, Orbital orb2) :
    _orb1(orb1), _orb2(orb2)
{
  if ( _orb2 < _orb1 ) std::swap(_orb1,_orb2);
  if ( _orb2.name().size() < _orb1.name().size() ) std::swap(_orb1,_orb2);
  if (_orb1.type()==Orbital::GenT) std::swap(_orb1,_orb2);
}

Orbital Kronecker::orb1() const
{   return _orb1;   }

Orbital Kronecker::orb2() const
{   return _orb2;   }

bool Kronecker::is_ordered(const Product< Orbital >& crobs, const Product< Orbital >& anobs) const
{
  bool orb1anni = ( crobs.find(_orb1) >= 0 );
  bool orb2anni = ( crobs.find(_orb2) >= 0 );
  bool orb1crea = ( anobs.find(_orb1) >= 0 );
  bool orb2crea = ( anobs.find(_orb2) >= 0 );
  if ( (orb1anni && orb2anni) || (orb1crea && orb2crea) ) {
    xout << "crobs: " << crobs << std::endl;
    xout << "anobs: " << anobs << std::endl;
    xout << "Kronecker: " << *this << std::endl;
    error("Type mismatch in crobs and anobs in term for Kronecker");
  }
  bool ordered = false;
  if ( orb1anni || orb2crea ) {
    ordered = false;
  } else if ( orb2anni || orb1crea ) {
    ordered = true;
  } else {
    xout << "crobs: " << crobs << std::endl;
    xout << "anobs: " << anobs << std::endl;
    xout << "Kronecker: " << *this << std::endl;
    error("Cannot guess the order of orbitals in Kronecker");
  }
  return ordered;
}

bool Kronecker::operator < (Kronecker const & k) const
{
    if ( _orb1<k._orb1 )
        return true;
    if ( k._orb1<_orb1 )
        return false;
    return _orb2<k._orb2;
}

Return Kronecker::replace(Orbital orb1, Orbital orb2, bool smart)
{
  Return rpl;
  rpl += _orb1.replace(orb1,orb2,smart);
  rpl += _orb2.replace(orb1,orb2,smart);
  return rpl;
}
Return Kronecker::replace(Spin spin1, Spin spin2, bool smart)
{
  Return rpl;
  rpl += _orb1.replace(spin1,spin2,smart);
  rpl += _orb2.replace(spin1,spin2,smart);
  return rpl;
}

std::ostream & operator << (std::ostream & o, Kronecker const & k)
{
    o << "\\delta_{" << k.orb1() << k.orb2() << "}";
    MyOut::pcurout->lenbuf += 1+2/MyOut::pcurout->wsi;
    return o;
}

