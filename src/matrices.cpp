#include "matrices.h"

Product< Orbital > Ops::genprodorb(short int exccl, const Orbital& occ, const Orbital& virt)
{
  Product<Orbital> porbs;
  std::string excl;
  Orbital orb;

  for (short i=0; i<exccl ; ++i)
  {  
    if (i>0) excl=num2str(i,std::dec);
    orb=Orbital(virt.name()+excl,virt.spin());
    porbs*=orb;
    orb=Orbital(occ.name()+excl,occ.spin());
    porbs*=orb;
  }
  return porbs;
}

Matrices::Matrices() // : _type(Interm)
{_antisymform=false;}
Matrices::Matrices(Ops::Type t, Product< Orbital > p, std::string name, Matrices::Spinsym matspinsym)
{
  _type=t;
  _orbs=p;
  if (t==Ops::Fock)
    _name="Fock";
  else if (t==Ops::FluctP)
    _name="FluctP";
  else if (t==Ops::XPert)
    _name="XPert";
  else
    _name=name;
  _matspinsym=matspinsym;
  if (t==Ops::FluctP)
    _antisymform=true;
  else
    _antisymform=false;
}
Ops::Type Matrices::type() const
{ return _type; }
Product< Orbital > Matrices::orbitals() const
{ return _orbs; }
std::string Matrices::name() const
{ return _name; }
bool Matrices::antisymform() const
{ return _antisymform; }

void Matrices::replace(Orbital orb1, Orbital orb2)
{
  for ( unsigned int i=0; i<_orbs.size(); ++i )
  {
    if (_orbs[i]==orb1) _orbs[i]=orb2;
  }
}
bool Matrices::expandantisym(bool firstpart)
{
  if (_antisymform)
  {
    if (_type!=Ops::FluctP)
      error("Can not expand antisymmetrical non-integral","Matrices::expandantisym");
    if (_orbs[0].spin()==Orbital::No)
      error("Can not expand antisymmetrical integral in space orbitals","Matrices::expandantisym");
    _antisymform=false;
    if (!firstpart)
    { // (PQ|RS) -> (PS|RQ)
      std::swap(_orbs[1],_orbs[3]);
    }
  }
  else
    return false;
  return true;
}
bool Matrices::operator<(const Matrices& t) const
{
  if ( _type < t._type ) return true;
  if ( t._type < _type ) return false;
  if ( _name < t._name ) return true;
  if ( t._name < _name ) return false;
  if (_type == Ops::FluctP)
  { // electron-symmetry
    if ( *this == t ) return false; // the Matrices are the same
  }
  return _orbs < t._orbs;
}
bool Matrices::operator==(const Matrices& t) const
{
  if ( _type != t._type || _name != t._name ) return false;
  if ( _orbs.size() != t._orbs.size() ) return false;
  if ( _orbs == t._orbs ) return true;
  if (_type == Ops::FluctP)
  { // electron-symmetry
    if (_orbs.subprod(0,1) == t._orbs.subprod(2,3) && _orbs.subprod(2,3) == t._orbs.subprod(0,1) ) return true;
  }
//  else if (_type == Ops::Exc || _type == Ops::Deexc)
  else if (InSet(_type, Ops::Exc,Ops::Deexc))
  { // electron-symmetry
 
   // for (unsigned int i=0; i<_orbs.size()/2; i++)
   // {
   //   
   //   for (unsigned int j=i; j<t._orbs.size()/2; j++)
   //   {
   //     if (this->spinsym(i*2) != t.spinsym(j*2)) break;
   //     
   //   }
   // }
  }
  return false;
}
void Matrices::reset_vertices()
{ _indx=-1; }

bool Matrices::vertices(long int ipos, Matrices& mat, long int ipos1, unsigned int indx)
{
  // compare types, excitation classes and index of matrix 
  if (_type != mat._type || _name != mat._name || _orbs.size() != mat._orbs.size() || _indx != mat._indx) return false;
  // in the case of external indices orbitals should match exactly
  //if ((_type == Ops::Exc0 || _type == Ops::Deexc0)&&_orbs[ipos]!=mat._orbs.at(ipos)) return false;
  // compare spin symmetries
  if (spinsym(ipos) != mat.spinsym(ipos1)) return false;
  // set index of matrix (if not set)
  if (_indx < 0) _indx=mat._indx=indx;
  return true;
}
void Matrices::set_connect(TCon2 connected2 )
{ _connected2=connected2; }
void Matrices::add_connect(long int con)
{ _connected2.insert(con); }
TCon2 Matrices::connected2() const
{ return _connected2; }

void Matrices::setkind(short int exccl, short int intlines, short int intvirt)
{
  _exccl=exccl;
  _intlines=intlines;
  _intvirt=intvirt;
}
Orbital Matrices::orbel(const Orbital& orb)
{
  Orbital orb1;
  bool next=true;
  for (unsigned int i=0; i<_orbs.size(); ++i)
  {
    if (_orbs[i]==orb) 
    {
      if(next)
        orb1=_orbs[i+1];
      else
        orb1=_orbs[i-1];
      break;
    }
    next=!next;
  }
  return orb1;
}
Orbital Matrices::orbel(const long int& ipos)
{
  assert( 0 <= ipos && unsigned(ipos) < _orbs.size() );
//  if (ipos<0 || unsigned(ipos)>=_orbs.size()) return Orbital();
  return (ipos%2==0?_orbs[ipos+1]:_orbs[ipos-1]);
}

Matrices::Spinsym Matrices::spinsym(long int ipos)
{
  if (_matspinsym==Triplet && ipos-2 <0) //first electron is triplet
    return Triplet;
  return Singlet;
}
void Matrices::set_no_spin()
{
  for (unsigned int i=0; i<_orbs.size(); ++i)
    _orbs[i].setspin(Orbital::No);
}

std::ostream & operator << (std::ostream & o, Matrices const & mat)
{
  if (mat.type()==Ops::Fock)
  {
    o << "f_{" << mat.orbitals() << "}";
    MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
  }
  else if (mat.type()==Ops::FluctP)
  {
    o << "(" << mat.orbitals().subprod(0,1) << "|" << mat.orbitals().subprod(2,3) << ")";
    MyOut::pcurout->lenbuf += 7;
  }
  else if (mat.type()==Ops::XPert)
  {
    o << "x_{" << mat.orbitals() << "}";
    MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
  }
  else if (InSet(mat.type(), Ops::Exc,Ops::Deexc,Ops::Interm))
  {
    o <<mat.name() <<"_{" << mat.orbitals() << "}";
    MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
  }
  else if (mat.type()==Ops::Number)
  {
    o <<mat.name();
    MyOut::pcurout->lenbuf += 2;
  }
  //o <<"["<<mat.connected2()<<"]" ;
  return o;
}

Permut::Permut()
{}
Permut::Permut(Product< Orbital > p1, Product< Orbital > p2) 
{
  assert( p1.size() == p2.size() );
  for ( lui io = 0; io < p1.size(); ++io ){
    assert( p1[io].type() == p2[io].type() );
    _orbs[ p1[io] ] = p2[io];
  }
}
Permut::Permut(Orbital o1, Orbital o2)
{
  assert( o1.type() == o2.type() );
  _orbs[o1] = o2;
}
Permut& Permut::operator+=(const Permut& p)
{
  for ( TPerMap::const_iterator pit = p._orbs.begin(); pit != p._orbs.end(); ++pit ){
    TPerMap::iterator it = _orbs.find(pit->first);
    if ( it == _orbs.end()){
      // the orbital not there yet
      _orbs[pit->first] = pit->second;
    } else { 
      // have to be completely the same
      assert( it->second == pit->second );
    }
  }
  return *this;
}
Permut& Permut::operator*=(const Permut& p)
{
  TPerMap orbs = _orbs;
  _orbs.clear();
  for ( TPerMap::const_iterator pit = p._orbs.begin(); pit != p._orbs.end(); ++pit ){
    TPerMap::iterator it = orbs.find(pit->second);
    if ( it == orbs.end()){
      // the orbital not there yet
      _orbs[pit->first] = pit->second;
    } else { 
      // replace it
      if ( pit->first != it->second ) // don't add permutation P(p,p)
        _orbs[pit->first] = it->second;
      orbs.erase(it);
    }
  }
  // add the rest
  _orbs.insert(orbs.begin(),orbs.end());
  return *this;
}

Product< Orbital > Permut::orbsfrom() const
{ 
  Product< Orbital > orbs;
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it )
    orbs.push_back(it->first);
  return orbs; 
  
}
Product< Orbital > Permut::orbsto() const
{ 
  Product< Orbital > orbs;
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it )
    orbs.push_back(it->second);
  return orbs;
}
bool Permut::operator<(const Permut& p) const
{
  if (_orbs.size()<p._orbs.size()) return true;
  if (p._orbs.size()<_orbs.size()) return false;
//  if ((*this)==p) return false;
  TPerMap::const_iterator pit = p._orbs.begin();
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it, ++pit ){
    if ( it->first < pit->first ) return true;
    if ( pit->first < it->first ) return false;
    if ( it->second < pit->second ) return true;
    if ( pit->second < it->second ) return false;
  }
  return false;
}
bool Permut::operator==(const Permut& p) const
{
  TPerMap::const_iterator pit = p._orbs.begin();
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it, ++pit )
    if ( it->first != pit->first || it->second != pit->second ) return false;
  return true;
}


std::ostream & operator << (std::ostream & o, Permut const & p)
{
  if (p.orbsfrom().size()>0) {
    o << "\\Perm{" << p.orbsfrom() << "," << p.orbsto() << "}";
    MyOut::pcurout->lenbuf += 3+(p.orbsfrom().size()+p.orbsto().size())/MyOut::pcurout->wsi;
    //o << "[" << MyOut::pcurout->lenbuf << "]";
  } else 
    o << "1";
  return o;
}

