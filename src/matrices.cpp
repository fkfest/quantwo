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
{
  _antisymform=false;
  _type = Ops::None;
}
Matrices::Matrices(Ops::Type t, Product< Orbital > p, short npairs, std::string name, Matrices::Spinsym matspinsym, bool antisymW)
{
  _type=t;
  _orbs=p;
  _npairs=npairs;
  assert(_npairs <= int(_orbs.size()/2));
  switch (t) {
    case Ops::Fock:
      _name = "Fock";
      break;
    case Ops::FluctP:
      _name = "FluctP";
      break;
    case Ops::XPert:
      _name = "XPert";
      break;
    case Ops::DensM:
      _name = "Gamma";
      break;
    default:
      _name=name;
  }
  _matspinsym=matspinsym;
  if (t==Ops::FluctP)
    _antisymform=antisymW;
  else if ( t==Ops::Exc && Input::iPars["prog"]["quan3"] > 0 ) {
    // make amplitudes antisymmetrical (now works only for doubles!)
    assert(_orbs.size() == 4);
    _antisymform=true;
  } else
    _antisymform=false;
}
Matrices::Matrices(const Kronecker& d)
{
  _type = Ops::Delta;
  _orbs.push_back(d.orb1());
  _orbs.push_back(d.orb2());
  _npairs = 1;
  _name = "delta";
  _matspinsym = Singlet;
  _antisymform = false;
}

Ops::Type Matrices::type() const
{ return _type; }
const Product< Orbital >& Matrices::orbitals() const
{ return _orbs; }

std::string Matrices::name() const
{ return _name; }
bool Matrices::antisymform() const
{ return _antisymform; }

Return Matrices::replace(Orbital orb1, Orbital orb2, bool smart)
{
  Return rpl;
  for ( unsigned int i=0; i<_orbs.size(); ++i )
  {
    rpl += _orbs[i].replace(orb1,orb2,smart);
  }
  return rpl;
}
bool Matrices::expandantisym(bool firstpart)
{
  if (_antisymform) {
    // works atm for doubles only!
    assert(_orbs.size() == 4);
    if (_type!=Ops::FluctP && !( _type==Ops::Exc && Input::iPars["prog"]["quan3"] > 0 ))
      error("Can not expand antisymmetrical non-integral","Matrices::expandantisym");
    if (_orbs[0].spin().type()==Spin::No)
      error("Can not expand antisymmetrical integral in space orbitals","Matrices::expandantisym");
    _antisymform=false;
    if (!firstpart) { // (PQ|RS) -> (PS|RQ)
      std::swap(_orbs[1],_orbs[3]);
    }
  }
  else
    return false;
  return true;
}
bool Matrices::nonsingldm() const
{
  if ( _type != Ops::DensM ) return false;
  assert( _orbs.size()%2 == 0 );
  assert( _orbs.size() == _cranorder.size() );
  bool dmsort = (Input::iPars["prog"]["dmsort"] > 0);
  if (dmsort) {
    // a^\dg(1)  a^\dg(2) ... a(2) a(1)
    for ( uint i = 0, j = _orbs.size()-1; i < _orbs.size()/2; ++i, --j ){
      if (_orbs[i].spin() != _orbs[j].spin() || _cranorder[i] != SQOpT::Creator || _cranorder[j] != SQOpT::Annihilator ) return true;
    }
  } else {
    // a^\dg(1) a(1) a^\dg(2) a(2)...
    for ( uint i = 0; i < _orbs.size()/2; i+=2 ){
      if (_orbs[i].spin() != _orbs[i+1].spin() || _cranorder[i] != SQOpT::Creator || _cranorder[i+1] != SQOpT::Annihilator ) return true;
    }
  }
  return false;
}
bool Matrices::operator<(const Matrices& t) const
{
  if ( _type < t._type ) return true;
  if ( t._type < _type ) return false;
  if ( _name < t._name ) return true;
  if ( t._name < _name ) return false;
  if (_type == Ops::FluctP) { // electron-symmetry
    if ( *this == t ) return false; // the Matrices are the same
  }
  return _orbs < t._orbs;
}
bool Matrices::operator==(const Matrices& t) const
{
  if ( _type != t._type || _name != t._name ) return false;
  if ( _orbs.size() != t._orbs.size() ) return false;
  if ( _orbs == t._orbs ) return true;
  if (_type == Ops::FluctP) { // electron-symmetry
    if (_orbs.subprod(0,1) == t._orbs.subprod(2,3) && _orbs.subprod(2,3) == t._orbs.subprod(0,1) ) return true;
  }// else if (InSet(_type, Ops::Exc,Ops::Deexc)) { // electron-symmetry
 
   // for (unsigned int i=0; i<_orbs.size()/2; i++) {
   //   for (unsigned int j=i; j<t._orbs.size()/2; j++) {
   //     if (this->spinsym(i*2) != t.spinsym(j*2)) break;
   //   }
   // }
  //}
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

void Matrices::set_conlines(Product< ConLine > conlines)
{ _conlines = conlines; }
void Matrices::set_conline(lui iorb, lui imat, lui idx)
{
  assert( iorb < _orbs.size() );
  if ( _conlines.size() == 0 )// first resize _conlines
    _conlines.resize(_orbs.size());
  _conlines[iorb] = ConLine(imat,idx);
}
void Matrices::add_conline(lui imat, lui idx)
{ _conlines.push_back(ConLine(imat,idx)); }
const Product< ConLine >& Matrices::conlines() const
{ return _conlines; }
const ConLine& Matrices::conline(lui iorb) const
{
  assert( iorb < _conlines.size() );
  return _conlines[iorb]; 
}
bool Matrices::is0() const
{
  if (_type == Ops::DensM){
    if ( _orbs.size()%2 != 0 ) return true;
    // only active orbitals!
    _foreach_cauto(Product<Orbital>,ito,_orbs){
      if (ito->type() != Orbital::Act) return true;
    }
    // number of creators == annihilators
    int cran = 0;
    _foreach_cauto(Product<SQOpT::Gender>,itca,_cranorder){
      if (*itca == SQOpT::Creator) 
        ++cran;
      else 
        --cran;
    }
    if (cran != 0 ) return true;
  }
  return false;
}

void Matrices::set_cran(const Product< SQOpT::Gender >& cran)
{
  assert( _type == Ops::DensM );
  assert( cran.size() == _orbs.size() );
  _cranorder = cran;
}

void Matrices::setkind(short int exccl, short int intlines, short int intvirt)
{
  _exccl=exccl;
  _intlines=intlines;
  _intvirt=intvirt;
}
long int Matrices::iorbel(lui ipos)
{ 
  assert( ipos < _orbs.size() );
  lui ipos1;
  if ( _type == Ops::DensM && Input::iPars["prog"]["dmsort"] > 0 )
    // different order of electrons: 1,2,...2, 1
    ipos1 = _orbs.size()-ipos-1;
  else
    ipos1 = ipos%2==0?ipos+1:ipos-1;
  
  if ( (short)ipos1 >= 2*_npairs )
    // one of the non-conserved electrons
    return -1;
  return ipos1;
}

Orbital Matrices::orbel(const Orbital& orb)
{
  long int ipos = _orbs.find(orb);
  if ( ipos >= 0 )
    return orbel(ipos);
  else
    return Orbital();
}
Orbital Matrices::orbel(const long int& ipos)
{
  long ipos1 = iorbel(ipos);
  if (ipos1 >= 0 )
    return _orbs[ipos1]; 
  else 
    return Orbital();
}

Matrices::Spinsym Matrices::spinsym(long int ipos)
{
  if (_matspinsym==Triplet && ipos-2 <0) //first electron is triplet
    return Triplet;
  return Singlet;
}
void Matrices::set_no_spin()
{
  Spin nospin(Spin::No);
  for (unsigned int i=0; i<_orbs.size(); ++i)
    _orbs[i].setspin(nospin);
}
void Matrices::permute(const Permut& p)
{
  for (unsigned int i=0; i<_orbs.size(); ++i){
    _orbs[i] = p.permutorb(_orbs[i]);
  }
}
std::string Matrices::plainname() const
{
  const std::string& plainsymbols = Input::sPars["syntax"]["plainsymb"];
  std::string plainnam;
  for ( uint i = 0; i < _name.size(); ++i ){
    if ( plainsymbols.find(_name[i]) != std::string::npos ){
      plainnam += _name[i];
    }
  }
  return plainnam;
}

std::ostream & operator << (std::ostream & o, Matrices const & mat)
{
  short clean = Input::iPars["output"]["clean"];
  std::string param("");
  if ( clean <= 0 ) param = "\\"+Input::sPars["command"]["parameter"] + " ";
  switch ( mat.type() ){
    case Ops::Fock:
      o << param << "f_{" << mat.orbitals() << "}";
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    case Ops::FluctP:
      if ( clean > 1 )
        o << "(" << mat.orbitals().subprod(0,1) << "|" << mat.orbitals().subprod(2,3) << ")";
      else 
        o << param << "\\" << Input::sPars["command"]["integral"] << "{" 
                   << mat.orbitals().subprod(0,1) << "}{" << mat.orbitals().subprod(2,3) << "}";
      MyOut::pcurout->lenbuf += 7;
      break;
    case Ops::XPert:
      o << param << "X_{" << mat.orbitals() << "}";
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    case Ops::DensM: {
      Product<Orbital> orbs(mat.orbitals());
      std::string name("\\"+Input::sPars["command"]["densmat"]);
      if ( mat.nonsingldm()){
        const Product<SQOpT::Gender> cran = mat.get_cran();
        assert(orbs.size() == cran.size());
        o << name << "^{}_{";
        for ( uint i = 0; i < orbs.size(); ++i ){
          o << orbs[i];
          if ( cran[i] == SQOpT::Creator ) o << "^{\\dg}";
        }
        o << "}";
      } else {
        std::ostringstream oss;
        bool dmsort = (Input::iPars["prog"]["dmsort"] > 0);
        // occ.indices
        if (dmsort) {
          for ( uint i = 0; i < orbs.size()/2; ++i ){
            oss << orbs[i];
          }
        } else {
          for ( uint i = 0; i < orbs.size(); i += 2 ){
            oss << orbs[i];
          }
        }
        IL::add2name(name,oss.str(),true);
        oss.str("");
        //virt. indices
        if (dmsort) {
          for ( uint i = orbs.size()-1; i >= orbs.size()/2 ; --i ){
            oss << orbs[i];
          }
        } else {
          for ( uint i = 1; i < orbs.size(); i += 2 ){
            oss << orbs[i];
          }
        }
        IL::add2name(name,oss.str(),false);
        o << param << name;
        MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      }
      break;
    }
    case Ops::Delta:
      o << param << "\\delta_{" << mat.orbitals() << "}";
      break;
    case Ops::Exc:
    case Ops::Deexc:
    case Ops::Interm: {
      Product<Orbital> orbs(mat.orbitals());
      std::string name(mat.name());
      std::ostringstream oss;
      // occ.indices
      for ( uint i = 1; i < orbs.size(); i += 2 ){
        oss << orbs[i];
      }
      IL::add2name(name,oss.str(),true);
      oss.str("");
      //virt. indices
      for ( uint i = 0; i < orbs.size(); i += 2 ){
        oss << orbs[i];
      }
      IL::add2name(name,oss.str(),false);
      o << param << name;
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    }
    case Ops::Number:
      o << param << mat.name();
      MyOut::pcurout->lenbuf += 2;
      break;
    case Ops::None:
    case Ops::Exc0:
    case Ops::Deexc0:
      break;
  }
  //o <<"["<<mat.connected2()<<"]" ;
  return o;
}

Permut::Permut()
{ dummy = 1;}
Permut::Permut(Product< Orbital > p1, Product< Orbital > p2) 
{
  dummy = 1;
  assert( p1.size() == p2.size() );
  for ( lui io = 0; io < p1.size(); ++io ){
    assert( p1[io].type() == p2[io].type() );
    _orbs[ p1[io] ] = p2[io];
  }
}
Permut::Permut(Orbital o1, Orbital o2)
{
  dummy = 1;
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
Permut& Permut::operator/=(const Permut& p)
{
  // generate a reverse map
  TPerMap orbs;
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it ){
    orbs[it->second] = it->first;
  }
  _orbs.clear();
  for ( TPerMap::const_iterator pit = p._orbs.begin(); pit != p._orbs.end(); ++pit ){
    TPerMap::iterator it = orbs.find(pit->second);
    if ( it == orbs.end()){
      // the orbital is not there
      _orbs[pit->second] = pit->first;
    } else {
      if ( pit->first != it->second ) // don't add permutation P(p,p)
        _orbs[it->second] = pit->first;
      orbs.erase(it);
    }
  }
  // add the rest
  for ( TPerMap::const_iterator it = orbs.begin(); it != orbs.end(); ++it ){
    if ( _orbs.find(it->second) != _orbs.end() ) error("Mismatch in permutaions!");
    _orbs[it->second] = it->first;
  }
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
Orbital Permut::permutorb(const Orbital& orb) const
{
  TPerMap::const_iterator it = _orbs.find(orb);
  if ( it == _orbs.end()){
    // the orbital is not there
    return orb;
  } else {
    return it->second;
  }
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
    o << "\\" << Input::sPars["command"]["permutation"] << "{" << p.orbsfrom() << "}{" << p.orbsto() << "}";
    MyOut::pcurout->lenbuf += 3+(p.orbsfrom().size()+p.orbsto().size())/MyOut::pcurout->wsi;
    //o << "[" << MyOut::pcurout->lenbuf << "]";
  } else 
    o << "1";
  return o;
}
