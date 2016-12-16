#include "matrix.h"

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

Matrix::Matrix() // : _type(Interm)
{
  create_Matrix(Ops::None,0,0,0,"",Singlet,false);
}
Matrix::Matrix(Ops::Type t, const Product< Orbital >& p, uint npairs, short lmel, short pmsym,
                   std::string name, Matrix::Spinsym matspinsym, bool antisymW)
{
  _orbs = p;
  create_Matrix(t,npairs,lmel,pmsym,name,matspinsym,antisymW);
}
Matrix::Matrix(Ops::Type t, const Product< Orbital >& pcrea, const Product< Orbital >& panni, 
               uint npairs, short int lmel, short pmsym, std::string name, 
               Matrix::Spinsym matspinsym, bool antisymW)
{
  assert( pcrea.size() == panni.size()+lmel );
  _orbs.reserve(pcrea.size()+panni.size());
  if ( t == Ops::DensM && Input::iPars["prog"]["dmsort"] > 0 ) {
    // different order of electrons: 1,2,...2, 1
    _foreach_cauto(Product<Orbital>,itorb,pcrea){
      _orbs.push_back(*itorb);
      _cranorder.push_back(SQOpT::Creator);
    }
    _foreach_crauto(Product<Orbital>,itorb,panni){
      _orbs.push_back(*itorb);
      _cranorder.push_back(SQOpT::Annihilator);
    }
    
  } else {
    uint 
      maxlen = std::max(pcrea.size(), panni.size());
    for ( uint iorb = 0; iorb < maxlen; ++iorb ){
      if ( iorb < pcrea.size() ) {
        _orbs.push_back(pcrea[iorb]);
        if (t == Ops::DensM) _cranorder.push_back(SQOpT::Creator);
      }
      if ( iorb < panni.size() ) {
        _orbs.push_back(panni[iorb]);
        if (t == Ops::DensM) _cranorder.push_back(SQOpT::Annihilator);
      }
    }
  }
  create_Matrix(t,npairs,lmel,pmsym,name,matspinsym,antisymW);
}

Matrix::Matrix(const Kronecker& d)
{
  _orbs.push_back(d.orb1());
  _orbs.push_back(d.orb2());
  create_Matrix(Ops::Delta,1,0,0,"",Singlet,false);
}
void Matrix::create_Matrix(Ops::Type t, uint npairs, short int lmel, short int pmsym,
                           std::string name, Matrix::Spinsym matspinsym, bool antisymW)
{
  _type=t;
  gen_name(name);
  _npairs=npairs;
  assert(_npairs <= _orbs.size()/2);
  _lmel = lmel;
  _pmsym = pmsym;
  assert(2*_npairs+std::abs(_lmel) == _orbs.size());
  assert((_orbs.size()-_lmel)%2 == 0 && (_orbs.size()+_lmel)%2 == 0 );
  _matspinsym=matspinsym;
  _internal = true;
  if (t==Ops::FluctP)
    _antisymform=antisymW;
  else if ( t==Ops::Exc && Input::iPars["prog"]["quan3"] > 0 ) {
    // make amplitudes antisymmetrical (now works only for doubles!)
    assert(_orbs.size() == 4);
    _antisymform=true;
  } else
    _antisymform=false;
  _exccl = _intlines = _intvirt = _orbtypeshash = 0;
}

void Matrix::gen_name(const std::string& name)
{
  std::string exc0 = Input::sPars["output"]["exc0"];
  switch (_type) {
    case Ops::Fock:
      _name = "Fock";
      break;
    case Ops::OneEl:
      _name = "OneEl";
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
    case Ops::Delta:
      _name = "delta";
    case Ops::Exc0:
    case Ops::Deexc0:
      if (exc0 != " "){
        _name = exc0;
        if ( _type == Ops::Deexc0 ) // add dagger
          IL::add2name(_name,Input::aPars["syntax"]["dg"].front());
        break;
      }
    default:
      _name=name;
  }
}
Ops::Type Matrix::type() const
{ return _type; }
const Product< Orbital >& Matrix::orbitals() const
{ return _orbs; }

std::string Matrix::name() const
{ return _name; }
bool Matrix::antisymform() const
{ return _antisymform; }
bool Matrix::is_internal(const TOrbSet& sumorbs)
{
  _internal = true;
  _foreach_cauto(Product<Orbital>, ito, _orbs){
    if ( sumorbs.count(*ito) == 0 ) {
      _internal = false;
      return false;
    }
  }
  return true;
}

void Matrix::combine(const Matrix& mat, const Set< uint >& dontorbs)
{
  if ( InSet(_type,Ops::Deexc0,Ops::Exc0) )
    _name = mat._name;
  else if ( !InSet(mat._type,Ops::Deexc0,Ops::Exc0) )
    IL::add2name(_name,mat._name);
  if ( dontorbs.size() == mat._orbs.size() )
    // don't have to unite anything
    return;
  if ( _matspinsym != mat._matspinsym || _matspinsym != Singlet )
    Error("Cannot combine non-Singlet matrices yet");
  if ( 2*_npairs != _orbs.size() || 2*mat._npairs != mat._orbs.size() )
    Error("Cannot combine non-conserving matrices yet");
  // combine orbitals (excluding those in dontorbs)
  for ( uint io = 0; io < mat._orbs.size(); ++io ){
    if ( dontorbs.count(io) == 0 ){
      _orbs.push_back(mat._orbs[io]);
      // make sure that we have both orbitals belonging to the same electron
      assert( dontorbs.count(mat.iorbel(io)) == 0 );
    }
  }
  _npairs = _orbs.size()/2;
}

Return Matrix::replace(Orbital orb1, Orbital orb2, bool smart)
{
  Return rpl;
  for ( unsigned int i=0; i<_orbs.size(); ++i )
  {
    rpl += _orbs[i].replace(orb1,orb2,smart);
  }
  return rpl;
}
Return Matrix::replace(Spin spin1, Spin spin2, bool smart)
{
  Return rpl;
  for ( unsigned int i=0; i<_orbs.size(); ++i )
  {
    rpl += _orbs[i].replace(spin1,spin2,smart);
  }
  return rpl;
}
bool Matrix::expandantisym(bool firstpart)
{
  if (_antisymform) {
    // works atm for doubles only!
    assert(_orbs.size() == 4);
    if (_type!=Ops::FluctP && !( _type==Ops::Exc && Input::iPars["prog"]["quan3"] > 0 ))
      error("Can not expand antisymmetrical non-integral","Matrix::expandantisym");
    if (_orbs[0].spin().type()==Spin::No)
      error("Can not expand antisymmetrical integral in space orbitals","Matrix::expandantisym");
    _antisymform=false;
    if (!firstpart) { // (PQ|RS) -> (PS|RQ)
      std::swap(_orbs[1],_orbs[3]);
    }
  }
  else
    return false;
  return true;
}
bool Matrix::nonsingldm() const
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
bool Matrix::operator<(const Matrix& t) const
{
  if ( _type < t._type ) return true;
  if ( t._type < _type ) return false;
  if ( _name < t._name ) return true;
  if ( t._name < _name ) return false;
  if ( _npairs < t._npairs ) return true;
  if ( t._npairs < _npairs ) return false;
  if (_type == Ops::FluctP) { // electron-symmetry
    if ( *this == t ) return false; // the Matrix are the same
  }
  return _orbs < t._orbs;
}
bool Matrix::operator==(const Matrix& t) const
{
  if ( _type != t._type || _name != t._name ) return false;
  if ( _orbs.size() != t._orbs.size() || _npairs != t._npairs || _lmel != t._lmel ) return false;
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
bool Matrix::equivalent(const Matrix& mat) const
{
  return ( _type == mat._type && _name == mat._name &&
           _npairs == mat._npairs && _lmel == mat._lmel && 
           _orbs.size() == mat._orbs.size() && _orbtypeshash == mat._orbtypeshash );
}
Equivalents Matrix::equivertices(uint offs) const
{
  Equivalents everts;
  // no symmetry in the residual
  if ( InSet(_type,Ops::Deexc0,Ops::Exc0) ) return everts;
  EquiVertices ev;
  // electrons (orbital pairs)
  for ( uint vert = 0; vert < _npairs; ++vert )
    if ( spinsym(vert*2) == Singlet )
      ev.add(vert+offs);
  if ( ev.size() > 1 ) everts.push_back(ev);
  ev.clear();
  // non-conserved electrons (single orbitals)
  for ( uint vert = _npairs; vert < nvertices(); ++vert )
    ev.add(vert+offs);
  if ( ev.size() > 1 ) everts.push_back(ev);
  return everts;
}
Product< Orbital > Matrix::crobs(bool anni) const
{
  Product<Orbital> orbs;
  // first indices - creators
  uint
    mult = 2, 
    offs = anni ? 1 : 0, 
    begin = 0, 
    end = _npairs,
    add = 1;
  if ( _type == Ops::DensM && Input::iPars["prog"]["dmsort"] > 0 ) {
    // different order of electrons: 1,2,...2, 1
    // note that here we exchange creators and annihilators in order to correspond to connections
    // i.e., for something like T^v_u \gamma^u_v we call "u" in gamma annihilator and "v" - creator!
    mult = 1;
    offs = 0;
    if ( !anni ){
      begin = 2*_npairs-1;
      end = _npairs-1;
      add = -1;
    }
  }
  for ( uint vert = begin; vert != end; vert+=add ){
    orbs.push_back(_orbs[mult*vert+offs]);
  }
  // non-conserved electrons
  if ( anni == (_lmel < 0) ){
    for ( uint vert = 2*_npairs; vert < _orbs.size(); ++vert ){
      orbs.push_back(_orbs[vert]);
    }
  }
  return orbs;
}

void Matrix::reset_vertices()
{ _indx=-1; }

bool Matrix::vertices(long int ipos, Matrix& mat, long int ipos1, unsigned int indx)
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
void Matrix::set_connect(TCon2 connected2 )
{ _connected2=connected2; }
void Matrix::add_connect(long int con)
{ _connected2.insert(con); }
TCon2 Matrix::connected2() const
{ return _connected2; }

void Matrix::set_conlines(Product< ConLine > conlines)
{ _conlines = conlines; }
void Matrix::set_conline(lui iorb, lui imat, lui idx)
{
  assert( iorb < _orbs.size() );
  if ( _conlines.size() == 0 )// first resize _conlines
    _conlines.resize(_orbs.size());
  _conlines[iorb] = ConLine(imat,idx);
}
void Matrix::add_conline(lui imat, lui idx)
{ _conlines.push_back(ConLine(imat,idx)); }
const Product< ConLine >& Matrix::conlines() const
{ return _conlines; }
const ConLine& Matrix::conline(lui iorb) const
{
  assert( iorb < _conlines.size() );
  return _conlines[iorb]; 
}
bool Matrix::is0() const
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

void Matrix::set_cran(const Product< SQOpT::Gender >& cran)
{
  assert( _type == Ops::DensM );
  assert( cran.size() == _orbs.size() );
  _cranorder = cran;
}

void Matrix::calc_orbtypeshash()
{
  _orbtypeshash = 0;
  std::vector<uint> ots; 
  _foreach_cauto(Product<Orbital>,ito,_orbs){
    uint itype = uint(ito->type());
    if (ots.size() < itype+1)
      ots.resize(itype+1);
    ots[itype] += 1;
  }
  lui off = 1;
  for ( uint io = 0; io < ots.size(); ++io ){
    _orbtypeshash += ots[io]*off;
    off *= _orbs.size();
  }
}
uint Matrix::diaglevel() const
{
  switch ( _type ){
    case Ops::Exc:
    case Ops::Exc0:
      return 0;
    case Ops::Deexc:
    case Ops::Deexc0:
      return 2;
    case Ops::Fock:
    case Ops::OneEl:
    case Ops::FluctP:
    case Ops::XPert:
    case Ops::DensM:
      return 1;
    default:
      Error("Unknown matrix in diagram");
      return 0;
  }
}

void Matrix::setkind(short int exccl, short int intlines, short int intvirt)
{
  _exccl=exccl;
  _intlines=intlines;
  _intvirt=intvirt;
  calc_orbtypeshash();
}
long int Matrix::iorbel(lui ipos) const
{ 
  assert( ipos < _orbs.size() );
  lui ipos1;
  if ( _type == Ops::DensM && Input::iPars["prog"]["dmsort"] > 0 )
    // different order of electrons: 1,2,...2, 1
    ipos1 = _orbs.size()-ipos-1;
  else
    ipos1 = ipos%2==0?ipos+1:ipos-1;
  
  if ( ipos1 >= 2*_npairs )
    // one of the non-conserved electrons
    return -1;
  return ipos1;
}

Orbital Matrix::orbel(const Orbital& orb) const
{
  long int ipos = _orbs.find(orb);
  if ( ipos >= 0 )
    return orbel(ipos);
  else
    return Orbital();
}
Orbital Matrix::orbel(const long int& ipos) const
{
  long ipos1 = iorbel(ipos);
  if (ipos1 >= 0 )
    return _orbs[ipos1]; 
  else 
    return Orbital();
}

Matrix::Spinsym Matrix::spinsym(long int ipos) const
{
  if (_matspinsym==Triplet && ipos-2 <0) //first electron is triplet
    return Triplet;
  return Singlet;
}
SQOpT::Gender Matrix::genderguess(uint ipos) const
{
  if (_cranorder.size() > 0) return _cranorder[ipos];
  // can't guess for non-conserved electrons, so return a placeholder
  if (ipos >= 2*_npairs) return SQOpT::Gen;
  // first creators, second annihilators
  return (ipos%2 == 0)? SQOpT::Creator : SQOpT::Annihilator;
}

void Matrix::set_no_spin()
{
  Spin nospin(Spin::No);
  for (unsigned int i=0; i<_orbs.size(); ++i)
    _orbs[i].setspin(nospin);
}
void Matrix::permute(const Permut& p)
{
  for (unsigned int i=0; i<_orbs.size(); ++i){
    _orbs[i] = p.permutorb(_orbs[i]);
  }
}
std::string Matrix::plainname() const
{
  const std::string& plainsymbols = Input::sPars["syntax"]["plainsymb"];
  std::string plainnam;
  for ( uint i = 0; i < _name.size(); ++i ){
    if ( plainsymbols.find(_name[i]) != std::string::npos ){
      plainnam += _name[i];
    }
  }
  if ( type() == Ops::FluctP ) {
    plainnam = integralnames();
  }
  return plainnam;
}
std::string Matrix::integralnames() const
{
  assert( _type == Ops::FluctP );
  std::string name = Input::sPars["syntax"]["coulombint"];
  for ( uint iorb = 0; iorb < _orbs.size(); ++iorb ) {
    assert( iorbel(iorb) >= 0 );
    if ( _orbs[iorb].type() != _orbs[iorbel(iorb)].type() ) {
      name = Input::sPars["syntax"]["exchangeint"];
      break;
    }
  }
  return name;
}

std::ostream & operator << (std::ostream & o, Matrix const & mat)
{
  short clean = Input::iPars["output"]["clean"];
  std::string exc0 = Input::sPars["output"]["exc0"];
  std::string tensor("");
  if ( clean <= 0 ) tensor = "\\"+Input::sPars["command"]["tensor"] + " ";
  switch ( mat.type() ){
    case Ops::Fock:
      o << tensor << "f_{" << mat.orbitals() << "}";
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    case Ops::OneEl:
      o << tensor << "h_{" << mat.orbitals() << "}";
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    case Ops::FluctP:
      if ( clean > 1 )
        o << "(" << mat.orbitals().subprod(0,1) << "|" << mat.orbitals().subprod(2,3) << ")";
      else 
        o << tensor << "\\" << Input::sPars["command"]["integral"] << "{" 
                   << mat.orbitals().subprod(0,1) << "}{" << mat.orbitals().subprod(2,3) << "}";
      MyOut::pcurout->lenbuf += 7;
      break;
    case Ops::XPert:
      o << tensor << "X_{" << mat.orbitals() << "}";
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    case Ops::Overlap:
      o << tensor << "S_{" << mat.orbitals() << "}";
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
        IL::add2name(name,oss.str(),true,false);
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
        IL::add2name(name,oss.str(),false,false);
        o << tensor << name;
        MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      }
      break;
    }
    case Ops::Delta:
      o << tensor << "\\delta_{" << mat.orbitals() << "}";
      break;
    case Ops::Exc0:
    case Ops::Deexc0:
      if (exc0 == " ") break;
    case Ops::Exc:
    case Ops::Deexc:
    case Ops::Interm: {
      std::string name(mat.name());
      std::ostringstream oss;
      // occ.indices
      oss << mat.crobs(true);
      IL::add2name(name,oss.str(),true,false);
      oss.str("");
      //virt. indices
      oss << mat.crobs(false);
      IL::add2name(name,oss.str(),false,false);
      o << tensor << name;
      MyOut::pcurout->lenbuf += 1+mat.orbitals().size()/MyOut::pcurout->wsi;
      break;
    }
    case Ops::Number:
      o << tensor << mat.name();
      MyOut::pcurout->lenbuf += 2;
      break;
    case Ops::None:
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

std::ostream & operator << (std::ostream & o, const EquiVertices& ev){
  _foreach_cauto(EquiVertices,it,ev){
    o << "(";
    _foreach_cauto(JointVertices,jv,*it){
      if (jv != it->begin())
        o << " ";
      o << *jv ;
    }
    o << ")";
  }
  return o;
}
std::ostream & operator << (std::ostream & o, const Equivalents& eqv){
  _foreach_cauto(Equivalents,it,eqv)
    o << "[" << *it << "]";
  return o;
}
