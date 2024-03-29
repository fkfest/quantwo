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

std::vector<uint> Matrix::calc_virtelvec(const Product<Orbital>& inorbs) const
{
  std::vector<uint> virtelvec (_npairs);
  for(uint i = 0; i < virtelvec.size(); i++){
    if( inorbs[i*2].type() == Orbital::Virt ) virtelvec[i]++;
    if( inorbs[(i*2)+1].type() == Orbital::Virt ) virtelvec[i]++;
  }
  return virtelvec;
}

std::vector<uint> Matrix::calc_nalpha(const Product<Orbital>& inorbs) const
{
  std::vector<uint> alphavec (inorbs.size()/2,0);
  for(uint i = 0; i < alphavec.size(); i++){
    if( inorbs[i*2].spin() == Spin::Up ) alphavec[i]++;
    if( inorbs[(i*2)+1].spin() == Spin::Up ) alphavec[i]++;
  }
  return alphavec;
}

void Matrix::itforder()
{
  Product<uint> permute(0,2,1,3);
  Product<uint> permute12(2,3,0,1);
  if (this->orbitals().size() == 4){
    uint virtel1, virtel2, occel;
    virtel1 = virtel2 = occel = 0;
    for(uint i = 0; i < 2; i++){
      if (this->orbitals()[i].type() == 2)
        ++virtel1;
      else
        ++occel;
    }
    for(uint i = 2; i < 4; i++){
      if (this->orbitals()[i].type() == 2)
        ++virtel2;
      else
        ++occel;
    }
    if (occel == 4 || virtel1 + virtel2 == 4) //(ij|kl) -> I1324[ikjl] or (ab|cd) -> I1324[acbd]
      this->set_orbs(this->get_orbs().refpro(permute));
    else if(virtel2 > virtel1|| ((virtel1 == virtel2) && this->orbitals()[0].type() == 1 && this->orbitals()[2].type() == 2))
      this->set_orbs(this->get_orbs().refpro(permute12));
  }
}

void Matrix::elemcoorder()
{
  Product<uint> permute(0,2,1,3);
  Product<uint> permute12(2,3,0,1);
  std::vector<uint> alphavec = calc_nalpha(this->orbitals());
  if ( this->type() == Ops::FluctP ){
    std::vector<uint> virtelvec = calc_virtelvec(this->orbitals());
    if ( alphavec[0] < alphavec[1] ){
      this->set_orbs(this->get_orbs().refpro(permute12));
    }
    else if ( alphavec[0] == alphavec[1] ){
      if ( this->get_orbs()[0].spin() == Spin::Down ){
        this->set_orbs(this->get_orbs().refpro(permute12));
      }
      else if ( virtelvec[0] < virtelvec[1] ){
        this->set_orbs(this->get_orbs().refpro(permute12));
      }
      else if ( virtelvec[0] == virtelvec[1] && this->get_orbs()[0].type() == Orbital::Occ ){
        this->set_orbs(this->get_orbs().refpro(permute12));
      }
    }
  }
  if ( this->get_orbs().size() == 4 ) this->set_orbs(this->get_orbs().refpro(permute));
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
    for (const auto& orb: pcrea){
      _orbs.push_back(orb);
      _cranorder.push_back(SQOpT::Creator);
    }
    _foreach_crauto(itorb,panni){
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

Matrix::Matrix(const Kronecker& d, bool ordered)
{
  if ( ordered ){
    _orbs.push_back(d.orb1());
    _orbs.push_back(d.orb2());
  } else {
    _orbs.push_back(d.orb2());
    _orbs.push_back(d.orb1());
  }
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
  if ( name == "U" && npairs == 3 )
    _threeelectronint = true;
  else
    _threeelectronint = false;
}

void Matrix::gen_name(const std::string& name)
{
  std::string exc0 = Input::sPars["output"]["exc0"];
  // replace default name
  if ( name == "T" ) {
    switch (_type) {
      case Ops::Fock:
        _name = "f";
        break;
      case Ops::OneEl:
        _name = "h";
        break;
      case Ops::FluctP:
        _name = "W";
        break;
      case Ops::XPert:
        _name = "X";
        break;
      case Ops::Overlap:
        _name = "S";
        break;
      case Ops::DensM:
        _name = "\\gamma";
        break;
      case Ops::Delta:
        _name = "\\delta";
        // fall through
      case Ops::Exc0:
        // fall through
      case Ops::Deexc0:
        if (exc0 != " "){
          _name = exc0;
          if ( _type == Ops::Deexc0 ) // add dagger
            IL::add2name(_name,Input::aPars["syntax"]["dg"].front());
          break;
        }
        // fall through
      default:
        _name=name;
    }
  } else {
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
  for (const auto& orb: _orbs){
    if ( sumorbs.count(orb) == 0 ) {
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
  if(orb2.spin().type() == Spin::Up || orb2.spin().type() == Spin::Down){
    for ( unsigned int i=0; i<_intorbs.size(); ++i ){
      _intorbs[i].replace(orb1,orb2,smart);
    }
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
  // is it an ordered density matrix?
  bool dmo = (_type == Ops::DensM && Input::iPars["prog"]["dmsort"] > 0);
  EquiVertices ev;
  OrbitalTypes curorbts;
  // electrons (orbital pairs)
  uint nextvert = 0;
  if ( nextvert < _npairs ) do {
    uint startvert = nextvert;
    nextvert = 0;
    curorbts = orbtypes4vertex(startvert,dmo);
    for ( uint vert = startvert; vert < _npairs; ++vert ) {
      if ( curorbts == orbtypes4vertex(vert,dmo) ){
        if ( dmo || spinsym(vert*2) == Singlet )
          ev.add(vert+offs);
      } else if ( nextvert == 0 ){
        nextvert = vert;
      }
    }
    if ( ev.size() > 1 ) everts.push_back(ev);
    ev.clear();
  } while (nextvert > 0);
  // non-conserved electrons (single orbitals)
  nextvert = _npairs;
  if ( nextvert < nvertices()) do {
    uint startvert = nextvert;
    nextvert = 0;
    curorbts = orbtypes4vertex(startvert,dmo);
    for ( uint vert = startvert; vert < nvertices(); ++vert ) {
      if ( curorbts == orbtypes4vertex(vert,dmo) ){
        ev.add(vert+offs);
      } else if ( nextvert == 0 ){
        nextvert = vert;
      }
    }
    if ( ev.size() > 1 ) everts.push_back(ev);
    ev.clear();
  } while (nextvert > 0);

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
  if ( _type == Ops::DensM ) {
    if ( Input::iPars["prog"]["dmsort"] > 0 ) {
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
    } else {
      // use order stored in _cranorder
      for ( uint iorb = 0; iorb < 2*_npairs; ++iorb ){
        if ( (_cranorder[iorb] == SQOpT::Creator) == anni )
          orbs.push_back(_orbs[iorb]);
      }
    }
//   } else if ( _type == Ops::Delta ) {
//     // exchange creators and annihilators in order to correspond to connections (cf. density matrices)
//     offs = anni ? 0 : 1;
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
OrbitalTypes Matrix::orbtypes4vertex(uint vertex, bool dmo) const
{
  OrbitalTypes orbts;
  if ( vertex < _npairs ){
    if ( dmo ) {
      // orbitals4vertices: 1 2 3 3 2 1
      orbts.push_back(_orbs[vertex].type());
      orbts.push_back(_orbs[_orbs.size()-1-vertex].type());
    } else {
      // orbitals4vertices: 1 1 2 2 3 3
      orbts.push_back(_orbs[2*vertex].type());
      orbts.push_back(_orbs[2*vertex+1].type());
    }
  } else {
    // non-conserved electrons
    assert( vertex < _npairs + std::abs(_lmel) );
    orbts.push_back(_orbs[_npairs+vertex].type());
  }
  return orbts;
}

void Matrix::reset_vertices()
{ _indx=-1; }

bool Matrix::samespin() const
{
  for( Product<Orbital>::const_iterator it = _orbs.begin(); it != std::prev(_orbs.end()); it++ ){
    if( it->spin() != std::next(it)->spin() ) return false;
  };
  return true;
}

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
    for (const auto& orb: _orbs){
      if (orb.type() != Orbital::Act) return true;
    }
    // number of creators == annihilators
    int cran = 0;
    for (const auto& ca:_cranorder){
      if (ca == SQOpT::Creator)
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

void Matrix::set_orbs( Product<Orbital>& crobs, Product<Orbital>& anobs )
{
  for( uint i = 0; i != _npairs; i++ ){
    _orbs[2*i] = crobs[i];
    _orbs[2*i+1] = anobs[i];
  }
}

void Matrix::calc_orbtypeshash()
{
  _orbtypeshash = 0;
  std::vector<uint> ots;
  for (const auto& orb: _orbs){
    uint itype = uint(orb.type());
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
long int Matrix::iorbel(lui ipos, bool physicist) const
{
  assert( ipos < _orbs.size() );
  lui ipos1;
  if (physicist){
    assert(_orbs.size() % 2 == 0);
    uint step = _orbs.size()/2;
    if( ipos + 1 > step )
      ipos1 = ipos - 2;
    else
      ipos1 = ipos + 2;
    return ipos1;
  }
  else{
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
void Matrix::set_no_el()
{
  for (unsigned int i=0; i<_orbs.size(); ++i)
    _orbs[i].setel(0);
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
  if ( type() == Ops::FluctP || type() == Ops::Fock || _threeelectronint ) {
    if ( (type() == Ops::Fock && Input::iPars["prog"]["algo"] == 2) || (Input::iPars["prog"]["algo"] == 2 && _threeelectronint) ) return plainnam;
    else return integralnames();
  }
  return plainnam;
}

std::string Matrix::integralnames() const
{
  if ( Input::iPars["prog"]["algo"] == 1 ){//ITF code
    return itfintegralnames();
  }
  else if ( Input::iPars["prog"]["algo"] == 2 ){//ElemCo code
    return elemcointegralnames();
  }
  else{
    error("Unknown algorithm in integralnames");
    return "1";
  }
}

std::string Matrix::itfintegralnames() const
{
  if (_type == Ops::Fock )
  {
    return itf1eint();
  }
  else if (_orbs.size() == 4) {
    return itf2eint();
  }
  else if (_orbs.size() == 6) {
    return itf3eint();
  }
  else
  {
    error("Unimplemented ITF name required.");
    __builtin_unreachable();
  }
}

std::string Matrix::itf1eint() const{
  assert(_orbs.size() == 2);
  assert(_type == Ops::Fock);
  if( this->_orbs[0].type() == Orbital::Occ && this->_orbs[1].type() == Orbital::Virt ) return "f21";
  else return "f";
}

std::string Matrix::itf2eint() const{
  assert(_orbs.size() == 4);
  assert(_type == Ops::FluctP);
  Product<Orbital> orbs = this->orbitals();
  uint virtel1, virtel2;
  virtel1 = virtel2 = 0;
  if (orbs[0].type() == Orbital::Virt) ++virtel1;
  if (orbs[1].type() == Orbital::Virt) ++virtel1;
  if (orbs[2].type() == Orbital::Virt) ++virtel2;
  if (orbs[3].type() == Orbital::Virt) ++virtel2;
  if(virtel2 > virtel1 || ((virtel1 == virtel2) && orbs[0].type() == Orbital::Occ && orbs[2].type() == Orbital::Virt)){
    Orbital orb1 = orbs[0];
    Orbital orb2 = orbs[1];
    orbs.eraseelem(1);
    orbs.eraseelem(0);
    orbs *= orb1;
    orbs *= orb2;
  }
  std::string name = Input::sPars["syntax"]["coulombint"]; //(ab|cd) (ij|kl) (ab|kl)
  if ( orbs[0].type() == Orbital::Occ ){ //(ij|kl)
//     if ((orbs[0].spin() == Spin::Up || orbs[0].spin() == Spin::Down) && orbs[0].spin() != orbs[2].spin()) name = "dI1234";
//     else name = "dI1324";}
    name = "dI1324";}
  else if (orbs[2].type() == Orbital::Occ) // (ab|kl)
    name = "dI1234";
  else{ //(ab|cd)
//     if((orbs[0].spin() == Spin::Up || orbs[0].spin() == Spin::Down) && orbs[0].spin() != orbs[2].spin()) name = "dI1234";
//     else name = "dI1324";
    name = "dI1324";
  }
  for ( uint iorb = 0; iorb < orbs.size(); ++iorb ) {
    assert( iorbel(iorb) >= 0 );
    if ( orbs[iorb].type() != orbs[iorbel(iorb)].type() ) {
      name = Input::sPars["syntax"]["exchangeint"];
      if (orbs[0].type() == Orbital::Virt ){ //(aj|kd) (aj|cl) (aj|kl) (ab|cl) (ab|kd)
        if (orbs[1].type() == Orbital::Occ ){ //(aj|kd) (aj|cl) (aj|kl)
          if (orbs[2].type() == Orbital::Virt){
            name = "dI1324";} //(aj|cl);
          else if (orbs[3].type() == Orbital::Virt){ //(aj|kd)
            if( orbs[0].spin() == Spin::Down && orbs[2].spin() == Spin::Up ) name = "dI3124";
            else name = "dI1342";
          }
        else{ //(aj|kl)
          if(orbs[0].spin() == Spin::Down && orbs[2].spin() == Spin::Up)
            name = "dI1423";
          else
            name = "dI1234";
        }
      }
      else{//(ab|cl) (ab|kd)
            if (orbs[2].type() == Orbital::Virt) //(ab|cl)
              if(orbs[0].spin() == Spin::Down && orbs[2].spin() == Spin::Up) // (AB|cl)
                name = "dI2314";
              else
                name = "dI1234";
            else if(orbs[0].spin() == Spin::Down && orbs[2].spin() == Spin::Up)
              name = "dI2341";
            else name ="dI1243";
          }
      }
      else{//(ib|kd) (ib|kl)
        if(orbs[3].type() == Orbital::Virt) //(ib|kd)
          name = "I3142";
        else if(orbs[0].spin() == Spin::Down && orbs[2].spin() == Spin::Up) //(IB|kl) which will be canonicalized for ITF to [BklI]
          name = "dI4123";
        else name ="dI2134";
      }
//       break;
    }
  }
  return name;
}

std::string Matrix::itf3eint() const{
  assert( _orbs.size() == 6 ); //three electron integrals
  std::string name;
  std::vector<uint> virtelvec = calc_virtelvec(_orbs);
  Product<uint> ref, refref;
  Product<Orbital> orbs = _orbs;
  ref.identity(_npairs);
  // for(uint i = 0; i < _npairs; i++){
  //   if( _orbs[i*2].type() == Orbital::Virt ) virtelvec[i]++;
  //   if( _orbs[(i*2)+1].type() == Orbital::Virt ) virtelvec[i]++;
  // }
  InsertionSortD(&virtelvec[0],&ref[0],virtelvec.size());
  for(Product<uint>::iterator it = ref.begin(); it != ref.end(); it++){
    refref.push_back(*it*2);
    refref.push_back(((*it*2)+1));
  }
  orbs = orbs.refpro(refref);
  std::fill(virtelvec.begin(), virtelvec.end(), 0);
  for(uint i = 0; i < _npairs; i++){
    if( orbs[i*2].type() == Orbital::Virt ) virtelvec[i]++;
    if( orbs[(i*2)+1].type() == Orbital::Virt ) virtelvec[i]++;
  }
  
  if( virtelvec[0] == virtelvec[1] && virtelvec[2] == 0 ){
    if( virtelvec[0] == 1 ){
      if( orbs[0].type() == Orbital::Virt ){
        if( orbs[2].type() == Orbital::Virt ) name = "dI132456";
        else name = "dI134256";
      }
      else name = "dI314256";
    }
    else if( virtelvec[0] == 2 ){
      name = "dI123456";
    }
  }
  else if( virtelvec[0] == 2 && virtelvec[1] == virtelvec[2] ){
    if( virtelvec[1] == 1 ){
      if(orbs[2].type() == Orbital::Occ) name = "dI125364";
      else name = "dI123564";
    }
    if( virtelvec[1] == 0 )
      name = "dI123456";
  }
  else{
    name = "dummy";}
  return name;
}

std::string Matrix::elemcointegralnames() const
{
  assert(_orbs.size() == 4);
  return elemco2eint();
}

std::string Matrix::elemco2eint() const
{
  assert( _type == Ops::FluctP );
  std::string name="dummy";
  bool exchange = false;
  for ( uint iorb = 0; iorb < this->orbitals().size(); ++iorb ) {
    assert( iorbel(iorb,true) >= 0 );
    if ( this->orbitals()[iorb].type() != this->orbitals()[iorbel(iorb,true)].type() )
      exchange = true;
  }
  const Product<Orbital>& orbs = this->orbitals();
  if ( !exchange ){ // <ik|jl> <ak|bl> <ac|bd> <ib|kd>
    if ( orbs[0].type() == Orbital::Occ ){ //<ik|jl> <ib|kd>
      if ( orbs[1].type() == Orbital::Occ ){
        name = "d_oooo";
        if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oOoO";
        else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OOOO";
      }
      else{ //<
        name = "d_ovov";
        if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oVoV";
        else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OVOV";
      }
    }
    else if (orbs[1].type() == Orbital::Occ){ // <ak|bl>
      name = "d_vovo";
      if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vOvO";
      else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VOVO";
    }
    else{ //<ac|bd>
      name = "d_vvvv";
      if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vVvV";
      else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VVVV";
    }
  }
  else{ //exchange integral
    if (orbs[0].type() == Orbital::Virt ){ // <vo|ov> <vv|oo> <vo|oo> <vv|vo> <vo|vv> <vv|ov>
      if (orbs[2].type() == Orbital::Occ ){ // <vv|oo> <vo|ov> <vo|oo> <vv|ov>
        if (orbs[1].type() == Orbital::Virt){ // <vv|oo> <vv|ov>
          if (orbs[3].type() == Orbital::Virt){
            name = "d_vvov";
            if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vVoV";
            else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VVOV";
          }
          else{
            name = "d_vvoo";
            if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vVoO";
            else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VVOO";
          }
        }
        else if (orbs[3].type() == Orbital::Virt){ //<vo|ov>
          name = "d_voov";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vOoV";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VOOV";
        }
        else{ //<vo|oo>
          name = "d_vooo";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vOoO";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VOOO";
        }
      }
      else{//(vv|vo) (vv|ov)
        if (orbs[1].type() == Orbital::Virt){ //<ac|bl>
          name = "d_vvvo";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vVvO";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VVVO";
        }
        else{
          name = "d_vovv";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_vOvV";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_VOVV";
        }
      }
    }
    else{//<oo|vv> <ov|vo> <ov|vv> <oo|vo> <oo|ov> <ov|oo>
      if(orbs[3].type() == Orbital::Virt){ // <oo|vv> <ov|vv> <oo|ov>
        if(orbs[1].type() == Orbital::Occ){// <oo|vv> <oo|ov> 
          if(orbs[2].type() == Orbital::Occ){ 
            name = "d_ooov";
            if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oOoV";
            else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OOOV";
          }
          else{
          name = "oovv";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "oOvV";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "OOVV";
          }
        }
        else{// <ov|vv>
          name = "d_ovvv";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oVvV";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OVVV";
        }
      }
      else{ //<ov|vo> <oo|vo> <ov|oo>
        if(orbs[1].type() == Orbital::Virt){//<ov|vo> <ov|oo>
          if(orbs[2].type() == Orbital::Virt){
            name = "d_ovvo";
            if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oVvO";
            else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OVVO";
          }
          else{
            name = "d_ovoo";
            if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oVoO";
            else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OVOO";
          }
        }
        else{
          name = "d_oovo";
          if( orbs[0].spin() == Spin::Up && orbs[1].spin() == Spin::Down ) name = "d_oOvO";
          else if( orbs[0].spin() == Spin::Down && orbs[1].spin() == Spin::Down ) name = "d_OOVO";
        }
      }
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
    case Ops::OneEl:
    case Ops::XPert:
    case Ops::Overlap:
      o << tensor << mat.name() << "_{" << mat.orbitals() << "}";
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
      // fall through
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

Return::Vals Permut::replace(const Orbital& orb1, const Orbital& orb2){
  for ( TPerMap::const_iterator it = _orbs.begin(); it != _orbs.end(); ++it ){
    if( it->second == orb1 ){
      Orbital orbfrom = it->first;
      _orbs.erase(it);
      _orbs[orbfrom] = orb2;
    }
  }
  TPerMap::const_iterator it = _orbs.find(orb1);
  if ( it == _orbs.end()){return Return::Delete;}
  Orbital orbto;
  orbto = it->second;
  _orbs.erase(it);
  _orbs[orb2] = orbto;
  return Return::Done;
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

Array<Product<Orbital>> Matrix::elecorbs(){
  Array<Product<Orbital>> outorbs;
  Product<Orbital> elecorbs;
  for( size_t i=0; i < _orbs.size(); i+=2 ){
    elecorbs.clear();
    for( size_t j=i; j<i+2; j++ ){
      elecorbs.push_back(_orbs[j]);
    }
    outorbs.push_back(elecorbs);
  }
  return outorbs;
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
  for (const auto& ver: ev){
    o << "(";
    for (const auto& jv: ver){
      o << " " << jv;
    }
    o << ")";
  }
  return o;
}
std::ostream & operator << (std::ostream & o, const Equivalents& eqv){
  for (const auto& e: eqv)
    o << "[" << e << "]";
  return o;
}
