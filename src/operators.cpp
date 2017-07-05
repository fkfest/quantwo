#include "operators.h"

/* SQOp::SQOp(std::string orb)
{
  _gender=(isupper((char)orb[0]) ? SQOpT::Creator : SQOpT::Annihilator);
  orb[0] = std::tolower(orb[0]);
  _orb = orb;
} */
SQOp::SQOp(SQOpT::Gender gender, Orbital orb)
{
  _gender=gender;
  _orb = orb;
  if (_orb.type()==Orbital::Occ && _gender==SQOpT::Creator) {
    _genderPH = SQOpT::Annihilator;
  } else if (_orb.type()==Orbital::Occ && _gender==SQOpT::Annihilator) {
    _genderPH = SQOpT::Creator;
  } else if (_orb.type()==Orbital::GenT || _orb.type()==Orbital::Act) {
    _genderPH = SQOpT::Gen;
  } else {
    _genderPH = _gender;
  }
}
Orbital SQOp::orb() const
{ return _orb; }
Return SQOp::replace(Orbital orb1, Orbital orb2, bool smart)
{
  return _orb.replace(orb1,orb2,smart);
}
Return SQOp::replace(Spin spin1, Spin spin2, bool smart)
{
  return _orb.replace(spin1,spin2,smart);
}

std::ostream & operator << (std::ostream & o, SQOp const & op)
{
  o << "\\op{" << op.orb() << "}";
  if ( op.gender()==SQOpT::Creator )
    o << "^\\dg";
  
  return o;
}

Oper::Oper()
{  
  _prefac=1; 
  _type=Ops::None;
  p_Term = 0;
}

Oper::Oper(Ops::Type type, bool antisym, Term* pTerm, const std::vector<OrbitalTypes>& orbtypes)
{
  assert( InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  std::string name;
  _type=type;
  p_Term = pTerm;
  if (type == Ops::Fock )
    name="f";
  if (type == Ops::OneEl )
    name="h";
  else if (type == Ops::FluctP )
    name="W";
  else
    name="T";
  create_Oper(name,antisym,orbtypes);
}
Oper::Oper(Ops::Type type, short int exccl, std::string name, int lm, int pmsym, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = pTerm;
  Orbital orb0, orb1;
  if (pTerm) {
    orb0 = pTerm->freeorbname(Orbital::Virt);
    orb1 = pTerm->freeorbname(Orbital::Occ);
  } else {
    orb0 = Orbital(std::string("A"));
    orb1 = Orbital(std::string("I"));
  }
  create_Oper(exccl,orb1,orb0,name,lm,pmsym);
}
Oper::Oper(Ops::Type type, short int exccl, Orbital occ, Orbital virt, std::string name, 
           int lm, int pmsym, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = pTerm;
  create_Oper(exccl,occ,virt,name,lm,pmsym);
//   xout << *this << std::endl;
}
Oper::Oper(Ops::Type type, short int exccl, const std::map< Orbital::Type, Orbital >& orbnames, 
           const std::vector<OrbitalTypes>& orbtypes, std::string name, int lm, int pmsym, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  assert( orbtypes.size() == 2 );
  assert( int(orbtypes[0].size()) == exccl );
  assert( int(orbtypes[1].size()) == exccl+lm );
  _type=type;
  p_Term = pTerm;
  create_Oper(exccl,orbnames,orbtypes,name,lm,pmsym);
}

Oper::Oper(Ops::Type type, short int exccl, const Product< Orbital >& occs, const Product< Orbital >& virts, 
           std::string name, int lm, int pmsym, Term* pTerm)
{
  assert( occs.size() + lm == virts.size() );
  assert( occs.size() == uint(exccl) );
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = pTerm;
  create_Oper(occs,virts,name,pmsym);
}
Oper::Oper(Ops::Type type, const Product< Orbital >& orbs, std::string name, int lm, int pmsym)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = 0;
  create_Oper(orbs,name,lm,pmsym);
}

Oper::Oper(Ops::Type type, short int exccl, const std::vector<OrbitalTypes>& orbtypes, 
           std::string name, int lm, int pmsym, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  assert( orbtypes.size() == 2 );
  assert( int(orbtypes[0].size()) == exccl );
  assert( int(orbtypes[1].size()) == exccl+lm );
  
  _type=type;
  p_Term = pTerm;
  std::map<Orbital::Type,Orbital> orbnames;
  if ( p_Term ){
    for ( uint i=0; i<orbtypes.size(); ++i ){
      _foreach_cauto(OrbitalTypes,iot,orbtypes[i]){
        if ( orbnames.count(*iot) == 0 ){
          // new orb type
          orbnames[*iot] = p_Term->freeorbname(*iot);
        }
      }
    }
  } else {
    TsPar& orbs = Input::sPars["syntax"];
    orbnames[Orbital::Occ] = Orbital(std::string(1,char(std::toupper(orbs["occorb"][0]))));
    orbnames[Orbital::Virt] = Orbital(std::string(1,char(std::toupper(orbs["virorb"][0]))));
    orbnames[Orbital::Act] = Orbital(std::string(1,char(std::toupper(orbs["actorb"][0]))));
    orbnames[Orbital::GenT] = Orbital(std::string(1,char(std::toupper(orbs["genorb"][0]))));
  }
  create_Oper(exccl,orbnames,orbtypes,name,lm,pmsym);
}

void Oper::create_Oper(const std::string& name, bool antisym, const std::vector<OrbitalTypes>& orbtypes)
{
  assert( InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  Product<Orbital> porbs;
  Matrix::Spinsym spinsym = Matrix::Singlet;
  Electron el = 1;
  if (p_Term) el = p_Term->nextelectron();
  // operators with general indices
  std::string orbname("P");
  if (orbtypes.size() > 0) {
    if (orbtypes.size() != 2 || orbtypes[0].size() == 0 || orbtypes[1].size() == 0) 
      error("Wrong type definition of Hamiltonian part","Oper::create_Oper");
    // specified type
    if (p_Term)
      orbname = p_Term->freeorbname(orbtypes[1][0],true).name();
    else
      error("Explicit types in Hamiltonian outside of term not implemented","Oper::create_Oper");
    orbname[0] = std::toupper(orbname[0]);
  }
  Orbital orb(orbname,el);
  _SQprod*=SQOp(SQOpT::Creator,orb);
  porbs*=orb;
  _orbs.insert(orb);
  _sumorbs.insert(orb);
  orbname = "Q";
  if (orbtypes.size() > 0) {
    // specified type
    if (p_Term)
      orbname = p_Term->freeorbname(orbtypes[0][0]).name();
    else
      error("Explicit types in Hamiltonian outside of term not implemented","Oper::create_Oper");
    orbname[0] = std::toupper(orbname[0]);
  }
  orb=Orbital(orbname,el);//same electron as in P
  porbs*=orb;
  _orbs.insert(orb);
  _sumorbs.insert(orb);
  _prefac=1;
  if ( InSet(_type, Ops::Fock,Ops::OneEl,Ops::XPert) ) {
    _SQprod*=SQOp(SQOpT::Annihilator,orb);
  } else {
   // we use chemical notation (PQ|RS) P^\dg R^\dg S Q
    ++el;
    if (p_Term) el = p_Term->nextelectron();
    orbname = "R";
    if (orbtypes.size() > 0) {
      if ( orbtypes[0].size() != 2 || orbtypes[1].size() != 2  ) 
        error("Wrong type definition of Hamiltonian part","Oper::create_Oper");
      // specified type
      if (p_Term) orbname = p_Term->freeorbname(orbtypes[1][1]).name();
      orbname[0] = std::toupper(orbname[0]);
    }
    orb=Orbital(orbname,el);
    _SQprod*=SQOp(SQOpT::Creator,orb);
    porbs*=orb;
    _orbs.insert(orb);
    _sumorbs.insert(orb);
    orbname = "S";
    if (orbtypes.size() > 0) {
      // specified type
      if (p_Term) orbname = p_Term->freeorbname(orbtypes[0][1]).name();
      orbname[0] = std::toupper(orbname[0]);
    }
    orb=Orbital(orbname,el);//same electron as in R
    _SQprod*=SQOp(SQOpT::Annihilator,orb);
    porbs*=orb;
    _orbs.insert(orb);
    _sumorbs.insert(orb);
    orb=porbs[1];
    _SQprod*=SQOp(SQOpT::Annihilator,orb);
    _prefac /= 4;
  }
  short npairs = porbs.size()/2;
  _mat=Matrix(_type,porbs,npairs,0,0,name,spinsym,antisym);
  // needed for expanding general normal ordered operators later
  if ( Input::iPars["prog"]["contrexcop"] > 1 ) move_SQprod();
}
void Oper::create_Oper(short int const & exccl,Orbital const & occ, Orbital const & virt, 
                       std::string const & name, int lm, int pmsym)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  std::string excl;
  Product<Orbital> occs, virts; 
  short 
    nocc = exccl,
    nvirt = exccl+lm,
    nmax = std::max(nocc,nvirt);
  for (short i=0; i<nmax ; ++i) {  
    if (i>0) excl=num2str(i,std::dec);
    if ( i < nocc )
      occs *= Orbital(occ.name()+excl,occ.spin());
    if ( i < nvirt )
      virts *= Orbital(virt.name()+excl,virt.spin());
  }
  create_Oper(occs,virts,name,pmsym);
}
void Oper::create_Oper(const short int& exccl, const std::map< Orbital::Type, Orbital >& orbnames, 
                       const std::vector<OrbitalTypes>& orbtypes, const std::string& name, int lm, int pmsym)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  std::string excl;
  Product<Orbital> occs, virts; 
  short 
    nocc = exccl,
    nvirt = exccl+lm;
  const OrbitalTypes & occtype(orbtypes[0]);
  const OrbitalTypes & virtype(orbtypes[1]);
  std::map< Orbital::Type, uint > num4type;
  assert(int(occtype.size()) == nocc);
  assert(int(virtype.size()) == nvirt);
  for (short i = 0; i < nocc; ++i) {  
    uint i4t = num4type[occtype[i]];
    ++num4type[occtype[i]];
    if (i4t > 0) {
      excl = num2str(i4t,std::dec);
    } else {
      excl = "";
    }
    const Orbital& occ(orbnames.at(occtype[i]));
    occs *= Orbital(occ.name()+excl,occ.spin());
  }
  for (short i = 0; i < nvirt; ++i) {  
    uint i4t = num4type[virtype[i]];
    ++num4type[virtype[i]];
    if (i4t > 0) {
      excl=num2str(i4t,std::dec);
    } else {
      excl = "";
    }
    const Orbital& virt(orbnames.at(virtype[i]));
    virts *= Orbital(virt.name()+excl,virt.spin());
  }
  create_Oper(occs,virts,name,pmsym);
}

void Oper::create_Oper(const Product< Orbital >& occs, const Product< Orbital >& virts, 
                       const std::string& name, int pmsym)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  Matrix::Spinsym spinsym = Matrix::Singlet;
  Product<Orbital> porbs;
  // excitation and deexcitation operators
  const Product<Orbital> 
    * p_orb0 = &virts,
    * p_orb1 = &occs;
  if (InSet(_type, Ops::Deexc,Ops::Deexc0)) std::swap(p_orb0,p_orb1);
  short 
    ncrea = p_orb0->size(),
    nanni = p_orb1->size(),
    nmax = std::max(ncrea,nanni);
  Electron el;
  for (unsigned short i = 0; i < nmax; ++i) {  
    Spin spin(Spin::No);
    bool setspindiff = (i == nmax-1 && spinsym != Matrix::Singlet);
    if ( p_Term ) 
      el = p_Term->nextelectron();
    else
      el = i+1;
    if ( i < ncrea ) {
      Orbital orb = (*p_orb0)[i];
      spin = orb.spin();
      if (setspindiff){ // triplet part -- use spin-difference
        assert(spin.type() != Spin::No);
        spin.settype(Spin::GenD);
      }
      spin.setel(el);
      orb.setspin(spin);
      porbs *= orb;
    }
    if ( i < nanni ) {
      Orbital orb = (*p_orb1)[i];
      spin = orb.spin();
      if (setspindiff){ // triplet part -- use spin-difference
        assert(spin.type() != Spin::No);
        spin.settype(Spin::GenD);
      }
      spin.setel(el);
      orb.setspin(spin);
      porbs *= orb;
    }
  }
  create_Oper(porbs,name,p_orb0->size()-p_orb1->size(),pmsym);
}

void Oper::create_Oper(const Product< Orbital >& orbs, const std::string& name, int lm, int pmsym)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  bool noprefac = (Input::iPars["prog"]["nobrafac"]) && InSet(_type, Ops::Exc0,Ops::Deexc0);
  int contrexcop = Input::iPars["prog"]["contrexcop"];
  Matrix::Spinsym spinsym = Matrix::Singlet;
  Product<SQOp> anniSQprod;
  // excitation and deexcitation operators
  assert( (orbs.size()-lm)%2 == 0 && (orbs.size()+lm)%2 == 0 );
  short 
    ncrea = (orbs.size()+lm)/2,
    nanni = (orbs.size()-lm)/2,
    nmax = std::max(ncrea,nanni),
    npairs = std::min(ncrea,nanni);
  // symmetry (hash of orbital-type -> number of such electrons) (for prefactor calculation)
  std::map<uint,uint> sym;
  std::map<uint,uint> symel;
  uint hash;
  _prefac = 1;
  Product<Orbital>::const_iterator itorb = orbs.begin();
  for (unsigned short i = 0; i < nmax; ++i) {  
    Electron el = 0;
    uint elhash = 0;
    if ( i < ncrea ) {
      _SQprod*=SQOp(SQOpT::Creator, *itorb);
      _orbs.insert(*itorb);
      if (!InSet(_type, Ops::Exc0,Ops::Deexc0)) 
        _sumorbs.insert(*itorb);
      el = itorb->spin().el();
      // add this type of electron-orbital
      hash = uint(itorb->type()) + 2*Orbital::MaxType*itorb->spin().spinhash();
      sym[hash] += 1;
      // hash for (any) electron
      elhash += uint(itorb->type()) + Orbital::MaxType*Orbital::MaxType*itorb->spin().spinhash(false);
      ++itorb;
    }
    if ( i < nanni ) {
      if (contrexcop > 0)
        anniSQprod *= SQOp(SQOpT::Annihilator, *itorb);
      else
        _SQprod *= SQOp(SQOpT::Annihilator, *itorb);
      _orbs.insert(*itorb);
      if (!InSet(_type, Ops::Exc0,Ops::Deexc0)) 
        _sumorbs.insert(*itorb);
      if ( el != 0 && el != itorb->spin().el() )
        error("Different electrons in the product!","Oper::create_Oper");
      // add this type of electron-orbital
      hash = Orbital::MaxType + uint(itorb->type()) + 2*Orbital::MaxType*itorb->spin().spinhash();
      sym[hash] += 1;
      // hash for (any) electron
      elhash += Orbital::MaxType*uint(itorb->type()) + 
                Orbital::MaxType*Orbital::MaxType*Spin::MaxType*itorb->spin().spinhash(false);
      ++itorb;
    }
    symel[elhash] += 1;
  }
  if (contrexcop > 0) {
    // add annihilation operators in the opposite order 
    _foreach_crauto(Product<SQOp>,ita,anniSQprod){
      _SQprod *= *ita;
    }
  }
  if (!noprefac) {
    // prefactor
    std::map<uint,uint>::const_iterator is; 
    _foreach(is,sym){
      for (uint i = 0; i < is->second; ++i)
        _prefac *= i+1;
    }
    if (spinintegr){
      // take into account the indistinguishability of the electrons
      _foreach(is,symel){
        for (uint i = 0; i < is->second; ++i)
          _prefac *= i+1;
      }
    }
    _prefac = 1/_prefac;
  }
  _mat=Matrix(_type,orbs,npairs,lm,pmsym,name,spinsym);
  if ( contrexcop > 1 ) {
    // general normal order
    if ( ncrea != nanni ) 
      error("General normal ordering (prog,contrexcop>1) is not implemented for non-conserving case","Oper::create_Oper"); 
    _sumops += gennormord();
  } 
}

TermSum Oper::gennormord()
{
  assert( _SQprod.size()%2 == 0 );
  TermSum ret;
  uint ncrea = _SQprod.size()/2;
  ret += Term(_SQprod);
  Oper op;
  switch (ncrea) {
    case (0) :
      break;
    case (1) :
      op = *this;
      ret -= op.gennormordterm(Inds(0),Inds(0));
      break;
    case (2) :
      op = *this;
      ret -= op.gennormordterm(Inds(0),Inds(0));
      op = *this;
      ret -= op.gennormordterm(Inds(1),Inds(1));
      op = *this;
      ret -= op.gennormordterm(Inds(0),Inds(1));
      op = *this;
      ret -= op.gennormordterm(Inds(1),Inds(0));
      op = *this;
      ret -= op.gennormordterm(Inds(0,1),Inds(0,1));
      break;
    default :
      error("General normal ordering (prog,contrexcop>1) is not implemented for excitation class "+ncrea,"Oper::create_Oper");
  }
  _SQprod.clear();
  return ret;
}

TermSum Oper::gennormordterm(const Product< uint >& creas, const Product< uint >& annis)
{
  assert( creas.size() == annis.size() );
  assert( _SQprod.size() >= creas.size() );
  Product<Orbital> pcrea, panni;
  std::vector<uint> mask(_SQprod.size(),1);
  TFactor fac(1);
  for ( uint iop = 0; iop < creas.size(); ++iop ){
    uint 
      icr = creas[iop],
      ian0 = annis[iop];
    assert( icr < _SQprod.size() && ian0 < _SQprod.size() );
    uint 
      ian = _SQprod.size() - 1 - ian0,
      ianorig = _SQprod.size() - 1 - icr;
    // should be ordered
    assert( _SQprod[icr].orb().getel() == _SQprod[ianorig].orb().getel() );
    Spin
      spinold = _SQprod[ian].orb().spin(),
      spinnew = _SQprod[ianorig].orb().spin();
    if ( spinold != spinnew ) {
      _SQprod[ianorig].replace(spinnew,spinold,false);
      _SQprod[ian].replace(spinold,spinnew,false);
      fac /= -2;
    }
    assert( _SQprod[icr].gender() == SQOpT::Creator );
    pcrea.push_back(_SQprod[icr].orb());
    mask[icr] = 0;
    assert( _SQprod[ian].gender() == SQOpT::Annihilator );
    panni.push_back(_SQprod[ian].orb());
    mask[ian] = 0;
  }
  // copy remaining SQ-operators
  Product<SQOp> sqops;
  for ( uint iop = 0; iop < _SQprod.size(); ++iop ){
    if ( mask[iop] > 0 )
      sqops.push_back(_SQprod[iop]);
  }
  _SQprod = sqops;
  
  Matrix gamma(Ops::DensM,pcrea,panni,creas.size());
  
  TermSum ret;
  if ( !gamma.is0() ){
    ret += gennormord();
    Term trmgamma;
    trmgamma *= gamma;
    ret *= trmgamma;
    ret *= fac;
  }
  return ret; 
}

void Oper::move_SQprod()
{
  if ( _SQprod.size() > 0 ) {
    assert ( _sumops.size() == 0 ); // otherwise one should multiply _SQprod with _sumops...
    _sumops += Term(_SQprod);
    _SQprod.clear();
  }
}

Matrix Oper::mat() const
{ return _mat;}
Product< SQOp > Oper::SQprod() const
{ return _SQprod;}
TFactor Oper::prefac() const
{ return _prefac;}
const TOrbSet & Oper::orbs() const
{ return _orbs;}
TOrbSet Oper::sumorbs() const
{ return _sumorbs;}



std::ostream & operator << (std::ostream & o, Oper const & op)
{
  if ( _todouble(_abs(_abs(op.prefac()) - 1)) > MyOut::pcurout->small) o << op.prefac();
  if ( op.sumorbs().size() > 0 ) o << "\\sum_{" << op.sumorbs() << "}";
  o << op.mat() << op.SQprod();
  if ( op.sumops().size() > 0 ) {
    o << "(" << op.sumops() << ")";
  }
  return o;
}

