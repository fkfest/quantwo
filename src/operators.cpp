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
}
SQOpT::Gender SQOp::gender() const
{ return _gender; }
SQOpT::Gender SQOp::genderPH() const
{ if (_orb.type()==Orbital::Occ && _gender==SQOpT::Creator) return SQOpT::Annihilator;
  if (_orb.type()==Orbital::Occ && _gender==SQOpT::Annihilator) return SQOpT::Creator;
  if (_orb.type()==Orbital::GenT || _orb.type()==Orbital::Act) return SQOpT::Gen;
  return _gender; }
Orbital SQOp::orb() const
{ return _orb; }
bool SQOp::operator==(const SQOp& o) const
{ return _gender==o._gender && _orb==o._orb; }
bool SQOp::operator<(const SQOp& o) const
{
  if ( _gender<o._gender )
    return true;
  if (o._gender<_gender )
    return false;
  return _orb<o._orb;
}
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

Oper::Oper(Ops::Type type, bool antisym, Term* pTerm)
{
  assert( InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  std::string name;
  _type=type;
  p_Term = pTerm;
  if (type == Ops::Fock )
    name="F";
  if (type == Ops::OneEl )
    name="h";
  else if (type == Ops::FluctP )
    name="W";
  else
    name="X";
  create_Oper(name,antisym);
}
Oper::Oper(Ops::Type type, short int exccl, std::string name, int lm, Term* pTerm)
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
  create_Oper(exccl,orb1,orb0,name,lm);
}
Oper::Oper(Ops::Type type, short int exccl, Orbital occ, Orbital virt, std::string name, int lm, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = pTerm;
  create_Oper(exccl,occ,virt,name,lm);
//   xout << *this << std::endl;
}
Oper::Oper(Ops::Type type, short int exccl, const std::map< Orbital::Type, Orbital >& orbnames, 
           const std::vector<OrbitalTypes>& orbtypes, std::string name, int lm, Term* pTerm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  assert( orbtypes.size() == 2 );
  assert( int(orbtypes[0].size()) == exccl );
  assert( int(orbtypes[1].size()) == exccl+lm );
  _type=type;
  p_Term = pTerm;
  create_Oper(exccl,orbnames,orbtypes,name,lm);
}

Oper::Oper(Ops::Type type, short int exccl, const Product< Orbital >& occs, const Product< Orbital >& virts, 
           std::string name, int lm, Term* pTerm)
{
  assert( occs.size() + lm == virts.size() );
  assert( occs.size() == uint(exccl) );
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = pTerm;
  create_Oper(occs,virts,name);
}
Oper::Oper(Ops::Type type, const Product< Orbital >& orbs, std::string name, int lm)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  _type=type;
  p_Term = 0;
  create_Oper(orbs,name,lm);
}

Oper::Oper(Ops::Type type, short int exccl, const std::vector<OrbitalTypes>& orbtypes, 
           std::string name, int lm, Term* pTerm)
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
  create_Oper(exccl,orbnames,orbtypes,name,lm);
}

void Oper::create_Oper(const std::string& name, bool antisym)
{
  assert( InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  Product<Orbital> porbs;
  Matrix::Spinsym spinsym = Matrix::Singlet;
  Electron el = 1;
  if (p_Term) el = p_Term->nextelectron();
  // operators with general indices
  Orbital orb(std::string("P"),el);
  _SQprod*=SQOp(SQOpT::Creator,orb);
  porbs*=orb;
  _orbs.insert(orb);
  _sumorbs.insert(orb);
  orb=Orbital(std::string("Q"),el);//same electron as in P
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
    orb=Orbital(std::string("R"),el);
    _SQprod*=SQOp(SQOpT::Creator,orb);
    porbs*=orb;
    _orbs.insert(orb);
    _sumorbs.insert(orb);
    orb=Orbital(std::string("S"),el);//same electron as in R
    _SQprod*=SQOp(SQOpT::Annihilator,orb);
    porbs*=orb;
    _orbs.insert(orb);
    _sumorbs.insert(orb);
    orb=porbs[1];
    _SQprod*=SQOp(SQOpT::Annihilator,orb);
    _prefac /= 4;
  }
  short npairs = porbs.size()/2;
  _mat=Matrix(_type,porbs,npairs,0,name,spinsym,antisym);
}
void Oper::create_Oper(short int const & exccl,Orbital const & occ, Orbital const & virt, 
                       std::string const & name, int lm)
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
  create_Oper(occs,virts,name);
}
void Oper::create_Oper(const short int& exccl, const std::map< Orbital::Type, Orbital >& orbnames, 
                       const std::vector<OrbitalTypes>& orbtypes, const std::string& name, int lm)
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
  create_Oper(occs,virts,name);
}

void Oper::create_Oper(const Product< Orbital >& occs, const Product< Orbital >& virts, const std::string& name)
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
  create_Oper(porbs,name,p_orb0->size()-p_orb1->size());
}

void Oper::create_Oper(const Product< Orbital >& orbs, const std::string& name, int lm)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::OneEl,Ops::XPert) );
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  bool noprefac = (Input::iPars["prog"]["nobrafac"]) && InSet(_type, Ops::Exc0,Ops::Deexc0);
  bool contrexcop = Input::iPars["prog"]["contrexcop"];
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
      if (contrexcop)
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
  if (contrexcop) {
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
  _mat=Matrix(_type,orbs,npairs,lm,name,spinsym);

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
  if (op.sumorbs().size()>0) o <<"\\sum_{"<<op.sumorbs()<<"}";
  o <<op.mat() << op.SQprod();
  return o;
}

