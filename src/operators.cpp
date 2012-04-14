#include "operators.h"

/* SQOp::SQOp(std::string orb)
{
  _gender=(isupper((char)orb[0]) ? Creator : Annihilator);
  orb[0] = std::tolower(orb[0]);
  _orb = orb;
} */
SQOp::SQOp(Gender gender, Orbital orb)
{
  _gender=gender;
  _orb = orb;
}
SQOp::Gender SQOp::gender() const
{ return _gender; }
SQOp::Gender SQOp::genderPH() const
{ if (_orb.type()==Orbital::Occ && _gender==Creator) return Annihilator;
  if (_orb.type()==Orbital::Occ && _gender==Annihilator) return Creator;
  if (_orb.type()==Orbital::GenT) return Gen;
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
void SQOp::replace(Orbital orb1, Orbital orb2)
{
  if (_orb==orb1) _orb=orb2;
}

std::ostream & operator << (std::ostream & o, SQOp const & op)
{
  o << "\\op{" << op.orb() << "}";
  if ( op.gender()==SQOp::Creator )
    o << "^\\dg";
  
  return o;
}

Oper::Oper()
{  
  _prefac=1.0; 
  _type=Ops::None;
}

Oper::Oper(Ops::Type type)
{
  assert( InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  std::string name;
  _type=type;
  if (type == Ops::Fock )
    name="F";
  else if (type == Ops::FluctP )
    name="W";
  else
    name="X";
  create_Oper(name);
}
Oper::Oper(Ops::Type type, short int exccl, std::string name)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  _type=type;
  Orbital orb0(std::string("A"));
  Orbital orb1(std::string("I"));
  create_Oper(exccl,orb1,orb0,name);
}
Oper::Oper(Ops::Type type, short int exccl, 
           void * term, Orbital (*freeorb)(void * term, Orbital::Type type), std::string name)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  _type=type;
  Orbital orb0(freeorb(term,Orbital::Virt));
  Orbital orb1(freeorb(term,Orbital::Occ));
  create_Oper(exccl,orb1,orb0,name);
}
Oper::Oper(Ops::Type type, short int exccl, Orbital occ, Orbital virt, std::string name)
{
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  _type=type;
  create_Oper(exccl,occ,virt,name);
}
Oper::Oper(Ops::Type type, short int exccl, const Product< Orbital >& occs, const Product< Orbital >& virts, 
           std::string name)
{
  assert( occs.size() == virts.size() );
  assert( occs.size() == uint(exccl) );
  assert( !InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  _type=type;
  create_Oper(occs,virts,name);
}

void Oper::create_Oper(const std::string& name)
{
  assert( InSet(_type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  Product<Orbital> porbs;
  // operators with general indices
  Orbital orb(std::string("P"));
  _SQprod*=SQOp(SQOp::Creator,orb);
  porbs*=orb;
  _sumindx*=orb;
  orb=Orbital(std::string("Q"));
  porbs*=orb;
  _sumindx*=orb;
  if ( InSet(_type, Ops::Fock,Ops::XPert) ) {
    _SQprod*=SQOp(SQOp::Annihilator,orb);
    _prefac=1.0;
  } else {
   // we use chemical notation (PQ|RS) P^\dg R^\dg S Q
    orb=Orbital(std::string("R"));
    _SQprod*=SQOp(SQOp::Creator,orb);
    porbs*=orb;
    _sumindx*=orb;
    orb=Orbital(std::string("S"));
    _SQprod*=SQOp(SQOp::Annihilator,orb);
    porbs*=orb;
    _sumindx*=orb;
    orb=Orbital(std::string("Q"));
    _SQprod*=SQOp(SQOp::Annihilator,orb);
    _prefac=1.0/4.0;
  }
  _mat=Matrices(_type,porbs,name);
}
void Oper::create_Oper(short int const & exccl,Orbital const & occ, Orbital const & virt, std::string const & name)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  std::string excl;
  Product<Orbital> occs, virts; 
  for (short i=0; i<exccl ; ++i) {  
    if (i>0) excl=num2str(i,std::dec);
    occs *= Orbital(occ.name()+excl,occ.spin());
    virts *= Orbital(virt.name()+excl,virt.spin());
  }
  create_Oper(occs,virts,name);
}
void Oper::create_Oper(const Product< Orbital >& occs, const Product< Orbital >& virts, const std::string& name)
{
  assert( !InSet(_type, Ops::FluctP,Ops::Fock,Ops::XPert) );
  assert( occs.size() == virts.size() );
  Product<Orbital> porbs;
  // excitation and deexcitation operators
  const Product<Orbital> 
    * p_orb0 = &virts,
    * p_orb1 = &occs;
  if (InSet(_type, Ops::Deexc,Ops::Deexc0)) std::swap(p_orb0,p_orb1);
  _prefac = 1.0;
  for (unsigned short i = 0; i < p_orb0->size(); ++i) {  
    _SQprod*=SQOp(SQOp::Creator, (*p_orb0)[i]);
    porbs *= (*p_orb0)[i];
    _sumindx *= (*p_orb0)[i];
    if (InSet(_type, Ops::Exc0,Ops::Deexc0)) 
      _fakesumindx *= (*p_orb0)[i];
    _SQprod *= SQOp(SQOp::Annihilator, (*p_orb1)[i]);
    porbs *= (*p_orb1)[i];
    _sumindx *= (*p_orb1)[i];
    if (InSet(_type, Ops::Exc0,Ops::Deexc0)) 
      _fakesumindx *= (*p_orb1)[i];
    _prefac = _prefac/double(i+1);
  }
  _prefac=_prefac*_prefac;
  _mat=Matrices(_type,porbs,name);
}

Matrices Oper::mat() const
{ return _mat;}
Product< SQOp > Oper::SQprod() const
{ return _SQprod;}
double Oper::prefac() const
{ return _prefac;}
Product< Orbital > Oper::sumindx() const
{ return _sumindx;}
Product< Orbital > Oper::realsumindx() const
{
  Product< Orbital > realsum;
  for (unsigned int i=0; i<_sumindx.size(); i++ )  
    if (_fakesumindx.find(_sumindx[i])<0)
      realsum*=_sumindx[i];
  return realsum;
}



std::ostream & operator << (std::ostream & o, Oper const & op)
{
  if ( std::abs(std::abs(op.prefac()) - 1.0) > Numbers::small) o << op.prefac();
  if (op.realsumindx().size()>0) o <<"\\sum_{"<<op.realsumindx()<<"}";
  o <<op.mat() << op.SQprod();
  return o;
}

