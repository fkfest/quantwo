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
  std::string name;
  _type=type;
  if ( InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) )
  {
    Orbital orb0(std::string("P")); //dummy
    Orbital orb1(std::string("Q")); //dummy
    if (type == Ops::Fock )
      name="F";
    else if (type == Ops::FluctP )
      name="W";
    else
      name="X";
    create_Oper(0,orb1,orb0,name);
  }
  else
    error("Use Oper(exccl,Type) to construct Cluster operators","Oper::Oper");
}
Oper::Oper(Ops::Type type, short int exccl, std::string name)
{
  _type=type;
  if ( InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) )
    error("Use Oper(Type) to construct Hamiltonian","Oper::Oper");
  else
  {
    Orbital orb0(std::string("A"));
    Orbital orb1(std::string("I"));
    create_Oper(exccl,orb1,orb0,name);
  }
}
Oper::Oper(Ops::Type type, short int exccl, 
           void * term, Orbital (*freeorb)(void * term, Orbital::Type type), std::string name)
{
  _type=type;
  if ( InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) )
    error("Use Oper(Type) to construct Hamiltonian","Oper::Oper");
  else
  {
    Orbital orb0(freeorb(term,Orbital::Virt));
    Orbital orb1(freeorb(term,Orbital::Occ));
    create_Oper(exccl,orb1,orb0,name);
  }
}
Oper::Oper(Ops::Type type, short int exccl, Orbital occ, Orbital virt, std::string name)
{
  _type=type;
  if ( InSet(type, Ops::FluctP,Ops::Fock,Ops::XPert) )
    error("Use Oper(Type) to construct Hamiltonian","Oper::Oper");
  create_Oper(exccl,occ,virt,name);
}

void Oper::create_Oper(short int const & exccl,Orbital const & occ, Orbital const & virt, std::string const & name)
{
  Product<Orbital> porbs;
  std::string excl;
  if ( InSet(_type, Ops::FluctP,Ops::Fock,Ops::XPert) )
  { // operators with general indices
    Orbital orb(std::string("P"));
    _SQprod*=SQOp(SQOp::Creator,orb);
    porbs*=orb;
    _sumindx*=orb;
    orb=Orbital(std::string("Q"));
    porbs*=orb;
    _sumindx*=orb;
    if ( InSet(_type, Ops::Fock,Ops::XPert) )
    {
      _SQprod*=SQOp(SQOp::Annihilator,orb);
      _prefac=1.0;
    }
    else
    {
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
  }
  else
  { // excitation and deexcitation operators
    Orbital orb0=virt;
    Orbital orb1=occ;
  
    if (InSet(_type, Ops::Deexc,Ops::Deexc0)) std::swap(orb0,orb1);
    Orbital orb;
    _prefac=1.0;
    for (short i=0; i<exccl ; ++i)
    {  
      if (i>0) excl=num2str(i,std::dec);
      orb=Orbital(orb0.name()+excl,orb0.spin());
      _SQprod*=SQOp(SQOp::Creator,orb);
      //in the case of blank excitation and deexcitation operators use blank matrix
 //     if (type!=Ops::Exc0 && type!=Ops::Deexc0) 
 //     {
        porbs*=orb;
        _sumindx*=orb;
        if (InSet(_type, Ops::Exc0,Ops::Deexc0)) 
          _fakesumindx*=orb;
 //     }
      orb=Orbital(orb1.name()+excl,orb1.spin());
      _SQprod*=SQOp(SQOp::Annihilator,orb);
 //     if (type!=Ops::Exc0 && type!=Ops::Deexc0)
 //     {
        porbs*=orb;
        _sumindx*=orb;
        if (InSet(_type, Ops::Exc0,Ops::Deexc0)) 
          _fakesumindx*=orb;
        _prefac=_prefac/double(i+1);
 //     }
    }
    _prefac=_prefac*_prefac;
  }
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

