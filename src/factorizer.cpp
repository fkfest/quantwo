#include "factorizer.h"


Factorizer::Factorizer(const Sum< Term, TFactor >& s)
{
  std::map<Orbital,const SlotType*> slotorbs;
  std::map<Matrices,const Tensor*> tensormats;
  // create an expression from the sum
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    Factor fac = (Factor) i->second;
    Term term = i->first;
    // add slottypes
    uint iorb = 0;
    Orbital orb;
    while ( (orb = term.orb(iorb)) != Orbital() ){
      slotorbs[orb] = _expression.add(Translators::orb2slot(orb));
      ++iorb;
    }
    Product<Matrices>::const_iterator im;
    _foreach(im,term.mat()){
      tensormats[*im] = _expression.add(Translators::mat2tensor(*im,slotorbs));
    }
    Sum<Term,TFactor> sumt = term.resolve_permutations();
    for ( Sum<Term,TFactor>::const_iterator ist = sumt.begin();ist != sumt.end(); ++ist ) {
      Factor fact = (Factor) ist->second;
      fact *= fac;
      
      
      xout << fact << ist->first << std::endl;
    }
  }
  xout << _expression << std::endl;
}





SlotType Translators::orb2slot(const Orbital& orb)
{
  switch (orb.type()) {
    case Orbital::Occ:
      return SlotType(Input::iPars["fact"]["nocc"],SlotType::Occ);
    case Orbital::Virt:
      return SlotType(Input::iPars["fact"]["nvir"],SlotType::Virt);
    case Orbital::Act:
      return SlotType(Input::iPars["fact"]["nact"],SlotType::Act);
    case Orbital::GenT:
    {
      int ngen = Input::iPars["fact"]["nocc"]+Input::iPars["fact"]["nvir"];
      if ( Input::iPars["prog"]["multiref"] > 0 ) ngen += Input::iPars["fact"]["nact"];
      return SlotType(ngen,SlotType::GenT);
    }
    default:
      error("Unknown orbital type!","Translators::orb2slot");
  }
  return SlotType();
}

Tensor Translators::mat2tensor(const Matrices& mat, const std::map< Orbital, const SlotType* >& slotorbs)
{
  SlotTs sts;
  const Product<Orbital>& orbs = mat.orbitals();
  Product<Orbital>::const_iterator iorb;
  _foreach(iorb,orbs){
    assert( slotorbs.count(*iorb) > 0 );
    sts.push_back(slotorbs.at(*iorb));
  }
  // generate cuts??
  // and parents??
  
  return Tensor(sts,mat.plainname());
}
