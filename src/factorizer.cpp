#include "factorizer.h"
#include <bitset> // for test

Factorizer::Factorizer(const Sum< Term, TFactor >& s)
{
  std::map<Orbital,const SlotType*> slotorbs;
  std::map<Matrices,const Tensor*> tensormats;

///////////////////////////////////// TEST ///////////////////////////////
  std::bitset<8> bs;
  bs[1] = true;
  bs [4] = true;
  bs [6] = true;
  bs [7] = true;
  xout << "bitset: " << bs << " " << bs.count() << std::endl;
  xout << "bitset2: " << (bs >> 2) << " " << (bs >> 2).count() << std::endl;

///////////////////////////////////// END TEST ///////////////////////////////
  
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
      Diagram diag = Translators::term2diagram(ist->first,slotorbs);
      xout << diag; 
      diag.binarize(_expression);
      // contractions
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
  Slots slotorder;
  Canonicalize(sts,slotorder);
  // look for default cuts
  std::string name = mat.plainname();
  
  // and parents??
  
  return Tensor(sts,mat.plainname());
}

Diagram Translators::term2diagram(const Term& term, const std::map<Orbital,const SlotType*>& slotorbs)
{
  Diagram diag;

  // all orbitals in term:
  Array<Orbital> orbitals;
  Array<Orbital>::const_iterator itorb;
  Orbital orb;
  uint iorb = 0; 
  while ( (orb = term.orb(iorb)) != Orbital() ){
    orbitals.push_back(orb);
    ++iorb;
  }
  // put the orbitals to diag
  if ( orbitals.size() > MAXNINDICES )
    error("Too many indices in the term. Increase MAXNINDICES!","Translators::term2diagram");
  _foreach(itorb,orbitals){
    assert( slotorbs.count(*itorb) > 0 );
    diag._slottypes.push_back(slotorbs.at(*itorb));
  }
  Slots slotorder;
  // use canonical order - then the intermediates will be in canonical order automatically 
  Canonicalize(diag._slottypes,slotorder);
  orbitals = orbitals.refarr(slotorder); 

  xout << "orbitals: " << orbitals << std::endl;

  Product<Matrices>::const_iterator im;
  uint nbareops = 0;
  _foreach(im,term.mat()){
    const Product<Orbital>& orbs = im->orbitals();
    Product<Orbital>::const_iterator itorb;
    SlotTs sts;
    Connections con;
    Slots positions;
    _foreach(itorb,orbs){
      assert( slotorbs.count(*itorb) > 0 );
      sts.push_back(slotorbs.at(*itorb));
      int ipos = orbitals.find(*itorb);
      assert( ipos >= 0 );
      positions.push_back(ipos);
      con.bitmask[ipos] = true;
    }
    slotorder = Slots();
    Canonicalize(sts,slotorder);
    // reorder positions according to the canonical order
    positions = positions.refarr(slotorder);
    // set slotref for bitset
    con.slotref.resize(positions.size());
    for ( uint ist = 0; ist < positions.size(); ++ist ){
      assert( con.bitmask[positions[ist]] );
      uint icnt = (con.bitmask>>(positions[ist]+1)).count();
      assert( icnt < con.slotref.size() );
      con.slotref[icnt] = ist;
    }
    if ( im->type() == Ops::Exc0 || im->type() == Ops::Deexc0 ) {
      diag._tensors.push_front(DiagramTensor(con,im->plainname()));
      ++nbareops;
    } else
      diag._tensors.push_back(DiagramTensor(con,im->plainname()));
  }
  if ( nbareops > 1 ) error("we can handle only upto one bare operator yet...", "Translators::term2diagram");
  return diag;
}
