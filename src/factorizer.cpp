#include "factorizer.h"
#include <bitset> // for test

Factorizer::Factorizer(const TermSum& s)
{
  std::map<Orbital,const SlotType*> slotorbs;
  std::map<Matrix,const Tensor*> tensormats;

  // create an expression from the sum
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    #ifdef _RATIONAL
    Factor fac = boost::rational_cast<Factor>(i->second);
    #else
    Factor fac = (Factor) i->second;
    #endif
    Term term = i->first;
    // add slottypes
    uint iorb = 0;
    Orbital orb;
    while ( (orb = term.orb(iorb)) != Orbital() ){
      slotorbs[orb] = _expression.add(Translators::orb2slot(orb));
      ++iorb;
    }
//    for (const auto& m: term.mat()){
//      tensormats[m] = _expression.add(Translators::mat2tensor(m,slotorbs));
//    }
    TermSum sumt = term.resolve_permutations();
    for ( TermSum::const_iterator ist = sumt.begin();ist != sumt.end(); ++ist ) {
      #ifdef _RATIONAL
      Factor fact = boost::rational_cast<Factor>(ist->second);
      #else
      Factor fact = (Factor) ist->second;
      #endif
      fact *= fac;
      Diagram diag = Translators::term2diagram(ist->first,fact,slotorbs,_expression);
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

Tensor Translators::mat2tensor(const Matrix& mat, const std::map< Orbital, const SlotType* >& slotorbs)
{
  SlotTs sts;
  const Product<Orbital>& orbs = mat.orbitals();
  for (const auto& orb: orbs){
    assert( slotorbs.count(orb) > 0 );
    sts.push_back(slotorbs.at(orb));
  }
  Slots slotorder;
  Canonicalize(sts,slotorder);
  // look for default cuts
  std::string name = mat.plainname();
  
  // and parents??
  
  return Tensor(sts,name);
}

Diagram Translators::term2diagram(const Term& term, Factor fact, const std::map< Orbital, const SlotType* >& slotorbs, const Expression& expr)
{
  const std::string& resultt = Input::sPars["syntax"]["result"];
  Diagram diag;
  
  // all orbitals in term:
  Array<Orbital> orbitals;
  Orbital orb;
  uint iorb = 0; 
  while ( (orb = term.orb(iorb)) != Orbital() ){
    orbitals.push_back(orb);
    ++iorb;
  }
  // put the orbitals to diag
  if ( orbitals.size() > MAXNINDICES )
    error("Too many indices in the term. Increase MAXNINDICES!","Translators::term2diagram");
  for (const auto& orb: orbitals){
    assert( slotorbs.count(orb) > 0 );
    diag._slottypes.push_back(slotorbs.at(orb));
  }
  Slots slotorder;
  // use canonical order - then the intermediates will be in canonical order automatically 
  Canonicalize(diag._slottypes,slotorder);
  orbitals = orbitals.refarr(slotorder); 

  uint nbareops = 0;
  for (const auto& m: term.mat()){
    const Product<Orbital>& orbs = m.orbitals();
    SlotTs sts;
    Connections con;
    Slots positions;
    for (const auto& orb: orbs){
      assert( slotorbs.count(orb) > 0 );
      sts.push_back(slotorbs.at(orb));
      int ipos = orbitals.find(orb);
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
    if ( m.type() == Ops::Exc0 || m.type() == Ops::Deexc0 ) {
      Tensor ten(sts,resultt);
      const Tensor * pTen = expr.find(ten,false);
      // the first tensor is the result tensor (residuum)
      diag.add(DiagramTensor(con,m.plainname()),pTen,true);
      ++nbareops;
    } else
      diag.add(DiagramTensor(con,m.plainname()));
  }
  #ifdef _RATIONAL
  assert( std::abs(std::abs(boost::rational_cast<Factor>(term.prefac())) - 1) < Numbers::verysmall );
  #else
  assert( std::abs(std::abs(term.prefac()) - 1) < Numbers::verysmall );
  #endif
  diag._fac = fact;
  if ( nbareops > 1 ) error("we can handle only upto one bare operator yet...", "Translators::term2diagram");
  return diag;
}
