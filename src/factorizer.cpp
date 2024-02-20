#include "factorizer.h"
#include <bitset> // for test

Factorizer::Factorizer(const std::vector<TermSum>& s)
{
  std::map<Orbital,const SlotType*> slotorbs;
  std::map<Matrix,const Tensor*> tensormats;

  for ( std::vector<TermSum>::const_iterator it=s.begin(); it!=s.end(); ++it) {
    // create an expression from the sum
    for ( TermSum::const_iterator i=it->begin();i!=it->end(); ++i) {
      Factor fac = _todouble(i->second);
      Term term = i->first;
      // add slottypes
      uint iorb = 0;
      Orbital orb;
      while ( (orb = term.orb(iorb)) != Orbital() ){
        slotorbs[orb] = _expression.add(Translators::orb2slot(orb));
        ++iorb;
      }
      TermSum sumt = term.resolve_permutations();
      for ( TermSum::const_iterator ist = sumt.begin();ist != sumt.end(); ++ist ) {
        Factor fact = _todouble(ist->second);
        fact *= fac;
        Term term = ist->first;
        Diagram diag = Translators::term2diagram(term,fact,slotorbs,_expression);
        if ( Input::iPars["prog"]["algo"] == 1 ) diag.binarize(_expression);
      }
    }
  }
}

SlotType Translators::orb2slot(const Orbital& orb)
{
  bool explspin = Input::iPars["prog"]["explspin"];
  switch (orb.type()) {
    case Orbital::Occ:
      if (explspin)
        if( orb.spin().type() == Spin::Up )
          return SlotType(Input::iPars["fact"]["nocc"],SlotType::OccA);
        else if( orb.spin().type() == Spin::Down )
          return SlotType(Input::iPars["fact"]["nocc"],SlotType::OccB);
        else
          error("What am i doing here?","Translators::orb2slot");
      else
        return SlotType(Input::iPars["fact"]["nocc"],SlotType::Occ);
      break;
    case Orbital::Virt:
      if (explspin)
        if( orb.spin().type() == Spin::Up ){
          return SlotType(Input::iPars["fact"]["nvir"],SlotType::VirtA);}
        else if( orb.spin().type() == Spin::Down ){ 
          return SlotType(Input::iPars["fact"]["nvir"],SlotType::VirtB);}
        else{
          error("What am i doing here?","Translators::orb2slot");}
      else
        return SlotType(Input::iPars["fact"]["nvir"],SlotType::Virt);
      break;
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
  Canonicalize(sts,slotorder,mat._threeelectronint);
  // look for default cuts
  std::string name = mat.plainname();

  // and parents??

  return Tensor(sts,name);
}

Diagram Translators::term2diagram(Term& term, Factor fact, const std::map< Orbital, const SlotType* >& slotorbs, Expression& expr)
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
  for (auto& m: term.get_mat()){
    if ( Input::iPars["prog"]["algo"] == 1 ){//ITF code
      m.itforder();
    }
    else if ( Input::iPars["prog"]["algo"] == 2 ){//ElemCo code
      m.elemcoorder();
    }
    else
      error("Unknown algorithm in prog, algo!");
    SlotTs sts;
    Connections con;
    //positions of orbitals (electron-order) in the tensor relative to the reference string
    //e.g. (ki|bj) relative to abijk -> 4213
    Slots positions;
    for (const auto& orb: m.get_orbs()){
      assert( slotorbs.count(orb) > 0 );
      sts.push_back(slotorbs.at(orb));
      int ipos = orbitals.find(orb);
      assert( ipos >= 0 );
      positions.push_back(ipos);
      con.bitmask[ipos] = true;
    }
    // stores info about performed permutations
    // e.g. kibj -> bkij slotorder = 2013
    slotorder = Slots();
    if ( m.type() == Ops::Exc0 || m.type() == Ops::Deexc0 || m.type() == Ops::Exc || Input::iPars["prog"]["algo"] < 2 ){
      Canonicalize(sts,slotorder,m._threeelectronint);
      // reorder positions according to the canonical order
      positions = positions.refarr(slotorder);
    }
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
  assert( std::abs(std::abs(_todouble(term.prefac())) - 1) < Numbers::verysmall );
  diag._fac = fact;
  if ( nbareops > 1 ) error("we can handle only upto one bare operator yet...", "Translators::term2diagram");
  expr._diagrams.push_back(diag);
  return diag;
}
