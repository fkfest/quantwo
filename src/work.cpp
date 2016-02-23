#include "work.h"

Sum< Term, TFactor > Q2::reduceSum(Sum< Term, TFactor > s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  bool brill = ( Input::iPars["prog"]["brill"] > 0 );
  bool quan3 = ( Input::iPars["prog"]["quan3"] > 0 );
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  bool usefock = Input::iPars["prog"]["usefock"];
  usefock = usefock && (Input::iPars["prog"]["noorder"]>0);
  Sum<Term,TFactor> sum,sum1;
  Term term,term1;
  bool added;
  TFactor prefac;

  if (usefock){
    // replace h matrices by fock matrices
    say("Use Fock matrices...");
    s = OneEl2Fock(s);
  }

  say("Reduce sum of terms");
  _xout3(s << std::endl);

  say("Antisymmetry...");
  sum.clear();
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    // expand antisymmetrized integrals
    sum1 = term.expand_antisym();
    sum1 *= i->second;
    sum += sum1;
  }
  _xout3(sum << std::endl);
  
  say("Kroneckers...");
  sum = Kroneckers(sum);
  _xout3(sum << std::endl);
  
  // set remaining general indices to the occupied or active space
  say("Handle general indices...");
  sum = GeneralIndices(sum);
  sum = ZeroTerms(sum);
  _xout3(sum << std::endl);

  if (spinintegr){
    // bring all the density matrices into singlet-order
    say("Singlet order...");
    sum = SingletDM(sum);
    _xout3(sum << std::endl);
  }

  say("Connections...");
  s = sum;
  sum.clear();
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    // generate Kallay's "triplets of integers"
    term.matrixkind();
    // set connections "map" of matrices
    term.setmatconnections();
    // is the term properly connected?
    if (!term.properconnect()) {
      continue;
    }
    sum += std::make_pair(term,i->second);
  }
  _xout3(sum << std::endl);
  if (quan3) {
    say("count electrons (a posteriori)...");
    s = sum;
    sum.clear();
    for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
      term=i->first;
      term.deleteNoneMats();
      term.setmatconnections();
      
      prefac=i->second*term.prefac();
      // remove prefactors in terms
      term.reset_prefac();
      added=false;
      for ( Sum<Term,TFactor>::iterator k=sum.begin();k!=sum.end(); ++k) {
        Permut perm;
        term1=k->first;
        if (term.equal(term1,perm)) {
          sum.erase(k);
          term1+=std::make_pair<Permut,TFactor>(perm,prefac);
          //std::cout<<"term old" << term1 <<std::endl;
          if ( !term1.term_is_0(minfac) ) sum+=term1;
          added=true;
          break;
        }
      }
      if (!added) {
        term+=std::make_pair<Permut,TFactor>(Permut(),prefac);
        if ( !term.term_is_0(minfac) ) sum+=term;
        //std::cout<<"term new" << term <<std::endl;
      }
//      Termel termel(term);
//      sum += std::make_pair(term,i->second);
    }
//    _xout3(sum << std::endl);
    return sum;
  }
  say("Diagrams and Spin-integration...");
  s = sum;
  sum.clear();
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    // use Brilloin condition 
    if ( brill && term.brilloin() ) continue;
    // remove "None" matrices
    term.deleteNoneMats();
    term.setmatconnections();
    term.spinintegration(spinintegr);
    sum += std::make_pair(term,i->second);
  }
  _xout3(sum << std::endl);

  say("Equal terms and permutations...");
  //std::cout << "SUM:" << sum << std::endl;
  sum = EqualTerms(sum,minfac);
  say("Remove terms with small prefactors...");
  // now remove everything with small prefactor
  sum = SmallTerms(sum,minfac);
  
  return sum;
}
Sum< Term, TFactor > Q2::Kroneckers(Sum< Term, TFactor > s)
{
  Sum<Term,TFactor> sum;
  Term term;
  bool printed = false;
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    // remove Kroneckers
    term.reduceTerm();
    if (term.kProd().size()>0 && term.prefac() != 0) {
      if (!printed)
        say("Still some Kroneckers left. Reducing electrons according to them...","Q2::reduceSum");
      term.reduceElectronsInTerm();
      if (!printed)
        say("... and transforming them to matrices...","Q2::reduceSum");
      term.krons2mats();
      printed = true;
    }
    sum += std::make_pair(term,i->second);
  }
  return sum;
}
Sum< Term, TFactor > Q2::OneEl2Fock(Sum< Term, TFactor > s)
{
  bool multiref = (Input::iPars["prog"]["multiref"] > 0);
  Term term;
  Sum<Term,TFactor> sum,sum1;
  Sum< Term, TFactor >::const_iterator its;
  _foreach(its,s){
    term = its->first;
    sum1 = term.oneel2fock(multiref);
    sum1 *= its->second;
    sum += sum1;
  }
  return sum;
}

Sum< Term, TFactor > Q2::SingletDM(Sum< Term, TFactor > s)
{
  Term term;
  uint iter = 0;
  for ( iter = 0; has_nonsingldm(s) && iter < 1000; ++iter ){
    xout << "has nonsingldm" << std::endl;
    Sum<Term,TFactor> sum,sum1;
    for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
      term = i->first;
      sum1 = term.dm2singlet();
      sum1 *= i->second;
      sum += sum1;
    }
    s = Kroneckers(sum);
  }
  if (iter == 1000) error("Could not bring density matrices into singlet order in 1000 iterations!");
  return s;
}
bool Q2::has_nonsingldm(const Sum< Term, TFactor >& s)
{
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) 
    if (i->first.has_nonsingldm()) return true;
  return false;
}
Sum< Term, TFactor > Q2::GeneralIndices(Sum< Term, TFactor > s)
{
  uint iter = 0;
  for ( iter = 0; has_generalindices(s) && iter < 1000; ++iter ){
    Sum<Term,TFactor> sum;
    for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
      Term term = i->first;
      Sum<Term,TFactor> sum1 = term.removegeneralindices();
      sum1 *= i->second;
      sum += sum1;
    }
    s = sum;
  }
  if (iter == 1000) error("Could not remove all general indices in 1000 iterations!");
  return s;
}
bool Q2::has_generalindices(const Sum< Term, TFactor >& s)
{
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) 
    if (i->first.has_generalindices()) return true;
  return false;
}
Sum< Term, TFactor > Q2::ZeroTerms(Sum< Term, TFactor > s)
{
  Sum<Term,TFactor> sum;
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    if ( ! i->first.removeit()){
      sum += std::make_pair(i->first,i->second);
    }
  }
  return sum;
}

Sum< Term, TFactor > Q2::EqualTerms(Sum< Term, TFactor > s, double minfac)
{
  Sum<Term,TFactor> sum;
  Term term,term1;
  TFactor prefac;
  bool added;
  for ( Sum<Term,TFactor>::const_iterator j=s.begin();j!=s.end(); ++j) {
    term=j->first;
    prefac=j->second*term.prefac();
    // remove prefactors in terms
    term.reset_prefac();
    
//     xout << "Term: " << term << std::endl;
//     UniGraph ug(term);
//     xout << ug << std::endl;
    
    
    added=false;
    for ( Sum<Term,TFactor>::iterator k=sum.begin();k!=sum.end(); ++k) {
      Permut perm;
      term1=k->first;
      if (term.equal(term1,perm)) {
        sum.erase(k);
        term1+=std::make_pair<Permut,TFactor>(perm,prefac);
//         std::cout<<"term old" << term1 <<std::endl;
        if ( !term1.term_is_0(minfac) ) sum+=term1;
        added=true;
        break;
      }
    }
    if (!added) {
      term+=std::make_pair<Permut,TFactor>(Permut(),prefac);
      if ( !term.term_is_0(minfac) ) sum+=term;
//       std::cout<<"term new" << term <<std::endl;
    }
  }
  return sum;
}
Sum< Term, TFactor > Q2::SmallTerms(Sum< Term, TFactor > s, double minfac)
{
  Sum<Term,TFactor> sum;
  Term term,term1;
  for ( Sum<Term,TFactor>::const_iterator j=s.begin(); j!=s.end(); ++j) {
    if ( _todouble(_abs(j->second)) < minfac ) continue;
    term1 = term = j->first;
    if ( _todouble(_abs(term.prefac())) < minfac ) continue;
    term1.reset_prefac();
    const Sum<Permut,TFactor>& perms = term.perm();
    for ( Sum<Permut,TFactor>::const_iterator it = perms.begin(); it != perms.end(); ++it ){
      if ( _todouble(_abs(it->second)) < minfac ) continue;
      term1 += *it; 
    }
    sum += term1;
  }
  return sum;
}
Sum< Term, TFactor > Q2::ResolvePermutaions(Sum< Term, TFactor > s)
{
  Sum<Term,TFactor> sum;
  Term term,term1;
  for ( Sum<Term,TFactor>::const_iterator j=s.begin(); j!=s.end(); ++j) {
    term1 = term = j->first;
    const Sum<Permut,TFactor>& perms = term.perm();
    term.reset_prefac();
    for ( Sum<Permut,TFactor>::const_iterator it = perms.begin(); it != perms.end(); ++it ){
      term1 = term; 
      term1.permute(it->first);
      term1 *= it->second;
      sum += term1;
    }
  }
  return sum;
}

Sum< Term, TFactor > Q2::postaction(Sum< Term, TFactor > s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  std::string divide = Input::sPars["act"]["divide"];
  Sum<Term,TFactor> sum = s;
  if ( divide != "" && divide != "1" ){
    Finput div(true);
    div.addline(divide);
    div.analyzeq();
    Sum<Term,TFactor> sum1;
    Sum<Permut,TFactor> divperm, divpermadd;
    sum1 = div.sumterms();
    Matrices OneMat;
    for ( Sum<Term,TFactor>::const_iterator it=sum1.begin();it!=sum1.end(); ++it) {
      if ( it->first.mat().size() != 1 && !(it->first.mat().front() == OneMat) )
        error("Not an One-Matrix in divide!");
      divpermadd.clear();
      divpermadd += it->first.perm();
      divpermadd *= it->first.prefac();
      divpermadd *= it->second;
      divperm += divpermadd;
    }
    sum1.clear();
//    xout << "permutations: " << divperm << std::endl;
    Term term;
    for ( Sum<Term,TFactor>::const_iterator i=sum.begin();i!=sum.end(); ++i) {
      term=i->first;
      Sum<Permut,TFactor> perms = term.perm();
      perms /= divperm;
      term.setperm(perms);
      sum1 +=  std::make_pair(term,i->second);
    }
    sum = sum1;
    sum1.clear();
    sum = ResolvePermutaions(sum);
    say("Equal terms and permutations...");
    sum = EqualTerms(sum,minfac);
    say("Remove terms with small prefactors...");
    sum = SmallTerms(sum,minfac);
  }
  return sum;
}

Sum< Term, TFactor > Q2::normalOrderPH(Sum< Term, TFactor > s)
{
  Sum<Term,TFactor> sum,sum0;
  Term term;
  say("Normal ordering");
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    sum0 += term.normalOrderPH_fullyContractedOnly();
    sum0 *= i->second;
    sum += sum0;
    sum0=Sum<Term,TFactor>();
  }
  return sum;
}
Sum< Term, TFactor > Q2::wick(Sum< Term, TFactor > s)
{
  int iwick = Input::iPars["prog"]["wick"];
  bool genwick = (iwick > 1);
  int noorder = Input::iPars["prog"]["noorder"];
  if (!genwick && noorder > 0 ) error("Cannot have non-ordered Hamiltonian with wick<2. Either set noorder=0 or wick=2");
  Sum<Term,TFactor> sum,sum0;
  Term term;
  _xout3(s << std::endl);
  say("Wick's theorem");
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    sum0 += term.wickstheorem(genwick,noorder);
    sum0 *= i->second;
    sum += sum0;
    sum0=Sum<Term,TFactor>();
  }
  return sum;
}

void Q2::printdiags(Output* pout, Sum< Term, TFactor > s)
{
  say("Diagrams...");
  *(pout->pout) << " Diagrams: " << std::endl;
    *(pout->pout) << std::endl;
  for ( Sum<Term,TFactor>::const_iterator it = s.begin(); it != s.end(); ++it){
    it->first.printdiag(pout);
  }
}

void Q2::printalgo(std::ofstream& out, Sum< Term, TFactor > s)
{
  const std::string& resultt = Input::sPars["syntax"]["result"];
  say("Algorithm...");
  
  Factorizer fact(s);
  
  out << "algorithm..." << std::endl;
  // external indices
  // .....
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i) {
    out << resultt << "[abij]" << " += ";
    // call term ...
    out << std::endl;
  }
}

