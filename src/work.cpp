#include "work.h"

TermSum Q2::evalEq(Finput& finput)
{
//       Input::verbose = 2;
  _xout1(finput << std::endl);
  TermSum sum_finp(finput.sumterms());
  _xout2(" = " << sum_finp << std::endl);
  TermSum sum_NO;
  if ( Input::iPars["prog"]["wick"] == 0 )
    sum_NO = Q2::normalOrderPH(sum_finp);
  else
    sum_NO = Q2::wick(sum_finp);
  _xout2(" = " << sum_NO << std::endl);
  TermSum sum_final1(Q2::reduceSum(sum_NO)),
    sum_final(Q2::postaction(sum_final1));
  _xout1(" = " << sum_final << std::endl);
  finput.sumterms(sum_final);
  return finput.sumterms();
}

TermSum Q2::reduceSum(TermSum s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  bool brill = ( Input::iPars["prog"]["brill"] > 0 );
  bool quan3 = ( Input::iPars["prog"]["quan3"] > 0 );
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  bool usefock = Input::iPars["prog"]["usefock"];
  usefock = usefock && (Input::iPars["prog"]["noorder"]>0);
  // 13.12.2016: temporary hack until a proper insert of intermediate tensors is implemented
  bool replaceE0 = Input::iPars["prog"]["replacee0"];
  replaceE0 = replaceE0 && (Input::iPars["prog"]["noorder"]>0);
  bool timing = ( Input::iPars["prog"]["cpu"] > 0 );
  std::clock_t c_start=0;
  TermSum sum,sum1;
  Term term,term1;
  bool added;
  TFactor prefac;

  if (usefock){
    if (timing) c_start = std::clock();
    // replace h matrices by fock matrices
    say("Use Fock matrices...");
    s = OneEl2Fock(s);
    if (timing) _CPUtiming("",c_start,std::clock());
  }
  _xout3(s << std::endl);

  if (replaceE0){
    if (timing) c_start = std::clock();
    // replace E^{\snam{0}} by <0|\op F|0>
    say("Replace E^0...");
    s = ReplaceE0(s);
    if (timing) _CPUtiming("",c_start,std::clock());
  }

  say("Reduce sum of terms");
  _xout3(s << std::endl);

  say("Antisymmetry...");
  if (timing) c_start = std::clock();
  sum.clear();
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    // expand antisymmetrized integrals
    sum1 = term.expand_antisym();
    sum1 *= i->second;
    sum += sum1;
  }
  _xout3(sum << std::endl);
  if (timing) _CPUtiming("",c_start,std::clock());

  say("Kroneckers...");
  if (timing) c_start = std::clock();
  sum = Kroneckers(sum);
  _xout3(sum << std::endl);
  if (timing) _CPUtiming("",c_start,std::clock());

  // set remaining general indices to the occupied or active space
  say("Handle general indices...");
  if (timing) c_start = std::clock();
  sum = GeneralIndices(sum);
  sum = ZeroTerms(sum);
  _xout3(sum << std::endl);
  if (timing) _CPUtiming("",c_start,std::clock());

  if (spinintegr){
    // bring all the density matrices into singlet-order
    say("Singlet order...");
    if (timing) c_start = std::clock();
    sum = SingletDM(sum);
    _xout3(sum << std::endl);
    if (timing) _CPUtiming("",c_start,std::clock());
  }

  say("Connections...");
  if (timing) c_start = std::clock();
  s = sum;
  sum.clear();
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
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
  if (timing) _CPUtiming("",c_start,std::clock());
  if (quan3) {
    say("count electrons (a posteriori)...");
    if (timing) c_start = std::clock();
    s = sum;
    sum.clear();
    for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
      term=i->first;
      term.deleteNoneMats();
      term.setmatconnections();

      prefac=i->second*term.prefac();
      // remove prefactors in terms
      term.reset_prefac();
      added=false;
      for ( TermSum::iterator k=sum.begin();k!=sum.end(); ++k) {
        Permut perm;
        term1=k->first;
        if (term.equal(term1,perm)) {
          sum.erase(k);
          term1+=std::make_pair(perm,prefac);
          //std::cout<<"term old" << term1 <<std::endl;
          if ( !term1.term_is_0(minfac) ) sum+=term1;
          added=true;
          break;
        }
      }
      if (!added) {
        term+=std::make_pair(Permut(),prefac);
        if ( !term.term_is_0(minfac) ) sum+=term;
        //std::cout<<"term new" << term <<std::endl;
      }
//      Termel termel(term);
//      sum += std::make_pair(term,i->second);
    }
//    _xout3(sum << std::endl);
    if (timing) _CPUtiming("",c_start,std::clock());
    return sum;
  }
  say("Diagrams and Spin-integration...");
  if (timing) c_start = std::clock();
  s = sum;
  sum.clear();
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
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
  if (timing) _CPUtiming("",c_start,std::clock());

  say("Equal terms and permutations...");
  if (timing) c_start = std::clock();
  //std::cout << "SUM:" << sum << std::endl;
  sum = EqualTerms(sum,minfac);
  if (timing) _CPUtiming("",c_start,std::clock());
  say("Remove terms with small prefactors...");
  if (timing) c_start = std::clock();
  // now remove everything with small prefactor
  sum = SmallTerms(sum,minfac);
  if (timing) _CPUtiming("",c_start,std::clock());
  // put overlaps if needed
  sum = VirtSpace(sum);

  return sum;
}
TermSum Q2::Kroneckers(const TermSum& s)
{
  TermSum sum;
  Term term;
  bool printed = false;
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
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
TermSum Q2::OneEl2Fock(const TermSum& s)
{
  bool multiref = (Input::iPars["prog"]["multiref"] > 0);
  int iusefock = Input::iPars["prog"]["usefock"];
  std::string decoration = "";
  if (iusefock > 1 && multiref) {
    // use only closed-shell fock
    multiref = false;
    decoration = "c";
  }
  Term term;
  TermSum sum,sum1;
  TermSum::const_iterator its;
  _foreach(its,s){
    term = its->first;
    sum1 = term.oneel2fock(decoration,multiref);
    sum1 *= its->second;
    sum += sum1;
  }
  return sum;
}
TermSum Q2::ReplaceE0(const TermSum& s)
{
  bool multiref = (Input::iPars["prog"]["multiref"] > 0);
  Term term;
  TermSum sum,sum1;
  TermSum::const_iterator its;
  _foreach(its,s){
    term = its->first;
    sum1 = term.replaceE0(multiref);
    sum1 *= its->second;
    sum += sum1;
  }
  return sum;
}

TermSum Q2::SingletDM(TermSum s)
{
  Term term;
  uint iter = 0;
  for ( iter = 0; has_nonsingldm(s) && iter < 1000; ++iter ){
    xout << "has nonsingldm" << std::endl;
    TermSum sum,sum1;
    for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
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
bool Q2::has_nonsingldm(const TermSum& s)
{
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i)
    if (i->first.has_nonsingldm()) return true;
  return false;
}
TermSum Q2::GeneralIndices(TermSum s)
{
  uint iter = 0;
  for ( iter = 0; has_generalindices(s) && iter < 1000; ++iter ){
    TermSum sum;
    for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
      Term term = i->first;
      TermSum sum1 = term.removegeneralindices();
      sum1 *= i->second;
      sum += sum1;
    }
    s = sum;
  }
  if (iter == 1000) error("Could not remove all general indices in 1000 iterations!");
  return s;
}
bool Q2::has_generalindices(const TermSum& s)
{
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i)
    if (i->first.has_generalindices()) return true;
  return false;
}
TermSum Q2::ZeroTerms(const TermSum& s)
{
  TermSum sum;
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    if ( ! i->first.removeit()){
      sum += std::make_pair(i->first,i->second);
    }
  }
  return sum;
}

TermSum Q2::EqualTerms(const TermSum& s, double minfac)
{
  int eqway = Input::iPars["prog"]["eqway"];
  TermSum sum;
  BigArray<UniGraph> ugraphs;
  BigArray<Term> uterms, newterms;
  Term term,term1;
  TFactor prefac;
  bool added;
  for ( TermSum::const_iterator j=s.begin();j!=s.end(); ++j) {
    term=j->first;
    prefac=j->second*term.prefac();
    // remove prefactors in terms
    term.reset_prefac();
    if (eqway > 0) {
      uterms.push_back(term);
      Term & uterm = uterms.back();
      uterm.set_prefac(prefac);
//       bool print = ( uterm.mat().size() == 2);
//       if (print) xout << "Term: " << uterm << std::endl;
      UniGraph ug(uterm);
//       xout << ug << std::endl;
      ug.minimize();
      added=false;
      for ( uint igr = 0; igr < ugraphs.size(); ++igr ){
        if ( ug.is_equal(ugraphs[igr]) ) {
//           if (print) xout << "before add " << newterms[igr] << std::endl;
          newterms[igr] += ug.permutation(ugraphs[igr]);
//           if (print) xout << "after add " << newterms[igr] << std::endl;
          added=true;
          break;
        }
      }
      if (added) {
        // not needed anymore
        uterms.pop_back();
      } else {
        term = ug.gen_term();
//         if (print) xout << term << std::endl;
        newterms.push_back(term);
        ugraphs.push_back(ug);
      }
    } else {
      for ( const Matrix& mat: term.mat() ){
        if (mat.has_pmsym()) {
          xout << "Tensor " << mat << " has plus/minus symmetry" << std::endl;
          error("Use eqway>0 for tensors with plus/minus symmetry");
        }
      }
      added=false;
      for ( TermSum::iterator k=sum.begin();k!=sum.end(); ++k) {
        Permut perm;
        term1=k->first;
        if (term.equal(term1,perm)) {
          sum.erase(k);
          term1+=std::make_pair(perm,prefac);
//         std::cout<<"term old" << term1 <<std::endl;
          if ( !term1.term_is_0(minfac) ) sum+=term1;
          added=true;
          break;
        }
      }
      if (!added) {
        term+=std::make_pair(Permut(),prefac);
        if ( !term.term_is_0(minfac) ) sum+=term;
//       std::cout<<"term new" << term <<std::endl;
      }
    }
  }
  if (eqway > 0) {
    for ( Term& newterm: newterms ) {
      if ( !newterm.term_is_0(minfac) ) {
        // for sorting... Remove if too slow
        newterm.matrixkind();
        newterm.setmatconnections();
        sum += newterm;
      }
    }
  }
  return sum;
}
TermSum Q2::SmallTerms(const TermSum& s, double minfac)
{
  TermSum sum;
  Term term,term1;
  for ( TermSum::const_iterator j=s.begin(); j!=s.end(); ++j) {
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
TermSum Q2::VirtSpace(const TermSum& s)
{
  std::string virtsp = Input::sPars["prog"]["virtspace"];
  std::transform(virtsp.begin(), virtsp.end(), virtsp.begin(), toupper);
  if ( virtsp != "PAO" && virtsp != "PNO" ) {
    // no overlaps
    return s;
  }
  // add overlap matrices
  xout << "add overlaps" << std::endl;
  TermSum sum;
  for ( const auto& ts: s ) {
    Term term = ts.first;
    term.addoverlaps();
    sum += std::make_pair(term,ts.second);
  }
  return sum;
}

TermSum Q2::ResolvePermutaions(const TermSum& s)
{
  TermSum sum;
  Term term,term1;
  for ( TermSum::const_iterator j=s.begin(); j!=s.end(); ++j) {
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

TermSum Q2::postaction(const TermSum& s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  std::string divide = Input::sPars["act"]["divide"];
  TermSum sum = s;
  if ( divide != "" && divide != "1" ){
    Finput div(true);
    div.addline(divide);
    div.analyzeline();
    div.analyzeq();
    TermSum sum1;
    Sum<Permut,TFactor> divperm, divpermadd;
    sum1 = div.sumterms();
    Matrix OneMat;
    for ( TermSum::const_iterator it=sum1.begin();it!=sum1.end(); ++it) {
      if ( it->first.mat().size() != 1 && !(it->first.mat().front() == OneMat) )
        error("Not an One-Matrix in divide!");
      divpermadd.clear();
      divpermadd += it->first.perm();
      divpermadd *= it->first.prefac();
      divpermadd *= it->second;
      divperm += divpermadd;
    }
    sum1.clear();
    xout << "divide by permutations: " << divperm << std::endl;
    Term term;
    for ( TermSum::const_iterator i=sum.begin();i!=sum.end(); ++i) {
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

TermSum Q2::normalOrderPH(const TermSum& s)
{
  bool timing = ( Input::iPars["prog"]["cpu"] > 0 );
  std::clock_t c_start=0;
  TermSum sum,sum0;
  Term term;
  say("Normal ordering");
  if (timing) c_start = std::clock();
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    sum0 += term.normalOrderPH_fullyContractedOnly();
    sum0 *= i->second;
    sum += sum0;
    sum0=TermSum();
  }
  if (timing) _CPUtiming("",c_start,std::clock());
  return sum;
}
TermSum Q2::wick(const TermSum& s)
{
  int iwick = Input::iPars["prog"]["wick"];
  bool genwick = (iwick > 1);
  int noorder = Input::iPars["prog"]["noorder"];
  if (!genwick && noorder > 0 ) error("Cannot have non-ordered Hamiltonian with wick<2. Either set noorder=0 or wick=2");
  bool timing = ( Input::iPars["prog"]["cpu"] > 0 );
  std::clock_t c_start=0;
  TermSum sum,sum0;
  Term term;
  _xout3(s << std::endl);
  say("Wick's theorem");
  if (timing) c_start = std::clock();
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    term=i->first;
    sum0 += term.wickstheorem(genwick,noorder);
    sum0 *= i->second;
    sum += sum0;
    sum0=TermSum();
  }
  if (timing) _CPUtiming("",c_start,std::clock());
  return sum;
}

void Q2::printdiags(Output* pout, const TermSum& s)
{
  say("Diagrams...");
  *(pout->pout) << " Diagrams: " << std::endl;
    *(pout->pout) << std::endl;
  for ( TermSum::const_iterator it = s.begin(); it != s.end(); ++it){
    it->first.printdiag(pout);
  }
}

void Q2::printalgo(std::ofstream& out, const TermSum& s)
{
  const std::string& resultt = Input::sPars["syntax"]["result"];
  say("Algorithm...");

  Factorizer fact(s);

  out << "algorithm..." << std::endl;
  // external indices
  // .....
  for ( TermSum::const_iterator i=s.begin();i!=s.end(); ++i) {
    out << resultt << "[abij]" << " += ";
    // call term ...
    out << std::endl;
  }
}

