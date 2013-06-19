#include "work.h"

Sum< Term, TFactor > Q2::reduceSum(Sum< Term, TFactor > s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  bool brill = ( Input::iPars["prog"]["brill"] > 0 );
  bool quan3 = ( Input::iPars["prog"]["quan3"] > 0 );
  Sum<Term,TFactor> sum,sum1;
  Term term,term1;
  bool added;
  TFactor prefac;
  say("Reduce sum of terms");
  say("Kroneckers, Connections and Antisymmetry...");
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i)
  {
    term=i->first;
    // remove Kroneckers
    term.reduceTerm();
    if (term.kProd().size()>0)
      error("Still some Kroneckers left. We can not handle it yet!","Q2::reduceSum");
    // generate Kallay's "triplets of integers"
    term.matrixkind();
    // set connections "map" of matrices
    term.setmatconnections();
    // is the term properly connected?
    if (!term.properconnect()) {
      continue;
    }
    // expand antisymmetrized integrals
    sum1 = term.expand_antisym();
    sum1 *= i->second;
    sum += sum1;
  }

  _xout3(sum << std::endl);
  if (quan3) {
    say("count electrons (a posteriori)...");
    sum1.clear();
    Sum< Term, TFactor > sum2(sum);
    for ( Sum<Term,TFactor>::const_iterator i=sum2.begin();i!=sum2.end(); ++i) {
      term=i->first;
      term.deleteNoneMats();
      term.setmatconnections();
      
      prefac=i->second*term.prefac();
      // remove prefactors in terms
      term.reset_prefac();
      added=false;
      for ( Sum<Term,TFactor>::iterator k=sum1.begin();k!=sum1.end(); ++k) {
        Permut perm;
        term1=k->first;
        if (term.equal(term1,perm)) {
          sum1.erase(k);
          term1+=std::make_pair<Permut,TFactor>(perm,prefac);
          //std::cout<<"term old" << term1 <<std::endl;
          if ( !term1.term_is_0(minfac) ) sum1+=term1;
          added=true;
          break;
        }
      }
      if (!added) {
        term+=std::make_pair<Permut,TFactor>(Permut(),prefac);
        if ( !term.term_is_0(minfac) ) sum1+=term;
        //std::cout<<"term new" << term <<std::endl;
      }
      
      
//      Termel termel(term);
//      sum1 += std::make_pair(term,i->second);
    }
//    _xout3(sum1 << std::endl);
    return sum1;
  }
    
  say("Diagrams and Spin-integration...");
  sum1.clear();
  for ( Sum<Term,TFactor>::const_iterator i=sum.begin();i!=sum.end(); ++i) {
    term=i->first;
    // use Brilloin condition 
    if ( brill && term.brilloin() ) continue;
    // remove "None" matrices
    term.deleteNoneMats();
    term.setmatconnections();
    term.spinintegration(spinintegr);
    sum1 += std::make_pair(term,i->second);
  }

  _xout3(sum1 << std::endl);

  say("Equal terms and permutations...");
  //std::cout << "SUM:" << sum1 << std::endl;
  sum = EqualTerms(sum1,minfac);
  say("Remove terms with small prefactors...");
  // now remove everything with small prefactor
  sum1 = SmallTerms(sum,minfac);
  
  return sum1;
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
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i)
  {
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
  Sum<Term,TFactor> sum,sum0;
  Term term;
  say("Wick's theorem");
  for ( Sum<Term,TFactor>::const_iterator i=s.begin();i!=s.end(); ++i)
  {
    term=i->first;
    sum0 += term.wickstheorem();
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
    it->first.printdiag(pout,it->second);
  }
}
