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
    if (!term.properconnect()) continue;
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
  sum.clear();
  // sum up all equal terms
  for ( Sum<Term,TFactor>::const_iterator j=sum1.begin();j!=sum1.end(); ++j) {
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
  say("Remove terms with small prefactors...");
  // now remove everything with small prefactor
  sum1.clear();
  for ( Sum<Term,TFactor>::const_iterator j=sum.begin(); j!=sum.end(); ++j) {
    if ( _todouble(_abs(j->second)) < minfac ) continue;
    term1 = term = j->first;
    if ( _todouble(_abs(term.prefac())) < minfac ) continue;
    term1.reset_prefac();
    const Sum<Permut,TFactor>& perms = term.perm();
    for ( Sum<Permut,TFactor>::const_iterator it = perms.begin(); it != perms.end(); ++it ){
      if ( _todouble(_abs(it->second)) < minfac ) continue;
      term1 += *it; 
    }
    sum1 += term1;
  }
  
  return sum1;
}
Sum< Term, TFactor > Q2::postaction(Sum< Term, TFactor > s)
{
  std::string divide = Input::sPars["act"]["divide"];
  Sum<Term,TFactor> sum = s;
  if (divide != ""){
    Sum<Term,TFactor> sum1;
    Sum<Permut,TFactor> divperm;
    divperm += Permut();
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
