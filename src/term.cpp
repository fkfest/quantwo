#include "term.h"

Term::Term() : _prefac(1.0){_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd) :
    _opProd(opProd), _prefac(1.0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd) :
    _opProd(opProd), _kProd(kProd), _prefac(1.0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}
    
Term::Term(const Product< SQOp >& opProd, const Product< Kronecker >& kProd, 
           const Product< Matrices >& mat, const Product< Orbital >& sumindx, const Product< Orbital >& realsumindx, 
           const double& prefac, const std::vector< Product<long int> >& connections) :
    _opProd(opProd), _kProd(kProd), _mat(mat), _sumindx(sumindx), _realsumindx(realsumindx), 
    _prefac(prefac), _connections(connections) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term& Term::operator*=(const Oper& op)
{
  _opProd *= op.SQprod();
  _mat *= op.mat();
  _sumindx *= op.sumindx();
  _realsumindx *= op.realsumindx();
  _prefac *= op.prefac();
  return *this;
}
Term& Term::operator*=(const double& fac)
{
  _prefac*=fac;
  return *this;
}
Term& Term::operator+=(const Permut& perm)
{
  _perm+=perm;
  return *this;
}
Term& Term::operator+=(const std::pair< Permut, double >& p)
{
  _perm+=p;
  return *this;
}
void Term::addconnection(const Product< long int >& connections)
{
  for (unsigned long int i=0; i<connections.size(); i++)
  { // Test
    if (connections[i]==0) 
      error("Zero in connections-map","Term::addconnection");
    if ((unsigned long)abs(connections[i])>_mat.size())
      error("Value in connections-map is larger than number of matrices in term","Term::addconnection");
  }
  _connections.push_back(connections);
}
void Term::addsummation(const Orbital& orb, short int excl)
{
  Orbital orb1;
  std::string excls;
  for (short i=0; i<excl ; ++i)
  {  
    if (i>0) excls=num2str(i,std::dec);
    orb1=Orbital(orb.name()+excls,orb.spin());
    if (_realsumindx.find(orb1)>=0)
      say("Strange, real summation runs already over "+orb1.name());
    _realsumindx*=orb1;
  }
}
void Term::addmatrix(const Matrices& mat)
{
  _mat *= mat;
}
void Term::replacematrix(const Matrices& mat, long unsigned int ipos)
{
  if (ipos >= _mat.size())
    error("The position is outside of this term: "+ipos,"Term::replacematrix");
  _mat[ipos]=mat;
}


Product<SQOp> Term::opProd() const
{   return _opProd; }

Product<Kronecker>  Term::kProd() const
{   return _kProd;  }
double Term::prefac() const
{ return _prefac; }
void Term::reset_prefac()
{ 
  _prefac=1.0;
  _perm=Sum<Permut,double>();
}
Product<Matrices> Term::mat() const
{ return _mat; }
Product< Orbital > Term::sumindx() const
{ return _sumindx; }
Product< Orbital > Term::realsumindx() const
{ return _realsumindx; }
Product< Orbital > Term::extindx() const
{
  // generate Product of external-lines orbitals
  Product<Orbital> peo;
  long int ipos,ipos1;
  for (unsigned int i=0; i<_mat.size(); i++)
  {
    for (unsigned int j=0; j<_mat[i].orbitals().size(); j++)
    {
      ipos=_realsumindx.find(_mat[i].orbitals().at(j));
      ipos1=peo.find(_mat[i].orbitals().at(j));
      if (ipos<0 && ipos1<0) // new external line
        peo*=_mat[i].orbitals().at(j);
    }
  }
  peo.resort();
  return peo;
}
Sum< Permut, double > Term::perm() const
{ return _perm; }
std::vector< Product< long int > > Term::connections() const
{ return _connections; }

bool Term::term_is_0(double minfac) const
{ 
  bool 
    is0 = (_opProd.size()==0 && _kProd.size()==0 && _mat.size()==0) || std::abs(_prefac)<minfac,
    loop = !is0;
  for ( Sum<Permut,double>::const_iterator it = _perm.begin(); loop && it != _perm.end(); ++it )
    loop = is0 = std::abs(it->second) < minfac;
  return is0;
}
bool Term::term_is_valid()
{
  for (unsigned long int i=0; i<_realsumindx.size(); i++)
  {
    if(_sumindx.find(_realsumindx[i])<0)
    {
      std::cout << *this << std::endl;
      error("Problem with summation index! May be a summation over excitations in a term without this excitaions?");
    }
  }
  return true;
}

bool Term::operator < (Term const & t) const
{
    if ( _kProd<t._kProd ) return true;
    if ( t._kProd<_kProd ) return false;
    if ( _opProd<t._opProd ) return true;
    if ( _opProd>t._opProd ) return false;
    if ( _nloops<t._nloops ) return false;
    if ( t._nloops<_nloops ) return true;
    if ( _mat<t._mat ) return true;
    if ( _mat>t._mat ) return false;
    if ( _sumindx<t._sumindx ) return true;
    if ( _sumindx>t._sumindx ) return false;
    if ( _perm<t._perm) return false;
    if ( t._perm<_perm) return true;
    if (_connections < t._connections) return true;
    if (t._connections < _connections) return false;
    return _prefac<t._prefac;
return true;
}
bool Term::equal(Term& t, Permut& perm) 
{
  
  if (_mat.size() != t._mat.size() ||
      _sumindx.size() != t._sumindx.size() ||
      _realsumindx.size() != t._realsumindx.size() ||
      _nintloops != t._nintloops || _nocc != t._nocc) return false;
  // generate Product of all orbitals and external-lines orbitals
  List<Orbital> peo(extindx()), peot(t.extindx());
  List<Orbital> po(peo),pot(peot); // start with external lines!
  po*=_realsumindx; // internal indices
  pot*=t._realsumindx; // internal indices
  if (po.size() != pot.size()) return false;

  //std::cout <<"term " << *this << "   " << t << std::endl;
  //std::cout <<"po,pot " << po << "   " << pot << std::endl;
  //std::cout <<"nocc,nintloops " << _nocc << "   " << t._nocc <<" and "<< _nintloops << "   " << t._nintloops << std::endl;
  Product<Matrices> mat, tmat;
  List<Orbital> por, port;
  Permut perm1;
  Orbital orb,orb1,orbt,orb1t;
  long int ipos=0,ipost=0;
  List<Orbital>::iterator it, jt;
  unsigned int ithis=0,ithist=0,i;
  bool equal=false,first,exter,extert,exter1,extert1;
  for (i=0; i<_mat.size(); i++) {
    _mat[i].reset_vertices();
    t._mat[i].reset_vertices();
  }
  while (po.size()>0) {
    orb=po.front();
    exter=(peo.find(orb)>=0); //external orbital
    i=0;
    equal=false;
    while (i<pot.size() && pot.size()>0 && !equal) {
      it = pot.begin();
      for ( unsigned int j = 0; j < i; ++j, ++it ){}
      orbt = *it;
      ++i;
      if (orb.type()!=orbt.type()) {
        equal=false;
        continue;
      }
      extert=(peot.find(orbt)>=0);
      if (exter!=extert) {// one is external orbital and the other not
        equal=false;
        continue;
      }
      orb1=orb;
      orb1t=orbt;
      port=por=List<Orbital>();
      mat=_mat;
      tmat=t._mat;
      perm1=perm;
      equal=true;
      first=true;
      exter1=extert1=false;
      if (exter) { //add orbitals to the product for removing
        por*=orb1;
        port*=orb1t;
        if (orb1 != orb1t) { // external orbitals not match -> add permutation
          perm1 *= Permut(orb1,orb1t);
        }
      }
      do {
        // find orbital which corresponds to the same electron
        for (unsigned int j=0; j<mat.size(); j++) {
          // dont search in the same matrix and in "external" matrices
          if((first || j!=ithis)&&(mat[j].type()!=Ops::Deexc0 && mat[j].type()!=Ops::Exc0)) {
            ipos=mat[j].orbitals().find(orb1);
            if (ipos>=0) {
              orb1=mat[j].orbel(ipos);
              ithis=j;
              break;
            }
          }
        }
        for (unsigned int j=0; j<tmat.size(); j++) {
          // dont search in the same matrix
          if((first || j!=ithist)&&(tmat[j].type()!=Ops::Deexc0 && tmat[j].type()!=Ops::Exc0)) {
            ipost=tmat[j].orbitals().find(orb1t);
            if (ipost>=0) {
              orb1t=tmat[j].orbel(ipost);
              ithist=j;
              break;
            }
          }
        }

        if (orb1.type() != orb1t.type() || !mat[ithis].vertices(ipos,tmat[ithist],ipost,ithis))
          equal=false;
        else { 
          exter1= (peo.find(orb1)>=0);
          extert1= (peot.find(orb1t)>=0);
          if (exter1!=extert1) { // if one of indices is external and the other not
            equal=false;
            break;
          }
          if (((orb1==orb) != (orb1t==orbt))&&!exter1) { // in one of matrices we have a loop and in the other not, and the index is not external!
              equal=false;
              break;
          }
          // add orbitals to the product for removing
          por*=orb1;
          port*=orb1t;
          first=false;
          if (exter1 && orb1 != orb1t) { // external orbitals not match -> add permutation
            perm1 *= Permut(orb1,orb1t);
          }
        }
      }while (orb1!=orb && equal && !exter1);
      if(equal) {
        // set _mat and t._mat to mat and tmat
        _mat=mat;
        t._mat=tmat;
        //set permutator to perm1
        perm=perm1;
        // remove por and port from product of orbitals (we dont need to handle this orbital again!)
        for ( jt = por.begin(); jt != por.end(); ++jt) {
          it = std::find(po.begin(),po.end(),*jt);
          assert( it != po.end() );
          po.erase(it);
        }
        for ( jt = port.begin(); jt != port.end(); ++jt) {
          it = std::find(pot.begin(),pot.end(),*jt);
          assert( it != pot.end() );
          pot.erase(it);
        }
      }
    }
    if (!equal && i>=pot.size())
      break;
  }
  return equal;
}

bool Term::equal_old(Term& t, Permut& perm) 
{
  
  if (_mat.size() != t._mat.size() ||
      _sumindx.size() != t._sumindx.size() ||
//      _realsumindx.size() != t._realsumindx.size() ||
      _nintloops != t._nintloops || _nocc != t._nocc) return false;
  // generate Product of all orbitals and external-lines orbitals
  Product<Orbital> peo(extindx()), peot(t.extindx());
  Product<Orbital> po(peo),pot(peot); // start with external lines!
  po*=_realsumindx; // internal indices
  pot*=t._realsumindx; // internal indices
  if (po.size() != pot.size()) return false;

  //std::cout <<"term " << *this << "   " << t << std::endl;
  //std::cout <<"po,pot " << po << "   " << pot << std::endl;
  //std::cout <<"nocc,nintloops " << _nocc << "   " << t._nocc <<" and "<< _nintloops << "   " << t._nintloops << std::endl;
  Product<Matrices> mat, tmat;
  Product<Orbital> por, port;
  Permut perm1;
  Orbital orb,orb1,orbt,orb1t;
  long int ipos=0,ipost=0;
  unsigned int ithis=0,ithist=0,i;
  bool equal=false,first,found,foundt,exter,extert,exter1,extert1;
  for (i=0; i<_mat.size(); i++) 
  {
    _mat[i].reset_vertices();
    t._mat[i].reset_vertices();
  }
  while (po.size()>0)
  {
    orb=po[0];
    exter=(peo.find(orb)>=0); //external orbital
    i=0;
    equal=false;
    while (i<pot.size() && pot.size()>0 && !equal)
    {
      orbt=pot[i];
      ++i;
      if (orb.type()!=orbt.type())
      {
        equal=false;
        continue;
      }
      extert=(peot.find(orbt)>=0);
      if (exter!=extert) 
      {// one is external orbital and the other not
        equal=false;
        continue;
      }
      orb1=orb;
      orb1t=orbt;
      port=por=Product<Orbital>();
      mat=_mat;
      tmat=t._mat;
      perm1=perm;
      equal=true;
      first=true;
      exter1=extert1=false;
      if (exter)
      { //add orbitals to the product for removing
        por*=orb1;
        port*=orb1t;
        if (orb1 != orb1t)
        { // external orbitals not match -> add permutation
          perm1 *= Permut(orb1,orb1t);
        }
      }
      do
      {
        // find orbital which corresponds to the same electron
        found=false;
        foundt=false;
        for (unsigned int j=0; j<mat.size()&&!(found&&foundt); j++)
        {// dont search in the same matrix and in "external" matrices
          if(!found&&(first || j!=ithis)&&(mat[j].type()!=Ops::Deexc0 && mat[j].type()!=Ops::Exc0)) 
          {
            ipos=mat[j].orbitals().find(orb1);
            if (ipos>=0)
            {
              found=true;
              orb1=mat[j].orbel(ipos);
              ithis=j;
            }
          }
          if(!foundt&&(first || j!=ithist)&&(tmat[j].type()!=Ops::Deexc0 && tmat[j].type()!=Ops::Exc0)) // dont search in the same matrix
          {
            ipost=tmat[j].orbitals().find(orb1t);
            if (ipost>=0)
            {
              foundt=true;
              orb1t=tmat[j].orbel(ipost);
              ithist=j;
            }
          }
        }

        if (orb1.type() != orb1t.type() || !mat[ithis].vertices(ipos,tmat[ithist],ipost,ithis))
          equal=false;
        else
        { 
          exter1= (peo.find(orb1)>=0);
          extert1= (peot.find(orb1t)>=0);
          if (exter1!=extert1)
          { // if one of indices is external and the other not
            equal=false;
            break;
          }
          if (((orb1==orb) != (orb1t==orbt))&&!exter1)
          { // in one of matrices we have a loop and in the other not, and the index is not external!
              equal=false;
              break;
          }
          // add orbitals to the product for removing
          por*=orb1;
          port*=orb1t;
          first=false;
          if (exter1 && orb1 != orb1t)
          { // external orbitals not match -> add permutation
            perm1 *= Permut(orb1,orb1t);
          }
        }
      }while (orb1!=orb && equal && !exter1);
      if(equal)
      {
        // set _mat and t._mat to mat and tmat
        _mat=mat;
        t._mat=tmat;
        //set permutator to perm1
        perm=perm1;
        // remove por and port from product of orbitals (we dont need to handle this orbital again!)
        for ( unsigned int j=0; j<por.size(); j++)
        {
          ipos=po.find(por[j]);
          if (ipos<0)
            error("Something strange with orbitals","Term::operator==");
          else
            po.erase(po.begin()+ipos);
          ipos=pot.find(port[j]);
          if (ipos<0)
            error("Something strange with orbitals","Term::operator==");
          else
            pot.erase(pot.begin()+ipos);
        }
      }
    }
    if (!equal && i>=pot.size())
      break;
  }
  return equal;
}


std::ostream & operator << (std::ostream & o, Term const & t)
{
 // o << "CN:" ;
 // for (unsigned long int i=0; i<t.connections().size(); i++)
 //   o << t.connections().at(i)<<":";
  std::streampos ipos0=o.tellp();
  bool printed = false;
  if ( std::abs(std::abs(t.prefac()) - 1.0) > MyOut::pcurout->small){
    o << t.prefac();
    MyOut::pcurout->lenbuf += o.tellp()-ipos0;
    printed = true;
  }
  
  if ( t.perm().size() > 1 || t.perm().begin()->second < 0 ){
    MyOut::pcurout->lenbuf++ ; // for "("
    o << "(" << t.perm() << ")";
    MyOut::pcurout->lenbuf++ ; // for ")"
  } else if ( std::abs(std::abs(t.perm().begin()->second) - 1.0) > MyOut::pcurout->small ||
              !t.perm().begin()->first.is1() ){ // don't print permutation 1
    if ( printed ) o << "*";
    o << t.perm();
  }
                                    
  if (t.realsumindx().size()>0) 
  {
    o <<"\\sum_{"<<t.realsumindx()<<"}";
    MyOut::pcurout->lenbuf += std::max(3,int(t.realsumindx().size())/MyOut::pcurout->wsi);
    MyOut::pcurout->hbufline = MyOut::pcurout->hsum;// set height of line to the height of \sum
  }
  o << t.mat(); //lenbuf will be handled in Matrices
  o << t.kProd(); //lenbuf will be handled in Kronecker
  if ( t.kProd().size()>0 && t.opProd().size()>0 )
    o << "*";
  o << t.opProd();
  MyOut::pcurout->lenbuf += t.opProd().size()*2;
return o;
}

Sum<Term, double> Term::normalOrder() const
{   return normalOrder(false);  }

Sum<Term, double> Term::normalOrder_fullyContractedOnly() const
{   return normalOrder(true);   }


Sum<Term, double> Term::normalOrder(bool fullyContractedOnly) const
{
Sum<Term, double>  sum;
    for ( unsigned int i=0 ; i+1<_opProd.size() ; ++i ) // iterate over Product<SQOp> but the last
    {
        // check if two consecutive operators need reordering
        if ( _opProd[i].gender()==SQOp::Annihilator  &&  _opProd[i+1].gender()==SQOp::Creator )
        {   // yes
        
            // handle 1st term
        Product<SQOp> p(_opProd); // copy Product<SQOp>
            std::swap(p[i], p[i+1]); // swap operators
            // check if we need downward recursion
            if ( !fullyContractedOnly || p[_opProd.size()-1].gender()==SQOp::Creator )
                sum -= Term(p, _kProd, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrder(fullyContractedOnly);

            // handle 2nd term 
        Product<SQOp> q(p);   // copy Product<SQOp>
          q.erase(q.begin()+i,q.begin()+i+2);
        // q.erase(std::find(q.begin(), q.end(), p[i]));        // remove
        // q.erase(std::find(q.begin(), q.end(), p[i+1]));      // operators
        
        Product<Kronecker>  d(_kProd);
            d *= Kronecker(p[i].orb(), p[i+1].orb());           // and add kronecker
            // check if we need downward recursion
            if ( !fullyContractedOnly || q.size()==0 || q[q.size()-1].gender()==SQOp::Creator )
                sum += Term(q, d, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrder(fullyContractedOnly);
            return sum;
        }
    }

    sum += *this;   // add current term to the intermediate sum
    return sum;
}

Sum<Term, double> Term::normalOrderPH() const
{   return normalOrderPH(false);  }

Sum<Term, double> Term::normalOrderPH_fullyContractedOnly() const
{   return normalOrderPH(true);   }


Sum<Term, double> Term::normalOrderPH(bool fullyContractedOnly) const
{
Sum<Term, double>  sum;
bool creat2, annih1;
    for ( unsigned int i=0 ; i+1<_opProd.size() ; ++i ) // iterate over Product<SQOp> but the last
    {
        // check if two consecutive operators need reordering
        creat2=(( _opProd[i].genderPH()==SQOp::Annihilator || _opProd[i].genderPH()==SQOp::Gen )
              &&  _opProd[i+1].genderPH()==SQOp::Creator );
        annih1=(_opProd[i].genderPH()==SQOp::Annihilator 
              && ( _opProd[i+1].genderPH()==SQOp::Creator || _opProd[i+1].genderPH()==SQOp::Gen));
        if ( creat2 || annih1 )
        {   // yes
            // handle 1st term
          Product<SQOp> p(_opProd); // copy Product<SQOp>
          std::swap(p[i], p[i+1]); // swap operators
          // check if we need downward recursion
          if ( !fullyContractedOnly || (p[_opProd.size()-1].genderPH()==SQOp::Creator)
                || (p[0].genderPH()==SQOp::Annihilator))
            sum -= Term(p, _kProd, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrderPH(fullyContractedOnly); 
          
          if ( _opProd[i].orb().type()==_opProd[i+1].orb().type() || 
               ( ( _opProd[i].orb().type()==Orbital::GenT || _opProd[i+1].orb().type()==Orbital::GenT ) 
                 && _opProd[i].gender()!=_opProd[i+1].gender() ))
          {
            // handle 2nd term 
            Product<SQOp> q(p);   // copy Product<SQOp>
            q.erase(q.begin()+i,q.begin()+i+2);
        
            Product<Kronecker>  d(_kProd);
            d *= Kronecker(p[i].orb(), p[i+1].orb());           // and add kronecker
            // check if we need downward recursion
            if ( !fullyContractedOnly || q.size()==0 || q[q.size()-1].genderPH()==SQOp::Creator 
                || q[0].genderPH()==SQOp::Annihilator)
              sum += Term(q, d, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrderPH(fullyContractedOnly);
          }
          return sum;
        }
    }
    if ( !fullyContractedOnly || _opProd.size()==0 ) 
      sum += *this;   // add current term to the intermediate sum
    return sum;
}
Sum< Term, double > Term::wickstheorem() const
{
  // generate "matrix" of indices to SQops
  std::vector< std::vector< long int > > opers;
  std::vector< long int > opermat;
  unsigned int m=0;
  for (unsigned int i=0; i<_opProd.size(); i++)
  {
    if (m==_mat.size())
    { // all SQops, which are not in Matrices have to be added as individual vectors
      opermat.push_back(i);
      opers.push_back(opermat);
      opermat=std::vector< long int >();
      m=0; // search from begin
    }   
    else if (_mat[m].orbitals().find(_opProd[i].orb())>=0)
      opermat.push_back(i);
    else 
    {
      m++;
      i--;
      if (opermat.size()>0) opers.push_back(opermat);
      opermat=std::vector< long int >();
    }
  }
  if (opermat.size()>0) opers.push_back(opermat);
  opermat=std::vector< long int >();
  
//  for (unsigned int i=0; i<opers.size(); i++)
//  {
//    for (unsigned int j=0; j<opers[i].size(); j++)
//      std::cout << opers[i][j] << " " ;
//    std::cout << std::endl;
//  }
  return wick(opers,opermat);
}

Sum<Term, double> Term::wick(std::vector< std::vector< long int > >& opers, std::vector< long int >& krons) const
{
  Sum<Term, double>  sum;
  if (opers.size()==0)
  { // no SQoperators left
    Product<SQOp> p;
    // generate Kroneckers
    Product<Kronecker> d;
    for (unsigned int i=0; i<krons.size(); i+=2)
      d*=Kronecker(_opProd[krons[i]].orb(),_opProd[krons[i+1]].orb());
    sum += Term(p,d,_mat, _sumindx, _realsumindx, _prefac, _connections);
    return sum;
  }
  unsigned long int istart,sign;
  long int curr=opers[0][0];
  SQOp::Gender gencurr=_opProd[curr].gender();
  Orbital::Type orbtypecurr=_opProd[curr].orb().type();
  // remove first SQop-index
  if (opers[0].size()<2)
  {
    opers.erase(opers.begin());
    istart=0;
  }
  else
  {
    opers[0].erase(opers[0].begin());
    istart=1;
  }
  // add first index to krons1
  krons.push_back(curr);
  sign=0;
  for ( unsigned int i=0 ; i<opers.size() ; ++i ) // iterate over all Operators 
  {
    if (i<istart) //but the first
    { //count SQops for sign
      sign+=opers[i].size();
      continue;
    }
    for ( unsigned int j=0 ; j<opers[i].size() ; ++j ) // iterate over all SQop in Operator i
    {
      // check if the first operator and operator i would yield a Kronecker
      if (_opProd[opers[i][j]].gender()!=gencurr && 
          (_opProd[opers[i][j]].orb().type()==orbtypecurr||
           orbtypecurr==Orbital::GenT || _opProd[opers[i][j]].orb().type()==Orbital::GenT))
      {
        // copy opers and krons
        std::vector< std::vector< long int > > opers1(opers);
        std::vector< long int > krons1(krons);
        // remove SQop-index
        if (opers1[i].size()==1)
          opers1.erase(opers1.begin()+i);
        else
          opers1[i].erase(opers1[i].begin()+j);
        // add index to "Kronecker"
        krons1.push_back(opers[i][j]);
        //call wick recursivly
        if ((sign+j)%2==0)
          sum+=wick(opers1,krons1);
        else
          sum-=wick(opers1,krons1);
      }
    }
    sign+=opers[i].size();
  }
  return sum;
}
void Term::setmatconnections()
{
  for (unsigned int i=0; i<_mat.size(); i++)
    for (unsigned int j=0; j<_mat[i].orbitals().size(); j++)
      for(unsigned int k=i+1; k<_mat.size();k++)
        if (_mat[k].orbitals().find(_mat[i].orbitals().at(j))>=0)
        {
          _mat[i].add_connect(k+1);
          _mat[k].add_connect(i+1);
        }
}
void Term::reduceTerm()
{
  int ipos1,ipos2,ikpr;
  Product<Kronecker> kpr(_kProd); // copy _kProd
  ikpr=0;
  for ( unsigned int i=0 ; i<_kProd.size() ; ++i ) // iterate over Product<Kronecker>
  {
    if (_kProd[i].orb1().type()!=_kProd[i].orb2().type() 
        && _kProd[i].orb2().type()!=Orbital::GenT && _kProd[i].orb1().type()!=Orbital::GenT )
    {
      _prefac=0.0; // Kronecker between two orbitals of different type
      return; 
    }
    // search for orbitals in (real) summations
    ipos1=_realsumindx.find(_kProd[i].orb1());
    ipos2=_realsumindx.find(_kProd[i].orb2());
    if ( ipos2>=0 ) // found orb2
    {
      _realsumindx.erase(_realsumindx.begin()+ipos2); // delete summation over orb2
      ipos2=_sumindx.find(_kProd[i].orb2());
      if ( ipos2<0)
        error("Strange, orbital not found in _sumindx","Term::reduceTerm");
      _sumindx.erase(_sumindx.begin()+ipos2);
      Q2::replace(_opProd,_kProd[i].orb2(),_kProd[i].orb1());
      Q2::replace(_mat,_kProd[i].orb2(),_kProd[i].orb1());
      kpr.erase(kpr.begin()+ikpr); // delete the Kronecker
      Q2::replace(kpr,_kProd[i].orb2(),_kProd[i].orb1());
      ikpr--;
    }
    else if ( ipos1>=0 ) // found orb1
    {
      _realsumindx.erase(_realsumindx.begin()+ipos1); // delete summation over orb2
      ipos1=_sumindx.find(_kProd[i].orb1());
      if ( ipos1<0)
        error("Strange, orbital not found in _sumindx","Term::reduceTerm");
      _sumindx.erase(_sumindx.begin()+ipos1);
      Q2::replace(_opProd,_kProd[i].orb1(),_kProd[i].orb2());
      Q2::replace(_mat,_kProd[i].orb1(),_kProd[i].orb2());
      kpr.erase(kpr.begin()+ikpr); // delete the Kronecker
      Q2::replace(kpr,_kProd[i].orb1(),_kProd[i].orb2());
      ikpr--;
    }
    ikpr++;
  }
  _kProd=kpr;
}
void Term::matrixkind()
{
  short exccl,intlines=0,intvirt=0;
  int ipos;
  for (unsigned int i=0; i<_mat.size(); i++)
  {
    // excitation class of operator (= #electrons = #orbitals/2)
    exccl=_mat[i].orbitals().size()/2;
    for (unsigned int j=0; j<_mat[i].orbitals().size(); j++)
    {
      // internal lines (have to be sumed up - search in _sumindx)
      ipos=_sumindx.find(_mat[i].orbitals().at(j));
      if (ipos>=0) 
      {
        ++intlines;
        if (_mat[i].orbitals().at(j).type()==Orbital::Virt)
          ++intvirt; // internal particle line
      }
    }
    _mat[i].setkind(exccl,intlines,intvirt);
  }
}
bool Term::expandintegral(bool firstpart)
{
  for (unsigned int i=0; i<_mat.size(); i++)
  {
    if (_mat[i].expandantisym(firstpart)) return true;
  }
  return false;
}
bool Term::antisymmetrized()
{
  for (unsigned int i=0; i<_mat.size(); i++)
  {
    if( _mat[i].antisymform()) return true;
  }
  return false;
}
Sum< Term, double > Term::expand_antisym(bool spinintegr)
{
  Sum< Term, double > sum;
  bool expand;
  Term term(*this); // copy term
  if (term.expandintegral(true))
  {
    // handle (PQ|RS) part
    sum+=term.expand_antisym(spinintegr);
    //handle (PS|RQ) part
    term=*this;
    expand=term.expandintegral(false);
    if(!expand)
      error("strange second expansion!","Term::expand_antisym");
    sum-=term.expand_antisym(spinintegr);
  }
  else
  {
    term.spinintegration(spinintegr);
    sum+=term;
  }
  return sum;
}

void Term::spinintegration(bool notfake)
{
  // generate Product of all orbitals and external-lines orbitals
  List<Orbital> peo(extindx());
  List<Orbital> po(peo); // start with external lines!
  po*=_realsumindx; // internal indices
  _nocc=_nintloops=_nloops=0;
  
  Orbital orb,orb1;
  Matrices::Spinsym spinsym=Matrices::Singlet;
  unsigned int ithis=0;
  long int ipos;
  List<Orbital>::iterator it;
  bool first, samespinsym,internalloop;
  while (po.size()>0)
  {
    orb=po.front();
    orb1=orb;
    first=true;
    samespinsym=true;
    internalloop=true;
    do {
      // remove orb1 from product of orbitals (we dont need to handle this orbital again!)
      it = std::find(po.begin(),po.end(),orb1);
      assert( it != po.end() );
      po.erase(it);
      // count number of occupied orbitals (for comparison)
      if (orb1.type()== Orbital::Occ) ++_nocc;
      // find orbital which corresponds to the same electron
      for (unsigned int j=0; j<_mat.size(); j++) {
        if(!first && j==ithis) continue; // dont search in the same matrix
        ipos=_mat[j].orbitals().find(orb1);
        if (ipos>=0) {
          orb1=_mat[j].orbel(orb1);
          ithis=j;
          if (first) {
            first=false;
            spinsym=_mat[j].spinsym(ipos);
          } else 
            samespinsym = samespinsym&&(spinsym==_mat[j].spinsym(ipos));
          break;
        }
      }
      internalloop=internalloop && (peo.find(orb1)<0);
    }while (orb1!=orb && _sumindx.find(orb1)>=0);
    if (samespinsym)
    {
      if (notfake) _prefac*=2.0;
      // count number of loops
      ++_nloops;
      // count number of internal loops
      if (internalloop) ++_nintloops;
    }
    else
      _prefac=0.0;
  }
  if (notfake)
  {// set no spin
    for (unsigned int i=0; i<_mat.size(); i++)
    {
      _mat[i].set_no_spin();
      if (InSet(_mat[i].type(),Ops::Exc, Ops::Deexc, Ops::Exc0, Ops::Deexc0, Ops::Interm ))
        for (unsigned int j=0; j<_mat[i].orbitals().size()/2; j++)
          _prefac *= double(j+1); // the symmetry of closed shell cluster operators is lower 
    }
    for (unsigned int i=0; i< _sumindx.size(); i++)
      _sumindx[i].setspin(Orbital::No);
    for (unsigned int i=0; i< _realsumindx.size(); i++)
      _realsumindx[i].setspin(Orbital::No);
  }
}
bool Term::properconnect() const
{
  long int imat,ipos;
  unsigned int i,j,k;
  Product<long int> notfound,found;
  for (i=0; i<_connections.size(); i++)
  { 
    notfound=Product<long int> ();
    for (j=0; j<_connections[i].size(); j++)
      notfound*=abs(_connections[i][j]);
    j=0;
    found=Product<long int> ();
    found*=notfound[0];
    notfound.erase(notfound.begin());
    while (notfound.size()>0 && j<found.size())
    {
      imat=found[j]-1;
      for (k=0; k<_mat[imat].connected2().size(); k++)
      {
        ipos=notfound.find(_mat[imat].connected2().at(k));
        if (ipos>=0)
        {
          found*=_mat[imat].connected2().at(k);
          notfound.erase(notfound.begin()+ipos);
        }
      }
      ++j;
    }
    if (_connections[i][0]>0)
    {// have to be connected
      if (notfound.size()>0) return false;
    }
    else
    {// have to be disconnected
      if (notfound.size()==0) return false;
    }
  }
  return true;
}

Orbital Term::freeorbname(Orbital::Type type)
{
  Orbital::Spin spin=Orbital::GenS;
  const std::string * ip_orbs;
  std::string lastorb=_lastorb[type].name();
  unsigned long int indx;
  if (type==Orbital::Occ)
    ip_orbs=&Orbital::occ_def;
  else if (type==Orbital::Virt)
    ip_orbs=&Orbital::virt_def;
  else
    ip_orbs=&Orbital::gen_def;
  if (lastorb.size()==0) 
  {
    lastorb=ip_orbs->at(0);
    _lastorb[type]=Orbital(lastorb,type,spin);
    return _lastorb[type];
  }
  for ( int i=lastorb.size()-1; i>=0; --i)
  {
    indx=ip_orbs->find(lastorb[i]);
    if(indx==std::string::npos)
      error("Something wrong with orbitals","Term::freeorbname");
    else if (indx<ip_orbs->size()-1) // not the last possible orbital index
    {
      lastorb[i]=ip_orbs->at(indx+1);
      _lastorb[type]=Orbital(lastorb,type,spin);
      return _lastorb[type];
    }
    lastorb[i]=ip_orbs->at(0);//set to first letter
  }
  //add one orbital index more
  lastorb+=ip_orbs->at(0);
  _lastorb[type]=Orbital(lastorb,type,spin);
  return _lastorb[type];
}
Orbital Term::getfreeorbname(void* Obj, Orbital::Type type)
{
  // explicitly cast to a pointer to Term
  Term * myself = (Term *) Obj;
  // call member
  return myself->freeorbname(type);
}
void Term::set_lastorb(Orbital orb, bool onlylarger)
{ 
  if ( onlylarger && orb.comp_letname(_lastorb[orb.type()]) <= 0 ) return; // do nothing
  _lastorb[orb.type()]=orb; 
  
}


Sum<Term,double> Q2::reduceSum(Sum<Term,double> s)
{
  double minfac = Input::fPars["prog"]["minfac"];
  Sum<Term,double> sum,sum1;
  Term term,term1;
  bool added;
  double prefac;
  say("Reduce sum of terms");
  say("Kroneckers, Connections, and Spin-integration...");
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  for ( Sum<Term,double>::const_iterator i=s.begin();i!=s.end(); ++i)
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
    if (!term.properconnect())
      continue;
    // expand antisymmetrized integrals and do spin integration
    sum=term.expand_antisym(spinintegr);
//    sum=term.expand_antisym(false);
    sum*=i->second;
    sum1+=sum;
  }
  say("Equal terms and permutations...");
  //std::cout << "SUM:" << sum1 << std::endl;
  sum=Sum<Term,double>();
  // sum up all equal terms
  for ( Sum<Term,double>::const_iterator j=sum1.begin();j!=sum1.end(); ++j) {
    term=j->first;
    prefac=j->second*term.prefac();
    // remove prefactors in terms
    term.reset_prefac();
    added=false;
    for ( Sum<Term,double>::iterator k=sum.begin();k!=sum.end(); ++k) {
      Permut perm;
      term1=k->first;
      if (term.equal(term1,perm)) {
        sum.erase(k);
        term1+=std::make_pair<Permut,double>(perm,prefac);
        //std::cout<<"term old" << term1 <<std::endl;
        if ( !term1.term_is_0(minfac) ) sum+=term1;
        added=true;
        break;
      }
    }
    if (!added) {
      term+=std::make_pair<Permut,double>(Permut(),prefac);
      if ( !term.term_is_0(minfac) ) sum+=term;
      //std::cout<<"term new" << term <<std::endl;
    }
  } 
  say("Remove terms with small prefactors...");
  // now remove everything with small prefactor
  sum1 = Sum<Term,double>();
  for ( Sum<Term,double>::const_iterator j=sum.begin(); j!=sum.end(); ++j) {
    if ( std::abs(j->second) < minfac ) continue;
    term1 = term = j->first;
    if ( std::abs(term.prefac()) < minfac ) continue;
    term1.reset_prefac();
    const Sum<Permut,double>& perms = term.perm();
    for ( Sum<Permut,double>::const_iterator it = perms.begin(); it != perms.end(); ++it ){
      if ( std::abs(it->second) < minfac ) continue;
      term1 += *it; //std::make_pair<Permut,double>(it->first,it->second);
    }
    sum1 += term1;
  }
  
  return sum1;
}
Sum< Term, double > Q2::normalOrderPH(Sum< Term, double > s)
{
  Sum<Term,double> sum,sum0;
  Term term;
  say("Normal ordering");
  for ( Sum<Term,double>::const_iterator i=s.begin();i!=s.end(); ++i)
  {
    term=i->first;
    sum0 += term.normalOrderPH_fullyContractedOnly();
    sum0 *= i->second;
    sum += sum0;
    sum0=Sum<Term,double>();
  }
  return sum;
}
Sum< Term, double > Q2::wick(Sum< Term, double > s)
{
  Sum<Term,double> sum,sum0;
  Term term;
  say("Wick's theorem");
  for ( Sum<Term,double>::const_iterator i=s.begin();i!=s.end(); ++i)
  {
    term=i->first;
    sum0 += term.wickstheorem();
    sum0 *= i->second;
    sum += sum0;
    sum0=Sum<Term,double>();
  }
  return sum;
}

template <class T>
void Q2::replace(Product< T > &p, Orbital orb1, Orbital orb2)
{
  for ( unsigned int i=0; i<p.size(); ++i )
  {
    replace(p[i],orb1, orb2);
  }
}
void Q2::replace(SQOp &op, Orbital orb1, Orbital orb2)
{ op.replace(orb1,orb2); }
void Q2::replace(Orbital &orb, Orbital orb1, Orbital orb2)
{ if (orb == orb1) orb=orb2; }
void Q2::replace(Matrices &mat, Orbital orb1, Orbital orb2)
{ mat.replace(orb1,orb2); }
void Q2::replace(Kronecker &kron, Orbital orb1, Orbital orb2)
{ kron.replace(orb1,orb2); }

