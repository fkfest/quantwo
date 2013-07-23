#include "term.h"

Term::Term() : _prefac(1){_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd) :
    _opProd(opProd), _prefac(1) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd) :
    _opProd(opProd), _kProd(kProd), _prefac(1) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}
    
Term::Term(const Product< SQOp >& opProd, const Product< Kronecker >& kProd, 
           const Product< Matrices >& mat, const TOrbSet& sumindx, const TOrbSet& realsumindx, 
           const TFactor& prefac, const std::vector< Product<long int> >& connections) :
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
Term& Term::operator*=(const TFactor& fac)
{
  _prefac*=fac;
  return *this;
}
Term& Term::operator*=(const Permut& perm)
{
  if ( _perm.size() == 0 )
    _perm += perm;
  else {
    Sum<Permut,TFactor> perms;
    Permut perm1;
    for ( Sum<Permut,TFactor>::iterator ip = _perm.begin(); ip != _perm.end(); ++ip ){
      perm1 = ip->first;
      perm1 *= perm;
      perms += std::pair< Permut, TFactor >(perm1,ip->second);
    }
    _perm = perms;
  }
  return *this;
}
Term& Term::operator+=(const Permut& perm)
{
  _perm+=perm;
  return *this;
}
Term& Term::operator+=(const std::pair< Permut, TFactor >& p)
{
  _perm+=p;
  return *this;
}
void Term::addconnection(const Product< long int >& connections)
{
  for (lui i=0; i<connections.size(); i++)
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
  std::pair<TOrbSet::iterator,bool> ret;
  for (short i=0; i<excl ; ++i) {  
    if (i>0) excls=num2str(i,std::dec);
    orb1=Orbital(orb.name()+excls,orb.spin());
    ret = _realsumindx.insert(orb1);
    if (!ret.second)
      say("Strange, real summation runs already over "+orb1.name());
  }
}
void Term::addmatrix(const Matrices& mat)
{
  _mat *= mat;
}
void Term::replacematrix(const Matrices& mat, lui ipos)
{
  if (ipos >= _mat.size())
    error("The position is outside of this term: "+ipos,"Term::replacematrix");
  _mat[ipos]=mat;
}


Product<SQOp> Term::opProd() const
{   return _opProd; }

Product<Kronecker>  Term::kProd() const
{   return _kProd;  }
TFactor Term::prefac() const
{ return _prefac; }
void Term::reset_prefac()
{ 
  _prefac=1;
  _perm=Sum<Permut,TFactor>();
}
void Term::permute(const Permut& p)
{
  assert(_opProd.size() == 0 && _kProd.size() == 0 );
  for ( Product<Matrices>::iterator it = _mat.begin(); it != _mat.end(); ++it ){
    it->permute(p);
  }
}

Product<Matrices> Term::mat() const
{ return _mat; }
TOrbSet Term::sumindx() const
{ return _sumindx; }
TOrbSet Term::realsumindx() const
{ return _realsumindx; }
TOrbSet Term::extindx() const
{
  // generate Product of external-lines orbitals
  TOrbSet peo;
  for ( lui i = 0; i < _mat.size(); ++i ) {
    for ( lui j = 0; j < _mat[i].orbitals().size(); ++j) {
      const Orbital & orb = _mat[i].orbitals()[j];
      if ( _realsumindx.count(orb) == 0 )
        peo.insert(orb);
    }
  }
  return peo;
}
Sum< Permut, TFactor > Term::perm() const
{ return _perm; }
std::vector< Product< long int > > Term::connections() const
{ return _connections; }

bool Term::term_is_0(double minfac) const
{ 
  bool 
    is0 = (_opProd.size()==0 && _kProd.size()==0 && _mat.size()==0) || _todouble(_abs(_prefac))<minfac,
    loop = !is0;
  for ( Sum<Permut,TFactor>::const_iterator it = _perm.begin(); loop && it != _perm.end(); ++it )
    loop = is0 = _todouble(_abs(it->second)) < minfac;
  return is0;
}
bool Term::term_is_valid()
{
  for ( TOrbSet::const_iterator it = _realsumindx.begin(); it != _realsumindx.end(); ++it ) {
    if(_sumindx.count(*it) == 0) {
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
  TOrbSet peo(extindx()), peot(t.extindx());
  List<Orbital> po(peo),pot(peot); // start with external lines!
  po*=_realsumindx; // internal indices
  pot*=t._realsumindx; // internal indices
  if (po.size() != pot.size()) return false;

  Product<Matrices> mat, tmat;
  Product<Orbital> por, port;
  Permut perm1;
  Orbital orb,orb1,orbt,orb1t;
  long int ipos=0,ipost=0;
  List<Orbital>::iterator it;
  Product<Orbital>::iterator jt;
  unsigned int ithis=0,ithist=0,i;
  bool equal=false,exter,extert,exter1,extert1,loop,loopt;
  std::vector< unsigned int > ordmat, ordmatt, oordmat, oordmatt;
  for (i=0; i<_mat.size(); i++) {
    _mat[i].reset_vertices();
    t._mat[i].reset_vertices();
    if ( InSet(_mat[i].type(),Ops::Deexc0,Ops::Exc0) )
      ordmat.push_back(i);
    else
      oordmat.push_back(i);
    if ( InSet(t._mat[i].type(),Ops::Deexc0,Ops::Exc0) )
      ordmatt.push_back(i);
    else
      oordmatt.push_back(i);
  }
  // put external matrices first
  for (i=0; i<oordmat.size(); i++) 
    ordmat.push_back(oordmat[i]);
  for (i=0; i<oordmatt.size(); i++) 
    ordmatt.push_back(oordmatt[i]);
  loop = false;
  while (po.size()>0) {
    orb=po.front();
    exter=(peo.count(orb)); //external orbital
    it = pot.begin();
    i = 0;
    equal=false;
    while (it != pot.end() && pot.size()>0 && !equal) {
      orbt = *it;
      ++it;
      ++i;
      if (orb.type()!=orbt.type()) {
        equal=false;
        continue;
      }
      extert=(peot.count(orbt));
      if (exter!=extert) {// one is external orbital and the other not
        equal=false;
        continue;
      }
      orb1=orb;
      orb1t=orbt;
      port = por = Product<Orbital>();
      mat=_mat;
      tmat=t._mat;
      perm1=perm;
      equal=true;
      exter1=extert1=false;
      if (exter) { //add orbitals to the product for removing
        por*=orb1;
        port*=orb1t;
        if (orb1 != orb1t) { // external orbitals not match -> add permutation
          perm1 += Permut(orb1,orb1t);
        }
      }
      assert(ordmat.size() == mat.size());
      for (unsigned int jo=0; jo<ordmat.size(); jo++) {
        unsigned int j = ordmat[jo];
        ipos = mat[j].orbitals().find(orb1);
        if (ipos >= 0) {
          ithis=j;
          break;
        }
      }
      assert(ordmatt.size() == tmat.size());
      for (unsigned int jo=0; jo<ordmatt.size(); jo++) {
        unsigned int j = ordmatt[jo];
        ipost = tmat[j].orbitals().find(orb1t);
        if (ipost >= 0) {
          ithist=j;
          break;
        }
      }
      do {
        // find orbital which corresponds to the same electron
        const ConLine& cl = mat[ithis].conline(ipos);
        ithis = cl.imat;
        ipos = cl.idx;
        ipos = mat[ithis].iorbel(ipos);
        if ( ipos < 0 ) break;
        orb1 = mat[ithis].orbitals()[ipos];
       
        const ConLine& clt = tmat[ithist].conline(ipost);
        ithist = clt.imat;
        ipost = clt.idx;
        ipost = tmat[ithist].iorbel(ipost);
        if ( ipost < 0 ) break;
        orb1t = tmat[ithist].orbitals()[ipost];
        
        if (orb1.type() != orb1t.type() || !mat[ithis].vertices(ipos,tmat[ithist],ipost,ithis))
          equal=false;
        else { 
          exter1= (peo.count(orb1));
          extert1= (peot.count(orb1t));
          if (exter1!=extert1) { // if one of indices is external and the other not
            equal=false;
            break;
          }
          loop = (orb1==orb);
          loopt = (orb1t==orbt);
          if ((loop != loopt)&&!exter1) { // in one of matrices we have a loop and in the other not, and the index is not external!
              equal=false;
              break;
          }
          // add orbitals to the product for removing
          por*=orb1;
          port*=orb1t;
          if (exter1 && orb1 != orb1t) { // external orbitals not match -> add permutation
            perm1 += Permut(orb1,orb1t);
          }
        }
      } while (equal && !loop && !exter1);
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

std::ostream & operator << (std::ostream & o, Term const & t)
{
 // o << "CN:" ;
 // for (unsigned long int i=0; i<t.connections().size(); i++)
 //   o << t.connections().at(i)<<":";
  std::streampos ipos0=o.tellp();
  bool printed = false;
  if ( _todouble(_abs(_abs(t.prefac()) - 1)) > MyOut::pcurout->small){
    o << t.prefac();
    MyOut::pcurout->lenbuf += o.tellp()-ipos0;
    printed = true;
  }
  
  if ( t.perm().size() > 1 || t.perm().begin()->second < 0 ){
    MyOut::pcurout->lenbuf++ ; // for "("
    o << "(" << t.perm() << ")";
    MyOut::pcurout->lenbuf++ ; // for ")"
  } else if ( _todouble(_abs(_abs(t.perm().begin()->second) - 1)) > MyOut::pcurout->small ||
              !t.perm().begin()->first.is1() ){ // don't print permutation 1
    if ( printed ) o << "*";
    o << t.perm();
  }
                                    
  if (t.realsumindx().size()>0) 
  {
    o <<"\\" << Input::sPars["command"]["sum"] <<"_{"<<t.realsumindx()<<"}";
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

Sum< Term, TFactor > Term::normalOrder() const
{   return normalOrder(false);  }

Sum< Term, TFactor > Term::normalOrder_fullyContractedOnly() const
{   return normalOrder(true);   }


Sum< Term, TFactor > Term::normalOrder(bool fullyContractedOnly) const
{
    Sum<Term, TFactor>  sum;
    for ( unsigned int i=0 ; i+1<_opProd.size() ; ++i ) // iterate over Product<SQOp> but the last
    {
        // check if two consecutive operators need reordering
        if ( _opProd[i].gender()==SQOp::Annihilator  &&  _opProd[i+1].gender()==SQOp::Creator ) {   // yes
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

Sum< Term, TFactor > Term::normalOrderPH() const
{   return normalOrderPH(false);  }

Sum< Term, TFactor > Term::normalOrderPH_fullyContractedOnly() const
{   return normalOrderPH(true);   }


Sum< Term, TFactor > Term::normalOrderPH(bool fullyContractedOnly) const
{
  Sum<Term, TFactor>  sum;
  bool creat2, annih1;
  for ( unsigned int i=0 ; i+1<_opProd.size() ; ++i ) {// iterate over Product<SQOp> but the last
    // check if two consecutive operators need reordering
    creat2=(( _opProd[i].genderPH()==SQOp::Annihilator || _opProd[i].genderPH()==SQOp::Gen )
              &&  _opProd[i+1].genderPH()==SQOp::Creator );
    annih1=(_opProd[i].genderPH()==SQOp::Annihilator 
              && ( _opProd[i+1].genderPH()==SQOp::Creator || _opProd[i+1].genderPH()==SQOp::Gen));
    if ( creat2 || annih1 ) {   // yes
      // handle 1st term
      Product<SQOp> p(_opProd); // copy Product<SQOp>
      std::swap(p[i], p[i+1]); // swap operators
      // check if we need downward recursion
      if ( !fullyContractedOnly || (p[_opProd.size()-1].genderPH()==SQOp::Creator)
                || (p[0].genderPH()==SQOp::Annihilator))
        sum -= Term(p, _kProd, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrderPH(fullyContractedOnly); 
          
      if ( _opProd[i].orb().type()==_opProd[i+1].orb().type() || 
               ( ( _opProd[i].orb().type()==Orbital::GenT || _opProd[i+1].orb().type()==Orbital::GenT ) 
                 && _opProd[i].gender()!=_opProd[i+1].gender() )) {
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
Sum< Term, TFactor > Term::wickstheorem() const
{
  // generate "matrix" of indices to SQops
  std::vector< std::vector< long int > > opers;
  std::vector< long int > opermat;
  unsigned int m=0;
  int cran = 0;
  for (unsigned int i=0; i<_opProd.size(); i++) {
    // calculate #creators - #annihilators
    if ( _opProd[i].gender() == SQOp::Creator ){
      ++cran;
    } else if ( _opProd[i].gender() == SQOp::Annihilator ){
      --cran;
    } else
      assert(false);
  }
  if ( cran != 0 ) return Sum< Term, TFactor>();
  for (unsigned int i=0; i<_opProd.size(); i++) {
    if (m==_mat.size()) { // all SQops, which are not in Matrices have to be added as individual vectors
      opermat.push_back(i);
      opers.push_back(opermat);
      opermat=std::vector< long int >();
      m=0; // search from begin
    } else if (_mat[m].orbitals().find(_opProd[i].orb())>=0)
      opermat.push_back(i);
    else {
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

Sum< Term, TFactor > Term::wick(std::vector< std::vector< long int > >& opers, std::vector< long int >& krons) const
{
  Sum<Term, TFactor>  sum;
  if (opers.size()==0) { // no SQoperators left
    Product<SQOp> p;
    // generate Kroneckers
    Product<Kronecker> d;
    for (unsigned int i=0; i<krons.size(); i+=2)
      d*=Kronecker(_opProd[krons[i]].orb(),_opProd[krons[i+1]].orb());
    sum += Term(p,d,_mat, _sumindx, _realsumindx, _prefac, _connections);
    return sum;
  }
  lui istart,sign;
  long int curr=opers[0][0];
  SQOp::Gender gencurr=_opProd[curr].gender();
  Orbital::Type orbtypecurr=_opProd[curr].orb().type();
  // remove first SQop-index
  if (opers[0].size()<2) {
    opers.erase(opers.begin());
    istart=0;
  } else {
    opers[0].erase(opers[0].begin());
    istart=1;
  }
  // add first index to krons1
  krons.push_back(curr);
  sign=0;
  for ( unsigned int i=0 ; i<opers.size() ; ++i ) { // iterate over all Operators 
    if (i<istart) { //but the first
      //count SQops for sign
      sign+=opers[i].size();
      continue;
    }
    for ( unsigned int j=0 ; j<opers[i].size() ; ++j ) {// iterate over all SQop in Operator i
      // check if the first operator and operator i would yield a Kronecker
      if (_opProd[opers[i][j]].gender()!=gencurr && 
          (_opProd[opers[i][j]].orb().type()==orbtypecurr||
           orbtypecurr==Orbital::GenT || _opProd[opers[i][j]].orb().type()==Orbital::GenT)) {
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
  long int kj;
  for (lui i=0; i<_mat.size(); ++i)
    for (lui j=0; j<_mat[i].orbitals().size(); ++j){
      const Orbital& ijorb = _mat[i].orbitals()[j];
      for(lui k=i+1; k<_mat.size(); ++k)
        if ( (kj = _mat[k].orbitals().find(ijorb)) >= 0) {
          _mat[i].add_connect(k+1);
          _mat[k].add_connect(i+1);
          _mat[i].set_conline(j,k,kj);
          _mat[k].set_conline(kj,i,j);
        }
    }
}
void Term::reduceTerm()
{
  int ikpr;
  TOrbSet::iterator it1, it2;
  Product<Kronecker> kpr(_kProd); // copy _kProd
  ikpr=0;
  for ( unsigned int i=0 ; i<_kProd.size() ; ++i ) // iterate over Product<Kronecker>
  {
    if (_kProd[i].orb1().type()!=_kProd[i].orb2().type() 
        && _kProd[i].orb2().type()!=Orbital::GenT && _kProd[i].orb1().type()!=Orbital::GenT )
    {
      _prefac=0; // Kronecker between two orbitals of different type
      return; 
    }
    // search for orbitals in (real) summations
    it1 = _realsumindx.find(_kProd[i].orb1());
    it2 = _realsumindx.find(_kProd[i].orb2());
    if ( it2 != _realsumindx.end() ) { // found orb2
      _realsumindx.erase(it2); // delete summation over orb2
      it2 = _sumindx.find(_kProd[i].orb2());
      if ( it2 == _sumindx.end() )
        error("Strange, orbital not found in _sumindx","Term::reduceTerm");
      _sumindx.erase(it2);
      Q2::replace(_opProd,_kProd[i].orb2(),_kProd[i].orb1());
      Q2::replace(_mat,_kProd[i].orb2(),_kProd[i].orb1());
      kpr.erase(kpr.begin()+ikpr); // delete the Kronecker
      Q2::replace(kpr,_kProd[i].orb2(),_kProd[i].orb1());
      ikpr--;
    } else if ( it1 != _realsumindx.end() ) { // found orb1
      _realsumindx.erase(it1); // delete summation over orb2
      it1 = _sumindx.find(_kProd[i].orb1());
      if ( it1 == _sumindx.end() )
        error("Strange, orbital not found in _sumindx","Term::reduceTerm");
      _sumindx.erase(it1);
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
static bool matisnone(const Matrices& mat)
{ return (mat.type() == Ops::None); };
void Term::deleteNoneMats()
{ _mat.erase(std::remove_if(_mat.begin(),_mat.end(),matisnone),_mat.end()); }
bool Term::brilloin() const
{
  for ( Product<Matrices>::const_iterator it = _mat.begin(); it != _mat.end(); ++it )
    if ( it->type() == Ops::Fock && it->orbitals()[0].type() != it->orbitals()[1].type() )
      return true;
  return false;
}

void Term::matrixkind()
{
  short exccl,intlines=0,intvirt=0;
  for (unsigned int i=0; i<_mat.size(); i++)
  {
    // excitation class of operator (= #electrons = #orbitals/2)
    exccl=_mat[i].orbitals().size()/2;
    for (unsigned int j=0; j<_mat[i].orbitals().size(); j++) {
      // internal lines (have to be sumed up - search in _sumindx)
      if ( _sumindx.count(_mat[i].orbitals()[j])) {
        ++intlines;
        if (_mat[i].orbitals()[j].type()==Orbital::Virt)
          ++intvirt; // internal particle line
      }
    }
    _mat[i].setkind(exccl,intlines,intvirt);
  }
}
bool Term::expandintegral(bool firstpart)
{
  for (unsigned int i=0; i<_mat.size(); i++)
    if (_mat[i].expandantisym(firstpart)) return true;
  return false;
}
bool Term::antisymmetrized()
{
  for (unsigned int i=0; i<_mat.size(); i++)
    if( _mat[i].antisymform()) return true;
  return false;
}
Sum< Term, TFactor > Term::expand_antisym()
{
  Sum< Term, TFactor > sum;
  bool expand;
  Term term(*this); // copy term
  if (term.expandintegral(true)) {
    // handle (PQ|RS) part
    sum+=term.expand_antisym();
    //handle (PS|RQ) part
    term=*this;
    expand=term.expandintegral(false);
    if(!expand)
      error("strange second expansion!","Term::expand_antisym");
    sum-=term.expand_antisym();
  } else {
    sum+=term;
  }
  return sum;
}

void Term::spinintegration(bool notfake)
{
  // generate Product of all orbitals and external-lines orbitals
  TOrbSet peo(extindx()), peo1(peo);
  // internal indices
  TOrbSet po(_realsumindx); 
  _nocc=_nintloops=_nloops=0;
  
  bool nonconserve = true, already_done = false;
  Orbital orb,orb1;
  Matrices::Spinsym spinsym=Matrices::Singlet;
  lui ithis = 0;
  long int ipos, iorb = 0;
  TOrbSet::iterator it;
  bool samespinsym;
  while ( po.size() > 0 || peo.size() > 0 ) {
    // start with external lines!
    if ( peo.size() > 0 )
      orb = *peo.begin();
    else
      orb = *po.begin();
    // find the orbital
    for (unsigned int j=0; j<_mat.size(); j++) {
      ipos=_mat[j].orbitals().find(orb);
      if (ipos >= 0) {
        iorb = ipos;
        ithis=j;
        spinsym=_mat[j].spinsym(ipos);
        break;
      }
    } 
    orb1 = orb;
    samespinsym = true;
    do {
      // remove orb1 from product of orbitals (we dont need to handle this orbital again!)
      it = peo.find(orb1);
      if ( it != peo.end() )
        peo.erase(it);
      else if ( (it = po.find(orb1)) != po.end() )
        po.erase(it);
      else {
        assert( false );
      }
      // count number of occupied orbitals (for comparison)
      if (orb1.type()== Orbital::Occ) ++_nocc;
      // find orbital which corresponds to the same electron
      const ConLine& cl = _mat[ithis].conline(iorb);
      ithis = cl.imat;
      ipos = cl.idx;
      samespinsym = samespinsym && (spinsym == _mat[ithis].spinsym(ipos));
      iorb = _mat[ithis].iorbel(ipos);
      if ( iorb < 0 ) break;
      orb1 = _mat[ithis].orbitals()[iorb];
      if (nonconserve){
        // test whether still in sets
        already_done = ( peo.count(orb1) == 0 && po.count(orb1) == 0 );
      }
    } while (orb1!=orb && !already_done && _sumindx.count(orb1));
    if (samespinsym) {
      if (notfake) _prefac*=2;
      // count number of loops
      ++_nloops;
      // count number of internal loops
      if ( peo1.find(orb) == peo1.end() ) ++_nintloops;
    } else
      _prefac=0;
  }
  if (notfake) {// set no spin
    for (unsigned int i=0; i<_mat.size(); i++) {
      _mat[i].set_no_spin();
      if (InSet(_mat[i].type(),Ops::Exc, Ops::Deexc, Ops::Exc0, Ops::Deexc0, Ops::Interm ))
        for (unsigned int j=0; j<_mat[i].orbitals().size()/2; j++)
          _prefac *= j+1; // the symmetry of closed shell cluster operators is lower 
    }
    TOrbSet sumindx;
    for ( TOrbSet::iterator it = _sumindx.begin(); it != _sumindx.end(); ++it ){
      orb = *it;
      orb.setspin(Orbital::No);
      sumindx.insert(orb);
    }
    _sumindx = sumindx;
    sumindx.clear();
    for ( TOrbSet::iterator it =_realsumindx.begin(); it != _realsumindx.end(); ++it ){
      orb = *it;
      orb.setspin(Orbital::No);
      sumindx.insert(orb);
    }
    _realsumindx = sumindx;
  }
}

bool Term::properconnect() const
{
  long int imat;
  for (lui i=0; i<_connections.size(); i++) { 
    TCon2 notfound;
    Product<long int> found;
    for (lui j=0; j<_connections[i].size(); j++)
      notfound.insert(abs(_connections[i][j]));
    found.push_back(*notfound.begin());
    notfound.erase(notfound.begin());
    lui j = 0;
    while (notfound.size() > 0 && j < found.size()) {
      imat = found[j] - 1;
      const TCon2 & icon2 = _mat[imat].connected2();
      for ( TCon2::const_iterator kt = icon2.begin(); kt != icon2.end(); ++kt) {
        TCon2::iterator ipos = notfound.find(*kt);
        if (ipos != notfound.end()) {
          found.push_back(*kt);
          notfound.erase(ipos);
        }
      }
      ++j;
    }
    if (_connections[i][0] > 0) {// have to be connected
      if (notfound.size()>0) return false;
    } else {// have to be disconnected
      if (notfound.size()==0) return false;
    }
  }
  return true;
}
void Term::printdiag(Output* pout, TFactor fac) const
{
  TsPar & diag = Input::sPars["diag"];
  short verb = Input::iPars["diag"]["printnames"]; // print names of operators and indices...
  
  *(pout->pout) << diag["bdiag"] << std::endl;
  
// "numbers"
  Product<Matrices> nums;
  for ( Product<Matrices>::const_iterator it = _mat.begin(); it != _mat.end(); ++it)
    if ( it->type() == Ops::Number ){
      nums.push_back(*it);
    }
  bool printperm = ( _perm.size() > 1 || _perm.begin()->second < 0 || 
                    _todouble(_abs(_abs(_perm.begin()->second) - 1)) > MyOut::pcurout->small ||
                    !_perm.begin()->first.is1());
  if ( printperm || nums.size() > 0){
    *(pout->pout) << diag["text"] << "{0}{$";
    *(pout->pout) << nums;
    if ( printperm )
      *(pout->pout) << "(" << _perm << ")";
    *(pout->pout) << "$}" << std::endl;
  }
  
  // put tensors
  lui im = 0;
  for ( Product<Matrices>::const_iterator it = _mat.begin(); it != _mat.end(); ++it, ++im ){
    switch ( it->type() ){
      case Ops::Fock:
        if ( verb > 1 )
          *(pout->pout) << diag["fockn"] << "{t" << im << "}" << std::endl;
        else
          *(pout->pout) << diag["fock"] << "{t" << im << "}" << std::endl;
        break;
      case Ops::FluctP:
        if ( verb > 1 )
          *(pout->pout) << diag["flucpotn"] << "{t" << im << "1}{t" << im << "2}" << std::endl;
        else
          *(pout->pout) << diag["flucpot"] << "{t" << im << "1}{t" << im << "2}" << std::endl;
        break;
      case Ops::Deexc:
        if ( verb > 0 )
          *(pout->pout) << diag["dexop"] << "[$" << it->name() << "$]{}{" << it->exccl() << "}{t" << im << "}" << std::endl;
        else
          *(pout->pout) << diag["dexop"] << "{}{" << it->exccl() << "}{t" << im << "}" << std::endl;
        break;
      case Ops::Exc:
        if ( verb > 0 )
          *(pout->pout) << diag["exop"] << "[$" << it->name() << "$]{}{" << it->exccl() << "}{t" << im << "}" << std::endl;
        else
          *(pout->pout) << diag["exop"] << "{}{" << it->exccl() << "}{t" << im << "}" << std::endl;
        break;
      case Ops::Deexc0:
        *(pout->pout) << diag["bdexop"] << "{" << it->exccl() << "}{t" << im << "}" << std::endl;
        break;
      case Ops::Exc0:
        *(pout->pout) << diag["bexop"] << "{" << it->exccl() << "}{t" << im << "}" << std::endl;
        break;
      case Ops::XPert:
        error("Diagram for external perturbations is not possible yet!");
        break;
      case Ops::Interm:
        error("Diagram for Intermediates is not possible yet!");
        break;
      case Ops::Number:
        break;
      case Ops::None:
        error("Why is Ops::None still there??");
      default:
        xout << *it << std::endl;
        error("Diagram for this matrix is not possible yet!");
    }
  }
  // put connection lines
  TOrbSet orbs;
  im = 0;
  for ( Product<Matrices>::const_iterator it = _mat.begin(); it != _mat.end(); ++it, ++im ){
    assert( it->conlines().size() == it->orbitals().size() );
    for ( lui iorb = 0; iorb < it->conlines().size(); ++iorb ){
      if (orbs.insert(it->orbitals()[iorb]).second) {// is new
        const ConLine & cl = it->conline(iorb);
        if ( iorb%2 == 0 ){
          assert( (it->type() < _mat[cl.imat].type() && it->orbitals()[iorb].type() == Orbital::Virt) ||
                  (it->type() > _mat[cl.imat].type() && it->orbitals()[iorb].type() == Orbital::Occ) );
          assert( cl.idx%2 == 1 );
          if ( verb > 2 )
            *(pout->pout) << diag["conline"]  << "[$" << it->orbitals()[iorb] << "$]{t" << im << int(iorb/2+1) << "}{t" << cl.imat << int(cl.idx/2+1) << "}" << std::endl;
          else
            *(pout->pout) << diag["conline"] << "{t" << im << int(iorb/2+1) << "}{t" << cl.imat << int(cl.idx/2+1) << "}" << std::endl;
        } else {
          assert( (it->type() > _mat[cl.imat].type() && it->orbitals()[iorb].type() == Orbital::Virt) ||
                  (it->type() < _mat[cl.imat].type() && it->orbitals()[iorb].type() == Orbital::Occ) );
          assert( cl.idx%2 == 0 );
          if ( verb > 2 )
            *(pout->pout) << diag["conline"] << "[$" << it->orbitals()[iorb] << "$]{t" << cl.imat << int(cl.idx/2+1) << "}{t" << im << int(iorb/2+1) << "}" << std::endl;
          else
            *(pout->pout) << diag["conline"] << "{t" << cl.imat << int(cl.idx/2+1) << "}{t" << im << int(iorb/2+1) << "}" << std::endl;
        }
      }
    }
  }
  *(pout->pout) << diag["ediag"] << std::endl;
}

Orbital Term::freeorbname(Orbital::Type type)
{
  TsPar& orbs = Input::sPars["syntax"];
  Orbital::Spin spin=Orbital::GenS;
  const std::string * ip_orbs;
  std::string lastorb=_lastorb[type].name();
  unsigned long int indx;
  if (type==Orbital::Occ)
    ip_orbs = & orbs["occorb"];
  else if (type==Orbital::Virt)
    ip_orbs = & orbs["virorb"];
  else if (type==Orbital::Act)
    ip_orbs = & orbs["actorb"];
  else
    ip_orbs = & orbs["genorb"];
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

