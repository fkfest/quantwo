#include "term.h"

Term::Term() : _prefac(1), _lastel(0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd) :
    _opProd(opProd), _prefac(1), _lastel(0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

Term::Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd) :
    _opProd(opProd), _kProd(kProd), _prefac(1), _lastel(0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}
    
Term::Term(const Product< SQOp >& opProd, const Product< Kronecker >& kProd, 
           const Product< Matrices >& mat, const TOrbSet& sumindx, const TOrbSet& realsumindx, 
           const TFactor& prefac, const std::vector< Product<long int> >& connections) :
    _opProd(opProd), _kProd(kProd), _mat(mat), _sumindx(sumindx), _realsumindx(realsumindx), 
    _prefac(prefac), _connections(connections), _lastel(0) {_nloops=_nintloops=_nocc=0;_perm+=Permut();}

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
Term& Term::operator*=(const Sum< Permut, TFactor >& perm)
{
  if ( _perm.size() == 0 ){
    _perm = perm;
    return *this;
  }
  Sum<Permut,TFactor> perms;
  for ( Sum<Permut,TFactor>::const_iterator ip = perm.begin(); ip != perm.end(); ++ip){
    Sum<Permut,TFactor> perms1;
    Permut perm1;
    for ( Sum<Permut,TFactor>::iterator i_p = _perm.begin(); i_p != _perm.end(); ++i_p ){
      perm1 = i_p->first;
      perm1 *= ip->first;
      TFactor fac = ip->second * i_p->second;
      perms += std::pair< Permut, TFactor >(perm1,fac);
    }
    _perm = perms;
  }
  return *this;
}

Term& Term::operator*=(const Term& t)
{
  _opProd *= t._opProd;
  _kProd *= t._kProd;
  _mat *= t._mat;
  _sumindx *= t._sumindx;
  _realsumindx *= t._realsumindx;
  _prefac *= t._prefac;
  *this *= t._perm;
  return *this;
}
Term& Term::operator*=(const Matrices& mat)
{
  _mat *= mat;
  return *this;
}

Sum< Term, TFactor > Term::times(const Sum< Term, TFactor >& s) const
{
  Sum< Term, TFactor > sum;
  for ( Sum<Term,TFactor>::const_iterator is = s.begin(); is != s.end(); ++is) {
    Term tt(*this);
    tt *= is->first;
    tt *= is->second;
    sum += tt;
  }
  return sum;
}
Sum< Term, TFactor > Term::times(const Sum< Matrices, TFactor >& s) const
{
  Sum< Term, TFactor > sum;
  for ( Sum<Matrices,TFactor>::const_iterator is = s.begin(); is != s.end(); ++is) {
    Term tt(*this);
    tt *= is->first;
    tt *= is->second;
    sum += tt;
  }
  return sum;
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
bool Term::removeit() const
{
  for ( uint i = 0; i < _mat.size(); ++i ){
    if (_mat[i].is0()) return true;
  }
  return false;
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
        if ( _opProd[i].gender()==SQOpT::Annihilator  &&  _opProd[i+1].gender()==SQOpT::Creator ) {   // yes
          // handle 1st term
          Product<SQOp> p(_opProd); // copy Product<SQOp>
          std::swap(p[i], p[i+1]); // swap operators
          // check if we need downward recursion
          if ( !fullyContractedOnly || p[_opProd.size()-1].gender()==SQOpT::Creator )
            sum -= Term(p, _kProd, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrder(fullyContractedOnly);

            // handle 2nd term 
          Product<SQOp> q(p);   // copy Product<SQOp>
          q.erase(q.begin()+i,q.begin()+i+2);
          // q.erase(std::find(q.begin(), q.end(), p[i]));        // remove
          // q.erase(std::find(q.begin(), q.end(), p[i+1]));      // operators
        
          Product<Kronecker>  d(_kProd);
          d *= Kronecker(p[i].orb(), p[i+1].orb());           // and add kronecker
          // check if we need downward recursion
          if ( !fullyContractedOnly || q.size()==0 || q[q.size()-1].gender()==SQOpT::Creator )
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
    creat2=(( _opProd[i].genderPH()==SQOpT::Annihilator || _opProd[i].genderPH()==SQOpT::Gen )
              &&  _opProd[i+1].genderPH()==SQOpT::Creator );
    annih1=(_opProd[i].genderPH()==SQOpT::Annihilator 
              && ( _opProd[i+1].genderPH()==SQOpT::Creator || _opProd[i+1].genderPH()==SQOpT::Gen));
    if ( creat2 || annih1 ) {   // yes
      // handle 1st term
      Product<SQOp> p(_opProd); // copy Product<SQOp>
      std::swap(p[i], p[i+1]); // swap operators
      // check if we need downward recursion
      if ( !fullyContractedOnly || (p[_opProd.size()-1].genderPH()==SQOpT::Creator)
                || (p[0].genderPH()==SQOpT::Annihilator))
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
        if ( !fullyContractedOnly || q.size()==0 || q[q.size()-1].genderPH()==SQOpT::Creator 
                || q[0].genderPH()==SQOpT::Annihilator)
              sum += Term(q, d, _mat, _sumindx, _realsumindx, _prefac, _connections).normalOrderPH(fullyContractedOnly);
      }
      return sum;
    }
  }
  if ( !fullyContractedOnly || _opProd.size()==0 ) 
    sum += *this;   // add current term to the intermediate sum
  return sum;
}

Sum< Term, TFactor > Term::wickstheorem(bool genw, int noord) const
{
  // generate "matrix" of indices to SQops
  TWOps opers;
  TWMats opermat;
  unsigned int m=0;
  int cran = 0;
  for (unsigned int i=0; i<_opProd.size(); i++) {
    // calculate #creators - #annihilators
    if ( _opProd[i].gender() == SQOpT::Creator ){
      ++cran;
    } else if ( _opProd[i].gender() == SQOpT::Annihilator ){
      --cran;
    } else
      assert(false);
  }
  if ( cran != 0 ) return Sum< Term, TFactor>();
  for (unsigned int i=0; i<_opProd.size(); i++) {
    if (m==_mat.size()) { // all SQops, which are not in Matrices have to be added as individual vectors
      opermat.push_back(i);
      opers.push_back(opermat);
      opermat=TWMats();
      m=0; // search from begin
    } else if (_mat[m].orbitals().find(_opProd[i].orb())>=0) {
      opermat.push_back(i);
      if ( noord > 1 || ( noord == 1 && InSet(_mat[m].type(),Ops::Fock,Ops::FluctP,Ops::XPert) ) ){
        // not normal ordered - add every corresponding SQOp as an individual operator
        opers.push_back(opermat);
        opermat=TWMats();
      }
    } else {
      m++;
      i--;
      if (opermat.size()>0) opers.push_back(opermat);
      opermat=TWMats();
    }
  }
  if (opermat.size()>0) opers.push_back(opermat);
  opermat=TWMats();
  
//   for (TWOps::iterator iop = opers.begin(); iop != opers.end(); ++iop) {
//     for ( TWMats::iterator ijop = iop->begin(); ijop != iop->end(); ++ijop)
//       std::cout << *ijop << " " ;
//     std::cout << std::endl;
//   }
  if (genw){
    TWMats densmat;
    return genwick(opers,opermat,densmat);
  } else {
    return wick(opers,opermat);
  }
}

Sum< Term, TFactor > Term::wick(TWOps& opers, TWMats& krons) const
{
  Sum<Term, TFactor>  sum;
  if (opers.size()==0) { // no SQoperators left
    Product<SQOp> p;
    // generate Kroneckers
    Product<Kronecker> d;
    assert(krons.size()%2 == 0);
    for (TWMats::iterator kr = krons.begin(); kr != krons.end(); ++kr){
      TWMats::iterator kr0 = kr;
      ++kr;
      d*=Kronecker(_opProd[*kr0].orb(),_opProd[*kr].orb());
    }
    sum += Term(p,d,_mat, _sumindx, _realsumindx, _prefac, _connections);
    return sum;
  }
  lui istart,sign;
  TWOps::iterator ifirstop = opers.begin();
  int curr=*(ifirstop->begin());
  SQOpT::Gender gencurr=_opProd[curr].gender();
  Orbital::Type orbtypecurr=_opProd[curr].orb().type();
  bool orbcurrgen = (orbtypecurr==Orbital::GenT);
  // remove first SQop-index
  if (ifirstop->size()<2) {
    opers.erase(ifirstop);
    istart=0;
  } else {
    ifirstop->erase(ifirstop->begin());
    istart=1;
  }
  // add first index to krons1
  krons.push_back(curr);
  sign=0;
  uint ii = 0;
  for ( TWOps::iterator iop = opers.begin(); iop != opers.end(); ++iop, ++ii ) { // iterate over all Operators 
    if (ii < istart) { //but the first
      //count SQops for sign
      sign+=iop->size();
      continue;
    }
    uint jj = 0;
    for ( TWMats::iterator ijop = iop->begin(); ijop != iop->end(); ++ijop, ++jj ) {// iterate over all SQop in Operator i
      // check if the first operator and operator i would yield a Kronecker
      const SQOp& op = _opProd[*ijop];
      Orbital::Type oporbtype = op.orb().type();
      if (op.gender()!=gencurr && 
          (orbcurrgen || oporbtype == orbtypecurr|| oporbtype == Orbital::GenT)) {
        // copy opers and krons
        TWOps opers1(opers);
        TWMats krons1(krons);
        // remove SQop-index
        TWOps::iterator iop1 = opers1.begin();
        std::advance(iop1,ii);
        if (iop1->size() == 1)
          opers1.erase(iop1);
        else {
          TWMats::iterator ijop1 = iop1->begin();
          std::advance(ijop1,jj);
          iop1->erase(ijop1);
        }
        // add index to "Kronecker"
        krons1.push_back(*ijop);
        //call wick recursivly
        if ((sign+jj)%2 == 0)
          sum+=wick(opers1,krons1);
        else
          sum-=wick(opers1,krons1);
      }
    }
    sign+=iop->size();
  }
  return sum;
}

Sum< Term, TFactor > Term::genwick(Term::TWOps& opers, const Term::TWMats& krons, Term::TWMats densmat) const
{
  Sum<Term, TFactor>  sum;
  if (opers.size()==0) { // no SQoperators left
    Product<SQOp> p;
    // generate Kroneckers
    Product<Kronecker> d;
    assert(krons.size()%2 == 0);
    for (TWMats::const_iterator kr = krons.begin(); kr != krons.end(); ++kr){
      TWMats::const_iterator kr0 = kr;
      ++kr;
      d*=Kronecker(_opProd[*kr0].orb(),_opProd[*kr].orb());
    }
    // add density matrices
    Product<Orbital> dmorbs;
    Product<SQOpT::Gender> dmcran;
    assert(densmat.size()%2 == 0);
    for (TWMats::const_iterator dm = densmat.begin(); dm != densmat.end(); ++dm){
      dmorbs.push_back(_opProd[*dm].orb());
      dmcran.push_back(_opProd[*dm].gender());
    }
    if (dmorbs.size() > 0){
      Product<Matrices> mat(_mat);
      short npair = dmorbs.size()/2;
      mat *= Matrices(Ops::DensM,dmorbs,npair);
      mat.back().set_cran(dmcran);
      sum += Term(p,d,mat, _sumindx, _realsumindx, _prefac, _connections);
    } else {
      sum += Term(p,d,_mat, _sumindx, _realsumindx, _prefac, _connections);
    }
    return sum;
  }
  lui istart,sign;
  TWOps::iterator ifirstop = opers.begin();
  int curr=*(ifirstop->begin());
  SQOpT::Gender gencurr=_opProd[curr].gender();
  Orbital::Type orbtypecurr=_opProd[curr].orb().type();
  bool orbcurrgen = (orbtypecurr==Orbital::GenT);
  bool orbcurract = (orbtypecurr==Orbital::Act);
  // remove first SQop-index
  if (ifirstop->size()<2) {
    opers.erase(ifirstop);
    istart=0;
  } else {
    ifirstop->erase(ifirstop->begin());
    istart=1;
  }
  sign=0;
  uint ii = 0;
  for ( TWOps::iterator iop = opers.begin(); iop != opers.end(); ++iop, ++ii ) { // iterate over all Operators 
    if (ii < istart) { //but the first
      //count SQops for sign
      sign+=iop->size();
      continue;
    }
    uint jj = 0;
    for ( TWMats::iterator ijop = iop->begin(); ijop != iop->end(); ++ijop, ++jj ) {// iterate over all SQop in Operator i
      // check if the first operator and operator i would yield a Kronecker
      const SQOp& op = _opProd[*ijop];
      Orbital::Type oporbtype = op.orb().type();
      if (op.gender()!=gencurr && 
          (orbcurrgen || oporbtype == orbtypecurr|| oporbtype == Orbital::GenT)) {
        // copy opers and krons
        TWOps opers1(opers);
        TWMats krons1(krons);
//         TWMats densmat1(densmat);
        // remove SQop-index
        TWOps::iterator iop1 = opers1.begin();
        std::advance(iop1,ii);
        if (iop1->size() == 1)
          opers1.erase(iop1);
        else {
          TWMats::iterator ijop1 = iop1->begin();
          std::advance(ijop1,jj);
          iop1->erase(ijop1);
        }
        // add index to "Kronecker"
        krons1.push_back(curr);
        krons1.push_back(*ijop);
        //call wick recursivly
        if ((sign+jj)%2 == 0)
          sum+=genwick(opers1,krons1,densmat);
        else
          sum-=genwick(opers1,krons1,densmat);
      }
    }
    sign+=iop->size();
  }
  if (orbcurrgen || orbcurract) {
    // add density matrix
    densmat.push_front(curr);
    if ((sign)%2 == 0)
      sum+=genwick(opers,krons,densmat);
    else
      sum-=genwick(opers,krons,densmat);
  }
  return sum;
}
Sum< Term, TFactor > Term::dmwickstheorem(const Matrices& dm) const
{
  assert( dm.type() == Ops::DensM );
  bool dmsort = (Input::iPars["prog"]["dmsort"] > 0);
  // generate "matrix" of indices to SQops
  TWMats opers, krons;
  const Product<Orbital>& orbs(dm.orbitals());
  // sort now in the singlet order
  assert( orbs.size()%2 == 0 );
  {
    std::map<Electron,int> els;
    const Product<SQOpT::Gender>& cran = dm.get_cran();
    assert( cran.size() == orbs.size() );
    for ( uint i = 0; i < orbs.size(); ++i ){
      opers.push_back(i);
      // every electron has to be presented by one creator and one annihilator!
      // otherwise the term is zero!
      if ( cran[i] == SQOpT::Creator ){
        els[orbs[i].spin().el()] += 1;
      } else if ( cran[i] == SQOpT::Annihilator ){
        els[orbs[i].spin().el()] -= 1;
      }
    }
    std::map<Electron,int>::const_iterator itel;
    _foreach(itel,els){
      if ( itel->second != 0 ) return Sum< Term, TFactor >();
    }
  }
  // only dmsort version is implemented yet...
  assert(dmsort);
  assert(dm.get_cran().size() == dm.orbitals().size());
  return dmwick(opers,krons,dm);
}
Sum< Term, TFactor > Term::dmwick(Term::TWMats& opers, const Term::TWMats& krons, const Matrices& dm) const
{
  Sum<Term, TFactor>  sum;
  const Product<Orbital>& orbs(dm.orbitals());
  const Product<SQOpT::Gender> cranorder(dm.get_cran());
  TWMats::iterator itelc = opers.begin(), itela = opers.end();
  --itela;
  uint nel =  opers.size()/2;
  uint sign = 0;
  for ( uint iel = 0; iel < nel; ++iel, ++itelc, --itela ){
    // find the first creator operator
    TWMats::iterator itelcr;
    uint cur = iel;
    for ( itelcr = itelc; itelcr != opers.end() && (cranorder[*itelcr] != SQOpT::Creator); ++itelcr, ++cur ){}
    if ( itelcr == opers.end() ) break; // try next time after resolving of kroneckers
//    assert( itelcr != opers.end() );
    
    for ( ; itelcr != itelc; --itelcr ){
      // move creator to position itelc
      TWMats::iterator it = itelcr;
      --it;
      --cur;
      std::swap(*it,*itelcr);
      if ( cranorder[*itelcr] == SQOpT::Annihilator ){
        TWMats opers1(opers);
        TWMats krons1(krons);
        // remove SQop-index
        TWMats::iterator iop1 = opers1.begin();
        std::advance(iop1,cur);
        opers1.erase(iop1);
        iop1 = opers1.begin();
        std::advance(iop1,cur);
        opers1.erase(iop1);
        // add index to "Kronecker"
        krons1.push_back(*it);
        krons1.push_back(*itelcr);
        //call wick recursivly
        if ((sign)%2 == 0)
          sum+=dmwick(opers1,krons1,dm);
        else
          sum-=dmwick(opers1,krons1,dm);
      }
      ++sign;
    }
    const Spin& spin( orbs[*itelc].spin() );
    // find annihilator of the same electron (from end)
    TWMats::iterator itelan;
    cur = opers.size() - 1 - iel;
    for ( itelan = itela; (cranorder[*itelan] != SQOpT::Annihilator || orbs[*itelan].spin() != spin ) && itelan != opers.begin(); --itelan, --cur ){}
    if ( itelan == opers.begin() ) break; // try next time after resolving of kroneckers
//    assert( itelan != opers.begin() );
    for ( ; itelan != itela; ++itelan, ++cur ){
      // move creator to position itelc
      TWMats::iterator it = itelan;
      ++it;
      std::swap(*it,*itelan);
      if ( cranorder[*itelan] == SQOpT::Creator ){
        TWMats opers1(opers);
        TWMats krons1(krons);
        // remove SQop-index
        TWMats::iterator iop1 = opers1.begin();
        std::advance(iop1,cur);
        opers1.erase(iop1);
        iop1 = opers1.begin();
        std::advance(iop1,cur);
        opers1.erase(iop1);
        // add index to "Kronecker"
        krons1.push_back(*it);
        krons1.push_back(*itelan);
        //call wick recursivly
        if ((sign)%2 == 0)
          sum+=dmwick(opers1,krons1,dm);
        else
          sum-=dmwick(opers1,krons1,dm);
      }
      ++sign;
    }
  }
  
  Product<SQOp> p(_opProd);
  // generate Kroneckers
  Product<Kronecker> d(_kProd);
  assert(krons.size()%2 == 0);
  for (TWMats::const_iterator kr = krons.begin(); kr != krons.end(); ++kr){
    TWMats::const_iterator kr0 = kr;
    ++kr;
    d *= Kronecker(orbs[*kr0],orbs[*kr]);
  }
  // add density matrix
  Product<Matrices> mat(_mat);
  if ( opers.size() > 0 ){
    Product<Orbital> dmorbs;
    Product<SQOpT::Gender> dmcran;
    for (TWMats::const_iterator dm = opers.begin(); dm != opers.end(); ++dm){
      dmorbs.push_back(orbs[*dm]);
      dmcran.push_back(cranorder[*dm]);
    }
    short npair = dmorbs.size()/2;
    mat *= Matrices(Ops::DensM,dmorbs,npair);
    mat.back().set_cran(dmcran);
//    if( !mat.back().nonsingldm() ) xout << "nonsingl" << mat << std::endl;
  }
  if ((sign)%2 == 0)
    sum += Term(p,d,mat, _sumindx, _realsumindx, _prefac, _connections);
  else
    sum -= Term(p,d,mat, _sumindx, _realsumindx, _prefac, _connections);
  return sum;
}

void Term::setmatconnections()
{
  long int kj;
  for (lui i=0; i<_mat.size(); ++i)
    for (lui j=0; j<_mat[i].orbitals().size(); ++j){
      const Orbital& ijorb = _mat[i].orbitals()[j];
      if ( (kj = _mat[i].orbitals().find(ijorb,j+1)) >= 0 ) {
        // self-connection (non-normal-ordered hamiltonian)
        assert(Input::iPars["prog"]["noorder"] > 0);
        _mat[i].add_connect(i+1);
//        _mat[i].add_connect(i+1);
        _mat[i].set_conline(j,i,kj);
        _mat[i].set_conline(kj,i,j);
        continue;
      } 
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
  TOrbSet::iterator it1, it2;
  for ( int i = 0 ; i < int(_kProd.size()) ; ++i ) // iterate over Product<Kronecker>
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
    Orbital 
      orb1 = _kProd[i].orb1(),
      orb2 = _kProd[i].orb2();
    bool insum = false;
    if ( it2 != _realsumindx.end() ) { // found orb2
      insum = true;
    } else if ( it1 != _realsumindx.end() ) { // found orb1
      insum = true;
      std::swap(it1,it2);
      std::swap(orb1,orb2);
    }
    if ( insum ) { // found in sum
      _realsumindx.erase(it2); // delete summation over orb2
      it2 = _sumindx.find(orb2);
      if ( it2 == _sumindx.end() )
        error("Strange, orbital not found in _sumindx","Term::reduceTerm");
      _sumindx.erase(it2);
      _kProd.erase(_kProd.begin()+i); // delete the Kronecker
      --i;
      this->replace(orb2,orb1);
    }
  }
}
void Term::krons2mats()
{
  Product<Kronecker>::const_iterator ik;
  _foreach(ik,_kProd){
    _mat.push_back(Matrices(*ik));
  }
  _kProd = Product<Kronecker>();
}

void Term::replace(Orbital orb1, Orbital orb2, bool smart)
{
  Q2::replace(_opProd,orb1,orb2,smart);
  Q2::replace(_mat,orb1,orb2,smart);
  Q2::replace(_kProd,orb1,orb2,smart);
  Q2::replace(_realsumindx,orb1,orb2,smart);
  Q2::replace(_sumindx,orb1,orb2,smart);
}

static bool matisnone(const Matrices& mat)
{ return (mat.type() == Ops::None); };
void Term::deleteNoneMats()
{ _mat.erase(std::remove_if(_mat.begin(),_mat.end(),matisnone),_mat.end()); }
bool Term::brilloin() const
{
  for ( Product<Matrices>::const_iterator it = _mat.begin(); it != _mat.end(); ++it )
    if ( it->type() == Ops::Fock && 
         (( it->orbitals()[0].type() == Orbital::Occ && it->orbitals()[1].type() == Orbital::Virt ) ||
          ( it->orbitals()[1].type() == Orbital::Occ && it->orbitals()[0].type() == Orbital::Virt )) )
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
    if (_mat[i].expandantisym(firstpart)) {
      if ( !firstpart ){
        // electrons have to be swaped in the matrix and all corresponding kroneckers
        assert(_mat[i].orbitals().size() == 4);
        Orbital 
          orb1 = _mat[i].orbitals()[1],
          orb3 = _mat[i].orbitals()[3],
          orbn1(orb1), orbn3(orb3);
        orbn1.setspin(orb3.spin());
        orbn3.setspin(orb1.spin());
        this->replace(orb1,orbn1,false);
        this->replace(orb3,orbn3,false);
      }
      return true;
    }
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
bool Term::has_nonsingldm() const
{
  for ( uint i = 0; i < _mat.size(); ++i )
    if (_mat[i].nonsingldm()) return true;
  return false;
}
Sum< Term, TFactor > Term::dm2singlet()
{
  Sum< Term, TFactor > sum;
  for ( uint i = 0; i < _mat.size(); ++i ){
    if (_mat[i].nonsingldm()){
      Matrices mat(_mat[i]);
      // bring into singlet order
      _mat.erase(_mat.begin()+i);
      sum = this->dmwickstheorem(mat);
      // can transform only one matrix per call
      return sum;
    }
  }
  sum += *this;
  return sum;
}
bool Term::has_generalindices() const
{
  TOrbSet::const_iterator it;
  _foreach(it,_sumindx)
    if (it->type() == Orbital::GenT) return true;
  return false;
}
Sum< Term, TFactor > Term::removegeneralindices()
{
  bool active = (Input::iPars["prog"]["multiref"] > 0);
  Sum< Term, TFactor > sum;
  this->set_lastorbs();
  Term tt(*this);
  TOrbSet::iterator it;
  _foreach(it,_sumindx){
    if ( it->type() == Orbital::GenT ){
      // replace
      Orbital orb = tt.freeorbname(Orbital::Occ);
      orb.setspin(it->spin());
      tt.replace(*it,orb);
      if (active){
        sum += tt;
        tt = *this;
        orb = tt.freeorbname(Orbital::Act);
        orb.setspin(it->spin());
        tt.replace(*it,orb);
        sum += tt;
        // have to repeat it
        return sum;
      }
    }
  }
  sum += tt;
  return sum;
}

void Term::spinintegration(bool notfake)
{
  // generate Product of all orbitals and external-lines orbitals
  TOrbSet peo(extindx()), peo1(peo);
  // internal indices
  TOrbSet po(_realsumindx); 
  _nocc=_nintloops=_nloops=0;
  Spin nospin(Spin::No);
  bool nonconserve = true, already_done = false;
  Orbital orb,orb1;
  Matrices::Spinsym spinsym=Matrices::Singlet;
  lui ithis = 0;
  long int ipos, iorb = 0;
  TOrbSet::iterator it;
  bool samespinsym, withdensmat;
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
    withdensmat = false;
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
      withdensmat = withdensmat || (_mat[ithis].type() == Ops::DensM);
      iorb = _mat[ithis].iorbel(ipos);
      if ( iorb < 0 ) break;
      orb1 = _mat[ithis].orbitals()[iorb];
      if (nonconserve){
        // test whether still in sets
        already_done = ( peo.count(orb1) == 0 && po.count(orb1) == 0 );
      }
    } while (orb1!=orb && !already_done && _sumindx.count(orb1));
    if (samespinsym) {
      if (notfake && !withdensmat) _prefac*=2;
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
    }
    TOrbSet sumindx;
    for ( TOrbSet::iterator it = _sumindx.begin(); it != _sumindx.end(); ++it ){
      orb = *it;
      orb.setspin(nospin);
      sumindx.insert(orb);
    }
    _sumindx = sumindx;
    sumindx.clear();
    for ( TOrbSet::iterator it =_realsumindx.begin(); it != _realsumindx.end(); ++it ){
      orb = *it;
      orb.setspin(nospin);
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
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spin = Spin::Gen;
  if (spinintegr) spin = Spin::GenS;
  const std::string * ip_orbs;
  std::string lastorb=_lastorb[type].letname();
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
void Term::set_lastorbs()
{
  TOrbSet::const_iterator it;
  _foreach(it,_sumindx){
    if (_lastorb[it->type()].name().size() == 0 || _lastorb[it->type()] < *it) 
      _lastorb[it->type()] = *it;
  }
}

Electron Term::nextelectron()
{
  ++_lastel;
  return _lastel;
}
Electron Term::getnextelectron(void* Obj)
{
  // explicitly cast to a pointer to Term
  Term * myself = (Term *) Obj;
  // call member
  return myself->nextelectron();
}
void Term::set_lastel(Electron el, bool onlylarger)
{
  if ( onlylarger && el < _lastel ) return;
  _lastel = el;
}


template <class T, class Q>
Return Q2::replace(Product< T >& p, Q orb1, Q orb2, bool smart)
{
  Return rpl;
  for ( unsigned int i=0; i<p.size(); ++i )
  {
    rpl += replace(p[i],orb1, orb2, smart);
  }
  return rpl;
}
template <class T, class Q>
Return Q2::replace(Set< T >& p, Q orb1, Q orb2, bool smart)
{
  Return rpl;
  Set<T> tmp;
  std::pair<typename Set<T>::iterator,bool> ret;
  for ( typename Set<T>::iterator it = p.begin(); it != p.end(); ++it )
  { 
    T t(*it);
    rpl += replace(t,orb1, orb2,smart);
    ret = tmp.insert(t);
    if (!ret.second)
      say("Strange, summation runs already over "+t.name());
  }
  p = tmp;
  return rpl;
}
template <class T>
Return Q2::replace(SQOp& op, T orb1, T orb2, bool smart)
{ return op.replace(orb1,orb2,smart); }
template <class T>
Return Q2::replace(T &orb, T orb1, T orb2, bool smart)
{ return orb.replace(orb1,orb2,smart);
  //if (orb == orb1) orb=orb2; 
}
template <class T>
Return Q2::replace(Matrices &mat, T orb1, T orb2, bool smart)
{ return mat.replace(orb1,orb2,smart); }
template <class T>
Return Q2::replace(Kronecker& kron, T orb1, T orb2, bool smart)
{ return kron.replace(orb1,orb2,smart); }

