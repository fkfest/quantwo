#include "term.h"

Term::Term() : _prefac(1), _lastel(0), _matconnectionsset(false) 
{
  _nloops=_nintloops=_nocc=0;
  _perm+=Permut();
}

Term::Term(Product<SQOp> const & opProd) :
    _opProd(opProd), _prefac(1), _lastel(0), _matconnectionsset(false) 
{
  _nloops=_nintloops=_nocc=0;
  _perm+=Permut();
}

Term::Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd) :
    _opProd(opProd), _kProd(kProd), _prefac(1), _lastel(0), _matconnectionsset(false) 
{
  _nloops=_nintloops=_nocc=0;
  _perm+=Permut();
}
    
Term::Term(const Product< SQOp >& opProd, const Product< Kronecker >& kProd, 
           const Product< Matrix >& mat, const TOrbSet& orbs, const TOrbSet& sumorbs, 
           const TFactor& prefac, const ConnectionsMap& connections) :
    _opProd(opProd), _kProd(kProd), _mat(mat), _orbs(orbs), _sumorbs(sumorbs), 
    _prefac(prefac), _connections(connections), _lastel(0), _matconnectionsset(false) 
{
  _nloops=_nintloops=_nocc=0;
  _perm+=Permut();
}

Term& Term::operator*=(const Oper& op)
{
  _opProd *= op.SQprod();
  _mat *= op.mat();
  _orbs *= op.orbs();
  _sumorbs *= op.sumorbs();
  _prefac *= op.prefac();
  if ( op.sumops().size() > 0 )
    _termsfacs *= op.sumops();
  assert( _opProd.size() == 0 || _termsfacs.size() == 0 );
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
#ifndef NDEBUG
  // all orbitals should differ!
  _foreach_cauto( TOrbSet, itorb, _orbs ){
    assert( t._orbs.find(*itorb) == t._orbs.end() );
  }
#endif
  _opProd *= t._opProd;
  _kProd *= t._kProd;
  _mat *= t._mat;
  _orbs *= t._orbs;
  _sumorbs *= t._sumorbs;
  _prefac *= t._prefac;
  *this *= t._perm;
  _termsfacs *= t._termsfacs;
  assert( _opProd.size() == 0 || _termsfacs.size() == 0 );
  return *this;
}
Term& Term::operator*=(const Matrix& mat)
{
  _mat *= mat;
  return *this;
}

TermSum Term::times(const TermSum& s) const
{
  TermSum sum;
  for ( TermSum::const_iterator is = s.begin(); is != s.end(); ++is) {
    Term tt(*this);
    tt *= is->first;
    tt *= is->second;
    sum += tt;
  }
  return sum;
}
TermSum Term::times(const Sum< Matrix, TFactor >& s) const
{
  TermSum sum;
  for ( Sum<Matrix,TFactor>::const_iterator is = s.begin(); is != s.end(); ++is) {
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
    ret = _sumorbs.insert(orb1);
    if (!ret.second)
      say("Strange, summation runs already over "+orb1.name());
  }
}
void Term::addsummation (const Product<Orbital> & orbs)
{
  std::pair<TOrbSet::iterator,bool> ret;
  for ( uint i = 0; i < orbs.size(); ++i ){
    ret = _sumorbs.insert(orbs[i]);
    if (!ret.second)
      say("Strange, summation runs already over "+orbs[i].name());
  }
}
void Term::addsummation(const Orbital& orb)
{
  std::pair<TOrbSet::iterator,bool> ret;
  ret = _sumorbs.insert(orb);
  if (!ret.second)
    say("Strange, summation runs already over "+orb.name());
}
void Term::addorb(const Orbital& orb)
{
  std::pair<TOrbSet::iterator,bool> ret;
  ret = _orbs.insert(orb);
  if (!ret.second)
    say("Strange, orbital is already there "+orb.name());
}

void Term::replacematrix(const Matrix& mat, lui ipos)
{
  if (ipos >= _mat.size())
    error("The position is outside of this term: "+any2str(ipos),"Term::replacematrix");
  _mat[ipos]=mat;
}

TermSum Term::expandtermsfacs()
{
  TermSum ts, tts;
  Product<TermSum> trmsfacs= _termsfacs;
  _termsfacs.clear();
  ts += *this;
  _foreach_cauto( Product<TermSum>, itfs, trmsfacs ){
    _foreach_cauto( TermSum, its, ts ){
      TermSum trms = its->first.times(*itfs);
      trms *= its->second;
      tts += trms;
    }
    ts = tts;
    tts.clear();
  }
#ifndef NDEBUG
  // all termsfacs in ts should be empty
  _foreach_cauto( TermSum, its, ts ){
    assert( its->first._termsfacs.size() == 0 );
  }
#endif
  return ts;
}

void Term::addoverlaps()
{
  TOrbSet orbs_done;
  //loop over matrices and check connections
  //put overlaps between amplitudes, and between amplitudes and external lines
  assert( _matconnectionsset );
  lui nmats = _mat.size();
  for ( lui im = 0; im < nmats; ++im ){
    if ( InSet(_mat[im].type(),Ops::Exc,Ops::Deexc,Ops::Exc0,Ops::Deexc0,Ops::DensM) ) {
      for ( uint iorb = 0;  iorb < _mat[im].orbitals().size(); ++iorb ) {
        const Orbital & orb = _mat[im].orbitals()[iorb];
        if ( orb.type() == Orbital::Virt ) {
          auto ret = orbs_done.insert(orb);
          if ( !ret.second ) {
            // already handled
            continue;
          }
          const ConLine& cl = _mat[im].conline(iorb);
          // check whether _mat[cl.imat] is excitation or deexcitation,
          // and an external line
          bool exc =false, exter = false;
          switch (_mat[cl.imat].type()) {
            case Ops::DensM:
            case Ops::Exc:
              exc = true;
              break;
            case Ops::Deexc0:
              exter = true;
              break;
            case Ops::Exc0:
              exc = true;
              exter = true;
              break;
            case Ops::Deexc:
              break;
            default:
              // overlap is not needed
              continue;
          }
          
          Orbital neworb(orb);
          do {
            neworb.add_prime();
            ret = _orbs.insert(neworb);
          } while ( !ret.second );
          _sumorbs.insert(neworb);
          orbs_done.insert(neworb);
          if ( exter ){
            // insert overlap, new orbital in _mat[im]
            _mat[im].set_orb(neworb,iorb);
          } else {
            // insert overlap, new orbital in _mat[cl.imat]
            _mat[cl.imat].set_orb(neworb,cl.idx);
          }
          Product<Orbital> oorbs;
          oorbs.resize(2);
          if ( exc ){
            oorbs[0] = orb;
            oorbs[1] = _mat[cl.imat].orbitals()[cl.idx];
          } else {
            oorbs[0] = _mat[cl.imat].orbitals()[cl.idx];
            oorbs[1] = orb;
          }
          _mat *= Matrix(Ops::Overlap,oorbs,1,0,0,"S");
        }
      }
    }
  }
  setmatconnections();
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
  for ( Product<Matrix>::iterator it = _mat.begin(); it != _mat.end(); ++it ){
    it->permute(p);
  }
}
TermSum Term::resolve_permutations() const
{
  TermSum sum;
  Sum<Permut,TFactor> empty_perm;
  for ( Sum<Permut,TFactor>::const_iterator itp = _perm.begin(); itp != _perm.end(); ++itp ){
    Term term = *this;
    term.setperm(empty_perm);
    TFactor fac = itp->second;
    term.permute(itp->first);
    sum += std::make_pair(term,fac);
  }
  return sum;
}

const Product< Matrix >& Term::mat() const
{ return _mat; }
const TOrbSet& Term::orbs() const
{ return _orbs; }
const TOrbSet& Term::sumorbs() const
{ return _sumorbs; }
TOrbSet Term::extindx() const
{
  // generate Product of external-lines orbitals
  TOrbSet peo;
  for ( lui i = 0; i < _mat.size(); ++i ) {
    for ( lui j = 0; j < _mat[i].orbitals().size(); ++j) {
      const Orbital & orb = _mat[i].orbitals()[j];
      if ( _sumorbs.count(orb) == 0 )
        peo.insert(orb);
    }
  }
  return peo;
}
TOrbSet Term::extcreaindx() const
{
  TOrbSet peo;
  _foreach_cauto(Product<Matrix>,im,_mat){
    if (!InSet(im->type(),Ops::Exc0,Ops::Deexc0)) continue;
    const Product<Orbital> orbs = im->orbitals();
    for ( uint i = 0; i < orbs.size(); ++i ){
      SQOpT::Gender gen = im->genderguess(i);
//       assert( gen != SQOpT::Gen );
      if ( gen == SQOpT::Creator && _sumorbs.count(orbs[i]) == 0 )
        peo.insert(orbs[i]);
    }
  }
  return peo;
}

Sum< Permut, TFactor > Term::perm() const
{ return _perm; }
ConnectionsMap Term::connections() const
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
  // orbitals
  Product<Orbital> from, to;
  for ( TOrbSet::const_iterator it = _sumorbs.begin(); it != _sumorbs.end(); ++it ) {
    if(_orbs.count(*it) == 0) {
      bool found = false;
      if ( it->getel() == 0 ) {
        // electron is not set
        _foreach_cauto(TOrbSet,itorb,_orbs){
          if ( it->equal(*itorb) ){
            from.push_back(*it);
            to.push_back(*itorb);
            found = true;
            break;
          }
        }
      }
      if (!found) {
        std::cout << *this << std::endl;
        error("Problem with summation index! May be a summation over excitations in a term without this excitaions?");
      }
    }
  }
  for ( uint i = 0; i < from.size(); ++i ) 
    this->replace(from[i],to[i],false);
  // check and set external matrices
  _foreach_auto( Product<Matrix>, itm, _mat ){
    itm->is_internal(_sumorbs) ;
    // external matrices only of Deexc0 or Exc0 type! 
    assert ( itm->internal() || InSet(itm->type(),Ops::Deexc0,Ops::Exc0) ); 
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
    if ( _orbs<t._orbs ) return true;
    if ( _orbs>t._orbs ) return false;
    if ( _perm<t._perm) return false;
    if ( t._perm<_perm) return true;
    if ( _termsfacs<t._termsfacs) return true;
    if ( t._termsfacs<_termsfacs) return false;
    if (_connections < t._connections) return true;
    if (t._connections < _connections) return false;
    return _prefac<t._prefac;
return true;
}
bool Term::equal(Term& t, Permut& perm) 
{
  
  if (_mat.size() != t._mat.size() ||
      _orbs.size() != t._orbs.size() ||
      _sumorbs.size() != t._sumorbs.size() ||
      _nintloops != t._nintloops || _nocc != t._nocc) return false;
//   xout << "compare " << *this << "  and  " << t << std::endl;
  // generate Product of all orbitals and external-lines orbitals
  TOrbSet peo(extindx()), peot(t.extindx()), 
          // external creator lines
          pceo(extcreaindx()), pceot(t.extcreaindx());
  List<Orbital> po(peo),pot(peot); // start with external lines!
  po*=_sumorbs; // internal indices
  pot*=t._sumorbs; // internal indices
  if (po.size() != pot.size()) return false;
  Product<Matrix> mat, tmat;
  Product<Orbital> por, port;
  Permut perm1;
  Orbital orb,orb1,orbt,orb1t;
  long int ipos=0,ipost=0;
  List<Orbital>::iterator it;
  Product<Orbital>::iterator jt;
  unsigned int ithis=0,ithist=0,i;
  bool equal=false,exter,extert,exter1,extert1,loop,loopt,
       // is it an external creator line?
       extcr,extcrt;
  Order ordmat, ordmatt, oordmat, oordmatt;
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
    extcr=(pceo.count(orb)); // external orbital from a creator operator
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
      extcrt=(pceot.count(orbt));
      if (extcr != extcrt){
        // one comes from external creator operator and the other not
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
          perm1 += Permut(orb1t,orb1);
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
      SQOpT::Gender genorb = mat[ithis].genderguess(ipos);
      assert(ordmatt.size() == tmat.size());
      ipost = -1;
      for (unsigned int jo=0; jo<ordmatt.size(); jo++) {
        unsigned int j = ordmatt[jo];
        ++ipost;
        ipost = tmat[j].orbitals().find(orb1t,ipost);
        if (ipost >= 0 ) {
          if ( !exter && tmat[j].genderguess(ipost) != genorb) { 
            // orbitals come from different gendered operators
            // exclude external indices from the gender-check because of the ambiguity in particle-nonconserved matrices
            --jo; // check same matrix again
            continue;
          }
          ithist=j;
          break;
        }
        ipost = -1;
      }
      if ( ipost < 0 ) {
        // not found
        equal=false;
        continue;
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
        
        if (orb1.type() != orb1t.type() || !mat[ithis].vertices(ipos,tmat[ithist],ipost,ithis) ||
            mat[ithis].genderguess(ipos) != tmat[ithist].genderguess(ipost) )
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
            perm1 += Permut(orb1t,orb1);
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
                                    
  if (t.sumorbs().size()>0) 
  {
    o <<"\\" << Input::sPars["command"]["sum"] <<"_{"<<t.sumorbs()<<"}";
    MyOut::pcurout->lenbuf += std::max(3,int(t.sumorbs().size())/MyOut::pcurout->wsi);
    MyOut::pcurout->hbufline = MyOut::pcurout->hsum;// set height of line to the height of \sum
  }
  o << t.mat(); //lenbuf will be handled in Matrix
  o << t.kProd(); //lenbuf will be handled in Kronecker
  if ( t.kProd().size()>0 && t.opProd().size()>0 )
    o << "*";
  o << t.opProd();
  MyOut::pcurout->lenbuf += t.opProd().size()*2;
  if ( t.termsfacs().size() > 0 ) {
    Product<TermSum> tts = t.termsfacs();
    _foreach_cauto( Product<TermSum>, itts, tts ){
      MyOut::pcurout->lenbuf++ ; // for "("
      o << "(" << *itts << ")";
      MyOut::pcurout->lenbuf++ ; // for ")"
    }
  }
return o;
}

TermSum Term::normalOrder() const
{   return normalOrder(false);  }

TermSum Term::normalOrder_fullyContractedOnly() const
{   return normalOrder(true);   }


TermSum Term::normalOrder(bool fullyContractedOnly) const
{
    TermSum  sum;
    for ( unsigned int i=0 ; i+1<_opProd.size() ; ++i ) // iterate over Product<SQOp> but the last
    {
        // check if two consecutive operators need reordering
        if ( _opProd[i].gender()==SQOpT::Annihilator  &&  _opProd[i+1].gender()==SQOpT::Creator ) {   // yes
          // handle 1st term
          Product<SQOp> p(_opProd); // copy Product<SQOp>
          std::swap(p[i], p[i+1]); // swap operators
          // check if we need downward recursion
          if ( !fullyContractedOnly || p[_opProd.size()-1].gender()==SQOpT::Creator )
            sum -= Term(p, _kProd, _mat, _orbs, _sumorbs, _prefac, _connections).normalOrder(fullyContractedOnly);

            // handle 2nd term 
          Product<SQOp> q(p);   // copy Product<SQOp>
          q.erase(q.begin()+i,q.begin()+i+2);
          // q.erase(std::find(q.begin(), q.end(), p[i]));        // remove
          // q.erase(std::find(q.begin(), q.end(), p[i+1]));      // operators
        
          Product<Kronecker>  d(_kProd);
          d *= Kronecker(p[i].orb(), p[i+1].orb());           // and add kronecker
          // check if we need downward recursion
          if ( !fullyContractedOnly || q.size()==0 || q[q.size()-1].gender()==SQOpT::Creator )
            sum += Term(q, d, _mat, _orbs, _sumorbs, _prefac, _connections).normalOrder(fullyContractedOnly);
          return sum;
        }
    }

    sum += *this;   // add current term to the intermediate sum
    return sum;
}

TermSum Term::normalOrderPH() const
{   return normalOrderPH(false);  }

TermSum Term::normalOrderPH_fullyContractedOnly() const
{   return normalOrderPH(true);   }


TermSum Term::normalOrderPH(bool fullyContractedOnly) const
{
  TermSum  sum;
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
        sum -= Term(p, _kProd, _mat, _orbs, _sumorbs, _prefac, _connections).normalOrderPH(fullyContractedOnly); 
          
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
              sum += Term(q, d, _mat, _orbs, _sumorbs, _prefac, _connections).normalOrderPH(fullyContractedOnly);
      }
      return sum;
    }
  }
  if ( !fullyContractedOnly || _opProd.size()==0 ) 
    sum += *this;   // add current term to the intermediate sum
  return sum;
}

TermSum Term::wickstheorem(bool genw, int noord) const
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
  if ( cran != 0 ) return TermSum();
  for (unsigned int i=0; i<_opProd.size(); i++) {
    if (m==_mat.size()) { // all SQops, which are not in Matrix have to be added as individual vectors
      opermat.push_back(i);
      opers.push_back(opermat);
      opermat=TWMats();
      m=0; // search from begin
    } else if (_mat[m].orbitals().find(_opProd[i].orb())>=0) {
      opermat.push_back(i);
      if ( noord > 1 || 
           ( noord == 1 && InSet(_mat[m].type(),Ops::Fock,Ops::OneEl,Ops::FluctP,Ops::XPert) ) ){
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

TermSum Term::wick(TWOps& opers, TWMats& krons) const
{
  TermSum  sum;
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
    sum += Term(p,d,_mat, _orbs, _sumorbs, _prefac, _connections);
    return sum;
  }
  lui istart,sign;
  TWOps::iterator ifirstop = opers.begin();
  int curr=*(ifirstop->begin());
  SQOpT::Gender gencurr=_opProd[curr].gender();
  Orbital::Type orbtypecurr=_opProd[curr].orb().type();
  bool orbcurrgen = (orbtypecurr==Orbital::GenT);
  if ( _opProd[curr].genderPH() == SQOpT::Creator ) {
    // quasi-Creator on the left --> this term is zero
    return sum;
  }
  // remove first SQop-index
  if (ifirstop->size() == 1) {
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
      if (op.gender()!=gencurr && op.genderPH() != SQOpT::Annihilator &&
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

TermSum Term::genwick(Term::TWOps& opers, const Term::TWMats& krons, Term::TWMats densmat) const
{
  TermSum  sum;
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
      Product<Matrix> mat(_mat);
      short npair = dmorbs.size()/2;
      mat *= Matrix(Ops::DensM,dmorbs,npair);
      mat.back().set_cran(dmcran);
      sum += Term(p,d,mat, _orbs, _sumorbs, _prefac, _connections);
    } else {
      sum += Term(p,d,_mat, _orbs, _sumorbs, _prefac, _connections);
    }
    return sum;
  }
  lui istart,sign;
  TWOps::iterator ifirstop = opers.begin();
  int curr = *(ifirstop->begin());
  SQOpT::Gender gencurr = _opProd[curr].gender();
  Orbital::Type orbtypecurr = _opProd[curr].orb().type();
  bool orbcurrgen = (orbtypecurr==Orbital::GenT);
  bool orbcurract = (orbtypecurr==Orbital::Act);
  if ( _opProd[curr].genderPH() == SQOpT::Creator ) {
    // quasi-Creator on the left --> this term is zero
    return sum;
  }
  // remove first SQop-index
  if (ifirstop->size() == 1) {
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
      if (op.gender()!=gencurr && op.genderPH() != SQOpT::Annihilator && 
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
TermSum Term::change2fock(uint imat, bool multiref) const
{
  assert( _mat[imat].type() == Ops::OneEl );
  TermSum sum;
  Term term(*this);
  // h_PQ = f_PQ - (PQ||KJ)\delta_KJ [-(PQ||TU)\gamma^T_U]
  Product<Orbital> orbs = _mat[imat].orbitals();
  TCon2 connected2 = _mat[imat].connected2();
  // f_PQ
  term._mat[imat] = Matrix(Ops::Fock,orbs,1);
  term._mat[imat].set_connect(connected2);
  sum += term;
  // (PQ||KJ)\delta_{KJ}
  term = *this;
  Orbital
    orb1 = term.freeorbname(Orbital::Occ),
    orb2 = term.freeorbname(Orbital::Occ);
  Electron el = term.nextelectron();
  orb1.setel(el);
  orb2.setel(el);
  orbs *= orb1;
  orbs *= orb2;
  term._mat[imat] = Matrix(Ops::FluctP,orbs,2);
  term._mat[imat].set_connect(connected2);
  term._kProd.push_back(Kronecker(orb1,orb2));
  term._sumorbs.insert(orb1);
  term._sumorbs.insert(orb2);
  term._orbs.insert(orb1);
  term._orbs.insert(orb2);
  sum -= term;
  if (multiref) {
    // (PQ||TU)\gamma^T_U
    term = *this;
    orbs = _mat[imat].orbitals();
    Orbital torb = term.freeorbname(Orbital::Act),
            uorb = term.freeorbname(Orbital::Act);
    Electron el = term.nextelectron();
    torb.setel(el);
    uorb.setel(el);
    orbs *= torb;
    orbs *= uorb;
    term._mat[imat] = Matrix(Ops::FluctP,orbs,2);
    term._mat[imat].set_connect(connected2);
    orbs.clear();
    orbs *= torb;
    orbs *= uorb;
    term._mat.push_back(Matrix(Ops::DensM,orbs,1));
    Product< SQOpT::Gender > cran;
    cran *= SQOpT::Creator;
    cran *= SQOpT::Annihilator;
    term._mat.back().set_cran(cran);
    term._sumorbs.insert(torb);
    term._orbs.insert(torb);
    term._sumorbs.insert(uorb);
    term._orbs.insert(uorb);
    sum -= term;
  }
  return sum;
}
TermSum Term::dmwickstheorem(const Matrix& dm) const
{
  assert( dm.type() == Ops::DensM );
  // generate "matrix" of indices to SQops
  TWMats opers, krons;
  const Product<Orbital>& orbs(dm.orbitals());
  // sort now in the singlet order
  assert( orbs.size()%2 == 0 );
  for ( uint i = 0; i < orbs.size(); ++i ){
    opers.push_back(i);
  }
  // only dmsort version is implemented yet...
  assert(Input::iPars["prog"]["dmsort"] > 0);
  assert(dm.get_cran().size() == dm.orbitals().size());
  return dmwick(opers,krons,dm);
}
TermSum Term::dmwick(Term::TWMats& opers, const Term::TWMats& krons, const Matrix& dm) const
{
  TermSum  sum;
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
  Product<Matrix> mat(_mat);
  if ( opers.size() > 0 ){
    Product<Orbital> dmorbs;
    Product<SQOpT::Gender> dmcran;
    for (TWMats::const_iterator dm = opers.begin(); dm != opers.end(); ++dm){
      dmorbs.push_back(orbs[*dm]);
      dmcran.push_back(cranorder[*dm]);
    }
    short npair = dmorbs.size()/2;
    mat *= Matrix(Ops::DensM,dmorbs,npair);
    mat.back().set_cran(dmcran);
//    if( !mat.back().nonsingldm() ) xout << "nonsingl" << mat << std::endl;
  }
  if ((sign)%2 == 0)
    sum += Term(p,d,mat, _orbs, _sumorbs, _prefac, _connections);
  else
    sum -= Term(p,d,mat, _orbs, _sumorbs, _prefac, _connections);
  return sum;
}

TermSum Term::replaceE0byfock(uint imat, bool multiref, bool replaceE0act) const
{
  assert( _mat[imat].name() == Input::sPars["hamilton"]["E0"] );
  TermSum sum;
  Term term(*this);
  // E^0 = \sum f_II [ + \sum f_TU \gamma^T_U]
  Product<Orbital> orbs;
  // f_II
  Orbital
    orb1 = term.freeorbname(Orbital::Occ);
  Electron el = term.nextelectron();
  orb1.setel(el);
  orbs *= orb1;
  orbs *= orb1;
  term._mat[imat] = Matrix(Ops::Fock,orbs,1);
  term._sumorbs.insert(orb1);
  term._orbs.insert(orb1);
  sum += term;
  if (multiref) {
    // f_TU \gamma^T_U
    term = *this;
    if (replaceE0act) {
      orbs.clear();
      Orbital torb = term.freeorbname(Orbital::Act),
              uorb = term.freeorbname(Orbital::Act);
      Electron el = term.nextelectron();
      torb.setel(el);
      uorb.setel(el);
      orbs *= torb;
      orbs *= uorb;
      term._mat[imat] = Matrix(Ops::Fock,orbs,1);
      term._mat.push_back(Matrix(Ops::DensM,orbs,1));
      Product< SQOpT::Gender > cran;
      cran *= SQOpT::Creator;
      cran *= SQOpT::Annihilator;
      term._mat.back().set_cran(cran);
      term._sumorbs.insert(torb);
      term._orbs.insert(torb);
      term._sumorbs.insert(uorb);
      term._orbs.insert(uorb);
    } else {
      std::string name = term._mat[imat].name();
      IL::add2name(name,"act",true,true);
      term._mat[imat].set_name(name);
    }
    sum += term;
  }
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
  _matconnectionsset = true;
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
    // search for orbitals in summations
    it1 = _sumorbs.find(_kProd[i].orb1());
    it2 = _sumorbs.find(_kProd[i].orb2());
    Orbital 
      orb1 = _kProd[i].orb1(),
      orb2 = _kProd[i].orb2();
    bool insum = false;
    if ( it2 != _sumorbs.end() ) { // found orb2
      insum = true;
    } else if ( it1 != _sumorbs.end() ) { // found orb1
      insum = true;
      std::swap(it1,it2);
      std::swap(orb1,orb2);
    }
    if ( insum ) { // found in sum
      _sumorbs.erase(it2); // delete summation over orb2
      it2 = _orbs.find(orb2);
      if ( it2 == _orbs.end() )
        error("Strange, orbital not found in _orbs","Term::reduceTerm");
      _orbs.erase(it2);
      _kProd.erase(_kProd.begin()+i); // delete the Kronecker
      --i;
      this->replace(orb2,orb1);
    }
  }
}
void Term::reduceElectronsInTerm()
{
  assert( _kProd.size() > 0 );
  TOrbSet::iterator it1, it2;
  _foreach_cauto(Product<Kronecker>,ik,_kProd){
    Spin 
      spin1 = ik->orb1().spin(),
      spin2 = ik->orb2().spin();
    it1 = Q2::findSpin<TOrbSet>(_sumorbs,spin1);
    it2 = Q2::findSpin<TOrbSet>(_sumorbs,spin2);
    if ( it1 != _sumorbs.end() ) { // found spin1
      std::swap(spin1,spin2);
    } 
    this->replace(spin2,spin1);
  }
}
void Term::krons2mats()
{
  // try to guess the proper order of the indices in the Kroneckers
  // creators and annihilators orbitals
  Product<Orbital> creators, annihilators;
  for ( const Matrix& mat: _mat ){
    creators *= mat.crobs();
    annihilators *= mat.crobs(true);
  }
  _foreach_cauto(Product<Kronecker>,ik,_kProd){
    _mat.push_back(Matrix(*ik,ik->is_ordered(creators,annihilators)));
  }
  _kProd = Product<Kronecker>();
}

void Term::replace(Orbital orb1, Orbital orb2, bool smart)
{
  Q2::replace(_opProd,orb1,orb2,smart);
  Q2::replace(_mat,orb1,orb2,smart);
  Q2::replace(_kProd,orb1,orb2,smart);
  Q2::replace(_sumorbs,orb1,orb2,smart);
  Q2::replace(_orbs,orb1,orb2,smart);
}
void Term::replace(Spin spin1, Spin spin2, bool smart)
{
  Q2::replace(_opProd,spin1,spin2,smart);
  Q2::replace(_mat,spin1,spin2,smart);
  Q2::replace(_kProd,spin1,spin2,smart);
  Q2::replace(_sumorbs,spin1,spin2,smart);
  Q2::replace(_orbs,spin1,spin2,smart);
}

static bool matisnone(const Matrix& mat)
{ return (mat.type() == Ops::None); }
void Term::combineMats(Matrix*& pMat, Product< Matrix >::iterator& it, const Set< uint >& norbs)
{
  if ( norbs.size() == it->orbitals().size() ){
    // all orbitals have to be removed
    it = _mat.erase(it);
  } else if ( pMat == 0 ){
    if ( norbs.size() > 0 ){
      // some of the orbitals have to be removed
      Matrix mat(it->type(),Product<Orbital>(),0);
      mat.combine(*it,norbs);
      *it = mat;
    }
    pMat = &(*it);
    ++it;
  } else {
    pMat->combine(*it,norbs);
    it = _mat.erase(it);
  }
}

void Term::deleteNoneMats(bool unite_exc0)
{ 
  _mat.erase(std::remove_if(_mat.begin(),_mat.end(),matisnone),_mat.end()); 
  if (unite_exc0) {
    Matrix 
      * pExc0 = 0,
      * pDexc0 = 0;
    Product<Matrix>::iterator it = _mat.begin();
    while ( it != _mat.end() ){
      if ( !InSet(it->type(),Ops::Deexc0,Ops::Exc0) ){
        ++it;
        continue;
      }
      if ( it->internal() ) {
        // the matrix is internal (comming from something like \sum_\mu T_\mu \tau_\mu)
        it = _mat.erase(it);
        continue;
      }
      Set<uint> norbs;
      const Product<Orbital> & orbs = it->orbitals();
      for ( uint io = 0; io < orbs.size(); ++io ){
        if ( _sumorbs.count(orbs[io]) ){
          // orbital is not really external
          norbs.insert(io);
        }
      }
      if ( it->type() == Ops::Deexc0 ) {
        combineMats(pDexc0,it,norbs);
      } else {
        combineMats(pExc0,it,norbs);
      }
    }
  }
  matrixkind();
}
bool Term::brilloin() const
{
  for ( Product<Matrix>::const_iterator it = _mat.begin(); it != _mat.end(); ++it )
    if ( it->type() == Ops::Fock && 
         (( it->orbitals()[0].type() == Orbital::Occ && it->orbitals()[1].type() == Orbital::Virt ) ||
          ( it->orbitals()[1].type() == Orbital::Occ && it->orbitals()[0].type() == Orbital::Virt )) )
      return true;
  return false;
}

void Term::matrixkind()
{
  short exccl,intlines=0,intvirt=0;
  for (unsigned int i=0; i<_mat.size(); i++) {
    // excitation class of operator (= #electrons = #orbitals/2)
    exccl=_mat[i].orbitals().size()/2;
    for (unsigned int j=0; j<_mat[i].orbitals().size(); j++) {
      // internal lines (have to be sumed up - search in _orbs)
      if ( _orbs.count(_mat[i].orbitals()[j])) {
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
        // electrons have to be swapped in the matrix and all corresponding kroneckers
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
TermSum Term::expand_antisym()
{
  TermSum sum;
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
TermSum Term::oneel2fock(bool multiref)
{
  TermSum sum;
  for ( uint i = 0; i < _mat.size(); ++i ){
    if (_mat[i].type() == Ops::OneEl){
      // replace "h" by "f - integrals"
      sum = this->change2fock(i,multiref);
      // can transform only one matrix per call
      return sum;
    }
  }
  sum += *this;
  return sum;
}
TermSum Term::replaceE0(bool multiref)
{
  std::string e0name = Input::sPars["hamilton"]["E0"];
  bool replaceE0act = (Input::iPars["prog"]["replacee0"] > 1);
  TermSum sum;
  for ( uint i = 0; i < _mat.size(); ++i ){
    if ( _mat[i].name() == e0name ){
      // replace "E^0" by "2\sum_i f_ii + \sum_tu f_tu \gamma^t_u"
      sum = this->replaceE0byfock(i,multiref,replaceE0act);
      // can transform only one matrix per call
      return sum;
    }
  }
  sum += *this;
  return sum;
}

bool Term::has_nonsingldm() const
{
  for ( uint i = 0; i < _mat.size(); ++i )
    if (_mat[i].nonsingldm()) return true;
  return false;
}
TermSum Term::dm2singlet()
{
  TermSum sum;
  for ( uint i = 0; i < _mat.size(); ++i ){
    if (_mat[i].nonsingldm()){
      this->dmelectrons(i);
      Matrix mat(_mat[i]);
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
bool Term::dmelectrons(uint imat)
{
  const Matrix & dm = _mat[imat];
  assert( dm.type() == Ops::DensM );
  const Product<Orbital>& orbs(dm.orbitals());
  assert( orbs.size()%2 == 0 );
  std::map<Electron,int> els;
  const Product<SQOpT::Gender>& cran = dm.get_cran();
  assert( cran.size() == orbs.size() );
  for ( uint i = 0; i < orbs.size(); ++i ){
    // every electron has to be presented by one creator and one annihilator!
    // otherwise, if just one creator-annihilator pair has different electrons - kronecker (singlet!)
    // call error if more than one...
    if ( cran[i] == SQOpT::Creator ){
      els[orbs[i].spin().el()] += 1;
    } else if ( cran[i] == SQOpT::Annihilator ){
      els[orbs[i].spin().el()] -= 1;
    }
  }
  bool replace_el = false;
  std::map<Electron,int>::const_iterator itel;
  Electron el1 = 0, el2 = 0;
  _foreach(itel,els){
    if ( itel->second > 0 && el1 == 0 ) {
      replace_el = true;
      el1 = itel->first;
    } else if ( itel->second < 0 && el2 == 0 ) {
      replace_el = true;
      el2 = itel->first;
    } else if ( itel->second != 0 ) {
      error("More than one annihilator-creator pair has different electrons","Term::dmwickstheorem");
//      return TermSum();
    }
  }
  if ( replace_el ) {
    assert( el1 > 0 && el2 > 0 );
    Spin
      spin1 = Spin(el1),
      spin2 = Spin(el2);
    TOrbSet::const_iterator
      it1 = Q2::findSpin<TOrbSet>(_sumorbs,spin1);
    if ( it1 != _sumorbs.end() ) { // found spin1
      std::swap(spin1,spin2);
    }
    warning("Assuming singlet excitations in density matrices " << *this);
    warning("Replace electrons" << spin1 << " " << spin2);
    this->replace(spin2,spin1);
  }
  return replace_el;
}
bool Term::has_generalindices() const
{
  _foreach_cauto(TOrbSet,it,_orbs)
    if (it->type() == Orbital::GenT) return true;
  return false;
}
TermSum Term::removegeneralindices()
{
  bool active = (Input::iPars["prog"]["multiref"] > 0);
  TermSum sum;
  this->set_lastorbs();
  Term tt(*this);
  _foreach_auto(TOrbSet,it,_orbs){
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
  TOrbSet peo(extindx());
  TOrbSet peo1(peo);
  // internal indices
  TOrbSet po(_sumorbs); 
  _nocc=_nintloops=_nloops=0;
  Spin nospin(Spin::No);
  bool nonconserve = true, already_done = false;
  Orbital orb,orb1;
  Matrix::Spinsym spinsym=Matrix::Singlet;
  lui ithis = 0;
  long int ipos, iorb = 0;
  TOrbSet::iterator it;
  bool samespinsym;
  uint number_of_densmat;
  bool dm_warning = false;
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
    number_of_densmat = 0;
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
      if ( _mat[ithis].type() == Ops::DensM ) ++number_of_densmat; 
      iorb = _mat[ithis].iorbel(ipos);
      if ( iorb < 0 ) break;
      orb1 = _mat[ithis].orbitals()[iorb];
      if (nonconserve){
        // test whether still in sets
        already_done = ( peo.count(orb1) == 0 && po.count(orb1) == 0 );
      }
    } while (orb1!=orb && !already_done && _orbs.count(orb1));
    if (samespinsym) {
      if ( notfake ) _prefac*=2;
      // reduce by a factor of two for every density matrix in the loop
      for (uint i = 0; i < number_of_densmat; ++i ){
        _prefac /= 2;
      }
      if ( number_of_densmat > 1 ) dm_warning = true;
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
    TOrbSet orbs;
    for ( TOrbSet::iterator it = _orbs.begin(); it != _orbs.end(); ++it ){
      orb = *it;
      orb.setspin(nospin);
      orbs.insert(orb);
    }
    _orbs = orbs;
    orbs.clear();
    for ( TOrbSet::iterator it =_sumorbs.begin(); it != _sumorbs.end(); ++it ){
      orb = *it;
      orb.setspin(nospin);
      orbs.insert(orb);
    }
    _sumorbs = orbs;
    if ( dm_warning ) {
      // It is probably relevant for explicitly inserted terms only (e.g. by creating fock from h)...
      warning("spin summation in " << *this << " relies on \\gamma^t\\alpha_u\\alpha = \\gamma^t\\beta_u\\beta = \\half \\gamma^t_u");
    }
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
void Term::printdiag(Output* pout) const
{
  TsPar & diag = Input::sPars["diag"];
  short verb = Input::iPars["diag"]["printnames"]; // print names of operators and indices...
  
  *(pout->pout) << diag["bdiag"] << std::endl;
  
// "numbers"
  Product<Matrix> nums;
  for ( Product<Matrix>::const_iterator it = _mat.begin(); it != _mat.end(); ++it)
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
  for ( Product<Matrix>::const_iterator it = _mat.begin(); it != _mat.end(); ++it, ++im ){
    switch ( it->type() ){
      case Ops::Fock:
        if ( verb > 1 )
          *(pout->pout) << diag["fockn"] << "{t" << im << "}" << std::endl;
        else
          *(pout->pout) << diag["fock"] << "{t" << im << "}" << std::endl;
        break;
      case Ops::OneEl:
        error("Line for one-El.operator not defined","printdiag");
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
  for ( Product<Matrix>::const_iterator it = _mat.begin(); it != _mat.end(); ++it, ++im ){
    assert( it->conlines().size() == it->orbitals().size() );
    for ( lui iorb = 0; iorb < it->conlines().size(); ++iorb ){
      if (orbs.insert(it->orbitals()[iorb]).second) {// is new
        const ConLine & cl = it->conline(iorb);
        if ( iorb%2 == 0 ){
          assert( (it->diaglevel() < _mat[cl.imat].diaglevel() && it->orbitals()[iorb].type() == Orbital::Virt) ||
                  (it->diaglevel() > _mat[cl.imat].diaglevel() && it->orbitals()[iorb].type() == Orbital::Occ) );
          assert( cl.idx%2 == 1 );
          if ( verb > 2 )
            *(pout->pout) << diag["conline"]  << "[$" << it->orbitals()[iorb] << "$]{t" << im << int(iorb/2+1) << "}{t" << cl.imat << int(cl.idx/2+1) << "}" << std::endl;
          else
            *(pout->pout) << diag["conline"] << "{t" << im << int(iorb/2+1) << "}{t" << cl.imat << int(cl.idx/2+1) << "}" << std::endl;
        } else {
          assert( (it->diaglevel() > _mat[cl.imat].diaglevel() && it->orbitals()[iorb].type() == Orbital::Virt) ||
                  (it->diaglevel() < _mat[cl.imat].diaglevel() && it->orbitals()[iorb].type() == Orbital::Occ) );
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

Orbital Term::freeorbname(Orbital::Type type, bool spinfree)
{
  TsPar& orbs = Input::sPars["syntax"];
  Spin::Type spin = Spin::Gen;
  if (spinfree) {
    spin = Spin::No;
  } else {
    bool spinintegr = Input::iPars["prog"]["spinintegr"];
    if (spinintegr) spin = Spin::GenS;
  }
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
  do {
    if (lastorb.size()==0) {
      lastorb=ip_orbs->at(0);
      _lastorb[type]=Orbital(lastorb,type,spin);
    } else {
      bool lastorbset = false;
      for ( int i=lastorb.size()-1; i>=0 && !lastorbset; --i) {
        indx=ip_orbs->find(lastorb[i]);
        if(indx==std::string::npos)
          error("Something wrong with orbitals","Term::freeorbname");
        else if (indx<ip_orbs->size()-1) {// not the last possible orbital index
          lastorb[i]=ip_orbs->at(indx+1);
          _lastorb[type]=Orbital(lastorb,type,spin);
          lastorbset = true;
        } else {
          lastorb[i]=ip_orbs->at(0);//set to first letter
        }
      }
      if (!lastorbset) {
        //add one orbital index more
        lastorb+=ip_orbs->at(0);
        _lastorb[type]=Orbital(lastorb,type,spin);
      }
    }
  } while (_lastorb[type].is_in_set(_orbs));
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
  _foreach_cauto(TOrbSet,it,_orbs){
    if (_lastorb[it->type()].name().size() == 0 || _lastorb[it->type()] < *it) 
      _lastorb[it->type()] = *it;
  }
}

Electron Term::nextelectron()
{
  if (_lastel == 0){
    // set to correct value (if there are already electrons )
    _foreach_cauto(TOrbSet,its,_orbs){
      Electron el = its->getel();
      if (el > _lastel)
        _lastel = el;
    }
  }
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

Orbital Term::orb(uint iorb) const
{
  if (iorb >= _orbs.size()) return Orbital();
  TOrbSet::const_iterator it = _orbs.begin();
  std::advance(it,iorb);
  return *it;
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
template <class T, class Q>
Return Q2::replace(T &orb, Q orb1, Q orb2, bool smart)
{ return orb.replace(orb1,orb2,smart);
  //if (orb == orb1) orb=orb2; 
}
template <class T>
Return Q2::replace(Matrix &mat, T orb1, T orb2, bool smart)
{ return mat.replace(orb1,orb2,smart); }
template <class T>
Return Q2::replace(Kronecker& kron, T orb1, T orb2, bool smart)
{ return kron.replace(orb1,orb2,smart); }
template < class T >
typename T::iterator Q2::findSpin(T orbs, const Spin& spin){
  typename T::iterator it;
  for (it=orbs.begin(); it != orbs.end() && it->spin() != spin; ++it);
  return it;
}

