#include "lexic.h"

Lelem::Lelem(std::string name, Lelem::Lex lex, Conn conn) : 
_name(name), _lex(lex), _conn(conn), _expandedbra(false){}
std::string Lelem::name() const
{ return _name; }
Lelem::Lex Lelem::lex() const
{ return _lex; }
Lelem::Conn Lelem::conn() const
{ return _conn; }
bool Lelem::expandedbra() const
{ return _expandedbra; }
Lelem Lelem::braexpanded() const
{ 
  Lelem result(*this);
  result._expandedbra=true;
  return result; 
}

bool Lelem::operator==(const Lelem& lel) const
{ return _lex==lel._lex && _name==lel._name; }


lui LEquation::closbrack(const Product< Lelem >& eqn, lui ipos) const
{
  lui i,ipos1=ipos;
  Lelem::Lex rk=Lelem::RPar,lk=eqn[ipos].lex();
  if (lk==Lelem::Bra)
    rk=Lelem::Ket;
  else if (lk==Lelem::LPar)
    rk=Lelem::RPar;
  else
    error("Not a bracket!","Lexic::closbrack"); 
  int nk=1;
  for ( i=ipos+1;i<eqn.size();++i)
  {
    if (eqn[i].lex()==lk) 
      ++nk; // count number of "("
    else if (eqn[i].lex()==rk) 
    {
      --nk; // count number of ")"
      if (nk==0) 
      {
        ipos1=i;
        break;
      }
    }
  }
  if ( nk != 0 ) 
    error("Number of brackets is incosistent: "+any2str(nk),"Lexic::closbrack"); 
  return ipos1;
}
lui LEquation::openbrack(const Product< Lelem >& eqn, lui ipos) const
{
  lui ipos1=ipos;
  Lelem::Lex lk=Lelem::LPar,rk=eqn[ipos].lex();
  if (rk==Lelem::Ket)
    lk=Lelem::Bra;
  else if (rk==Lelem::RPar)
    lk=Lelem::LPar;
  else
    error("Not a bracket!","Lexic::openbrack"); 
  int nk=-1;
  for ( long int i=ipos-1;i>=0;--i) {
    if (eqn[i].lex()==lk) {
      ++nk; // count number of "("
      if (nk==0) {
        ipos1=i;
        break;
      }
    }
    else if (eqn[i].lex()==rk) 
      --nk; // count number of ")"
  }
  if ( nk != 0 ) 
    error("Number of brackets is incosistent: "+any2str(nk),"Lexic::openbrack"); 
  return ipos1;
}
Product< long int > LEquation::addconnections(const Product< Lelem >& aterm, lui beg, lui end) const
{
  Product<long int> connection;
  if (aterm[end].conn()==Lelem::Normal) return connection;
  for (lui i=beg+1; i<end; i++)
    if (aterm[i].lex()==Lelem::Oper)
    {
      if(aterm[end].conn()==Lelem::Connect)
        connection *= i+1;
      else
        connection *= -(i+1);
    }
  return connection;
}

bool LEquation::extractit()
{
  // expand custom operators
  _eqn = expandnewops(_eqn);
  // expand parentheses
  _eqn=expandeqn(_eqn,_connections);
  // remove redundant connections
  for (lui i=0; i<_connections.size();i++)
    for (lui j=0; j<_connections[i].size();j++)
      if (InSet(_eqn[abs(_connections[i][j])-1].lex(), Lelem::Num,Lelem::Frac)) //connection to a number
        _connections[i].erase(_connections[i].begin()+j);
  for (lui i=0; i<_connections.size();i++)
    if (_connections[i].size()<2) //smaller than two elements "connected"
      _connections.erase(_connections.begin()+i);
  if (_connections.size()==0) return true;
  for (lui i=0; i<_connections.size()-1;i++)
    for (lui j=i+1; j<_connections.size();j++)
      if (_connections[i]==_connections[j]) //same connection
      {
        _connections.erase(_connections.begin()+j);
        --j;
      }
  for (lui k = 0; k<_connections.size();k++)
    _xout2("final Connection " << k << ": " << _connections[k] << std::endl);
  return true;
}
Product< Lelem > LEquation::expandnewops(const Product< Lelem >& eqn) const
{
  Product<Lelem> result;
  for (lui i = 0; i < eqn.size(); ++i ){
    result *= eqn[i];
    lui j = result.size()-1;
    do { 
      if ( result[j].lex() == Lelem::Oper ) {
        NewOpMap::const_iterator itnewop = _newops.find(result[j].name());
        if ( itnewop != _newops.end() ) {
          const Product<Lelem>& nop = itnewop->second;
          result.erase(result.begin()+j);
          result.insert(result.begin()+j,nop.begin(),nop.end());
        }
      }
      ++j;
    } while( j < result.size() && result.size() < lui(Input::iPars["prog"]["maxlels"]) );
    if ( result.size() >= lui(Input::iPars["prog"]["maxlels"]) )
      error("Equation is too long. Cyclic dependencies in \\newoperator?");
  }
  return result;
}

lui LEquation::elem(const Product< Lelem >& aterm, lui beg, bool bk) const
{
  lui i, end, nk=0;
  bool braket=false;
  
  for ( i=beg; i<aterm.size() ; i++ ) {
    if (aterm[i].lex()==Lelem::LPar) ++nk;
    if (aterm[i].lex()==Lelem::RPar) --nk;
    if (bk) {
      if (aterm[i].lex()==Lelem::Bra) braket=true;
      if (aterm[i].lex()==Lelem::Ket) braket=false;
    }
    if (InSet(aterm[i].lex(), Lelem::Plus,Lelem::Minus) && i!= beg && nk==0 && !braket ) 
      break; // end of term
  }
  if (nk!=0 || braket) error("Check input, Mismatch in parentheses or bra/ket","Lexic::term");
  if (i==beg || i==0 ) 
    end=0;
  else 
    end=i-1;
  return end;
}

Product< Lelem > LEquation::expandeqn(const Product< Lelem >& eqn, std::vector< Product<long int> > & connections) const
{
  Product<Lelem> result=eqn, res;
  Product<long int>  connect;
  std::vector< Product<long int> > con,con1;
  lui beg=0, end, i,j,conbeg,lastpos;
  
  while (!expanded(result)) {
    res=result;
    result=Product<Lelem>();
    beg=0;
    conbeg=0;
    con=connections;
    connections=std::vector< Product<long int> >();
    while ((end=term(res,beg))!=0) {
      std::vector< Product<long int> > con1;
      // shift connections
      for (i=conbeg;i<con.size();i++) {
        if (abs(con[i].front())>(long int)beg && abs(con[i].back())-2<(long int)end)
          for (j=0;j<con[i].size();j++)
            if (con[i][j]>0)
              connect*=con[i][j]-beg;
            else
              connect*=con[i][j]+beg;
        else
          break;
        con1.push_back(connect);
        connect=Product<long int>();
        ++conbeg;
      }
      // construct a term as a subvector and expand it:
      lastpos=result.size();
      result*=expandterm(res.subprod(beg,end),con1);
      // shift connections back
      for (i=0;i<con1.size();i++) {
        for (j=0;j<con1[i].size();j++)
          if (con1[i][j]>0)
            connect*=con1[i][j]+lastpos;
          else
            connect*=con1[i][j]-lastpos;
        connections.push_back(connect);
        connect=Product<long int>();
      }
      beg=end+1;
    }
  }
  return result;
}

Product< Lelem > LEquation::expandterm(const Product< Lelem >& aterm, std::vector< Product<long int> > & connections) const
{
  lui i;
//  std::cout << "Term: " << aterm << std::endl;
//  for (unsigned int k = 0; k<connections.size();k++)
//    std::cout << "in connection " << k << ": " << connections[k] << std::endl;
  for (i=0; i<aterm.size(); i++)
    if (aterm[i].lex()==Lelem::LPar || (aterm[i].lex()==Lelem::Bra && !aterm[i].expandedbra())) 
      return expandpar(aterm,i,connections);
  return aterm;
}
Product< Lelem > LEquation::expandpar(const Product< Lelem >& aterm, lui beg, 
                                   std::vector< Product< long int > >& connections) const
{ // e.g., aterm=-a(b+c)d
  lui i=0, end,ipos,start=0,ijcon,iposres,lenb,ipar;//lena,
  end=closbrack(aterm,beg);
  Product<Lelem> result, inpar=aterm.subprod(beg+1,end-1);
  std::vector< Product<long int> > con(connections);
  con.push_back(addconnections(aterm,beg,end));
  connections=std::vector< Product<long int> >();
  Product<long int> connect;
  bool coninterm,coninpar;
  if (aterm[beg].lex()==Lelem::Bra) 
    ipar=1;
  else
    ipar=0;
  // start of term (without sign)
  if (InSet(aterm[0].lex(), Lelem::Plus,Lelem::Minus)) start=1;
  while (i<inpar.size()) {
    ipos=elem(inpar,i);
    //add sign
    if (inpar[i].lex()==Lelem::Minus) {// minus in parenthesis
      if (aterm[0].lex()==Lelem::Minus) 
        result *= Lelem("",Lelem::Plus); // -1*-1 == +1
      else
        result *= inpar[i]; // +1*-1 == -1
      ++i;
    } else if (inpar[i].lex()==Lelem::Plus) {// plus in parenthesis
      if (aterm[0].lex()==Lelem::Minus) 
        result *= aterm[0]; // -1*+1 == -1
      else
        result *= inpar[i]; // +1*+1 == +1
      ++i;
    } else if (start>0) {// term has a sign, but no sign by expression in parenthesis
      result *= aterm[0];
    }
    iposres=result.size()-start;
    // add "a" (without sign)
    if (beg>start) result *= aterm.subprod(start,beg-1);
    if (aterm[beg].lex()==Lelem::Bra) result *= aterm[beg].braexpanded();
    // add "b" (or "c")
    result *= inpar.subprod(i,ipos);
    // add "d"
    if (aterm[end].lex()==Lelem::Ket) result *= aterm[end];
    if (end<aterm.size()-1) result *= aterm.subprod(end+1,aterm.size());
    // handle connections
    lenb=ipos-i+1;
    for (lui k=0;k<con.size();k++) {
      for (lui j=0;j<con[k].size();j++) {
        ijcon=abs(con[k][j]);
        coninterm=(ijcon<=beg)||(ijcon>end+1);// connection in a or d
        coninpar=((ijcon>i+beg+1)&&(ijcon<ipos+beg+3));// connection in b (or c)
        coninterm=coninterm||coninpar;
        if (coninterm) { // change ijcon
          if (coninpar)
            ijcon+=iposres-i-1+ipar;
          else if (ijcon<=beg)
            ijcon+=iposres;
          else if (ijcon>end+1)
            ijcon+=iposres+lenb-(end-beg+1)+ipar*2;
          if (con[k][j]>0)
            connect*=ijcon;
          else
            connect*=-ijcon;
        }
      }
      if (connect.size()>0) {
        connections.push_back(connect);
        connect=Product<long int>();
      }
    }
    i=ipos+1;
  }
  return result;
}
bool LEquation::expanded(const Product< Lelem >& eqn) const
{ 
  for (lui i=0; i<eqn.size(); i++) {
    if (eqn[i].lex()==Lelem::Bra && !eqn[i].expandedbra())
      return false;
    if (eqn[i].lex()==Lelem::LPar)
      return false;
  }
  return true;
}
bool LEquation::do_sumterms(bool excopsonly )
{
  lui beg=0;
  bool plus=true, bra=false, ket=false;
  if (!expanded(_eqn))
    error("Expand the lexic equation first!","Lexic::do_sumterms");
  Term term;
  reset_term(term);
  Product<long int> indxoperterm;
  for (lui i=0; i<_eqn.size(); i++) {
    if(InSet(_eqn[i].lex(), Lelem::Bra,Lelem::Ket)) { // handle bra/ket
      if (_eqn[i].lex()==Lelem::Bra) {
        if (bra) 
          error("Cannot handle two BRAs in one term yet...");
        else
          bra=true;
      } else if (ket)
        error("Cannot handle two KETs in one term yet...");
      else
        ket=true;
      term*=handle_braket(_eqn[i],term,excopsonly);
      indxoperterm.push_back(i+1);
    } else if (InSet(_eqn[i].lex(), Lelem::Minus,Lelem::Plus)) { // add current term and save the sign of the next term
      bra=ket=false; // reset bra and ket variables
      if (!excopsonly) {
        if(i>0) addterm(term,plus,beg,i-1,indxoperterm,excopsonly);
        plus=(_eqn[i].lex()==Lelem::Plus);
        beg=i+1;
        reset_term(term);
        indxoperterm.clear();
      }
    } else if (InSet(_eqn[i].lex(), Lelem::Frac,Lelem::Num)) { // add prefactor
      term*=handle_factor(_eqn[i]);
    } else if (_eqn[i].lex()==Lelem::Oper) { // handle Operator
      term*=handle_operator(_eqn[i],term,excopsonly);
      indxoperterm.push_back(i+1);
    } else if (_eqn[i].lex()==Lelem::Sum) { // handle \sum
      if (!excopsonly)
        _sumsterm *= _eqn[i];
    } else if (_eqn[i].lex()==Lelem::Param) { // handle Parameter
      if (!excopsonly)
        _paramterm*=_eqn[i];
    } else if (_eqn[i].lex()==Lelem::Perm) { // handle Permutation
      if (!excopsonly)
        term *= handle_permutation(_eqn[i]);
    } else if (_eqn[i].lex()==Lelem::Times) { // handle Multiplication
      // don't do anything
    } else if (_eqn[i].lex()==Lelem::Div) { // handle Division
      error("Sorry, cannot handle Division!","Lexic::do_sumterms");
    } else {
      xout << _eqn[i] << std::endl;
      error(_eqn[i].name()+" is not implemented yet...","Lexic::do_sumterms");
    }
  }
  // add last term
  if(_eqn.size()>0) addterm(term,plus,beg,_eqn.size()-1,indxoperterm,excopsonly);
  return true;
}
void LEquation::reset_term(Term& term) const
{
  term=Term();
  term.addmatrix(Matrices());
  if (_excops.size()>0) {
    for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
      Orbital::Type ot = static_cast<Orbital::Type>(it);
      _foreach_cauto( LExcitationMap,itex,_excops) {
        TOrb4Type::const_iterator ito4e = itex->second._orbs4excops.find(ot);
        if ( ito4e != itex->second._orbs4excops.end() ) {
          term.set_lastorb(ito4e->second);
          break;
        }
      }
    }
  }
}

void LEquation::addterm(Term& term, bool plus, lui beg, lui end, 
                     Product<long int > const & indxoperterm, bool excopsonly)
{
  double minfac = Input::fPars["prog"]["minfac"];
  if( excopsonly || term.term_is_0(minfac)) {
    //reset parameter-info
    handle_parameters(term,true);
    return; // dont add zero term
  }
  // handle sums in term
  for ( uint i = 0; i < _sumsterm.size(); ++i ) handle_sum(_sumsterm[i],term);
  // reset sums information
  _sumsterm=Product<Lelem>();
  // handle parameters
  handle_parameters(term);
//   if( term.term_is_0(minfac) ) {
//     //reset parameter-info
//     handle_parameters(term,true);
//     return; // dont add zero term
//   }
  //add connections to term
  Product<long int> connect;
  long int ipos;
  for (unsigned long int i=0;i<_connections.size();i++) {
    if (abs(_connections[i].front())>(long int)beg && abs(_connections[i].back())-2<(long int)end) {
      for (unsigned long int j=0;j<_connections[i].size();j++) {
        ipos=indxoperterm.find(abs(_connections[i][j]));
        if (ipos<0)
          error("Connected operator is not in indxoperterm","Lexic::addterm");
        if (_connections[i][j]>0)
          connect*=ipos+2;
        else
          connect*=-ipos-2;
      }
      _xout2("Connections in Term #" << _sumterms.size()+1 << ": " <<connect<<std::endl);
      term.addconnection(connect);
      connect=Product<long int>();
    }
  }
  // validate term
  term.term_is_valid();
  // add term
  if(plus) 
    _sumterms += term;
  else
    _sumterms -= term;
}

Oper LEquation::handle_braket(const Lelem& lel, Term& term, bool excopsonly)
{
  const TParArray& refs = Input::aPars["syntax"]["ref"];
  const TParArray& csfs = Input::aPars["syntax"]["csf"];
  std::string lelnam=lel.name();
  if (InSet(lelnam, refs))
    return Oper(); // Reference, blank operator
    
  lui 
    ibeg = 0,
    iend = IL::nextwordpos(lelnam,ibeg,false);
  if (InSet(lelnam.substr(ibeg,iend-ibeg), csfs)){
    return handle_explexcitation(term,lelnam.substr(iend),(lel.lex()==Lelem::Bra),excopsonly);
  } else {
    int lm = 0;
    return handle_excitation(term,lelnam,(lel.lex()==Lelem::Bra),lm,excopsonly);
  }
}
Oper LEquation::handle_explexcitation(Term& term, const std::string& name, bool dg, bool excopsonly, bool phi)
{
  lui ipos, ipos1;
  short excl;
  Product<Orbital> occs, virts;
  if ( name[0] != '_' && name[0] != '^' )
    error("Doesn't start with _ or ^ :"+name,"Lexic::handle_explexcitation");
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  lui up, down;
  down=IL::lexfind(name,"_");
  up=IL::lexfind(name,"^");
  if (phi){ // in phi we have an opposite order of indices!
    std::swap(up,down);
  }
  if (up!=std::string::npos && up!=name.size()-1){
    ipos = up;
    ipos = IL::skip(name,ipos,"{}_^ ");
    while ( (ipos < down) == (up < down) && (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
      if (ipos < down && ipos1 >= down) ipos1 = down;
      Orbital orb(IL::plainname(name.substr(ipos,ipos1-ipos)),spintype);
      occs *= orb;
      if ( orb.type() == Orbital::Virt ) 
        warning("Do you really want to have orbital " << orb << " as occupied?");
      ipos = IL::skip(name,ipos1,"{}_^ ");
    }
  }
  if (down!=std::string::npos && down!=name.size()-1){
    ipos = down;
    ipos = IL::skip(name,ipos,"{}_^ ");
    while ( (ipos < up) == (down < up) && (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
      if (ipos < up && ipos1 >= up) ipos1 = up;
      Orbital orb(IL::plainname(name.substr(ipos,ipos1-ipos)),spintype);
      virts *= orb;
      if ( orb.type() == Orbital::Occ ) warning("Do you really want to have orbital " << orb << " as virtual?");
      ipos = IL::skip(name,ipos1,"{}_^ ");
    }
  }
  
  excl = occs.size();
  // set lastorb (if smaller)
  for ( uint i = 0; i < occs.size(); ++i )
    term.set_lastorb(Orbital(occs[i].letname(),spintype),true);
  for ( uint i = 0; i < virts.size(); ++i )
    term.set_lastorb(Orbital(virts[i].letname(),spintype),true);
  //make sure that we haven't used these orbital names already
  correct_orbs(term,occs);
  correct_orbs(term,virts);
  int lmelec = virts.size()-occs.size();
  if (excopsonly) return Oper();
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0,excl,occs,virts,"",lmelec,&term);  
  else
    return Oper(Ops::Exc0,excl,occs,virts,"",lmelec,&term);
}
void LEquation::correct_orbs(Term& term, const Product<Orbital>& orbs)
{
  if ( _excops.size() > 0 ){
    bool spinintegr = Input::iPars["prog"]["spinintegr"];
    Spin::Type spintype = Spin::Gen;
    if (spinintegr) spintype = Spin::GenS;
    //make sure that we haven't used these orbital names already
    for ( uint i = 0; i < orbs.size(); ++i ){
      _foreach_auto(LExcitationMap,itex,_excops){
        for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
          Orbital::Type ot = static_cast<Orbital::Type>(it);
          TOrb4Type::iterator ito4e = itex->second._orbs4excops.find(ot);
          if ( ito4e != itex->second._orbs4excops.end() &&
               ito4e->second == Orbital(orbs[i].letname(),spintype) )  // is there, change it
            ito4e->second = term.freeorbname(ot);
        }
      }
    }
  }
}

Oper LEquation::handle_excitation(Term& term, const std::string& name, bool dg, int lmel, bool excopsonly)
{
//   bool multiref = ( Input::iPars["prog"]["multiref"] > 0 );
  std::string newname, nameadd,namedown;
  short excl;
  bool newdg;
  int lmelec;
  std::vector< Product<Orbital::Type> > orbtypes;
  bool updown=handle_namupdown(newname,excl,nameadd,namedown,newdg,lmelec,orbtypes,name);
  dg = (dg != newdg);
  if (lmelec != 0 && lmel != 0 && lmelec != lmel )
    error("Mismatch in non-conserving class in "+name,"LEquation::handle_excitation");
  if (lmelec == 0 ) lmelec = lmel;
  // find excitation class
  if (!updown && excl == 0 && lmelec <= 0)
    error("No excitation class in "+name,"LEquation::handle_excitation");
  LExcitationMap::iterator itex = _excops.find(name);
  TOrb4Type orb4t;
  if ( itex == _excops.end() ) {// first run (don't distinguish terms)
    if (orbtypes.size() > 0){
      for ( uint i = 0; i < orbtypes.size(); ++i ){
        _foreach_cauto(Product<Orbital::Type>,iot,orbtypes[i]){
          if (orb4t.count(*iot) == 0)
            orb4t[*iot] = term.freeorbname(*iot);
        }
      }
    } else {
      orb4t[Orbital::Occ] = term.freeorbname(Orbital::Occ);
      orb4t[Orbital::Virt] = term.freeorbname(Orbital::Virt);
    }
    //TODO: implement Triplet!
    Matrices::Spinsym spinsym = Matrices::Singlet;
    _excops[name] = LExcitationInfo(orb4t,excl,spinsym);
  } else {// may be second run...
    orb4t = itex->second._orbs4excops;
    itex->second._posexcopsterm = term.mat().size(); //set position of this operator in the term
  }
  if (excopsonly) return Oper();
  // create \tau_{excl}
  if (orbtypes.size() == 0){
    if (dg)
      return Oper(Ops::Deexc0,excl,orb4t[Orbital::Occ],orb4t[Orbital::Virt],"",lmelec,&term);  
    else
      return Oper(Ops::Exc0,excl,orb4t[Orbital::Occ],orb4t[Orbital::Virt],"",lmelec,&term);
  } else {
    if (dg)
      return Oper(Ops::Deexc0,excl,orb4t,orbtypes,"",lmelec,&term);  
    else
      return Oper(Ops::Exc0,excl,orb4t,orbtypes,"",lmelec,&term);
  }
}

TFactor LEquation::handle_factor(const Lelem& lel) const
{
  TFactor fac;
  lui ipos=0, ipos1;
  std::string lelnam=lel.name();
  if (lel.lex()==Lelem::Num) {
    double facd;
    if(!str2num<double>(facd,lelnam,std::dec))
        error("Factor is not a number "+lelnam,"Lexic::handle_factor");
    if ( typeid(TFactor) != typeid(double) ){
      // try to bring it to integer over integer form
      // NOTE: won't work for 1e-2 etc.
      long int 
        denom = std::pow(10,lelnam.size()),
        nom = facd*denom;
      fac = nom;
      fac /= denom;
    } else 
      fac = facd;
  } else {
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos);
    if ( typeid(TFactor) == typeid(double) ){
      double facd;
      if(!str2num<double>(facd,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Numerator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Lexic::handle_factor");
      fac = facd;
    } else {
      // NOTE: won't work for non-integer nominators or denominators 
      long int nom;
      if(!str2num<long int>(nom,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Numerator is not an integer "+lelnam.substr(ipos,ipos1-ipos),"Lexic::handle_factor");
      fac = nom;
    }
    ipos=lelnam.find("/");
    if(ipos==std::string::npos)
      error("Something wrong with frac "+lelnam,"Lexic::handle_factor");
    ++ipos;
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos);
    if ( typeid(TFactor) == typeid(double) ){
      double fac1;
      if(!str2num<double>(fac1,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Denominator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Lexic::handle_factor");
      fac /= fac1;
    } else {
      // NOTE: won't work for non-integer nominators or denominators 
      long int denom;
      if(!str2num<long int>(denom,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Denominator is not an integer "+lelnam.substr(ipos,ipos1-ipos),"Lexic::handle_factor");
      fac /= denom;
    }
  }
  return fac;
}
Oper LEquation::handle_operator(const Lelem& lel, Term& term, bool excopsonly)
{
  const TParArray& bexcops = Input::aPars["syntax"]["bexcop"];
  TsPar& hms = Input::sPars["hamilton"];
  std::string lelnam=lel.name(),name,nameadd,namedown;
  short excl;
  bool dg;
  int lmelec;
  std::vector< Product<Orbital::Type> > orbtypes;
  bool updown=handle_namupdown(name,excl,nameadd,namedown,dg,lmelec,orbtypes,lelnam);
  // parts of Hamilton operator
  if ( InSet(name, hms)) {
    if (excopsonly) return Oper();
    if (updown)
      say("Sub- and superscripts in Hamiltonian will be ignored: "+lelnam);
    if ( name==hms["fock"] ) return Oper(Ops::Fock,true,&term);
    if ( name==hms["oneelop"] ) return Oper(Ops::OneEl,true,&term);
    if ( name==hms["flucpot"] ) return Oper(Ops::FluctP,true,&term);
    if ( name==hms["dflucpot"] ) return Oper(Ops::FluctP,false,&term);
    if ( name==hms["perturbation"] ) return Oper(Ops::XPert,true,&term);
  }
  // excitation class
  if (namedown=="")
    error("No excitation class in operator "+lelnam,"LEquation::handle_operator");
  IL::add2name(name,nameadd); // add nameadd to name (as superscript)
  if (InSet(name.substr(0,4), bexcops)) { // bare excitation operator
    assert(orbtypes.size() == 0);
    return handle_excitation(term,namedown,dg,lmelec,excopsonly);
  }
  if (excopsonly) return Oper();
  if (excl == 0 && lmelec <= 0)
    error("Excitation class in "+lelnam,"LEquation::handle_operator");
  if (orbtypes.size() == 0){
    if(dg)
      return Oper(Ops::Deexc,excl,name,lmelec,&term);
    else
      return Oper(Ops::Exc,excl,name,lmelec,&term);
  } else {
    if(dg)
      return Oper(Ops::Deexc,excl,orbtypes,name,lmelec,&term);
    else
      return Oper(Ops::Exc,excl,orbtypes,name,lmelec,&term);
  }
}
bool LEquation::handle_namupdown(std::string& name, short int& excl, std::string& nameup, std::string& namedown, bool& dg, 
                                int& lmel, std::vector< Product< Orbital::Type > >& orbtypes, const std::string& lelnam) const
{ 
  const TParArray& dgs = Input::aPars["syntax"]["dg"];
  const TParArray& lessmore = Input::aPars["syntax"]["lessmore"];
  bool foundupdown = false;
  lui up, down, iposnam, ipos, ipos1;
  down=IL::lexfind(lelnam,"_");
  up=IL::lexfind(lelnam,"^");
  // last position of name of operator
  iposnam=std::min(up,down);
  assert( iposnam > 0 );
  --iposnam;
  iposnam=std::min(std::size_t(iposnam),lelnam.size()-1);
  name=lelnam.substr(0,iposnam+1);
  std::string upname, downname;
  if (up!=std::string::npos && up!=lelnam.size()-1){
    upname = lelnam.substr(up);
  }  
  if (down!=std::string::npos && down!=lelnam.size()-1){
    downname = lelnam.substr(down);
  }  
  if (up < down){
    upname = upname.substr(0,down-up);
  } else if (down < up) {
    downname = downname.substr(0,up-down);
  }
  
  excl=0;
  nameup="";
  namedown="";
  dg=false;
  lmel=0;
  // less (negative number) or more (positive) electrons after the operator
  uint lmsize = 0;
  for ( TParArray::const_iterator itlm = lessmore.begin(); itlm != lessmore.end(); ++itlm ){
    assert( itlm == lessmore.begin() || lmsize == itlm->size() );
    lmsize = itlm->size();
  }
  if (up!=std::string::npos && up!=lelnam.size()-1){
    // handle superscript
    foundupdown = true;
    lui last = std::string::npos;
    ipos=1;
    ipos=IL::skip(upname,ipos,"{} ");
    while((ipos1=IL::nextwordpos(upname,ipos,false))!=ipos && ipos < last ) {
      std::string nampart(upname.substr(ipos,ipos1-ipos));
      if ( nameup.size() > 0 && upname[ipos]!='}' ) nameup += " ";
      if(InSet(nampart, dgs)) {
        dg=true;
        nameup += dgs.front();
      } else if (upname[ipos]!='}') {
        nameup += nampart;
        if (nampart.size() > lmsize ){
          // is it less/more?
          std::string lmstr(nampart.substr(0,lmsize));
          if ( InSet(lmstr,lessmore) && str2num<int>(lmel,nampart.substr(lmsize),std::dec)){
            // it is a non-conserving operator
            if ( lmstr == lessmore.front() ) lmel = -lmel;
          }
        }
      }
      ipos=ipos1;
    }
  }
  if (down!=std::string::npos && down!=lelnam.size()-1){
    // handle subscript
    foundupdown = true;
    ipos=1;
    ipos=IL::skip(downname,ipos,"{} ");
    ipos1=IL::nextwordpos(downname,ipos,false);
    namedown = downname.substr(ipos,ipos1-ipos);
    if ( str2num<short>(excl,namedown,std::dec) ) {
      if ( handle_orbtypes(orbtypes,downname.substr(ipos1))){
        assert(int(orbtypes[0].size()) == excl);
        assert(int(orbtypes[1].size()) == lmel+excl);
      }
    } else {
      // subscript is not an excitation class, probably rather something like \mu_2
      ipos1 = IL::closbrack(downname,1);
      namedown = downname.substr(ipos,ipos1-ipos);
    }
  }
  return foundupdown;
}
bool LEquation::handle_orbtypes(std::vector< Product< Orbital::Type > >& orbtypes, const std::string& string) const
{
  bool foundorbtypes = false;
  lui up, down, ipos, ipos1;
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  down=IL::lexfind(string,"_");
  up=IL::lexfind(string,"^");
  if (up!=std::string::npos && up!=string.size()-1){
    ipos = up;
    ipos = IL::skip(string,ipos,"{}_^ ");
    Product<Orbital::Type> occtypes;
    while ( (ipos < down) == (up < down) && (ipos1 = IL::nextwordpos(string,ipos,true,false)) != ipos ){//non greedy
      Orbital orb(IL::plainname(string.substr(ipos,ipos1-ipos)),spintype);
      occtypes *= orb.type();
      if ( orb.type() == Orbital::Virt ) xout << "WARNING: Do you really want to have orbital " << orb << " as occupied?" << std::endl;
      ipos = IL::skip(string,ipos1,"{}_^ ");
    }
    if ( occtypes.size() > 0 ) foundorbtypes = true;
    orbtypes.push_back(occtypes);
  }
  if (down!=std::string::npos && down!=string.size()-1){
    ipos = down;
    ipos = IL::skip(string,ipos,"{}_^ ");
    Product<Orbital::Type> virtypes;
    while ( (ipos < up) == (down < up) && (ipos1 = IL::nextwordpos(string,ipos,true,false)) != ipos ){//non greedy
      Orbital orb(IL::plainname(string.substr(ipos,ipos1-ipos)),spintype);
      virtypes *= orb.type();
      if ( orb.type() == Orbital::Occ ) xout << "WARNING: Do you really want to have orbital " << orb << " as virtual?" << std::endl;
      ipos = IL::skip(string,ipos1,"{}_^ ");
    }
    if ( virtypes.size() > 0 ) foundorbtypes = true;
    orbtypes.push_back(virtypes);
  }
  if (!foundorbtypes) orbtypes.clear();
  return foundorbtypes;
}

void LEquation::handle_sum(const Lelem& lel, Term& term) const
{
  lui ipos, ipos1, up, down;
  long int iposexcn;
  std::string lelnam=lel.name(),name;
  down=IL::lexfind(lelnam,"_");
  up=IL::lexfind(lelnam,"^");
  if (up!=std::string::npos)
    error("Sum from-to is not implemented yet: "+lelnam);
  if (down==std::string::npos)
    error("Sum without summation indices: "+lelnam);
  ipos=down+1;
  while (ipos<lelnam.size()) {
    ipos=IL::skip(lelnam,ipos,"{}, ");
    if (ipos==lelnam.size()) break;
    ipos1=IL::nextwordpos(lelnam,ipos);
    name=lelnam.substr(ipos,ipos1-ipos);
    LExcitationMap::const_iterator itex = _excops.find(name);
    if (itex != _excops.end()) {
      iposexcn=itex->second._posexcopsterm;
      if (iposexcn>=0) {
        const Product<Matrices>& mats = term.mat();
        assert( uint(iposexcn) < mats.size() );
        term.addsummation(mats[iposexcn].orbitals());
      } else
        say("Sum is not present in this term: "+lelnam);
    }
    else
      say("No excitation operator which would correspond to summation index "+name);
    ipos=ipos1;
  }
}
Permut LEquation::handle_permutation(const Lelem& lel) const
{
  Product<Orbital> orbs1, orbs2;
  lui ipos=0, ipos1;
  std::string name, lelnam=lel.name();
  ipos=IL::skip(lelnam,ipos,"{} ");
  ipos1=IL::nextwordpos(lelnam,ipos);
  name = lelnam.substr(ipos,ipos1);
  ipos = 0;
  while ( (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
    orbs1 *= Orbital(IL::plainname(name.substr(ipos,ipos1-ipos)));
    ipos = IL::skip(name,ipos1,"{}_^ ");
  }
  ipos=lelnam.find("/");
  if(ipos==std::string::npos)
    error("Something wrong with frac "+lelnam,"Lexic::handle_factor");
  ++ipos;
  ipos=IL::skip(lelnam,ipos,"{} ");
  ipos1=IL::nextwordpos(lelnam,ipos);
  name = lelnam.substr(ipos,ipos1);
  ipos = 0;
  while ( (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
    orbs2 *= Orbital(IL::plainname(name.substr(ipos,ipos1-ipos)));
    ipos = IL::skip(name,ipos1,"{}_^ ");
  }
  return Permut(orbs1,orbs2);
}

void LEquation::handle_parameters(Term& term, bool excopsonly)
{
  if (!excopsonly) {// handle saved parameters
    lui ipos, ipos1, iposnam, up, down;
    int iposexcn;
    std::string lelnam,name,nameadd,excn;
    for (unsigned int i=0; i<_paramterm.size(); i++) {
      lelnam=_paramterm[i].name();
      down=IL::lexfind(lelnam,"_");
      up=IL::lexfind(lelnam,"^");
      // last position of name of parameter
      iposnam=std::min(up,down)-1;
      name=lelnam.substr(0,iposnam+1);
      if (up==std::string::npos || up==lelnam.size()-1) {
        // no superscript
      } else {
        ipos=up+1;
        ipos=IL::skip(lelnam,ipos,"{} ");
        while((ipos1=IL::nextwordpos(lelnam,ipos))!=ipos) {
          if (lelnam[ipos]!='}')
            nameadd+=lelnam.substr(ipos,ipos1-ipos);
          ipos=ipos1;
        }
      } 
      IL::add2name(name,nameadd); // add nameadd to name (as superscript)
      
      // handle subscript
      if (down==std::string::npos || down==lelnam.size()-1) { // no subscript, parameter is a "number"
        term.addmatrix(Matrices(Ops::Number,Product<Orbital>(),0,name));
      } else {
        ipos=down+1;
        ipos=IL::skip(lelnam,ipos,"{} ");
        ipos1=IL::nextwordpos(lelnam,ipos);
        excn=lelnam.substr(ipos,ipos1-ipos);
        LExcitationMap::const_iterator itex = _excops.find(excn);
        if ( itex != _excops.end() ) {
          iposexcn = itex->second._posexcopsterm;
          if (iposexcn>=0) {
            Matrices mat(Ops::Interm,
                         term.mat()[iposexcn].orbitals(),
                         itex->second._exccls, name,itex->second._spinsymexcs);
            term.replacematrix(mat,iposexcn);
          } else
            say("Parameter is not present in this term: "+lelnam);
        } else
        // TODO : add parameters with explicit excitations
          error("Unknown excitation in parameter "+excn);
      }
    }
  }
  // reset all information
  _foreach_auto(LExcitationMap,itex,_excops) {
    itex->second.reset_term_info();
  }
  _paramterm=Product<Lelem>();
}


std::ostream& operator<<(std::ostream& o, const Lelem& lel)
{
  TsPar& commands = Input::sPars["command"];
  if (lel.lex()==Lelem::Bra) 
    o << "< ";
  else if (lel.lex()==Lelem::Ket)
    o << "| ";
  else if (lel.lex()==Lelem::LPar)
    o << "(";
  else if (lel.lex()==Lelem::RPar)
    o << ")";
  else if (lel.lex()==Lelem::Oper)
    o << "\\"<< commands["operator"]<<" ";
  else if (lel.lex()==Lelem::Param)
    o << "\\"<< commands["parameter"]<<" ";
  else if (lel.lex()==Lelem::Perm)
    o << "\\"<< commands["permutation"];
  else if (lel.lex()==Lelem::Plus)
    o << "+";
  else if (lel.lex()==Lelem::Minus)
    o << "-";
  else if (lel.lex()==Lelem::Times)
    o << "*";
  else if (lel.lex()==Lelem::Div)
    o << "/";
  else if (lel.lex()==Lelem::Sum)
    o << "\\"<< commands["sum"];
  
  o << lel.name() << " ";
  if (lel.lex()==Lelem::Bra) 
    o << "| ";
  else if (lel.lex()==Lelem::Ket)
    o << "> ";
  return o;
}

std::ostream& operator<<(std::ostream& o, const LEquation& inp)
{
  for (unsigned long int i=0; i<inp.eqn().size(); i++)
    o << inp.eqn().at(i);
  return o;
}
