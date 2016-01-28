#include "equation.h"

Product< Orbital > LExcitationInfo::orbitals(bool dg) const
{
  if ( dg ) {
    // dagger
    Product<Orbital> orbs;
    assert( (_orbs.size()-std::abs(_lmel))%2 == 0 );
    uint npairs = (_orbs.size()-std::abs(_lmel))/2;
    for ( uint i = 0; i < npairs; ++i ){
      orbs.push_back(_orbs[2*i+1]);
      orbs.push_back(_orbs[2*i]);
    }
    for ( uint i = 2*npairs; i < _orbs.size(); ++i )
      orbs.push_back(_orbs[i]);

    return orbs;
  } else {
    return _orbs;
  }
}

LExcitationMap::iterator LExcitationMap::get_add(const std::string& name, int lmel )
{
  LExcitationMap::iterator itex = this->find(name);
  if ( itex != this->end() ) return itex;
  // not there yet. add it. 
#define _LPN LParsedName
  LParsedName exc(name,_LPN::Lmel|_LPN::Dg|_LPN::Excl|_LPN::Orbtypes);
#undef _LPN
  short excl = exc.excl;
  std::vector<OrbitalTypes> orbtypes = exc.orbtypes;
  
  if ( lmel != 0 && exc.lmel != lmel )
    error("Mismatch in non-conserving class in "+name,"LEquation::handle_excitation");
  if (lmel == 0 ) lmel = exc.lmel;
  // find excitation class
  if (!exc.foundsscipt && excl == 0 && lmel <= 0)
    error("No excitation class in "+name,"LEquation::handle_excitation");
  if (orbtypes.empty()) { // create default orbtypes (occ(excl), virt(excl+lmel))
    orbtypes.push_back(OrbitalTypes(Orbital::Occ,excl));
    orbtypes.push_back(OrbitalTypes(Orbital::Virt,excl+lmel));
  } else {
    assert(orbtypes.size() == 2);
    if ( int(orbtypes[0].size()) != excl || int(orbtypes[1].size()) != excl+lmel ){
      error("Inconsistency in orbital types and excitation class","LExcitationMap::get_add");
    }
  }
  TOrb4Type orb4t;
  for ( uint i = 0; i < orbtypes.size(); ++i ){
    _foreach_cauto(OrbitalTypes,iot,orbtypes[i]){
      if (orb4t.count(*iot) == 0)
        orb4t[*iot] = _globalterm.freeorbname(*iot);
    }
  }
  //TODO: implement Triplet!
  Matrices::Spinsym spinsym = Matrices::Singlet;
  Ops::Type opstype;
  if ( exc.dg )
    opstype = Ops::Deexc0;
  else
    opstype = Ops::Exc0;
  Oper op(opstype,excl,orb4t,orbtypes,"",lmel,&_globalterm);
  
  return (this->insert(std::make_pair(name,LExcitationInfo(op.mat().orbitals(),lmel,spinsym)))).first;
}
void LExcitationMap::set_lastorbs(const Product< Orbital >& orbs, Spin::Type spintype)
{
  // set lastorb (if smaller)
  for ( uint i = 0; i < orbs.size(); ++i )
    _globalterm.set_lastorb(Orbital(orbs[i].letname(),spintype),true);
}

void LExcitationMap::correct_orbs(const Product< Orbital >& orbs)
{
  if ( this->size() == 0 ) return;
  //make sure that we haven't used these orbital names already
  for ( uint i = 0; i < orbs.size(); ++i ){
    std::string newname;
    _foreach_auto(LExcitationMap,itex,*this){
      _foreach_auto(Product<Orbital>,itorb,itex->second._orbs){
        if ( itorb->letname() == orbs[i].letname() ){
          if (newname.empty()) newname = _globalterm.freeorbname(itorb->type()).letname();
          itorb->replace_letname(newname);
        }
      }
    }
  }
}

LParsedName::LParsedName(const std::string& namein, uint try2set)
              : lmel(0),dg(false),excl(0)
{
  std::string upname, downname;
  foundsscipt = IL::nameupdown(name,upname,downname,namein);
  if ( try2set == Name ) return;
  if (!upname.empty()) 
    this->parse_superscript(upname,try2set);
  
  if (!downname.empty()) 
    this->parse_subscript(downname,try2set);
}

void LParsedName::parse_superscript(const std::string& up, uint try2set)
{
  const TParArray& dgs = Input::aPars["syntax"]["dg"];
  const TParArray& lessmore = Input::aPars["syntax"]["lessmore"];
  lui ipos, ipos1;
  uint lmsize = 0;
  // less (negative number) or more (positive) electrons after the operator
  for ( TParArray::const_iterator itlm = lessmore.begin(); itlm != lessmore.end(); ++itlm ){
    lmsize = itlm->size();
    if ( itlm != lessmore.begin() && lmsize != itlm->size() ) {
      error("String length of less and more differ (syntax,lessmore)!","LParsedName::parse_superscript");
    }
  }
  lui last = std::string::npos;
  ipos=1;
  ipos=IL::skip(up,ipos,"{} ");
  while((ipos1=IL::nextwordpos(up,ipos,false))!=ipos && ipos < last ) {
    std::string nampart(up.substr(ipos,ipos1-ipos));
    if ( nameadd.size() > 0 && up[ipos]!='}' ) nameadd += " ";
    if(InSet(nampart, dgs)) {
      dg=true;
      nameadd += dgs.front();
    } else if (up[ipos]!='}') {
      nameadd += nampart;
      if (try2set&Lmel && nampart.size() > lmsize ){
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
void LParsedName::parse_subscript(const std::string& down, uint try2set)
{
  lui ipos, ipos1;
  ipos=1;
  ipos=IL::skip(down,ipos,"{} ");
  ipos1=IL::nextwordpos(down,ipos,false);
  excitation = down.substr(ipos,ipos1-ipos);
  if ( try2set&Excl && str2num<short>(excl,excitation,std::dec) ) {
    if ( gen_orbtypes(down.substr(ipos1))){
      assert(int(orbtypes[0].size()) == excl);
      assert(int(orbtypes[1].size()) == lmel+excl);
    }
  } else {
    // subscript is not an excitation class, probably rather something like \mu_2
    ipos1 = IL::closbrack(down,1);
    excitation = down.substr(ipos,ipos1-ipos);
  }
}
bool LParsedName::gen_orbtypes(const std::string& string)
{
  std::string name, up, down;
  IL::nameupdown(name,up,down,string);
  orbtypes.push_back(OrbitalTypes(up,true));
  orbtypes.push_back(OrbitalTypes(down,false));
  if ( orbtypes[0].empty() && orbtypes[1].empty()) orbtypes.clear();
  return !orbtypes.empty();
}

bool LEquation::extractit()
{
  // expand custom operators
  _eqn = _eqn.expandnewops(_newops);
  // expand parentheses
  _eqn.expand(_connections);
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


bool LEquation::do_sumterms(bool excopsonly )
{
  lui beg=0;
  bool plus=true, bra=false, ket=false;
  if (!_eqn.expanded())
    error("Expand the lexic equation first!","Lexic::do_sumterms");
  Term term;
  reset_term(term);
  Product<long int> indxoperterm;
  for (lui i=0; i<_eqn.size(); i++) {
    const Lelem & lel = _eqn[i];
    Lelem::Lex lex = lel.lex();
    if(InSet(lex, Lelem::Bra,Lelem::Ket)) { // handle bra/ket
      if ( lex == Lelem::Bra ) {
        if (bra) 
          error("Cannot handle two BRAs in one term yet...");
        else
          bra = true;
      } else if (ket) {
        error("Cannot handle two KETs in one term yet...");
      } else
        ket = true;
      term *= handle_braket(lel,term,excopsonly);
      indxoperterm.push_back(i+1);
    } else if (InSet(lex, Lelem::Minus,Lelem::Plus)) { // add current term and save the sign of the next term
      bra = ket = false; // reset bra and ket variables
      if (!excopsonly) {
        if ( i > 0 ) addterm(term,plus,beg,i-1,indxoperterm,excopsonly);
        plus = ( lex == Lelem::Plus );
        beg = i+1;
        reset_term(term);
        indxoperterm.clear();
      }
    } else if (InSet(lex, Lelem::Frac,Lelem::Num)) { // add prefactor
      term *= handle_factor(lel);
    } else if (lex == Lelem::Oper) { // handle Operator
      term *= handle_operator(lel,term,excopsonly);
      indxoperterm.push_back(i+1);
    } else if (lex == Lelem::Sum) { // handle \sum
      if (!excopsonly)
        _sumsterm.add(lel);
    } else if (lex == Lelem::Param) { // handle Parameter
      if (!excopsonly)
        term.addmatrix(handle_parameter(lel));
      indxoperterm.push_back(i+1);
    } else if (lex == Lelem::Perm) { // handle Permutation
      if (!excopsonly)
        term *= handle_permutation(lel);
    } else if (lex == Lelem::Times) { // handle Multiplication
      // don't do anything
    } else if (lex == Lelem::Div) { // handle Division
      error("Sorry, cannot handle Division!","Lexic::do_sumterms");
    } else {
      xout << lel << std::endl;
      error(lel.name()+" is not implemented yet...","Lexic::do_sumterms");
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
    term.copy_lastorbs(_excops.orbsterm());
  }
}

void LEquation::addterm(Term& term, bool plus, lui beg, lui end, 
                     Product<long int > const & indxoperterm, bool excopsonly)
{
  double minfac = Input::fPars["prog"]["minfac"];
  if( excopsonly || term.term_is_0(minfac)) return; // dont add zero term
  // handle sums in term
  for ( uint i = 0; i < _sumsterm.size(); ++i ) handle_sum(_sumsterm[i],term);
  // reset sums information
  _sumsterm=LelString();
  //add connections to term
  Product<long int> connect;
  long int ipos;
  for (unsigned long int i=0;i<_connections.size();i++) {
    if (abs(_connections[i].front())>(long int)beg && abs(_connections[i].back())-2<(long int)end) {
      for (unsigned long int j=0;j<_connections[i].size();j++) {
        ipos=indxoperterm.find(abs(_connections[i][j]));
        if (ipos<0)
          error("Connected operator is not in indxoperterm","LEquation::addterm");
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
    return handle_excitation(lelnam,(lel.lex()==Lelem::Bra),lm,excopsonly);
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
  if (excopsonly) {
    // update orbitals in excops
    _excops.set_lastorbs(occs,spintype);
    _excops.set_lastorbs(virts,spintype);
    //make sure that we haven't used these orbital names already
    _excops.correct_orbs(occs);
    _excops.correct_orbs(virts);
    return Oper();
  }
  excl = occs.size();
  // set lastorb (if smaller)
  for ( uint i = 0; i < occs.size(); ++i )
    term.set_lastorb(Orbital(occs[i].letname(),spintype),true);
  for ( uint i = 0; i < virts.size(); ++i )
    term.set_lastorb(Orbital(virts[i].letname(),spintype),true);
  int lmelec = virts.size()-occs.size();
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0,excl,occs,virts,"",lmelec,&term);  
  else
    return Oper(Ops::Exc0,excl,occs,virts,"",lmelec,&term);
}

Oper LEquation::handle_excitation(const std::string& name, bool dg, int lmel, bool excopsonly)
{
  LExcitationMap::iterator itex = _excops.get_add(name,lmel);
  if (excopsonly) return Oper();
  const LExcitationInfo & info = itex->second;
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0, info.orbitals(dg),"",info.lmel(dg));  
  else
    return Oper(Ops::Exc0,info.orbitals(dg),"",info.lmel(dg));
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
  LParsedName op(lel.name(),LParsedName::Name);
  // bare excitation operator 
  bool bare_excop = (InSet(op.name, bexcops));
#define _LPN LParsedName
  uint try2set = _LPN::Lmel|_LPN::Dg;
  if ( bare_excop ){
    try2set |= _LPN::Orbs|_LPN::Excitation;
  } else {
    try2set |= _LPN::Nameadd|_LPN::Excl|_LPN::Orbtypes;
  }
#undef _LPN
  op = LParsedName(lel.name(),try2set);
  std::string name = op.name;
  int lmelec = op.lmel;
  
  // parts of Hamilton operator
  if ( InSet(name, hms)) {
    if (excopsonly) return Oper();
    if (op.foundsscipt)
      say("Sub- and superscripts in Hamiltonian will be ignored: "+lel.name());
    if ( name==hms["fock"] ) return Oper(Ops::Fock,true,&term);
    if ( name==hms["oneelop"] ) return Oper(Ops::OneEl,true,&term);
    if ( name==hms["flucpot"] ) return Oper(Ops::FluctP,true,&term);
    if ( name==hms["dflucpot"] ) return Oper(Ops::FluctP,false,&term);
    if ( name==hms["perturbation"] ) return Oper(Ops::XPert,true,&term);
  }
  // excitation class
  if (op.excitation == "")
    error("No excitation class in operator "+lel.name(),"LEquation::handle_operator");
  IL::add2name(name,op.nameadd); // add nameadd to name (as superscript)
  if (bare_excop) { // bare excitation operator
    return handle_excitation(op.excitation,op.dg,lmelec,excopsonly);
  }
  if (excopsonly) return Oper();
  if (op.excl == 0 && lmelec <= 0)
    error("Excitation class in "+lel.name(),"LEquation::handle_operator");
  if (op.orbtypes.size() == 0){
    if(op.dg)
      return Oper(Ops::Deexc,op.excl,name,lmelec,&term);
    else
      return Oper(Ops::Exc,op.excl,name,lmelec,&term);
  } else {
    if(op.dg)
      return Oper(Ops::Deexc,op.excl,op.orbtypes,name,lmelec,&term);
    else
      return Oper(Ops::Exc,op.excl,op.orbtypes,name,lmelec,&term);
  }
}

void LEquation::handle_sum(const Lelem& lel, Term& term) const
{
  lui ipos, ipos1, up, down;
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
      term.addsummation(itex->second.orbitals());
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
Matrices LEquation::handle_parameter(const Lelem& lel)
{
#define _LPN LParsedName
  LParsedName op(lel.name(),_LPN::Lmel|_LPN::Dg|_LPN::Orbs|_LPN::Excitation|_LPN::Nameadd);
#undef _LPN
  std::string name = op.name;
  IL::add2name(name,op.nameadd); // add nameadd to name (as superscript)
  if ( op.excitation.empty() ) // no subscript, parameter is a "number"
    return Matrices(Ops::Number,Product<Orbital>(),0,name);
  // TODO check for orbitals here
  //
  LExcitationMap::const_iterator itex = _excops.get_add(op.excitation,op.lmel);
  
  return Matrices(Ops::Interm,itex->second.orbitals(op.dg),itex->second.exccls(op.dg), 
                  name,itex->second.spinsymexcs());
}


std::ostream& operator<<(std::ostream& o, const LEquation& inp)
{
  for (unsigned long int i=0; i<inp.eqn().size(); i++)
    o << inp.eqn().at(i);
  return o;
}