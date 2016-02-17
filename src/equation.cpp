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
    Error("Mismatch in non-conserving class in "+name);
  if (lmel == 0 ) lmel = exc.lmel;
  // find excitation class
  if (!exc.foundsscipt && excl == 0 && lmel <= 0)
    Error("No excitation class in "+name);
  if (orbtypes.empty()) { // create default orbtypes (occ(excl), virt(excl+lmel))
    orbtypes.push_back(OrbitalTypes(Orbital::Occ,excl));
    orbtypes.push_back(OrbitalTypes(Orbital::Virt,excl+lmel));
  } else {
    assert(orbtypes.size() == 2);
    if ( int(orbtypes[0].size()) != excl || int(orbtypes[1].size()) != excl+lmel ){
      Error("Inconsistency in orbital types and excitation class");
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

LParsedName::LParsedName(const std::string& namein, uint try2set, bool strict)
              : lmel(0),dg(false),excl(-1)
{
  std::string upname, downname;
  foundsscipt = IL::nameupdown(name,upname,downname,namein);
  if ( try2set == Name ) return;
 
  const TParArray& csfs = Input::aPars["syntax"]["csf"];
  // we swap sub- and superscript indices in Phi
  bool phi = InSet(name,csfs);
  
  if (!upname.empty()) 
    this->parse_superscript(upname,try2set);
  
  if (!downname.empty()) 
    this->parse_subscript(downname,try2set,strict);

  if ( found_orbs() && phi ){
    // swap for \phi 
    Product<Orbital> tmp(occ);
    occ = virt;
    virt = tmp;
  }
  if ( found_orbs() ) {
    assert( excl == -1 );
    excl = occ.size();
    lmel = virt.size() - occ.size();
  }
  // few checks
  if (!orbtypes.empty() &&
      ( int(orbtypes[0].size()) != excl ||
        int(orbtypes[1].size()) != lmel+excl )) {
    error(namein+": Inconsistency in the number of orbital types and the excitation class!", 
          "LParsedName");
  }
  _foreach_cauto(Product<Orbital>,itorb,occ){
    if (itorb->type() == Orbital::Virt && strict )
      warning("Do you really want to have orbital " << *itorb << " as occupied?");
  }
  _foreach_cauto(Product<Orbital>,itorb,virt){
    if (itorb->type() == Orbital::Occ && strict )
      warning("Do you really want to have orbital " << *itorb << " as virtual?");
  }
}

void LParsedName::parse_superscript(const std::string& up, uint try2set)
{
  const TParArray& dgs = Input::aPars["syntax"]["dg"];
  const TParArray& lessmore = Input::aPars["syntax"]["lessmore"];
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  const std::string& supername = Input::sPars["command"]["supername"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  lui ipos, ipos1;
  ipos=0;
  while((ipos1=IL::nextwordpos(up,ipos,true,false))!=ipos && ipos < up.size() ) {
    std::string word(up.substr(ipos,ipos1-ipos));
    IL::delbrack(word);
    lui 
      iposw = 0, 
      iposw1 = IL::nextwordpos(word,iposw,false,true);
    std::string mainpart(word.substr(iposw,iposw1-iposw));
    if ( nameadd.size() > 0 && up[ipos]!='}' ) nameadd += " ";
    if( try2set&Dg && InSet(word, dgs) ) {
      dg=true;
      nameadd += dgs.front();
    
    } else if ( try2set&Lmel && InSet(word,lessmore) ){
      // it is less/more
      ipos=ipos1;
      ipos1=IL::nextwordpos(up,ipos,false);
      std::string lmstr = up.substr(ipos,ipos1-ipos);
      nameadd += word+lmstr;
      if ( str2num<int>(lmel,lmstr,std::dec) ){
        // it is a non-conserving operator
        if ( word == lessmore.front() ) lmel = -lmel;
      } else {
        Error("Number of non-conserved electrons not recognized in "+up);
      }
    } else if ( try2set&Orbs && mainpart != "\\"+supername && !found_excitation() ){
      Orbital orb(IL::plainname(word),spintype);
      occ.push_back(orb);
      
    } else if ( try2set&Nameadd ){
      if ( mainpart == "\\"+supername && mainpart == word){
        // \snam{sfd}: add next word
        ipos1 = IL::nextwordpos(up,ipos1,true,false);
        word = up.substr(ipos,ipos1-ipos);
      }
      nameadd += word;
    } else {
      Error("Unknown part "+word+" in superscript "+up);
    }
    ipos=ipos1;
  }
}
void LParsedName::parse_subscript(const std::string& down, uint try2set, bool strict)
{
  const TParArray& excits = Input::aPars["syntax"]["excitation"];
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  lui ipos, ipos1;
  ipos = 0;
//   xout << "down: " << down << std::endl;
  while( (ipos1=IL::nextwordpos(down,ipos,true,false))!=ipos && ipos < down.size() ) {
    std::string word(down.substr(ipos,ipos1-ipos));
    IL::delbrack(word);
    lui 
      iposw = 0, 
      iposw1 = IL::nextwordpos(word,iposw,false,true);
    std::string mainpart(word.substr(iposw,iposw1-iposw));
    IL::delbrack(mainpart);
//     xout << "word: " << word << " mainpart " << mainpart << std::endl;
    short exclass;
    
    if ( try2set&Excl && str2num<short>(exclass,mainpart,std::dec) ) {
      // excitation class
      if (strict && (found_excitation() || found_orbs())) Error("Two excitations in "+down+
        "\n Use {} if there is more than one digit in the excitation class.");
      excl = exclass;
      if( gen_orbtypes(word.substr(iposw1)) && !(try2set&Orbtypes) )
        Error("Orbtypes present although not asked for in "+down);
      
    } else if ( try2set&Excitation && InSet(mainpart,excits) ){
      // something like \mu_2
      if (strict && (found_excitation() || found_orbs())) Error("Two excitations in "+down);
      excitation = word; 
      
    } else if ( try2set&Orbs ){
      if (strict && found_excitation()) Error("Excitations and orbitals at the same time in "+down);
      Orbital orb(IL::plainname(word),spintype);
      virt.push_back(orb);
      
    } else {
      Error("Unknown part "+word+" in subscript "+down);
    }
    ipos = ipos1;
    ipos = IL::skip(down,ipos,"} ");
  }
  
  // check
  if ( excl > 0 && !excitation.empty() ) 
    error("Found excitation class and explicit excitation in "+down,"LParsedName::parse_subscript"); 
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
        term.addsummation(handle_sum(lel));
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
void LEquation::correct_orbs(Term& term, const Product< Orbital >& occs, 
                             const Product< Orbital >& virts, Spin::Type spintype, bool excopsonly)
{
  if (excopsonly) {
    // update orbitals in excops
    _excops.set_lastorbs(occs,spintype);
    _excops.set_lastorbs(virts,spintype);
    //make sure that we haven't used these orbital names already
    _excops.correct_orbs(occs);
    _excops.correct_orbs(virts);
  } else {
    // set lastorb (if smaller)
    _foreach_cauto( Product<Orbital>,itorb,occs)
      term.set_lastorb(Orbital(itorb->letname(),spintype),true);
    _foreach_cauto( Product<Orbital>,itorb,occs)
      term.set_lastorb(Orbital(itorb->letname(),spintype),true);
  }
}

Oper LEquation::handle_braket(const Lelem& lel, Term& term, bool excopsonly)
{
  const TParArray& refs = Input::aPars["syntax"]["ref"];
  std::string lelnam=lel.name();
  if (InSet(lelnam, refs))
    return Oper(); // Reference, blank operator
  return handle_excitation(term,lelnam,(lel.lex()==Lelem::Bra),0,excopsonly);
}
Oper LEquation::handle_excitation(Term& term, const std::string& name, 
                                  bool dg, int lmel, bool excopsonly)
{
  const TParArray& excits = Input::aPars["syntax"]["excitation"];
  bool spinintegr = Input::iPars["prog"]["spinintegr"];
  Spin::Type spintype = Spin::Gen;
  if (spinintegr) spintype = Spin::GenS;
  LParsedName op(name,LParsedName::Name);
  std::string excit;
  if ( InSet(op.name,excits) ){
    // name is already the excitation-name
    excit = name;
  } else {
#define _LPN LParsedName
    op = LParsedName(name,_LPN::Orbs|_LPN::Nameadd|_LPN::Excitation|_LPN::Dg|_LPN::Lmel);
#undef _LPN
    if (op.found_orbs()){
      correct_orbs(term,op.occ,op.virt,spintype,excopsonly);
    } else if (op.found_excitation()){
      excit = op.excitation;
    } else {
      Error("Neither orbitals nor excitations in the excitation operator "+name);
    }
    if (op.lmel != lmel) {
      if (lmel == 0) {
        lmel = op.lmel;
      } else {
        Error("Inconsistency in number of non-conserved electrons in "+name);
      }
    }
    dg = ( dg != op.dg );
  }
  LExcitationMap::iterator itex;
  if ( !excit.empty() ){ // add this \mu_i
    itex = _excops.get_add(excit,lmel);
  }
  // create \tau_{excl}
  if (excopsonly) {
    return Oper();
  } else if (dg && op.found_orbs()) {
    return Oper(Ops::Deexc0,op.excl,op.occ,op.virt,"",lmel,&term);
  } else if (dg) {
    return Oper(Ops::Deexc0, itex->second.orbitals(dg),"",itex->second.lmel(dg));
  } else if (op.found_orbs()) {
    return Oper(Ops::Exc0,op.excl,op.occ,op.virt,"",lmel,&term);
  } else {
    return Oper(Ops::Exc0,itex->second.orbitals(dg),"",itex->second.lmel(dg));
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
  LParsedName op(lel.name(),LParsedName::Name);
  if (InSet(op.name, bexcops)) { // bare excitation operator
    return handle_excitation(term,lel.name(),false,0,excopsonly);
  }
#define _LPN LParsedName
  op = LParsedName(lel.name(),_LPN::Lmel|_LPN::Dg|_LPN::Nameadd|_LPN::Excl|_LPN::Orbtypes);
#undef _LPN
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
  if (op.excl < 0)
    Error("No excitation class in operator "+lel.name());
  IL::add2name(name,op.nameadd); // add nameadd to name (as superscript)
  if (excopsonly) return Oper();
  if (op.excl == 0 && lmelec <= 0)
    Error("Excitation class in "+lel.name());
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

Product<Orbital> LEquation::handle_sum(const Lelem& lel)
{
#define _LPN LParsedName
  LParsedName op(lel.name(),_LPN::Orbs|_LPN::Excitation|_LPN::Nameadd,false);
#undef _LPN
  if (!op.nameadd.empty())
    error("Sum from-to is not implemented yet: "+lel.name());
  if (op.excitation.empty() && !op.found_orbs())
    error("Sum without summation indices: "+lel.name());
  const std::string& excs = op.excitation;
  Product<Orbital> orbs(op.orbs());
  // iterate through excitations
  lui ipos = 0;
  while (ipos < excs.size()) {
    ipos = IL::skip(excs,ipos,"{}, ");
    if (ipos == excs.size()) break;
    lui ipos1 = IL::nextwordpos(excs,ipos);
    std::string name = excs.substr(ipos,ipos1-ipos);
    LExcitationMap::const_iterator itex = _excops.get_add(name);
    if (itex != _excops.end()) {
      orbs *= itex->second.orbitals();
//      term.addsummation(itex->second.orbitals());
    }
    ipos=ipos1;
  }
  return orbs;
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
  if ( op.excitation.empty() && !op.found_orbs() ) // no subscript, parameter is a "number"
    return Matrices(Ops::Number,Product<Orbital>(),0,name);
  
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