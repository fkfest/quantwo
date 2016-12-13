#include "finput.h"

Finput::Finput(bool eq) : 
_eq(eq){}

Finput::Finput(std::string paramspath) :
_eq(false)
{
  InitInpars(paramspath);
}

void Finput::InitInpars(std::string paramspath)
{
  std::string finp_file(paramspath+"params.reg");
  std::ifstream finp;
  finp.open(finp_file.c_str());
  if (finp.is_open()) {
    std::string 
      line, type, set, name;
    while (finp.good()) {
      std::getline (finp,line);
      if ( !line.empty()){
        _xout3(line << std::endl);
        set = IL::key(line,"set");
        if ( set.size() == 0 ) continue;
        type = IL::key(line,"type");
        name = IL::key(line,"name");
        if ( type.size() == 0 )
          error("no type is given");
        else if ( type[0] == 's' ){
          Input::sPars[set][name] = IL::key(line,"value");
        } else if ( type[0] == 'i' ){
          int x;
          if (!str2num<int>(x,IL::key(line,"value"),std::dec))
            error("integer parameter is not integer :"+line);
          Input::iPars[set][name] = x;
        } else if ( type[0] == 'f' ){
          double x; 
          if (!str2num<double>(x,IL::key(line,"value"),std::dec))
            error("float parameter is not float :"+line); 
          Input::fPars[set][name] = x;
        } else if ( type[0] == 'a' ){
          Input::aPars[set][name] = IL::parray(IL::key(line,"value"));
        } else
          error("unknown type in params.reg");
//       finput+=line;
      }
    }
  }
  else
    error(finp_file+": Bad input-parameters file!");
  finp.close();
}

bool Finput::addline(const std::string& line)
{
  const TParArray& beqs = Input::aPars["syntax"]["beq"];
  const TParArray& eeqs = Input::aPars["syntax"]["eeq"];
  const TParArray& newcs = Input::aPars["syntax"]["newcommand"];
  const TParArray& newops = Input::aPars["syntax"]["newoperator"];
  const TParArray& comments = Input::aPars["syntax"]["comment"];
  short iprint = Input::iPars["output"]["level"];
  if ( Input::iPars["prog"]["noorder"] > 0 ) {
    // default non-normal-ordered Hamiltonian: \op H = \op h + \op W
    TsPar& newops = Input::sPars["newoperator"];
    newops["H"] = "\\op h + \\op W";
  }
  if ( iprint > 0 && _eq )
    _ineq.push_back(line);
  if ( iprint > 1 && !_eq )
    _inlines.push_back(line);
  bool neweq = false;
  lui ipos=0, ipend;
  // and skip " " on begin
  ipos = IL::skip(line,ipos," ");
  // line without front-spaces and comments
  std::string linesp;
  // remove comments
  for (unsigned long int i=ipos; i<line.size(); i++) {
    if(InSet(line.substr(i,1),comments)) // comment
      break;
    linesp += line[i];
  }
  // remove spaces at the end
  ipend = IL::skipr(linesp,linesp.size()," ");
  linesp = linesp.substr(0,ipend);
  ipos = 0;
  ipend = IL::nextwordpos(linesp,ipos,false);
  if (ipos == ipend){ // empty line
  } else if (InSet(linesp, beqs)) {// begin equation
    _input="";
    if ( iprint > 1 && !_eq )
      _inlines.pop_back();
    _eq=true;
    _eqns.push_back(LEquation());
    analyzenewops();
  } else if (InSet(linesp, eeqs)) {// end equation
    if ( iprint > 0 && _eq )
      _ineq.pop_back();
    _eq=false;
    analyzeline();
    _input="";
    neweq = true;
  } else if (InSet(linesp.substr(ipos,ipend-ipos), newcs)) {// newcommand
    ipos = IL::addnewcom(linesp,ipend);
  } else if (InSet(linesp.substr(ipos,ipend-ipos), newops)) {// newoperator
    ipos = IL::addnewcom(linesp,ipend,"newoperator");
  } else if (!_eq){
    IL::changePars(linesp,ipos);
  } else {
    _input += linesp;
    _input += " "; // add space for separation
  }
  return neweq;
}
std::string Finput::input() const
{ return _input; }
LelString Finput::eqn() const
{ 
  if ( _eqns.size() > 0 )
    return _eqns.back().eqn();
  else
    return LelString();
}
TermSum Finput::sumterms() const
{ 
  if ( _eqns.size() > 0 )
    return _eqns.back().sumterms();
  else
    return TermSum();
}
void Finput::sumterms( const TermSum& ts )
{ 
  if ( _eqns.size() > 0 )
    return _eqns.back().sumterms(ts);
  else
    error("No equations yet!");
}

bool Finput::analyzeq()
{
  assert( _eqns.size() > 0 );
  insert_tensors();
  _eqns.back().extractit();
  _eqns.back().do_sumterms(true);
  _lhs.push_back( _eqns.back().do_sumterms() );
  return true;
}
bool Finput::analyzeline()
{
  assert( _eqns.size() > 0 );
  LEquation& eqn = _eqns.back();
  char ch;
  lui i=0, ipos, ipos1;
  Lelem::Conn conn=Lelem::Normal;
  while (i<_input.size()) {
    ch=_input[i];
    if (ch=='<') { // bra
      ipos=_input.find('|',i);
      if (ipos==std::string::npos)
        error("Can not find bra","Finput::analyzeline");
      else {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        eqn += Lelem(_input.substr(i,ipos1-i),Lelem::Bra);
      }
      i=ipos;
    } else if (ch=='|') { // ket
      ipos=_input.find('>',i);
      if (ipos==std::string::npos)
        error("Can not find ket","Finput::analyzeline");
      else {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        //connections
        if (_input.substr(ipos+1,2)=="_C") conn=Lelem::Connect;
        else if (_input.substr(ipos+1,2)=="_D") conn=Lelem::Disconnect;
        else conn=Lelem::Normal;
        eqn += Lelem(_input.substr(i,ipos1-i),Lelem::Ket,conn);
        if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
          ipos+=2;
      }
      i=ipos;
    } else if (ch=='\\') { // command
      ++i;
      ipos=analyzecommand(i);
      i=ipos-1;
    } else if (ch=='=') { // assignment
      eqn += Lelem("",Lelem::Equal);
    } else if (ch=='+') { // plus
      eqn += Lelem("",Lelem::Plus);
    } else if (ch=='-') { // minus
      eqn += Lelem("",Lelem::Minus);
    } else if (ch=='/') { // Division
      eqn += Lelem("",Lelem::Div);
    } else if (ch=='*') { // times
      eqn += Lelem("",Lelem::Times);
    } else if (ch=='(') { // left parenthesis
      eqn += Lelem("",Lelem::LPar);
    } else if (ch==')') { // right parenthesis
      //connections
      if (_input.substr(i+1,2)=="_C") conn=Lelem::Connect;
      else if (_input.substr(i+1,2)=="_D") conn=Lelem::Disconnect;
      else conn=Lelem::Normal;
      eqn += Lelem("",Lelem::RPar,conn);
      if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
        i+=2;
    } else if (ch=='[') { // left commutator bracket
      eqn += Lelem("",Lelem::LCom);
    } else if (ch==']') { // right commutator bracket
      eqn += Lelem("",Lelem::RCom);
    } else if (ch==',') { // comma from the commutator
      eqn += Lelem("",Lelem::Comma);
    } else if (isdigit(ch)) { // number
      ipos=IL::nextwordpos(_input,i);
      eqn += Lelem(_input.substr(i,ipos-i),Lelem::Num);
      i=ipos-1;
    } else if (InSet(ch, '&',' ')) { // do nothing
    } else
      say("Character "+std::string(1,ch)+" is ignored");
    ++i;
  }
  return true;
}
lui Finput::analyzecommand(lui ipos)
{
  assert( _eqns.size() > 0 );
  LEquation& eqn = _eqns.back();
  const TParArray& skipops = Input::aPars["syntax"]["skipop"];
  TsPar& commands = Input::sPars["command"];
  // custom commands
  TsPar& newcom = Input::sPars["newcommand"];
  lui ipos1 = IL::nextwordpos(_input,ipos);
  std::string str = _input.substr(ipos,ipos1-ipos);
  ipos = ipos1;
  lui 
    ipos2, ipos3;
  if (str==commands["operator"]) { // operators
    ipos1=IL::nextwordpos(_input,ipos);
    eqn += Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Oper);
  } else if (str==commands["tensor"]) { // tensor
    ipos1=IL::nextwordpos(_input,ipos);
    std::string name = _input.substr(ipos,ipos1-ipos);
    if ( name == "\\"+commands["integral"] ){ // two sets of indices {pq}{rs}
       ipos1=IL::nextwordpos(_input,ipos1);
       ipos1=IL::nextwordpos(_input,ipos1);
    } 
    eqn += Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Tensor);
  } else if (str==commands["fraction"]) { // fraction
    ipos1=IL::nextwordpos(_input,ipos);
    ipos2=ipos1;
    ipos3=IL::nextwordpos(_input,ipos2);
    eqn += Lelem(_input.substr(ipos,ipos1-ipos)+"/"+_input.substr(ipos2,ipos3-ipos2),Lelem::Frac);
    ipos1=ipos3;
  } else if (str.substr(0,commands["sum"].size())==commands["sum"]) { // sum
    eqn += Lelem(str.substr(commands["sum"].size()),Lelem::Sum);
  } else if (str==commands["permutation"]) { // permutation
    ipos1=IL::nextwordpos(_input,ipos);
    ipos2=ipos1;
    ipos3=IL::nextwordpos(_input,ipos2);
    eqn += Lelem(_input.substr(ipos,ipos1-ipos)+"/"+_input.substr(ipos2,ipos3-ipos2),Lelem::Perm);
    ipos1=ipos3;
  } else if (InSet(str, skipops)){//,"left","right","lk","rk","\\"))
  } else if (newcom.count(str)){// custom command
    // replace and go back
    assert( ipos > str.size() );
    ipos1 = _input.rfind('\\',ipos-str.size());
    assert( ipos1 != std::string::npos );
    _input.replace(ipos1,ipos-ipos1,newcom[str]);
  } else 
    error("Unknown command in equation! "+str,"Finput::analyzecommand");
  return ipos1;
}
void Finput::analyzenewops()
{
  assert( _input.size() == 0 );
  assert( _eqns.size() > 0 );
  LEquation& eqn = _eqns.back();
  TsPar& newops = Input::sPars["newoperator"];
  eqn = LEquation();
  for ( TsPar::iterator iop = newops.begin(); iop != newops.end(); ++iop ){
    _input = iop->second;
    eqn +=  Lelem("",Lelem::LPar);
    analyzeline();
    eqn +=  Lelem("",Lelem::RPar);
    eqn.addnewop(iop->first,eqn.eqn());
    eqn.reseteq();
  }
  _input = "";
}

void Finput::insert_tensors()
{
  // at least the current eqn should be there
  assert( _eqns.size() > 0 );
  // lhs from the current eqn is not there yet
  assert( _eqns.size() == _lhs.size()+1 );
  for ( uint i = 0; i < _lhs.size(); ++i ){
    
  }
}

std::ostream& operator<<(std::ostream& o, const Finput& inp)
{
  for (unsigned long int i=0; i<inp.eqn().size(); i++)
    o << inp.eqn().at(i);
  return o;
}

