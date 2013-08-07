#include "finput.h"

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
    analyzenewops();
  } else if (InSet(linesp, eeqs)) {// end equation
    if ( iprint > 0 && _eq )
      _ineq.pop_back();
    _eq=false;
    analyzeq();
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
Product< Lelem > Finput::eqn() const
{ return _eqn.eqn(); }
Sum< Term, TFactor > Finput::sumterms() const
{ return _eqn.sumterms(); }

bool Finput::analyzeq()
{
  analyzeit();
  _eqn.extractit();
  _eqn.do_sumterms(true);
  _eqn.do_sumterms();
  return true;
}
bool Finput::analyzeit()
{
  char ch;
  lui i=0, ipos, ipos1;
  Lelem::Conn conn=Lelem::Normal;
  while (i<_input.size()) {
    ch=_input[i];
    if (ch=='<') { // bra
      ipos=_input.find('|',i);
      if (ipos==std::string::npos)
        error("Can not find bra","Finput::analyzeit");
      else {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        _eqn += Lelem(_input.substr(i,ipos1-i),Lelem::Bra);
      }
      i=ipos;
    } else if (ch=='|') { // ket
      ipos=_input.find('>',i);
      if (ipos==std::string::npos)
        error("Can not find ket","Finput::analyzeit");
      else {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        //connections
        if (_input.substr(ipos+1,2)=="_C") conn=Lelem::Connect;
        else if (_input.substr(ipos+1,2)=="_D") conn=Lelem::Disconnect;
        else conn=Lelem::Normal;
        _eqn += Lelem(_input.substr(i,ipos1-i),Lelem::Ket,conn);
        if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
          ipos+=2;
      }
      i=ipos;
    } else if (ch=='\\') { // command
      ++i;
      ipos=IL::nextwordpos(_input,i);
      ipos=analyzecommand(_input.substr(i,ipos-i),ipos);
      i=ipos-1;
    } else if (ch=='+') { // plus
      _eqn += Lelem("",Lelem::Plus);
    } else if (ch=='-') { // minus
      _eqn += Lelem("",Lelem::Minus);
    } else if (ch=='/') { // Division
      _eqn += Lelem("",Lelem::Div);
    } else if (ch=='*') { // times
      _eqn += Lelem("",Lelem::Times);
    } else if (InSet(ch, '(','[')) { // left parenthesis
      _eqn += Lelem("",Lelem::LPar);
    } else if (InSet(ch, ')',']')) { // right parenthesis
      //connections
      if (_input.substr(i+1,2)=="_C") conn=Lelem::Connect;
      else if (_input.substr(i+1,2)=="_D") conn=Lelem::Disconnect;
      else conn=Lelem::Normal;
      _eqn += Lelem("",Lelem::RPar,conn);
      if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
        i+=2;
    } else if (isdigit(ch)) { // number
      ipos=IL::nextwordpos(_input,i);
      _eqn += Lelem(_input.substr(i,ipos-i),Lelem::Num);
      i=ipos-1;
    } else if (InSet(ch, '&',' ')) { // do nothing
    } else
    //  std::cout << "Character " << ch << " is ignored" << std::endl;
      say("Character "+std::string(1,ch)+" is ignored");
    ++i;
  }
  return true;
}
lui Finput::analyzecommand(const std::string& str, lui ipos)
{
  const TParArray& skipops = Input::aPars["syntax"]["skipop"];
  TsPar& commands = Input::sPars["command"];
  // custom commands
  TsPar& newcom = Input::sPars["newcommand"];
  lui 
    ipos1=ipos, 
    ipos2, ipos3;
  if (str==commands["operator"]) { // operators
    ipos1=IL::nextwordpos(_input,ipos);
    _eqn += Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Oper);
  } else if (str==commands["parameter"]) { // parameters
    ipos1=IL::nextwordpos(_input,ipos);
    _eqn += Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Par);
  } else if (str==commands["fraction"]) { // fraction
    ipos1=IL::nextwordpos(_input,ipos);
    ipos2=ipos1;
    ipos3=IL::nextwordpos(_input,ipos2);
    _eqn += Lelem(_input.substr(ipos,ipos1-ipos)+"/"+_input.substr(ipos2,ipos3-ipos2),Lelem::Frac);
    ipos1=ipos3;
  } else if (str.substr(0,commands["sum"].size())==commands["sum"]) { // sum
    _eqn += Lelem(str.substr(commands["sum"].size()),Lelem::Sum);
  } else if (str==commands["permutation"]) { // permutation
    ipos1=IL::nextwordpos(_input,ipos);
    ipos2=ipos1;
    ipos3=IL::nextwordpos(_input,ipos2);
    _eqn += Lelem(_input.substr(ipos,ipos1-ipos)+"/"+_input.substr(ipos2,ipos3-ipos2),Lelem::Perm);
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
  TsPar& newops = Input::sPars["newoperator"];
  _eqn = Equation();
  for ( TsPar::iterator iop = newops.begin(); iop != newops.end(); ++iop ){
    _input = iop->second;
    _eqn +=  Lelem("",Lelem::LPar);
    analyzeit();
    _eqn +=  Lelem("",Lelem::RPar);
    _eqn.addnewop(iop->first,_eqn.eqn());
    _eqn.reseteq();
  }
  _input = "";
}

lui Equation::closbrack(const Product< Lelem >& eqn, lui ipos)
{
  lui i,ipos1=ipos;
  Lelem::Lex rk=Lelem::RPar,lk=eqn[ipos].lex();
  if (lk==Lelem::Bra)
    rk=Lelem::Ket;
  else if (lk==Lelem::LPar)
    rk=Lelem::RPar;
  else
    error("Not a bracket!","Finput::closbrack"); 
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
    error("Number of brackets is incosistent: "+nk,"Finput::closbrack"); 
  return ipos1;
}
lui Equation::openbrack(const Product< Lelem >& eqn, lui ipos)
{
  lui ipos1=ipos;
  Lelem::Lex lk=Lelem::LPar,rk=eqn[ipos].lex();
  if (rk==Lelem::Ket)
    lk=Lelem::Bra;
  else if (rk==Lelem::RPar)
    lk=Lelem::LPar;
  else
    error("Not a bracket!","Finput::openbrack"); 
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
    error("Number of brackets is incosistent: "+nk,"Finput::openbrack"); 
  return ipos1;
}
Product< long int > Equation::addconnections(const Product< Lelem >& aterm, lui beg, lui end)
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

bool Equation::extractit()
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
Product< Lelem > Equation::expandnewops(const Product< Lelem >& eqn)
{
  Product<Lelem> result;
  for (lui i = 0; i < eqn.size(); ++i ){
    result *= eqn[i];
    lui j = result.size()-1;
    do { 
      if ( result[j].lex() == Lelem::Oper && _newops.count(result[j].name()) ) {
        const Product<Lelem>& nop = _newops[result[j].name()];
        result.erase(result.begin()+j);
        result.insert(result.begin()+j,nop.begin(),nop.end());
      }
      ++j;
    } while( j < result.size() && result.size() < lui(Input::iPars["prog"]["maxlels"]) );
    if ( result.size() >= lui(Input::iPars["prog"]["maxlels"]) )
      error("Equation is too long. Cyclic dependencies in \\newoperator?");
  }
  return result;
}

lui Equation::elem(const Product< Lelem >& aterm, lui beg, bool bk)
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
  if (nk!=0 || braket) error("Check input, Mismatch in parentheses or bra/ket","Finput::term");
  if (i==beg || i==0 ) 
    end=0;
  else 
    end=i-1;
  return end;
}

lui Equation::term(const Product< Lelem >& eqn, lui beg)
{ return elem(eqn,beg,true); }

Product< Lelem > Equation::expandeqn(const Product< Lelem >& eqn, std::vector< Product<long int> > & connections)
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

Product< Lelem > Equation::expandterm(const Product< Lelem >& aterm, std::vector< Product<long int> > & connections)
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
Product< Lelem > Equation::expandpar(const Product< Lelem >& aterm, lui beg, 
                                   std::vector< Product< long int > >& connections)
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
bool Equation::expanded(const Product< Lelem >& eqn)
{ 
  for (lui i=0; i<eqn.size(); i++) {
    if (eqn[i].lex()==Lelem::Bra && !eqn[i].expandedbra())
      return false;
    if (eqn[i].lex()==Lelem::LPar)
      return false;
  }
  return true;
}
bool Equation::do_sumterms(bool excopsonly )
{
  lui i,beg=0,nterm=0;
  bool plus=true, ok=true, bra=false, ket=false;
  if (!expanded(_eqn))
    error("Expand finput first!","Finput::do_sumterms");
  Term term;
  if (_excops.size()>0) { // set last orbital of the term to the last orbital of pure excitation operators (has to be set before)
    for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
      Orbital::Type ot = static_cast<Orbital::Type>(it);
      for ( int i = _orbs4excops.size()-1; i >= 0; --i )
        if (_orbs4excops[i].count(ot))
          term.set_lastorb(_orbs4excops[i][ot]);
    }
  }
  term.addmatrix(Matrices());
  Product<long int> indxoperterm;
  for (i=0; i<_eqn.size(); i++) {
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
      term*=handle_braket(_eqn[i],term);
      indxoperterm *= i+1;
    } else if (InSet(_eqn[i].lex(), Lelem::Minus,Lelem::Plus)) { // add current term and save the sign of the next term
      bra=ket=false; // reset bra and ket variables
      if (!excopsonly) {
        if(i>0) addterm(term,plus,beg,i-1,indxoperterm,nterm,excopsonly);
        plus=(_eqn[i].lex()==Lelem::Plus);
        beg=i+1;
        term=Term();
        term.addmatrix(Matrices());
        if (_excops.size()>0) {
          for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
            Orbital::Type ot = static_cast<Orbital::Type>(it);
            for ( int i = _orbs4excops.size()-1; i >= 0; --i )
              if (_orbs4excops[i].count(ot))
                term.set_lastorb(_orbs4excops[i][ot]);
          }
        }
        indxoperterm=Product<long int>();
      }
    } else if (InSet(_eqn[i].lex(), Lelem::Frac,Lelem::Num)) { // add prefactor
      term*=handle_factor(_eqn[i]);
    } else if (_eqn[i].lex()==Lelem::Oper) { // handle Operator
      term*=handle_operator(_eqn[i],term,excopsonly);
      indxoperterm *= i+1;
    } else if (_eqn[i].lex()==Lelem::Sum) { // handle \sum
      if (!excopsonly)
        handle_sum(_eqn[i],term);
    } else if (_eqn[i].lex()==Lelem::Par) { // handle Parameter
      if (!excopsonly)
        _paramterm*=_eqn[i];
    } else if (_eqn[i].lex()==Lelem::Perm) { // handle Permutation
      if (!excopsonly)
        term *= handle_permutation(_eqn[i]);
    } else if (_eqn[i].lex()==Lelem::Times) { // handle Multiplication
      // don't do anything
    } else if (_eqn[i].lex()==Lelem::Div) { // handle Division
      error("Sorry, cannot handle Division!","Finput::do_sumterms");
    } else {
      xout << _eqn[i] << std::endl;
      error(_eqn[i].name()+" is not implemented yet...","Finput::do_sumterms");
    }
  }
  // add last term
  if(_eqn.size()>0) addterm(term,plus,beg,_eqn.size()-1,indxoperterm,nterm,excopsonly);
  return ok;
}
void Equation::addterm(Term& term, bool plus, lui beg, lui end, 
                     Product<long int > const & indxoperterm, lui & nterm, bool excopsonly)
{
  double minfac = Input::fPars["prog"]["minfac"];
  if( excopsonly || term.term_is_0(minfac)) {
    //reset parameter-info
    handle_parameters(term,true);
    return; // dont add zero term
  }
  ++nterm;
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
          error("Connected operator is not in indxoperterm","Finput::addterm");
        if (_connections[i][j]>0)
          connect*=ipos+2;
        else
          connect*=-ipos-2;
      }
      _xout2("Connections in Term #" << nterm << ": " <<connect<<std::endl);
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

Oper Equation::handle_braket(const Lelem& lel, Term& term)
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
    return handle_explexcitation(term,lelnam.substr(iend),(lel.lex()==Lelem::Bra));
  } else {
    return handle_excitation(term,lelnam,(lel.lex()==Lelem::Bra));
  }
}
Oper Equation::handle_explexcitation(Term& term, const std::string& name, bool dg)
{
  lui ipos, ipos1;
  short excl;
  Product<Orbital> occs, virts;
  if ( name[0] != '_' && name[0] != '^' )
    error("Doesn't start with _ or ^ :"+name,"Finput::handle_explexcitation");
  ipos = 0;
  ipos = IL::skip(name,ipos,"{}_^ ");
  while ( (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
    // TODO change it proper to singlet and triplet operators, handle sub and superscripts!!
    Orbital orb(IL::plainname(name.substr(ipos,ipos1-ipos)),Spin::GenS);
    if ( orb.type() == Orbital::Occ )
      occs *= orb;
    else if ( orb.type() == Orbital::Virt )
      virts *= orb;
    else
      error("general orbitals in excitation operators are not supported! "+name);
    ipos = IL::skip(name,ipos1,"{}_^ ");
  }
//   if ( occs.size() != virts.size() )
//     error("can handle particle-number conserving operators only! "+name);
  excl = occs.size();
  // set lastorb (if smaller)
  for ( uint i = 0; i < occs.size(); ++i )
    term.set_lastorb(Orbital(occs[i].letname(),Spin::GenS),true);
  for ( uint i = 0; i < virts.size(); ++i )
    term.set_lastorb(Orbital(virts[i].letname(),Spin::GenS),true);
  if ( _orbs4excops.size() > 0 ){
    //make sure that we haven't use these orbital names already
    for ( uint i = 0; i < occs.size(); ++i ){
      for ( uint iex = 0; iex < _orbs4excops.size(); ++iex ){
        for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
          Orbital::Type ot = static_cast<Orbital::Type>(it);
          if ( _orbs4excops[iex][ot] == Orbital(occs[i].letname(),Spin::GenS) )  // is there, change it
            _orbs4excops[iex][ot] = term.freeorbname(ot);
        }
      }
    }
    for ( uint i = 0; i < virts.size(); ++i ){
      for ( uint iex = 0; iex < _orbs4excops.size(); ++iex )
        for ( uint it = Orbital::Occ; it < Orbital::MaxType; ++it ){
          Orbital::Type ot = static_cast<Orbital::Type>(it);
          if ( _orbs4excops[iex][ot] == Orbital(virts[i].letname(),Spin::GenS) )  // is there, change it
            _orbs4excops[iex][ot] = term.freeorbname(ot);
        }
    }
  }
  int lmelec = virts.size()-occs.size();
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0,excl,occs,virts,"",lmelec,&term);  
  else
    return Oper(Ops::Exc0,excl,occs,virts,"",lmelec,&term);
}

Oper Equation::handle_excitation(Term& term, const std::string& name, bool dg, int lmel)
{
//   bool multiref = ( Input::iPars["prog"]["multiref"] > 0 );
  long int ipos2;
  std::string newname, nameadd,namedown;
  short excl;
  bool newdg;
  int lmelec;
  std::vector< Product<Orbital::Type> > orbtypes;
  bool updown=handle_namupdown(newname,excl,nameadd,namedown,newdg,lmelec,orbtypes,name);
  dg = (dg != newdg);
  if (lmelec != 0 && lmel != 0 && lmelec != lmel )
    error("Mismatch in non-conserving class in "+name,"Equation::handle_excitation");
  if (lmelec == 0 ) lmelec = lmel;
  // find excitation class
  if (!updown && excl == 0 && lmelec <= 0)
    error("No excitation class in "+name,"Equation::handle_excitation");
  ipos2=_excops.find(name);
  TOrb4Type orb4t;
  if (ipos2<0) {// first run (don't distinguish terms)
    _excops*=name;
    Product<Orbital::Type>::const_iterator iot;
    if (orbtypes.size() > 0){
      for ( uint i = 0; i < orbtypes.size(); ++i ){
        _foreach(iot,orbtypes[i]){
          if (orb4t.count(*iot) == 0)
            orb4t[*iot] = term.freeorbname(*iot);
        }
      }
    } else {
      orb4t[Orbital::Occ] = term.freeorbname(Orbital::Occ);
      orb4t[Orbital::Virt] = term.freeorbname(Orbital::Virt);
    }
    _orbs4excops.push_back(orb4t);
    _exccls*=excl;
    _spinsymexcs*=Matrices::Singlet; //TODO: implement Triplet!
    _posexcopsterm*=-1; // initialize
  } else {// may be second run...
    orb4t = _orbs4excops[ipos2];
    _posexcopsterm[ipos2]=term.mat().size(); //set position of this operator in the term
  }
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

TFactor Equation::handle_factor(const Lelem& lel)
{
  TFactor fac;
  lui ipos=0, ipos1;
  std::string lelnam=lel.name();
  if (lel.lex()==Lelem::Num) {
    double facd;
    if(!str2num<double>(facd,lelnam,std::dec))
        error("Factor is not a number "+lelnam,"Finput::handle_factor");
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
        error("Numerator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
      fac = facd;
    } else {
      // NOTE: won't work for non-integer nominators or denominators 
      long int nom;
      if(!str2num<long int>(nom,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Numerator is not an integer "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
      fac = nom;
    }
    ipos=lelnam.find("/");
    if(ipos==std::string::npos)
      error("Something wrong with frac "+lelnam,"Finput::handle_factor");
    ++ipos;
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos);
    if ( typeid(TFactor) == typeid(double) ){
      double fac1;
      if(!str2num<double>(fac1,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Denominator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
      fac /= fac1;
    } else {
      // NOTE: won't work for non-integer nominators or denominators 
      long int denom;
      if(!str2num<long int>(denom,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Denominator is not an integer "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
      fac /= denom;
    }
  }
  return fac;
}
Oper Equation::handle_operator(const Lelem& lel, Term& term, bool excopsonly)
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
    if ( name==hms["flucpot"] ) return Oper(Ops::FluctP,true,&term);
    if ( name==hms["dflucpot"] ) return Oper(Ops::FluctP,false,&term);
    if ( name==hms["perturbation"] ) return Oper(Ops::XPert,true,&term);
  }
  // excitation class
  if (namedown=="")
    error("No excitation class in operator "+lelnam,"Equation::handle_operator");
  IL::add2name(name,nameadd); // add nameadd to name (as superscript)
  if (InSet(name.substr(0,4), bexcops)) { // bare excitation operator
    assert(orbtypes.size() == 0);
    return handle_excitation(term,namedown,dg,lmelec);
  }
  if (excopsonly) return Oper();
  if (excl == 0 && lmelec <= 0)
    error("Excitation class in "+lelnam,"Equation::handle_operator");
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
bool Equation::handle_namupdown(std::string& name, short int& excl, std::string& nameup, std::string& namedown, bool& dg, 
                                int& lmel, std::vector< Product< Orbital::Type > >& orbtypes, const std::string& lelnam)
{ 
  const TParArray& dgs = Input::aPars["syntax"]["dg"];
  const TParArray& lessmore = Input::aPars["syntax"]["lessmore"];
  bool foundupdown = false;
  lui up, down, iposnam, ipos, ipos1;
  down=IL::lexfind(lelnam,"_");
  up=IL::lexfind(lelnam,"^");
  // last position of name of operator
  iposnam=std::min(up,down)-1;
  iposnam=std::min(iposnam,lelnam.size()-1);
  name=lelnam.substr(0,iposnam+1);
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
    ipos=up+1;
    ipos=IL::skip(lelnam,ipos,"{} ");
    while((ipos1=IL::nextwordpos(lelnam,ipos))!=ipos) {
      std::string nampart(lelnam.substr(ipos,ipos1-ipos));
      if(InSet(nampart, dgs)) {
        dg=true;
        nameup += dgs.front();
      } else if (lelnam[ipos]!='}') {
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
    ipos=down+1;
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos,false);
    namedown = lelnam.substr(ipos,ipos1-ipos);
    str2num<short>(excl,namedown,std::dec);
    if ( handle_orbtypes(orbtypes,lelnam.substr(ipos1))){
      assert(int(orbtypes[0].size()) == excl);
      assert(int(orbtypes[1].size()) == lmel+excl);
    }
  }
  return foundupdown;
}
bool Equation::handle_orbtypes(std::vector< Product< Orbital::Type > >& orbtypes, const std::string& string)
{
  bool foundorbtypes = false;
  lui up, down, ipos, ipos1;
  down=IL::lexfind(string,"_");
  up=IL::lexfind(string,"^");
  if (down!=std::string::npos && down!=string.size()-1){
    ipos = down;
    ipos = IL::skip(string,ipos,"{}_^ ");
    Product<Orbital::Type> occtypes;
    while ( (ipos < up) == (down < up) && (ipos1 = IL::nextwordpos(string,ipos,true,false)) != ipos ){//non greedy
      Orbital orb(IL::plainname(string.substr(ipos,ipos1-ipos)),Spin::GenS);
      occtypes *= orb.type();
      ipos = IL::skip(string,ipos1,"{}_^ ");
    }
    if ( occtypes.size() > 0 ) foundorbtypes = true;
    orbtypes.push_back(occtypes);
  }
  if (up!=std::string::npos && up!=string.size()-1){
    ipos = up;
    ipos = IL::skip(string,ipos,"{}_^ ");
    Product<Orbital::Type> virtypes;
    while ( (ipos < down) == (up < down) && (ipos1 = IL::nextwordpos(string,ipos,true,false)) != ipos ){//non greedy
      Orbital orb(IL::plainname(string.substr(ipos,ipos1-ipos)),Spin::GenS);
      virtypes *= orb.type();
      ipos = IL::skip(string,ipos1,"{}_^ ");
    }
    if ( virtypes.size() > 0 ) foundorbtypes = true;
    orbtypes.push_back(virtypes);
  }
  if (!foundorbtypes) orbtypes.clear();
  return foundorbtypes;
}

void Equation::handle_sum(const Lelem& lel, Term& term)
{
  lui ipos, ipos1, up, down;
  long int iposnam;
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
    iposnam=_excops.find(name);
    if (iposnam >= 0) {
      term.addsummation(_orbs4excops[iposnam][Orbital::Occ],_exccls[iposnam]);
      term.addsummation(_orbs4excops[iposnam][Orbital::Virt],_exccls[iposnam]);
    }
    else
      say("No excitation operator which would correspond to summation index "+name);
    ipos=ipos1;
  }
}
Permut Equation::handle_permutation(const Lelem& lel)
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
    error("Something wrong with frac "+lelnam,"Finput::handle_factor");
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

void Equation::handle_parameters(Term& term, bool excopsonly)
{
  if (!excopsonly) {// handle saved parameters
    lui ipos, ipos1, iposnam, up, down;
    int iposexcn,indxexcn;
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
        indxexcn=_excops.find(excn);
        if (indxexcn>=0) {
          iposexcn=_posexcopsterm[indxexcn];
          if (iposexcn>=0) {
            Matrices mat(Ops::Interm,
                         Ops::genprodorb(_exccls[indxexcn],_orbs4excops[indxexcn][Orbital::Occ],_orbs4excops[indxexcn][Orbital::Virt]),
                         _exccls[indxexcn], name,_spinsymexcs[indxexcn]);
            term.replacematrix(mat,iposexcn);
          } else
            say("Parameter is not present in this term: "+lelnam);
        } else
        // TODO : add parameters with explicit excitations
          error("Unknown excitation in parameter"+excn);
      }
    }
  }
  // reset all information
  _posexcopsterm.assign(_excops.size(),-1);
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
  else if (lel.lex()==Lelem::Par)
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

std::ostream& operator<<(std::ostream& o, const Finput& inp)
{
  for (unsigned long int i=0; i<inp.eqn().size(); i++)
    o << inp.eqn().at(i);
  return o;
}

