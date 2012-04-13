#include "finput.h"

const std::string IL::separator=" \\{}()[]+-*/&<>|";
const std::string IL::gluer="^_";
const std::string IL::brackets="(){}[]<>";
// define possible operator names (the order is very important sometimes!)
const std::string Finput::commands[5]={ "op", "prm", "sum", "frac", "half"};
const std::string Finput::beqs[2]={"\\beq", "\\begin(equation)"};
const std::string Finput::eeqs[2]={"\\eeq", "\\end(equation)"};
const std::string Finput::comments[2]={"%", "\\comment"};
const std::string Finput::skipops[5]={"left","right","lk","rk","\\"};
const std::string Finput::refs[3]={"0","\\Phi_0", "HF"};
const std::string Finput::csfs[2]={"\\Phi", "\\phi"};
const std::string Finput::dgs[2]={"\\dg","\\dagger"};
const std::string Finput::bexcops[2]={"\\tau","\\Tau"};
const char Finput::hms[4]={'H', 'F', 'W', 'X'}; 

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

std::string IL::key(const std::string& line, const std::string& keyword)
{
  lui ipos; 
  while (true) {
    ipos = line.find(keyword);
    if ( ipos != std::string::npos ){
      ipos += keyword.size();
      ipos = skip(line,ipos," ");
      if ( line[ipos] == '=' )
	break;
    } else
      error(keyword+" not found");
  }
  ++ipos;
  lui ipend = endword(line,ipos);
  return line.substr(ipos,ipend-ipos);
}
lui IL::skip(const std::string& str, const lui& ipos, const std::string& what)
{
  lui ires=ipos;
  while (ires<str.size()&& what.find(str[ires])!=std::string::npos )
    ++ires;
  return ires;
}
lui IL::endword(const std::string& line, lui& ipos)
{
  assert(ipos < line.size());
  ipos = skip(line,ipos," ");
  char end = ',';
  if ( line[ipos] == '"' ){
    end = '"'; ++ipos;
  }else if ( line[ipos] == '{' ){
    end = '}'; ++ipos;
  }
  lui ipend;
  for ( ipend = ipos+1; ipend < line.size() && line[ipend] != end; ++ipend ){}
  return ipend;
}
lui IL::closbrack(const std::string& str, lui ipos)
{
  lui i=brackets.find(str[ipos]),ipos1=ipos;
  if (i==std::string::npos) 
    error(str[ipos]+"is not a bracket!","IL::closbrack"); 
  char lk(brackets[i]), rk(brackets[i+1]); // left and right brackets
  int nk=1;
  for ( i=ipos+1;i<str.size();++i)
  {
    if (str[i]==lk) 
      ++nk; // count number of "("
    else if (str[i]==rk) 
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
    error("Number of brackets is incosistent: "+nk,"IL::closbrack"); 
  return ipos1;
}

lui IL::nextwordpos(const std::string& str, lui& ipos, bool glue, bool greedy)
{
  lui nwpos;
  ipos=IL::skip(str,ipos," "); // remove spaces
  // e.g. from ab_{cd}ef
  if ( greedy && glue ){
    // ab_{cd}
    for (nwpos=ipos ;nwpos < str.size() && (separator.find(str[nwpos])==std::string::npos || 
                    nwpos==ipos || gluer.find(str[nwpos-1])!=std::string::npos); ++nwpos) 
      if (char(str[nwpos])=='{') nwpos=closbrack(str,nwpos);
  } else if ( greedy && !glue ) {
    // ab
    for (nwpos=ipos ;nwpos < str.size() && ((separator.find(str[nwpos])==std::string::npos && 
                    gluer.find(str[nwpos])==std::string::npos) || nwpos==ipos); ++nwpos) 
      if (char(str[nwpos])=='{') nwpos=closbrack(str,nwpos);
  } else if ( !greedy && glue ){
    // a (and b_{cd})
    for (nwpos=ipos ;nwpos < str.size() && (gluer.find(str[nwpos])!=std::string::npos ||
                   nwpos==ipos || gluer.find(str[nwpos-1])!=std::string::npos); ++nwpos) 
      if (char(str[nwpos])=='{') nwpos=closbrack(str,nwpos);
  } else {
    // a (and b and {cd} )
    for (nwpos=ipos ;nwpos < str.size() && nwpos==ipos; ++nwpos) //stupid construction... :)
      if (char(str[nwpos])=='{') nwpos=closbrack(str,nwpos);
  }
  return nwpos;
}
inline std::string IL::plainname(std::string name)
{
  std::string pname;
  for ( uint i = 0; i < name.size(); ++i )
    if ( !InSet(name[i],"_^{}") )
      pname += name[i];
  return pname;
}


Finput::Finput() : 
_eq(false){}

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
	std::cout << line << std::endl;
	type = IL::key(line,"type");
	set = IL::key(line,"set");
	name = IL::key(line,"name");
	
	if ( type.size() == 0 )
	  error("no type is given");
	else if ( type[0] == 's' ){
	  Input::sInpars[name] = IL::key(line,"value");
	} else if ( type[0] == 'i' ){
	  int x;
	  str2num<int>(x,IL::key(line,"value"),std::dec);
	  Input::iInpars[name] = x;
	} else if ( type[0] == 'f' ){
	  double x; 
	  str2num<double>(x,IL::key(line,"value"),std::dec); 
	  Input::fInpars[name] = x;
	} else if ( type[0] == 'a' ){
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

Finput& Finput::operator+=(const std::string& line)
{
  unsigned long int ipos=0;
  ipos=IL::skip(line,ipos," "); // skip space characters
  if (InSet(line.substr(ipos), beqs))
  {// begin equation
    _input="";
    _eq=true;
  }
  else if (InSet(line.substr(ipos), eeqs))
  {// end equation
    _eq=false;
    analyzeit();
    extractit();
    do_sumterms(true);
    do_sumterms();
    _input="";
  }
  else
    for (unsigned long int i=ipos; i<line.size(); i++)
    {
      if(InSet(line.substr(i,1),comments)) // comment
        break;
      _input+=line[i];
    }
  _input+=" "; // add space for separation
  return *this;
}
std::string Finput::input() const
{ return _input; }
Product< Lelem > Finput::eqn() const
{ return _eqn; }
Sum< Term, double > Finput::sumterms() const
{ return _sumterms; }
bool Finput::analyzeit()
{
  char ch;
  unsigned long int i=0, ipos, ipos1;
  Lelem::Conn conn=Lelem::Normal;
  while (i<_input.size())
  {
    ch=_input[i];
    if (ch=='<')
    { // bra
      ipos=_input.find('|',i);
      if (ipos==std::string::npos)
        error("Can not find bra","Finput::analyzeit");
      else
      {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        _eqn*=Lelem(_input.substr(i,ipos1-i),Lelem::Bra);
      }
      i=ipos;
    }
    else if (ch=='|')
    { // ket
      ipos=_input.find('>',i);
      if (ipos==std::string::npos)
        error("Can not find ket","Finput::analyzeit");
      else
      {
        ++i;
        ipos1=IL::nextwordpos(_input,i);
        //connections
        if (_input.substr(ipos+1,2)=="_C") conn=Lelem::Connect;
        else if (_input.substr(ipos+1,2)=="_D") conn=Lelem::Disconnect;
        else conn=Lelem::Normal;
        _eqn*=Lelem(_input.substr(i,ipos1-i),Lelem::Ket,conn);
        if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
          ipos+=2;
      }
      i=ipos;
    }
    else if (ch=='\\')
    { // command
      ++i;
      ipos=IL::nextwordpos(_input,i);
      ipos=analyzecommand(_input.substr(i,ipos-i),ipos);
      i=ipos-1;
    }
    else if (ch=='+')
    { // plus
      _eqn*=Lelem("",Lelem::Plus);
    }
    else if (ch=='-')
    { // minus
      _eqn*=Lelem("",Lelem::Minus);
    }
    else if (ch=='/')
    { // Division
      _eqn*=Lelem("",Lelem::Div);
    }
    else if (ch=='*')
    { // times
      _eqn*=Lelem("",Lelem::Times);
    }
    else if (InSet(ch, '(','['))
    { // left parenthesis
      _eqn*=Lelem("",Lelem::LPar);
    }
    else if (InSet(ch, ')',']'))
    { // right parenthesis
      //connections
      if (_input.substr(i+1,2)=="_C") conn=Lelem::Connect;
      else if (_input.substr(i+1,2)=="_D") conn=Lelem::Disconnect;
      else conn=Lelem::Normal;
      _eqn*=Lelem("",Lelem::RPar,conn);
      if (InSet(conn, Lelem::Connect,Lelem::Disconnect))
        i+=2;
    }    
    else if (isdigit(ch))
    { // number
      ipos=IL::nextwordpos(_input,i);
      _eqn*=Lelem(_input.substr(i,ipos-i),Lelem::Num);
      i=ipos-1;
    }  
    else if (InSet(ch, '&',' '))
      ; // do nothing
    else
    //  std::cout << "Character " << ch << " is ignored" << std::endl;
      say("Character "+std::string(1,ch)+" is ignored");
    ++i;
  }
  return true;
}
unsigned long int Finput::analyzecommand(const std::string& str,unsigned long int ipos)
{
  unsigned long int ipos1=ipos, ipos2, ipos3;
  if (str==commands[0])//"op")
  { // operators
    ipos1=IL::nextwordpos(_input,ipos);
    _eqn*=Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Oper);
  }
  else if (str==commands[1])//"prm")
  { // parameters
    ipos1=IL::nextwordpos(_input,ipos);
    _eqn*=Lelem(_input.substr(ipos,ipos1-ipos),Lelem::Par);
  }
  else if (str.substr(0,3)==commands[2])//"sum")
  { // sum
    _eqn*=Lelem(str.substr(3),Lelem::Sum);
  }
  else if (str.substr(0,4)==commands[3])//"frac")
  { // fraction
    ipos1=IL::nextwordpos(_input,ipos);
    ipos2=ipos1;
    ipos3=IL::nextwordpos(_input,ipos2);
    _eqn*=Lelem(_input.substr(ipos,ipos1-ipos)+"/"+_input.substr(ipos2,ipos3-ipos2),Lelem::Frac);
    ipos1=ipos3;
  }
  else if (str==commands[4])//"half")
  { // a half
    _eqn*=Lelem("0.5",Lelem::Num);
  }
  else if (!InSet(str, skipops))//,"left","right","lk","rk","\\"))
    error("Unknown command in equation! "+str,"Finput::analyzecommand");
  return ipos1;
}
long unsigned int Finput::closbrack(const Product< Lelem >& eqn, long unsigned int ipos)
{
  unsigned long int i,ipos1=ipos;
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
long unsigned int Finput::openbrack(const Product< Lelem >& eqn, long unsigned int ipos)
{
  unsigned long int ipos1=ipos;
  Lelem::Lex lk=Lelem::LPar,rk=eqn[ipos].lex();
  if (rk==Lelem::Ket)
    lk=Lelem::Bra;
  else if (rk==Lelem::RPar)
    lk=Lelem::LPar;
  else
    error("Not a bracket!","Finput::openbrack"); 
  int nk=-1;
  for ( long int i=ipos-1;i>=0;--i)
  {
    if (eqn[i].lex()==lk) 
    {
      ++nk; // count number of "("
      if (nk==0) 
      {
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
Product< long int > Finput::addconnections(const Product< Lelem >& aterm,long unsigned int beg, unsigned long int end)
{
  Product<long int> connection;
  if (aterm[end].conn()==Lelem::Normal) return connection;
  for (unsigned long int i=beg+1; i<end; i++)
    if (aterm[i].lex()==Lelem::Oper)
    {
      if(aterm[end].conn()==Lelem::Connect)
        connection *= i+1;
      else
        connection *= -(i+1);
    }
  return connection;
}

bool Finput::extractit()
{
  // expand Hamiltonians
  _eqn=expandH(_eqn);
  // expand parentheses
  _eqn=expandeqn(_eqn,_connections);
  // remove redundant connections
  for (unsigned long int i=0; i<_connections.size();i++)
    for (unsigned long int j=0; j<_connections[i].size();j++)
      if (InSet(_eqn[abs(_connections[i][j])-1].lex(), Lelem::Num,Lelem::Frac)) //connection to a number
        _connections[i].erase(_connections[i].begin()+j);
  for (unsigned long int i=0; i<_connections.size();i++)
    if (_connections[i].size()<2) //smaller than two elements "connected"
      _connections.erase(_connections.begin()+i);
  if (_connections.size()==0) return true;
  for (unsigned long int i=0; i<_connections.size()-1;i++)
    for (unsigned long int j=i+1; j<_connections.size();j++)
      if (_connections[i]==_connections[j]) //same connection
      {
        _connections.erase(_connections.begin()+j);
        --j;
      }
  for (unsigned int k = 0; k<_connections.size();k++)
    std::cout << "final Connection " << k << ": " << _connections[k] << std::endl;
  return true;
}
Product<Lelem> Finput::expandH(Product<Lelem> const & eqn)
{
  Product<Lelem> result;
  for (unsigned long int i=0; i<eqn.size(); i++ )
  {
    if (eqn[i].lex()==Lelem::Oper && eqn[i].name().at(0)==hms[0])//H
    {// (\op F + \op W)
      result *= Lelem("",Lelem::LPar);
      result *= Lelem(hms[1]+eqn[i].name().substr(1),Lelem::Oper);//F
      result *= Lelem("",Lelem::Plus);
      result *= Lelem(hms[2]+eqn[i].name().substr(1),Lelem::Oper);//W
      result *= Lelem("",Lelem::RPar);
    }
    else  
      result *= eqn[i];
  }
  return result;
}

long unsigned int Finput::elem(const Product< Lelem >& aterm, long unsigned int beg, bool bk)
{
  unsigned long int i, end, nk=0;
  bool braket=false;
  
  for ( i=beg; i<aterm.size() ; i++ )
  {
    if (aterm[i].lex()==Lelem::LPar) ++nk;
    if (aterm[i].lex()==Lelem::RPar) --nk;
    if (bk)
    {
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

long unsigned int Finput::term(const Product< Lelem >& eqn, long unsigned int beg)
{ return elem(eqn,beg,true); }

Product< Lelem > Finput::expandeqn(const Product< Lelem >& eqn, std::vector< Product<long int> > & connections)
{
  Product<Lelem> result=eqn, res;
  Product<long int>  connect;
  std::vector< Product<long int> > con,con1;
  unsigned long int beg=0, end, i,j,conbeg,lastpos;
  
  while (!expanded(result))
  {
    res=result;
    result=Product<Lelem>();
    beg=0;
    conbeg=0;
    con=connections;
    connections=std::vector< Product<long int> >();
    while ((end=term(res,beg))!=0)
    {
      std::vector< Product<long int> > con1;
      // shift connections
      for (i=conbeg;i<con.size();i++)
      {
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
      for (i=0;i<con1.size();i++)
      {
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

Product< Lelem > Finput::expandterm(const Product< Lelem >& aterm, std::vector< Product<long int> > & connections)
{
  unsigned long int i;
//  std::cout << "Term: " << aterm << std::endl;
//  for (unsigned int k = 0; k<connections.size();k++)
//    std::cout << "in connection " << k << ": " << connections[k] << std::endl;
  for (i=0; i<aterm.size(); i++)
    if (aterm[i].lex()==Lelem::LPar || (aterm[i].lex()==Lelem::Bra && !aterm[i].expandedbra())) 
      return expandpar(aterm,i,connections);
  return aterm;
}
Product< Lelem > Finput::expandpar(const Product< Lelem >& aterm, long unsigned int beg, 
                                   std::vector< Product< long int > >& connections)
{ // e.g., aterm=-a(b+c)d
  unsigned long int i=0, end,ipos,start=0,ijcon,iposres,lena,lenb,ipar;
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
  while (i<inpar.size())
  {
    ipos=elem(inpar,i);
    //add sign
    if (inpar[i].lex()==Lelem::Minus) // minus in parenthesis
    {
      if (aterm[0].lex()==Lelem::Minus) 
        result *= Lelem("",Lelem::Plus); // -1*-1 == +1
      else
        result *= inpar[i]; // +1*-1 == -1
      ++i;
    }
    else if (inpar[i].lex()==Lelem::Plus) // plus in parenthesis
    {
      if (aterm[0].lex()==Lelem::Minus) 
        result *= aterm[0]; // -1*+1 == -1
      else
        result *= inpar[i]; // +1*+1 == +1
      ++i;
    }
    else if (start>0) // term has a sign, but no sign by expression in parenthesis
    {
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
    lena=beg-start;
    lenb=ipos-i+1;
    for (unsigned long int k=0;k<con.size();k++)
    {
      for (unsigned long int j=0;j<con[k].size();j++)
      {
        ijcon=abs(con[k][j]);
        coninterm=(ijcon<=beg)||(ijcon>end+1);// connection in a or d
        coninpar=((ijcon>i+beg+1)&&(ijcon<ipos+beg+3));// connection in b (or c)
        coninterm=coninterm||coninpar;
        if (coninterm)
        { // change ijcon
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
      if (connect.size()>0)
      {
        connections.push_back(connect);
        connect=Product<long int>();
      }
    }
    i=ipos+1;
  }
  return result;
}
bool Finput::expanded(const Product< Lelem >& eqn)
{ 
  for (unsigned long int i=0; i<eqn.size(); i++)
  {
    if (eqn[i].lex()==Lelem::Bra && !eqn[i].expandedbra())
      return false;
    if (eqn[i].lex()==Lelem::LPar)
      return false;
  }
  return true;
}
bool Finput::do_sumterms(bool excopsonly )
{
  unsigned long int i,beg=0,nterm=0;
  bool plus=true, ok=true, bra=false, ket=false;
  if (!expanded(_eqn))
    error("Expand finput first!","Finput::do_sumterms");
  Term term;
  if (_excops.size()>0)
  { // set last orbital of the term to the last orbital of pure excitation operators (has to be set before)
    term.set_lastorb(_occexcops.back());
    term.set_lastorb(_virexcops.back());
  }
  Product<long int> indxoperterm;
  for (i=0; i<_eqn.size(); i++)
  {
    if(InSet(_eqn[i].lex(), Lelem::Bra,Lelem::Ket))
    { // handle bra/ket
      if (_eqn[i].lex()==Lelem::Bra)
      {
        if (bra) 
          error("Can not handle two BRAs in one term yet...");
        else
          bra=true;
      }
      else if (ket)
        error("Can not handle two KETs in one term yet...");
      else
        ket=true;
      term*=handle_braket(_eqn[i],term);
      indxoperterm *= i+1;
    }
    else if (InSet(_eqn[i].lex(), Lelem::Minus,Lelem::Plus))
    { // add current term and save the sign of the next term
      bra=ket=false; // reset bra and ket variables
      if (!excopsonly)
      {
        if(i>0) addterm(term,plus,beg,i-1,indxoperterm,nterm,excopsonly);
        plus=(_eqn[i].lex()==Lelem::Plus);
        beg=i+1;
        term=Term();
        if (_excops.size()>0)
        {
          term.set_lastorb(_occexcops.back());
          term.set_lastorb(_virexcops.back());
        }
        indxoperterm=Product<long int>();
      }
    }
    else if (InSet(_eqn[i].lex(), Lelem::Frac,Lelem::Num))
    { // add prefactor
      term*=handle_factor(_eqn[i]);
    }
    else if (_eqn[i].lex()==Lelem::Oper)
    { // handle Operator
      term*=handle_operator(_eqn[i],term,excopsonly);
      indxoperterm *= i+1;
    }
    else if (_eqn[i].lex()==Lelem::Sum)
    { // handle \sum
      if (!excopsonly)
        handle_sum(_eqn[i],term);
    }
    else if (_eqn[i].lex()==Lelem::Par)
    { // handle Parameter
      if (!excopsonly)
        _paramterm*=_eqn[i];
    }
    else if (_eqn[i].lex()==Lelem::Div)
    { // handle Division
      error("Sorry, cannot handle Division!","Finput::do_sumterms");
    }
    else
    {
      error(_eqn[i].name()+" is not implemented yet...","Finput::do_sumterms");
    }
  }
  // add last term
  if(_eqn.size()>0) addterm(term,plus,beg,_eqn.size()-1,indxoperterm,nterm,excopsonly);
  return ok;
}
void Finput::addterm(Term& term, bool plus, long unsigned int beg, long unsigned int end, 
                     Product<long int > const & indxoperterm, unsigned long int & nterm, bool excopsonly)
{
  if(excopsonly || term.term_is_0())
  {
    //reset parameter-info
    handle_parameters(term,true);
    return; // dont add zero term
  }
  ++nterm;
  // handle parameters
  handle_parameters(term);
  //add connections to term
  Product<long int> connect;
  long int ipos;
  for (unsigned long int i=0;i<_connections.size();i++)
  {
    if (abs(_connections[i].front())>(long int)beg && abs(_connections[i].back())-2<(long int)end)
    {
      for (unsigned long int j=0;j<_connections[i].size();j++)
      {
        ipos=indxoperterm.find(abs(_connections[i][j]));
        if (ipos<0)
          error("Connected operator is not in indxoperterm","Finput::addterm");
        if (_connections[i][j]>0)
          connect*=ipos+1;
        else
          connect*=-ipos-1;
      }
      std::cout << "Connections in Term #" << nterm << ": " <<connect<<std::endl;
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

Oper Finput::handle_braket(const Lelem& lel, Term& term)
{
  std::string lelnam=lel.name();
  if (InSet(lelnam, refs))
    return Oper(); // Reference, blank operator
    
  lui 
    ibeg = 0,
    iend = IL::nextwordpos(lelnam,ibeg,false);
  if (InSet(lelnam.substr(ibeg,iend-ibeg), csfs)){
    return handle_explexcitation(lelnam.substr(iend),(lel.lex()==Lelem::Bra),term);
  } else
    return handle_excitation(lelnam,(lel.lex()==Lelem::Bra),term);
}
Oper Finput::handle_explexcitation(const std::string& name, bool dg, Term& term)
{
  lui ipos, ipos1;
  long int ipos2;
  short excl;
  Product<Orbital> occs, virts;
  if ( name[0] != '_' && name[0] != '^' )
    error("Doesn't start with _ or ^ :"+name,"Finput::handle_explexcitation");
  ipos = 0;
  ipos = IL::skip(name,ipos,"{}_^ ");
  while ( (ipos1 = IL::nextwordpos(name,ipos,true,false)) != ipos ){//non greedy
    Orbital orb(IL::plainname(name.substr(ipos,ipos1-ipos)),Orbital::GenS);
    if ( orb.type() == Orbital::Occ )
      occs *= orb;
    else if ( orb.type() == Orbital::Virt )
      virts *= orb;
    else
      error("general orbitals in excitation operators are not supported! "+name);
    ipos = IL::skip(name,ipos1,"{}_^ ");
  }
  if ( occs.size() != virts.size() )
    error("can handle particle-number conserving operators only! "+name);
  excl = occs.size();
  // set lastorb (if smaller)
  for ( uint i = 0; i < occs.size(); ++i )
    term.set_lastorb(Orbital(occs[i].letname(),Orbital::GenS),true);
  for ( uint i = 0; i < virts.size(); ++i )
    term.set_lastorb(Orbital(virts[i].letname(),Orbital::GenS),true);
  if ( _occexcops.size() > 0 ){
    //make sure that we haven't use these orbital names already
    for ( uint i = 0; i < occs.size(); ++i ){
      ipos2=_occexcops.find( Orbital(occs[i].letname(),Orbital::GenS));
      if (ipos2>=0)  // is there, change it
        _occexcops[ipos2] = term.freeorbname(Orbital::Occ);
    }
  }
  if ( _virexcops.size() > 0 ){
    for ( uint i = 0; i < virts.size(); ++i ){
      ipos2=_virexcops.find(Orbital(virts[i].letname(),Orbital::GenS));
      if (ipos2>=0)  // is there, change it
        _virexcops[ipos2] = term.freeorbname(Orbital::Virt);
    }
  }
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0,excl,occs,virts);  
  else
    return Oper(Ops::Exc0,excl,occs,virts);
}

Oper Finput::handle_excitation(const std::string& name, bool dg, Term& term)
{
  unsigned long int ipos,ipos1;
  long int ipos2;
  short excl;
//  if (name.compare(0,3,"\\mu") != 0
//      && name.compare(0,3,"\\nu") != 0
//      && name.compare(0,4,"\\rho") != 0)
//    error("Unknown excitation "+name,"Finput::handle_excitation");
  // find excitation class
  ipos=name.find("_");
  if (ipos==std::string::npos || ipos==name.size()-1)
    error("No excitation class in excitation "+name,"Finput::handle_excitation");
  else
  {
    ++ipos;
    ipos=IL::skip(name,ipos,"{} ");
    ipos1=IL::nextwordpos(name,ipos);
    if(!str2num<short>(excl,name.substr(ipos,ipos1-ipos),std::dec))
      error("Excitation class "+name.substr(ipos,ipos1-ipos),"Finput::handle_excitation");
  }
  ipos2=_excops.find(name);
  Orbital occ,virt;
  if (ipos2<0) {// first run (don't distinguish terms)
    _excops*=name;
    occ=term.freeorbname(Orbital::Occ);
    virt=term.freeorbname(Orbital::Virt);
    _occexcops*=occ;
    _virexcops*=virt;
    _exccls*=excl;
    _spinsymexcs*=Matrices::Singlet; //TODO: implement Triplet!
    _posexcopsterm*=-1; // initialize
  } else {// may be second run...
    occ=_occexcops[ipos2];
    virt=_virexcops[ipos2];
    _posexcopsterm[ipos2]=term.mat().size(); //set position of this operator in the term
  }
  // create \tau_{excl}
  if (dg)
    return Oper(Ops::Deexc0,excl,occ,virt);  
  else
    return Oper(Ops::Exc0,excl,occ,virt);
}

double Finput::handle_factor(const Lelem& lel)
{
  double fac,fac1;
  unsigned long int ipos=0, ipos1;
  std::string lelnam=lel.name();
  if (lel.lex()==Lelem::Num)
  {
    if(!str2num<double>(fac,lelnam,std::dec))
      error("Factor is not a number "+lelnam,"Finput::handle_factor");
  }
  else
  {
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos);
    if(!str2num<double>(fac,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Numerator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
    ipos=lelnam.find("/");
    if(ipos==std::string::npos)
      error("Something wrong with frac "+lelnam,"Finput::handle_factor");
    ++ipos;
    ipos=IL::skip(lelnam,ipos,"{} ");
    ipos1=IL::nextwordpos(lelnam,ipos);
    if(!str2num<double>(fac1,lelnam.substr(ipos,ipos1-ipos),std::dec))
        error("Denominator is not a number "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_factor");
    fac=fac/fac1;
  }
  return fac;
}
Oper Finput::handle_operator(const Lelem& lel, Term& term, bool excopsonly)
{
  unsigned long int ipos, ipos1, iposnam, up, down;
  std::string lelnam=lel.name(),name,nameadd;
  short excl;
  bool dg=false;
  // parts of Hamilton operator
  down=lelnam.find("_");
  up=lelnam.find("^");
  if ( InSet(lelnam[0], hms))
  {
    if (excopsonly) return Oper();
    if (up!=std::string::npos || down!=std::string::npos)
      say("Sub- and superscripts in Hamiltonian will be ignored: "+lelnam);
    if ( lelnam[0]==hms[1] ) return Oper(Ops::Fock);
    if ( lelnam[0]==hms[2] ) return Oper(Ops::FluctP);
    if ( lelnam[0]==hms[3] ) return Oper(Ops::XPert);
  }
  // last position of name of operator
  iposnam=std::min(up,down)-1;
  // excitation class
  if (down==std::string::npos || down==lelnam.size()-1)
    error("No excitation class in operator "+lelnam,"Finput::handle_operator");
  name=lelnam.substr(0,iposnam+1);
  if (up==std::string::npos || up==lelnam.size()-1)
    ; // no dagger
  else
  {
    ipos=up+1;
    ipos=IL::skip(lelnam,ipos,"{} ");
    while((ipos1=IL::nextwordpos(lelnam,ipos))!=ipos)
    {
      if(InSet(lelnam.substr(ipos,ipos1-ipos), dgs))
      {
        dg=true;
        nameadd+=dgs[0];
      }
      else if (lelnam[ipos]!='}')
        nameadd+=lelnam.substr(ipos,ipos1-ipos);
      ipos=ipos1;
    }
  } 
  add2name(name,nameadd); // add nameadd to name (as superscript)
  ipos=down+1;
  ipos=IL::skip(lelnam,ipos,"{} ");
  ipos1=IL::nextwordpos(lelnam,ipos);
  if (InSet(name.substr(0,4), bexcops))
  { // bare excitation operator
    return handle_excitation(lelnam.substr(ipos,ipos1-ipos),dg,term);
  }
  if (excopsonly) return Oper();
  if(!str2num<short>(excl,lelnam.substr(ipos,ipos1-ipos),std::dec))
    error("Excitation class "+lelnam.substr(ipos,ipos1-ipos),"Finput::handle_operator");
  if(dg)
    return Oper(Ops::Deexc,excl,(void*) &term, Term::getfreeorbname,name);
  else
    return Oper(Ops::Exc,excl,(void*) &term, Term::getfreeorbname,name);
}
void Finput::handle_sum(const Lelem& lel, Term& term)
{
  unsigned long int ipos, ipos1, up, down;
  long int iposnam;
  std::string lelnam=lel.name(),name;
  down=lelnam.find("_");
  up=lelnam.find("^");
  if (up!=std::string::npos)
    error("Sum from-to is not implemented yet: "+lelnam);
  if (down==std::string::npos)
    error("Sum without summation indices: "+lelnam);
  ipos=down+1;
  while (ipos<lelnam.size())
  {
    ipos=IL::skip(lelnam,ipos,"{}, ");
    if (ipos==lelnam.size()) break;
    ipos1=IL::nextwordpos(lelnam,ipos);
    name=lelnam.substr(ipos,ipos1-ipos);
    iposnam=_excops.find(name);
    if (iposnam >= 0)
    {
      term.addsummation(_occexcops[iposnam],_exccls[iposnam]);
      term.addsummation(_virexcops[iposnam],_exccls[iposnam]);
    }
    else
      say("No excitation operator which would correspond to summation index "+name);
    ipos=ipos1;
  }
}
void Finput::handle_parameters(Term& term, bool excopsonly)
{
  if (!excopsonly)
  {// handle saved parameters
    unsigned long int ipos, ipos1, iposnam, up, down;
    int iposexcn,indxexcn;
    std::string lelnam,name,nameadd,excn;
    for (unsigned int i=0; i<_paramterm.size(); i++)
    {
      lelnam=_paramterm[i].name();
      down=lelnam.find("_");
      up=lelnam.find("^");
      // last position of name of parameter
      iposnam=std::min(up,down)-1;
      name=lelnam.substr(0,iposnam+1);
      if (up==std::string::npos || up==lelnam.size()-1)
        ; // no superscript
      else
      {
        ipos=up+1;
        ipos=IL::skip(lelnam,ipos,"{} ");
        while((ipos1=IL::nextwordpos(lelnam,ipos))!=ipos)
        {
          if (lelnam[ipos]!='}')
            nameadd+=lelnam.substr(ipos,ipos1-ipos);
          ipos=ipos1;
        }
      } 
      add2name(name,nameadd); // add nameadd to name (as superscript)
      // handle subscript
      if (down==std::string::npos || down==lelnam.size()-1)
      { // no subscript, parameter is a "number"
        term.addmatrix(Matrices(Ops::Number,Product<Orbital>(),name));
      }
      else
      {
        ipos=down+1;
        ipos=IL::skip(lelnam,ipos,"{} ");
        ipos1=IL::nextwordpos(lelnam,ipos);
        excn=lelnam.substr(ipos,ipos1-ipos);
        indxexcn=_excops.find(excn);
        if (indxexcn>=0)
        {
          iposexcn=_posexcopsterm[indxexcn];
          if (iposexcn>=0)
          {
            Matrices mat(Ops::Interm,
                         Ops::genprodorb(_exccls[indxexcn],_occexcops[indxexcn],_virexcops[indxexcn]),
                         name,_spinsymexcs[indxexcn]);
            term.replacematrix(mat,iposexcn);
          }
          else
            say("Parameter is not present in this term: "+lelnam);
        }
        else
          error("Unknown excitation in parameter"+excn);
      }
    }
  }
  // reset all information
  _posexcopsterm.assign(_excops.size(),-1);
  _paramterm=Product<Lelem>();
}
void Finput::add2name(std::string& name, const std::string& nameadd)
{
  unsigned long int ipos,ipos1;
  if (nameadd.size()>0)
  { // add nameadd to name (as superscript)
    ipos=name.rfind("^");
    if (ipos!=std::string::npos)
    { // there is already a superscript
      ++ipos;
      ipos1=IL::nextwordpos(name,ipos);
      name.insert(ipos,"{");
      name.insert(ipos1,nameadd+"}");
    }
    else
    { // no superscript yet
      name += "^{"+nameadd+"}";
    }
  }
}


std::ostream& operator<<(std::ostream& o, const Lelem& lel)
{
  if (lel.lex()==Lelem::Bra) 
    o << "< ";
  else if (lel.lex()==Lelem::Ket)
    o << "| ";
  else if (lel.lex()==Lelem::LPar)
    o << "(";
  else if (lel.lex()==Lelem::RPar)
    o << ")";
  else if (lel.lex()==Lelem::Oper)
    o << "\\"<< Finput::commands[0]<<" ";
  else if (lel.lex()==Lelem::Par)
    o << "\\"<< Finput::commands[1]<<" ";
  else if (lel.lex()==Lelem::Plus)
    o << "+";
  else if (lel.lex()==Lelem::Minus)
    o << "-";
  else if (lel.lex()==Lelem::Times)
    o << "*";
  else if (lel.lex()==Lelem::Div)
    o << "/";
  else if (lel.lex()==Lelem::Sum)
    o << "\\"<< Finput::commands[2];
  
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

