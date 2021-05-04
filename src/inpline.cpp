#include "inpline.h"

std::string IL::key(const std::string& line, lui& ipos, const std::string& keyword)
{
  do {
    ipos = line.find(keyword);
    if ( ipos != std::string::npos ){
      ipos += keyword.size();
      ipos = skip(line,ipos," ");
      if ( line[ipos] == '=' ){
        ++ipos;
        lui
          ipend = endword(line,ipos),
          ipbeg = ipos;
        ipos = skip(line,ipend," \"},");
        return line.substr(ipbeg,ipend-ipbeg);
      }
    }
  } while (ipos != std::string::npos);
  return std::string();
}
std::string IL::key(const std::string& line, const std::string& keyword)
{
  lui ipos = 0;
  return IL::key(line,ipos,keyword);
}

TParArray IL::parray(const std::string& str)
{
  std::string listsep = Input::sPars["syntax"]["listseparator"];
  TParArray res;
  lui
    ipos = 0,
    ipend = endword(str,ipos,listsep);

  while( ipend != ipos ){
    res.push_back(str.substr(ipos,ipend-ipos));
    ipos = ipend;
    if ( ipos < str.size() ) ++ipos;
    ipos = skip(str,ipos,listsep);
    ipend = endword(str,ipos,listsep);
  }
  return res;
}
lui IL::addnewcom(const std::string& str, lui ipos, std::string what)
{
  if ( str[ipos] != '{' )
    error("bad \\"+what+" "+str);
  lui ipend = closbrack(str,ipos);
  ++ipos;
  ipos=skip(str,ipos," \\"); //we don't have backslash in command-names
  if ( ipos == ipend )
    error("empty name in \\"+what+" "+str);
  std::string name(str.substr(ipos,skipr(str,ipend," ")-ipos));
  ipos = ipend+1;
  if ( ipos >= str.size() || str[ipos] != '{' )
    error("bad \\"+what+" "+str);
  ipend = closbrack(str,ipos);
  ++ipos;
  ipos=skip(str,ipos," ");
  if ( ipos == ipend )
    error("empty value in \\"+what+" "+str);
  std::string value(str.substr(ipos,skipr(str,ipend," ")-ipos));
  Input::sPars[what][name] = value;
  return ipend+1;
}
void IL::changePars(const std::string& str, lui ipos)
{
  bool stype, itype, ftype, atype;
  stype = itype = ftype = atype = false;
  lui ipend = IL::nextwordpos(str,ipos);
  if ( ipend == ipos ) return; // empty line
  std::string set = str.substr(ipos,ipend-ipos);
  // search set in parameters
  if ( Input::sPars.count(set) )
    stype = true;
  if ( Input::iPars.count(set) )
    itype = true;
  if ( Input::fPars.count(set) )
    ftype = true;
  if ( Input::aPars.count(set) )
    atype = true;
  if ( !(stype || itype || ftype || atype) ) {
    // print a warning and return
    _xout0("Unrecognized set: " << str << std::endl);
    return;
  }
  ipos = IL::skip(str,ipend," :,");
  // go through all names in the line
  while( (ipend = IL::nextwordpos(str,ipos)) != ipos ){
    // search name in parameter-set
    std::string name = str.substr(ipos,ipend-ipos);
    if ( stype && ( Input::sPars[set].count(name) || set.compare(0,3,"new") == 0 ) ){
      Input::sPars[set][name] = IL::key(str,ipos,name);
    } else if ( itype && Input::iPars[set].count(name)){ // change the parameter
      int x;
      if (!str2num<int>(x,IL::key(str,ipos,name),std::dec))
        error("integer input parameter is not integer in"+set+":"+name);
      Input::iPars[set][name] = x;
    } else if ( ftype && Input::fPars[set].count(name)){ // change the parameter
      double x;
      if (!str2num<double>(x,IL::key(str,ipos,name),std::dec))
        error("float input parameter is not float in "+set+":"+name);
      Input::fPars[set][name] = x;
    } else if ( atype && Input::aPars[set].count(name)){
      Input::aPars[set][name] = IL::parray(IL::key(str,ipos,name));
    } else {
      _xout0(str.substr(ipos) << std::endl);
      error(set+":"+name+" : unrecognized name");
    }
  }
}

lui IL::skip(const std::string& str, const lui& ipos, const std::string& what)
{
  lui ires=ipos;
  while (ires<str.size()&& what.find(str[ires])!=std::string::npos )
    ++ires;
  return ires;
}
lui IL::skipr(const std::string& str, const lui& ipos, const std::string& what)
{
  lui ires = std::min(std::size_t(ipos),str.size());
  while (ires > 0 && what.find(str[ires-1])!=std::string::npos )
    --ires;
  return ires;
}
void IL::delbrack(std::string& str, lui ipos, std::string brackets)
{
  lui endstr = str.size();
  for ( ; ipos < endstr; ++ipos, --endstr ) {
    //skip empty spaces
    ipos = skip(str,ipos," ");
    endstr = skipr(str,endstr," ");
    if ( brackets.find(str[ipos]) == std::string::npos ||
         closbrack(str,ipos)+1 != endstr ) {
      break;
    }
  }
  if ( str.size() > ipos ) str = str.substr(ipos,endstr-ipos);
}

lui IL::endword(const std::string& line, lui& ipos, std::string separ)
{
//  assert(ipos < line.size());
  ipos = skip(line,ipos," ");
  if ( separ.find(" ") < separ.size() ) {
    // if space is one of the separators we have to skip all separators!
    ipos = skip(line,ipos,separ);
  }
  if ( ipos >= line.size() ) return ipos;
  std::string end;
  if ( line[ipos] == '"' ){
    end = '"'; ++ipos;
  }else if ( line[ipos] == '{' ){
    end = '}'; ++ipos;
  }
  if ( end.empty() ) {
    lui iposres = std::string::npos;
    for (const auto& is: separ ){
      iposres = std::min(iposres, lexfind(line,std::string(1,is),ipos));
    }
    return iposres;
  }
  return lexfind(line,end,ipos);
}
lui IL::closbrack(const std::string& str, lui ipos)
{
  const std::string& brackets = Input::sPars["syntax"]["brackets"];
  lui i=brackets.find(str[ipos]),ipos1=ipos;
  if (i==std::string::npos)
    error(any2str(str[ipos])+"is not a bracket!","IL::closbrack");
  char lk(brackets[i]), rk(brackets[i+1]); // left and right brackets
  int nk=1;
  bool backslashed = false;
  for ( i=ipos+1;i<str.size();++i) {
    if ( backslashed ) {
      // don't count backslashed brackets
    } if (str[i]==lk) {
      ++nk; // count number of "("
    } else if (str[i]==rk) {
      --nk; // count number of ")"
      if (nk==0) {
        ipos1=i;
        break;
      }
    }
    backslashed = !backslashed && ( str[i] == '\\' );
  }
  if ( nk != 0 )
    error("Number of brackets is incosistent: "+any2str(nk),"IL::closbrack");
  return ipos1;
}

lui IL::nextwordpos(const std::string& str, lui& ipos, bool glue, bool greedy)
{
  const std::string& separator = Input::sPars["syntax"]["separator"];
  const std::string& gluer = Input::sPars["syntax"]["gluer"];
  lui nwpos;
  ipos=IL::skip(str,ipos," "); // remove spaces
  if ( ipos < str.size() && str[ipos] == '%' ) return ipos+1; // comment sign is one word

  bool
    glued = false,
    is_command = false,
    backslashed = false,
    breaknext = false;
  for (nwpos=ipos ;nwpos < str.size(); ++nwpos){
    bool
      is_separator = (separator.find(str[nwpos]) != std::string::npos),
      is_gluer = (gluer.find(str[nwpos]) != std::string::npos);
    if ( !glue && is_gluer ) {
      // gluer is a separator, too
      is_separator = true;
      is_gluer = false;
    }
    if ( is_gluer ) breaknext = false;
    if ( breaknext ) break;

    if ( nwpos == ipos || glued ) {
      if ( glued && is_gluer ) error("Two consecutive gluers in "+str,"IL::nextwordpos");
      if ( char(str[nwpos])=='{' ) {
        nwpos=closbrack(str,nwpos);
        breaknext = true;
      } else if ( char(str[nwpos])=='\\' ) {
        backslashed = true;
      } else if ( is_separator || !greedy ) {
        breaknext = true;
      }
    } else if ( backslashed ) {
      backslashed = false;
      if ( is_separator || is_gluer ) {
        breaknext = true;
      } else {
        is_command = true;
      }
    } else if ( is_separator ) {
      break;
    } else if ( is_command ) {
      if ( is_gluer ) {
        is_command = false;
      }
    } else if ( !greedy && !is_gluer ) {
      breaknext = true;
    }

    glued = is_gluer;
  }
  return nwpos;
}
std::string IL::plainname(std::string name)
{
  std::string pname;
  for ( uint i = 0; i < name.size(); ++i )
    if ( !InSet(name[i],"_^{}") )
      pname += name[i];
  return pname;
}
lui IL::lexfind(const std::string& str, const std::string& sstr, const lui& ipos)
{
  lui ipos1, ires(ipos), icur(0);
  int level = 0;
  bool found(false);
  for ( ipos1 = ipos; ipos1 < str.size() && !found; ++ipos1 ){
    if ( level == 0 ){
      // search on the current level only!
      if ( str[ipos1] == sstr[icur] ){
        if ( icur == 0 ) ires = ipos1;
        ++icur;
        if ( icur == sstr.size() ) found = true;
      } else if ( icur > 0 ){
        icur = 0;
        ipos1 = ires;
        continue;
      }
    }
    if ( str[ipos1] == '{' ) ++level;
    if ( str[ipos1] == '}' ) --level;
  }
  if (!found)
    ires = std::string::npos;
  return ires;
}
bool IL::nameupdown(std::string& name, std::string& nameup, std::string& namedown, const std::string& lelnam)
{
  if ( lelnam.empty() ) return false;
  lui down = IL::lexfind(lelnam,"_");
  lui up = IL::lexfind(lelnam,"^");
  // last position of name of operator
  lui iposnam = std::min(up,down);
  iposnam = std::min(std::size_t(iposnam),lelnam.size());
  name = lelnam.substr(0,iposnam);
  if (up < lelnam.size()-1){
    nameup = lelnam.substr(up);
  }
  if (down < lelnam.size()-1){
    namedown = lelnam.substr(down);
  }
  if (up < down){
    nameup = nameup.substr(0,down-up);
  } else if (down < up) {
    namedown = namedown.substr(0,up-down);
  }
  // get rid of ^ _ and left-right {}
  delbrack(nameup,1);
  delbrack(namedown,1);
  return !(nameup.empty() && namedown.empty());
}
void IL::add2name(std::string& name, const std::string& nameadd, bool superscript, bool snam)
{
  const std::string& supername = Input::sPars["command"]["supername"]; // environment for name addition
  unsigned long int ipos,ipos1;
  std::string ch("^");
  if (!superscript) ch = "_";
  if (nameadd.size()>0) { // add nameadd to name (as superscript)
    ipos=IL::lexfind(name,ch);
    if (ipos!=std::string::npos) { // there is already a superscript
      ++ipos;
      if ( snam && name.compare(ipos,supername.size()+1,"\\"+supername))
        ipos += supername.size()+2;
      ipos1=IL::nextwordpos(name,ipos,false);
      name.insert(ipos,"{");
      ++ipos1;
      name.insert(ipos1,nameadd+"}");
    } else { // no superscript yet
      if (snam)
        name += ch + "{\\"+ supername + "{" + nameadd + "}}";
      else
        name += ch + "{" + nameadd + "}";
    }
  }
}
