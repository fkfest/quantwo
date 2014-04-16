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
  TParArray res;
  lui 
    ipos = 0,
    ipend = endword(str,ipos);
  while( ipend != ipos ){
    res.push_back(str.substr(ipos,ipend-ipos));
    ipos = ipend+1;
    ipos = skip(str,ipos," ,");
    ipend = endword(str,ipos);
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
    if ( stype && Input::sPars[set].count(name) ){
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
  lui ires = std::min(ipos,str.size());
  while (ires > 0 && what.find(str[ires-1])!=std::string::npos )
    --ires;
  return ires;
}
lui IL::endword(const std::string& line, lui& ipos)
{
//  assert(ipos < line.size());
  ipos = skip(line,ipos," ");
  if ( ipos >= line.size() ) return ipos;
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
  const std::string& brackets = Input::sPars["syntax"]["brackets"];
  lui i=brackets.find(str[ipos]),ipos1=ipos;
  if (i==std::string::npos) 
    error(any2str(str[ipos])+"is not a bracket!","IL::closbrack"); 
  char lk(brackets[i]), rk(brackets[i+1]); // left and right brackets
  int nk=1;
  for ( i=ipos+1;i<str.size();++i) {
    if (str[i]==lk) 
      ++nk; // count number of "("
    else if (str[i]==rk) {
      --nk; // count number of ")"
      if (nk==0) {
        ipos1=i;
        break;
      }
    }
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
void IL::add2name(std::string& name, const std::string& nameadd, bool superscript)
{
  unsigned long int ipos,ipos1;
  std::string ch("^");
  if (!superscript) ch = "_";
  if (nameadd.size()>0) { // add nameadd to name (as superscript)
    ipos=IL::lexfind(name,ch);
    if (ipos!=std::string::npos) { // there is already a superscript
      ++ipos;
      ipos1=IL::nextwordpos(name,ipos,false);
      name.insert(ipos,"{");
      ++ipos1;
      name.insert(ipos1,nameadd+"}");
    } else { // no superscript yet
      name += ch + "{" + nameadd + "}";
    }
  }
}
