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


lui LelString::closbrack(lui ipos, Lelem::Lex find) const
{
  lui i,ipos1=ipos;
  Lelem::Lex rk=Lelem::RPar,lk=(*this)[ipos].lex();
  if (lk==Lelem::Bra)
    rk=Lelem::Ket;
  else if (lk==Lelem::LPar)
    rk=Lelem::RPar;
  else if (lk==Lelem::LCom)
    rk=Lelem::RCom;
  else
    error("Not a bracket!","Lexic::closbrack"); 
  Lelem::Lex what = rk;
  if ( find != Lelem::None ) what = find;
  int nk=1;
  for ( i=ipos+1;i<this->size();++i) {
    if ( nk == 1 && (*this)[i].lex()==what) {
      --nk;
      ipos1=i;
      break;
    }
    if ((*this)[i].lex()==lk) 
      ++nk; // count number of "("
    else if ((*this)[i].lex()==rk) {
      --nk; // count number of ")"
    }
  }
  if ( nk != 0 ) 
    error("Number of brackets is incosistent: "+any2str(nk),"Lexic::closbrack");
  if ( ipos1 == ipos )
    error("Not found in "+any2str((*this)[ipos]),"Lexic::closbrack");
  return ipos1;
}
lui LelString::openbrack(lui ipos) const
{
  lui ipos1=ipos;
  Lelem::Lex lk=Lelem::LPar,rk=(*this)[ipos].lex();
  if (rk==Lelem::Ket)
    lk=Lelem::Bra;
  else if (rk==Lelem::RPar)
    lk=Lelem::LPar;
  else if (rk==Lelem::RCom)
    lk=Lelem::LCom;
  else
    error("Not a bracket!","Lexic::openbrack"); 
  int nk=-1;
  for ( long int i=ipos-1;i>=0;--i) {
    if ((*this)[i].lex()==lk) {
      ++nk; // count number of "("
      if (nk==0) {
        ipos1=i;
        break;
      }
    }
    else if ((*this)[i].lex()==rk) 
      --nk; // count number of ")"
  }
  if ( nk != 0 ) 
    error("Number of brackets is incosistent: "+any2str(nk),"Lexic::openbrack"); 
  return ipos1;
}
Product< long int > LelString::addconnections(lui beg, lui end) const
{
  Product<long int> connection;
  if ((*this)[end].conn()==Lelem::Normal) return connection;
  for (lui i=beg+1; i<end; i++)
    if ((*this)[i].lex()==Lelem::Oper) {
      if((*this)[end].conn()==Lelem::Connect)
        connection *= i+1;
      else
        connection *= -(i+1);
    }
  return connection;
}
LelString LelString::expandnewops(const NewOpMap& newops) const
{
  LelString result;
  for (lui i = 0; i < this->size(); ++i ){
    result.add((*this)[i]);
    lui j = result.size()-1;
    do { 
      if ( result[j].lex() == Lelem::Oper ) {
        NewOpMap::const_iterator itnewop = newops.find(result[j].name());
        if ( itnewop != newops.end() ) {
          const LelString& nop = itnewop->second;
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
lui LelString::elem(lui beg, bool bk) const
{
  lui i, end, nk=0, ncom=0;
  bool braket=false;
  
  for ( i = beg; i < this->size() ; i++ ) {
    if ((*this)[i].lex()==Lelem::LPar) ++nk;
    if ((*this)[i].lex()==Lelem::RPar) --nk;
    if ((*this)[i].lex()==Lelem::LCom) ++ncom;
    if ((*this)[i].lex()==Lelem::RCom) --ncom;
    if (bk) {
      if ((*this)[i].lex()==Lelem::Bra) braket=true;
      if ((*this)[i].lex()==Lelem::Ket) braket=false;
    }
    if (InSet((*this)[i].lex(), Lelem::Plus,Lelem::Minus) && 
        i!= beg && nk==0 && ncom==0 && !braket ) 
      break; // end of term
  }
  if (nk!=0 || ncom!=0 || braket) error("Check input, Mismatch in parentheses or bra/ket","Lexic::elem");
  if (i==beg || i==0 ) 
    end=0;
  else 
    end=i-1;
  return end;
}
LelString LelString::expandcom(lui beg) const
{ // e.g., [a,b] --> (a)(b)-(b)(a)
  lui 
    end = this->closbrack(beg),
    comma = this->closbrack(beg, Lelem::Comma);
  LelString result, 
    beforecomma = this->substring(beg+1,comma-1),
    aftercomma = this->substring(comma+1,end-1);
  result.add(Lelem("",Lelem::LPar)); // (
  result.add(Lelem("",Lelem::LPar)); //  (
  result.add(beforecomma);           //   a 
  result.add(Lelem("",Lelem::RPar)); //  )
  result.add(Lelem("",Lelem::LPar)); //  (   
  result.add(aftercomma);            //   b
  result.add(Lelem("",Lelem::RPar)); //  )
  result.add(Lelem("",Lelem::Minus));//  - 
  result.add(Lelem("",Lelem::LPar)); //  (
  result.add(aftercomma);            //   b
  result.add(Lelem("",Lelem::RPar)); //  )
  result.add(Lelem("",Lelem::LPar)); //  (
  result.add(beforecomma);           //   a
  result.add(Lelem("",Lelem::RPar)); //  )
  result.add(Lelem("",Lelem::RPar)); // )
  return result; 
}
LelString LelString::expandpar(lui beg, ConnectionsMap& connections) const
{ // e.g., this=-a(b+c)d
  lui i=0, end,ipos,start=0,ijcon,iposres,lenb,ipar;//lena,
  end=this->closbrack(beg);
  LelString result, inpar=this->substring(beg+1,end-1);
  ConnectionsMap con(connections);
  con.push_back(this->addconnections(beg,end));
  connections = ConnectionsMap();
  Product<long int> connect;
  bool coninterm,coninpar;
  if ((*this)[beg].lex()==Lelem::Bra) 
    ipar=1;
  else
    ipar=0;
  // start of term (without sign)
  if (InSet(this->begin()->lex(), Lelem::Plus,Lelem::Minus)) start=1;
  const Lelem& el0 = this->front();
  while (i<inpar.size()) {
    ipos=inpar.elem(i);
    //add sign
    if (inpar[i].lex()==Lelem::Minus) {// minus in parenthesis
      if (el0.lex()==Lelem::Minus) 
        result.add(Lelem("",Lelem::Plus)); // -1*-1 == +1
      else
        result.add(inpar[i]); // +1*-1 == -1
      ++i;
    } else if (inpar[i].lex()==Lelem::Plus) {// plus in parenthesis
      if (el0.lex()==Lelem::Minus) 
        result.add(el0); // -1*+1 == -1
      else
        result.add(inpar[i]); // +1*+1 == +1
      ++i;
    } else if (start>0) {// term has a sign, but no sign by expression in parenthesis
      result.add(el0);
    }
    iposres=result.size()-start;
    // add "a" (without sign)
    if (beg>start) result.add(this->substring(start,beg-1));
    if ((*this)[beg].lex()==Lelem::Bra) result.add((*this)[beg].braexpanded());
    // add "b" (or "c")
    result.add(inpar.substring(i,ipos));
    // add "d"
    if ((*this)[end].lex()==Lelem::Ket) result.add((*this)[end]);
    if (end < this->size()-1) result.add(this->substring(end+1,this->size()));
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
LelString LelString::expandterm(ConnectionsMap& connections) const
{
//  std::cout << "Term: " << *this << std::endl;
//  for (unsigned int k = 0; k<connections.size();k++)
//    std::cout << "in connection " << k << ": " << connections[k] << std::endl;
  lui i = 0;
  _foreach_cauto(LelString,itl,*this){
    if (itl->lex()==Lelem::LPar || (itl->lex()==Lelem::Bra && !itl->expandedbra())) 
      return this->expandpar(i,connections);
    ++i;
  }
  return *this;
}
bool LelString::expanded() const
{ 
  _foreach_cauto(LelString,itl,*this) {
    if (itl->lex()==Lelem::Bra && !itl->expandedbra())
      return false;
    if (itl->lex()==Lelem::LPar)
      return false;
  }
  return true;
}
bool LelString::expanded_com() const
{ 
  _foreach_cauto(LelString,itl,*this) {
    if (itl->lex()==Lelem::LCom )
      return false;
  }
  return true;
}
void LelString::expand( ConnectionsMap& connections)
{
  LelString res;
  Product<long int>  connect;
  ConnectionsMap con,con1;
  lui beg=0, end, i,j,conbeg,lastpos;
  
  while (!expanded()) {
    res=*this;
    *this = LelString();
    beg=0;
    conbeg=0;
    con=connections;
    connections = ConnectionsMap();
    while ((end=res.term(beg))!=0) {
      ConnectionsMap con1;
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
      lastpos=this->size();
      this->add(res.substring(beg,end).expandterm(con1));
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
}
void LelString::expand_commutators()
{
  LelString res;
  lui beg=0, end;
  while (!expanded_com()) {
    res = *this;
    *this = LelString();
    for ( beg = 0; beg < res.size(); ++beg ) {
      if ( res[beg].lex() == Lelem::LCom ) {
        end=res.closbrack(beg);
        this->add(res.substring(beg,end).expandcom(0));
        beg=end;
      } else {
        this->add(res[beg]);
      }
    }
  }
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
  else if (lel.lex()==Lelem::LCom)
    o << "[";
  else if (lel.lex()==Lelem::RCom)
    o << "]";
  else if (lel.lex()==Lelem::Comma)
    o << ",";
  else if (lel.lex()==Lelem::Oper)
    o << "\\"<< commands["operator"]<<" ";
  else if (lel.lex()==Lelem::Tensor)
    o << "\\"<< commands["tensor"]<<" ";
  else if (lel.lex()==Lelem::Perm)
    o << "\\"<< commands["permutation"];
  else if (lel.lex()==Lelem::Equal)
    o << "=";
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
