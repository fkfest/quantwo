#include "globals.h"

Return& Return::operator+=(const Return& ret)
{
  switch (ret._val){
    case Done:
    case Delete:
    case Repeat:
      _val = _val | ret._val;
      break;
    case Change_sign:
      _val = _val ^ ret._val;
  }
  return *this;
}

Output::Output()
{
  pout = &std::cout;
  setvals();
}

Output::Output(std::ostream& o)
{
  pout=&o;
  setvals();
}
void Output::setvals()
{
  maxlenline = Input::iPars["output"]["maxlenline"];
  maxnlines = Input::iPars["output"]["maxnlines"];
  wsi = Input::iPars["output"]["widthsubidx"];
  hsum = Input::fPars["output"]["sumheight"];
  if (maxlenline == 0) maxlenline = 80;
  if (maxnlines == 0) maxnlines = 42;
  if (wsi == 0) wsi = 1;
  if (hsum < Numbers::small) hsum = 1.8;
  lenline=lenbuf=nlines=npages=0;
  hline=hbufline=1.0;
  inequation=false;
  small = std::pow(10,-pout->precision());
}

void Output::flushbuf(bool linebreak)
{
  uint offset = 60;
  std::list<char> pm = {'+', '-', '('};
  if(buf.str().size() > maxlenline+offset){
    std::streampos ipos0, ipos;
    ipos0 = 0;
    ipos = ipos0;
    for(char& c : buf.str()){
      if(ipos > maxlenline+offset && InSet(c,pm) && linebreak){
        *pout <<"\\nl"<<std::endl;
        lenline=0;
        nlines+=hline;
        hline=1.0;
        ipos = 0;
      }
      *pout << c;
      ipos += 1;
    }
  }
  else
    *pout << buf.str();
  buf.str("");
  hline=std::max(hline, hbufline);
  hbufline=1.0;
  lenline += lenbuf;
  lenbuf=0;
}
void Output::newpageeqn()
{
  *pout << " \\newpg" << std::endl;
  eeq(false);
  beq();
  nlines=0;
  ++npages;
  flushbuf();
}
void Output::newlineeqn()
{
  *pout <<"\\nl"<<std::endl;
  lenline=0;
  nlines+=hline;
  hline=1.0;
  flushbuf();
}
void Output::beq()
{
  *pout << "\\beq" << std::endl;
  *pout << "&&";
  inequation=true;
}
void Output::eeq(bool doflush)
{
  if (doflush) flushbuf();
  *pout << "\\eeq" << std::endl;
  inequation=false;
  lenline=0;
  hline=1.0;
}
bool Output::breaklongline()
{
  if (buf.str().size()==0) return false; // nothing in buffer
  if (inequation && lenline>0 && lenline+lenbuf > maxlenline)
    {// new line
      if (int(nlines)>maxnlines)
        newpageeqn();
      else
        newlineeqn();
      return true;
    }
  return false;
}

namespace math {
long int gcd(long int n1, long int n2) {
  if ( n2 == 0 || n2 == 1 || n1 == 1 ) return 1;
  n1 = std::abs(n1);
  n2 = std::abs(n2);
  // Euclidean algorithm
#if true
  long int tmp;
  while ( n2 != 0 ) {
    tmp = n2;
    n2 = n1 % n2;
    n1 = tmp;
  }
#else
  while ( n1 != n2 ) {
    if (n1 > n2)
      n1 -= n2;
    else
      n2 -= n1;
  }
#endif
  return n1;
}
TRational abs(const TRational& f){return TRational(std::abs(f.numerator()),f.denominator());}
double todouble(const TRational& f){return double(f.numerator())/double(f.denominator());}
}
TRational operator/(long int i, const TRational& f){return TRational(i*f.denominator(),f.numerator());}
std::ostream & operator << (std::ostream & o, TRational const & p){
  int digits = 0, number = std::max(std::abs(p.numerator()),std::abs(p.denominator()));
  if ( p.numerator() < 0 || p.denominator() < 0 ) ++digits;
  do { number /= 10; ++digits; } while(number);
  if ( p.denominator() != 1 )
    o << "\\frac{" << p.numerator() << "}{" << p.denominator() << "}" ;
  else
    o << p.numerator();
  MyOut::pcurout->lenbuf += digits;
  return o;
}
int Input::verbose = 0;
TsParSet Input::sPars;
TiParSet Input::iPars;
TfParSet Input::fPars;
TaParSet Input::aPars;
Output MyOut::defout;
Output * MyOut::pcurout = &MyOut::defout;
