#include "globals.h"

Output::Output()
{
  maxlenline = 0;
  maxnlines = 0;
  wsi = 0;
  hsum = 0.0;
  pout = &std::cout;
  setvals();
}

Output::Output(std::ostream& o)
{ 
  maxlenline = Input::iPars["output"]["maxlenline"];
  maxnlines = Input::iPars["output"]["maxnlines"];
  wsi = Input::iPars["output"]["widthsubidx"];
  hsum = Input::fPars["output"]["sumheight"];
  pout=&o;
  setvals();
}
void Output::setvals()
{
  if (maxlenline == 0) maxlenline = 80;
  if (maxnlines == 0) maxnlines = 42;
  if (wsi == 0) wsi = 1;
  if (hsum < Numbers::small) hsum = 1.8;
  lenline=lenbuf=nlines=npages=0; 
  hline=hbufline=1.0;
  inequation=false;
  small = std::pow(10,-pout->precision());
}

void Output::flushbuf()
{
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
  lenline=0;
  hline=1.0;
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

Output MyOut::defout;
Output * MyOut::pcurout = &MyOut::defout;
double Input::minfac=1e-10;
TsParSet Input::sPars;
TiParSet Input::iPars;
TfParSet Input::fPars;
TaParSet Input::aPars;
