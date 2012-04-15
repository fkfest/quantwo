#include "globals.h"

Output::Output()
{
  lenline=nlines=npages=0; 
  hbufline=1.0;
  inequation=false;
  pout = &std::cout;
  small = std::pow(10,-pout->precision());
}

Output::Output(std::ostream& o)
{ 
  lenline=lenbuf=nlines=npages=0; 
  hbufline=1.0;
  inequation=false;
  pout=&o;
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

CurOut* CurOut::pInstance(0);
CurOut* CurOut::Create(Output* pOut_)
{
  if ( pInstance ){ //destroy the previous instance
//    if (delout) delete pOut;
    delete pInstance;
  }
  pInstance = new CurOut(pOut_);
  return pInstance;
}
CurOut* CurOut::Instance()
{
  if (pInstance == 0)
    pInstance = new CurOut();
  return pInstance;
}

Output MyOut::defout;
Output * MyOut::pcurout = &MyOut::defout;
double Input::minfac=1e-10;
TsParSet Input::sPars;
TiParSet Input::iPars;
TfParSet Input::fPars;
TaParSet Input::aPars;
