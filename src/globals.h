#ifndef GLOBAL_H 
#define GLOBAL_H

#include <iostream>
#include <sstream>
#include <map>
#include <list>
#include <cmath>
#include "utilities.h"

typedef long unsigned int lui;
typedef std::list< std::string > TParArray;
typedef std::map< std::string, std::string > TsPar;
typedef std::map< std::string, int > TiPar;
typedef std::map< std::string, double > TfPar;
typedef std::map< std::string, TParArray > TaPar;

typedef std::map< std::string, TsPar > TsParSet;
typedef std::map< std::string, TiPar > TiParSet;
typedef std::map< std::string, TfPar > TfParSet;
typedef std::map< std::string, TaPar > TaParSet;


namespace Numbers
{
  static const double small=1.e-16;
}
// global variables and functions for output
class Output
{
  public:
  Output();
  Output(std::ostream & o);
  void setvals();
  // 1+small == 1 in output
  double small;
  // length of current line
  int lenline;
  // length of buffer
  int lenbuf;
  // number of lines on page
  double nlines;
  // number of pages
  int npages;
  // height of the current line
  double hline;
  // height of the current line in buffer
  double hbufline;
  // output
  std::ostream * pout;
  // buffer, will be flushed after newline or newpage
  std::ostringstream buf;
  // max length of line
  int maxlenline;
  // max number of lines on page
  int maxnlines;
  // how many super and subscript indices one needs to fill a position
  int wsi;
  // height of \sum
  double hsum;
  // start new page of equations
  void newpageeqn();
  // new line of equations
  void newlineeqn();
  // begin equation
  void beq();
  // end equation (+flushbuf)
  void eeq(bool doflush=true);
  // flash buffer
  void flushbuf();
  // if true: this is in an equation
  bool inequation;
  // break line if to long
  bool breaklongline();
};
namespace MyOut
{
  extern Output defout;
  extern Output * pcurout;
}

// global variables and functions from input
namespace Input
{
  // verbosity
  extern int verbose;
  // input-parameters
  //string parameters
  extern TsParSet sPars;
  //integer parameters
  extern TiParSet iPars;
  //float parameters
  extern TfParSet fPars;
  //array of string parameters
  extern TaParSet aPars;
}

#define xout std::cout
#define _xout(level,x) \
    { \
       if ( Input::verbose >= (level) ) \
       { xout << x; } \
    }
// output unless Input::verbose is below 0
#define _xout0(x) _xout(0,x)
// output unless Input::verbose is below 1
#define _xout1(x) _xout(1,x)
// output unless Input::verbose is below 2
#define _xout2(x) _xout(2,x)
// output unless Input::verbose is below 3
#define _xout3(x) _xout(3,x)

// class Inpars
// {
// public:
//   static Inpars * getInpars();
//   static void setInparsFileName(std::string filename){ m_filename = filename; };
// private:
//   Inpars();
//   static std::string m_filename;
//   static Inpars * m_configInstance;
//   // input-parameters
//   typedef std::map< std::string, std::string > TsInpars;
//   typedef std::map< std::string, int > TiInpars;
//   typedef std::map< std::string, double > TfInpars;
//   TsInpars sInpars;
//   TiInpars iInpars;
//   TfInpars fInpars;
// };
#endif
