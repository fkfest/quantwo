#ifndef GLOBAL_H 
#define GLOBAL_H

#include <iostream>
#include <sstream>
#include <map>

typedef long unsigned int lui;

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
  static const int maxlenline=80;
  // max number of lines on page
  static const int maxnlines=42;
  // how many super and subscript indices one needs to fill a position
  static const int wsi=1;
  // height of \sum
  static const double hsum=2;
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
  typedef std::map< std::string, std::string > TsInpars;
  typedef std::map< std::string, int > TiInpars;
  typedef std::map< std::string, double > TfInpars;
// global variables and functions from input
namespace Input
{
  // terms with prefactor smaller than this will be neglected
  extern double minfac;
  // input-parameters
  extern TsInpars sInpars;
  extern TiInpars iInpars;
  extern TfInpars fInpars;
}

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
