#ifndef GLOBAL_H 
#define GLOBAL_H

#include <iostream>
#include <sstream>

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
// global variables and functions from input
namespace Input
{
  // result in spin or space orbitals
  extern bool dospinintegr;
  // terms with prefactor smaller than this will be neglected
  extern double minfac;
}
 
#endif
