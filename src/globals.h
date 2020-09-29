#ifndef GLOBAL_H 
#define GLOBAL_H

#include <iostream>
#include <sstream>
#include <map>
#include <list>
#include <cmath>

typedef long unsigned int lui;
typedef unsigned int uint;
typedef std::list< std::string > TParArray;
typedef std::map< std::string, std::string > TsPar;
typedef std::map< std::string, int > TiPar;
typedef std::map< std::string, double > TfPar;
typedef std::map< std::string, TParArray > TaPar;

typedef std::map< std::string, TsPar > TsParSet;
typedef std::map< std::string, TiPar > TiParSet;
typedef std::map< std::string, TfPar > TfParSet;
typedef std::map< std::string, TaPar > TaParSet;

namespace math {
// greatest common divisor
long int gcd(long int n1, long int n2);
}
// rational numbers
class TRational 
{
public:
  TRational() {};
  TRational(long int n,long int d = 1) : _numerator(n), _denominator(d){normalize();};
//   TRational(const TRational& f) : _numerator(f._numerator),_denominator(f._denominator){};
  long int numerator() const { return _numerator; };
  long int denominator() const { return _denominator; };
  
  TRational & operator*=(const TRational& f){ _numerator*= f._numerator; _denominator*= f._denominator; normalize(); return *this;};
  TRational & operator/=(const TRational& f){ _numerator*= f._denominator; _denominator*= f._numerator; normalize(); return *this;};
  TRational & operator+=(const TRational& f){ _numerator*= f._denominator; _numerator+= f._numerator*_denominator; _denominator*= f._denominator; normalize(); return *this;};
  TRational & operator-=(const TRational& f){ _numerator*= f._denominator; _numerator-= f._numerator*_denominator; _denominator*= f._denominator; normalize(); return *this;};
  
  TRational operator*(const TRational& f) const { return TRational(_numerator*f._numerator,_denominator*f._denominator);};
  TRational operator/(const TRational& f) const { return TRational(_numerator*f._denominator,_denominator*f._numerator);};
  TRational operator+(const TRational& f) const { return TRational(_numerator*f._denominator+f._numerator*_denominator,_denominator*f._denominator);};
  TRational operator-(const TRational& f) const { return TRational(_numerator*f._denominator-f._numerator*_denominator,_denominator*f._denominator);};
  TRational operator-() const { return TRational(-_numerator,_denominator);};
  
  bool operator<(const TRational& f) const { return (_numerator*f._denominator < f._numerator*_denominator);};
  bool operator<=(const TRational& f) const { return (_numerator*f._denominator <= f._numerator*_denominator);};
  bool operator>(const TRational& f) const { return (_numerator*f._denominator > f._numerator*_denominator);};
  bool operator>=(const TRational& f) const { return (_numerator*f._denominator >= f._numerator*_denominator);};
  bool operator==(const TRational& f) const { return (_numerator*f._denominator == f._numerator*_denominator);};
  bool operator!=(const TRational& f) const { return (_numerator*f._denominator != f._numerator*_denominator);};
private:
  long int _numerator = 0, _denominator = 1;  
  void normalize() {long int div = math::gcd(_numerator,_denominator); 
                    if (_denominator < 0) div = -div; 
                    _numerator /= div; _denominator /= div; };
};
TRational operator/(long int i, const TRational& f);
namespace math{
TRational abs(const TRational& f);
double todouble(const TRational& f);
}
std::ostream & operator << (std::ostream & o, TRational const & p);

#ifdef _RATIONAL
typedef TRational TFactor;
#define _abs math::abs
#define _todouble math::todouble
#else
typedef double TFactor;
#define _abs std::abs
#define _todouble
#endif

namespace Numbers
{
  // "zero"
  static const double verysmall=1.e-18;
  static const double small=1.e-16;
  static const int big=1000;
}

class Return{
public:
  enum Vals{
    Done = 0,
    Delete = 1,
    Change_sign = 2,
    Repeat = 4
  };
  Return(Vals val = Done) : _val(val){};
  Return & operator +=(const Return& ret);
  char _val;
};

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

#define _foreach(It,Array) for ( (It) = (Array).begin(); (It) != (Array).end(); ++(It) )
#define _foreach_auto(It,Array) for ( auto It = (Array).begin(); It != (Array).end(); ++It )
#define _foreach_rauto(It,Array) for ( auto It = (Array).rbegin(); It != (Array).rend(); ++It )
#define _foreach_cauto(It,Array) for ( auto It = (Array).begin(); It != (Array).end(); ++It )
#define _foreach_crauto(It,Array) for ( auto It = (Array).rbegin(); It != (Array).rend(); ++It )

// print timing
#define _CPUtiming(what,start,end) xout << std::fixed << std::setprecision(2) << "CPU time " << what << 1000.0*(end-start)/CLOCKS_PER_SEC << " ms\n";

#endif
