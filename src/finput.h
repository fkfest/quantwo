#ifndef Finput_H
#define Finput_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>
#include "utilities.h"
#include "product.h"
#include "term.h"
#include "globals.h"
#include "inpline.h"
#include "lexic.h"

/*!
    Input analyzer 
*/
class Finput {
public:
  // constructor
  Finput ( bool eq = false );
  // constructor + init input-parameters
  Finput( std::string paramspath );
  // add string
  bool addline( const std::string& line );
  // analyze equation
  bool analyzeq();
  // get input
  std::string input() const;
  Product<Lelem> eqn() const;
  // get sum of terms
  Sum<Term,TFactor> sumterms() const;
  // analyze input
  bool analyzeit();
  // clear all arrays
  void clear() {_inlines.clear(); _ineq.clear(); _input.clear(); _eqn = Equation(); _eq = false;};
  // return input lines
  const std::vector<std::string> & inlines() const { return _inlines;};
  const std::vector<std::string> & ineq() const { return _ineq;};
  
private:
  // initialyse default input-parameters 
  void InitInpars(std::string paramspath);
  // analyse all \newop's and save as a map in _eqn.newops
  void analyzenewops();
  // analyze command coming after backslash
  lui analyzecommand(std::string const & str, lui ipos);
  // variables
  std::string _input;
  bool _eq;
  Equation _eqn;
  std::vector<std::string> _inlines;
  std::vector<std::string> _ineq;
};

std::ostream & operator << (std::ostream & o, Finput const & inp);

#endif

