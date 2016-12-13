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
#include "equation.h"

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
  // last equation
  LelString eqn() const;
  // get sum of terms
  TermSum sumterms() const;
  // set sum of terms
  void sumterms( const TermSum& ts );
  // analyze input line
  bool analyzeline();
  // try to insert tensors into the last lexic expression from the previously calculated expressions
  void insert_tensors();
  // clear all arrays
  void clear() {_inlines.clear(); _ineq.clear(); _input.clear(); _eq = false;};
  // return input lines
  const std::vector<std::string> & inlines() const { return _inlines;};
  const std::vector<std::string> & ineq() const { return _ineq;};
  
private:
  // initialyse default input-parameters 
  void InitInpars(std::string paramspath);
  // analyse all \newop's and save as a map in _eqns.back().newops
  void analyzenewops();
  // analyze command coming after backslash at position ipos-1
  lui analyzecommand(lui ipos);
  // variables
  std::string _input;
  bool _eq;
  // left-hand-sides
  std::vector<Matrix> _lhs;
  // equations
  std::vector<LEquation> _eqns;
  std::vector<std::string> _inlines;
  std::vector<std::string> _ineq;
};

std::ostream & operator << (std::ostream & o, Finput const & inp);

#endif

