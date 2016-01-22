#ifndef Inpline_H
#define Inpline_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>
#include "utilities.h"
#include "globals.h"


//! functions to analyze input line
namespace IL{
  // get key which corresponds to keyword in parameter registry
  std::string key(const std::string& line, lui& ipos, const std::string& keyword);
  // get key which corresponds to keyword in parameter registry
  std::string key(const std::string& line, const std::string& keyword);
  // generate TParArray from string of parameters (e.g. "\dg","\dagger" -> \dg and \dagger )
  TParArray parray(const std::string& str);
  // add new command from newcommand or new operator from newoperator
  lui addnewcom(const std::string& str, lui ipos, std::string what = "newcommand");
  // change default parameters (set:name=value)
  void changePars(const std::string& str, lui ipos);
  // skip all characters in str beginning from ipos, which are present in what
  lui skip(const std::string& str, const lui& ipos, const std::string& what);
  // skip all characters in str to the left from ipos, which are present in what
  // ipos and result are end()-like, i.e. one-based
  lui skipr(const std::string& str, const lui& ipos, const std::string& what);
  // end of word (may be in " )
  lui endword(const std::string& line, lui& ipos);
  // get positions of next word (begin and end)
  // use gluer if glue is true, otherwise a gluer is a separator too
  // greedy: only separators will separate words, otherwise only gluer and {} can glue symbols together
  lui nextwordpos(const std::string& str, lui& ipos, bool glue = true, bool greedy = true);
  // get position of the closing bracket (which corresponds to the bracket on ipos)
  lui closbrack(const std::string& str, lui ipos);
  // generate plain name out of latex name, i.e. skip "^_{}" (e.g., a_{15}^2 -> a152)
  std::string plainname(std::string name);
  // find position of a substring sstr on the current level of the string str (i.e. don't search inside of {})
  lui lexfind(const std::string& str, const std::string& sstr, const lui& ipos = 0);
  // divide lelnam into name, superscript and subscript. returns true if up or down is there
  bool nameupdown(std::string& name, std::string& nameup, std::string& namedown, const std::string& lelnam);
  // add nameadd to name (as superscript if superscript=true) (inside \snam{} if snam is true)
  void add2name(std::string & name, std::string const & nameadd, bool superscript = true, bool snam = true);
}

#endif
