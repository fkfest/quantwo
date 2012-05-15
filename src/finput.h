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

/*!
    Lexic elements
*/
class Lelem {
  public:
  // enumerate lexic
  enum Lex {Bra, Ket, LPar, RPar, Oper, Par, Num, Frac, Plus, Minus, Times, Div, Sum };
  // enumerate types expressions in parantheses (Normal, Connected, Disconnected, ...)
  enum Conn {Normal, Connect, Disconnect }; 
  // constructor from name and Lex
  Lelem (std::string name, Lex lex, Conn conn=Normal);
  std::string name() const;
  Lex lex() const;
  Conn conn() const;
  // return if the bra was already expanded
  bool expandedbra() const;
  // return same element with a status "bra is expanded"
  Lelem braexpanded() const;
  // check equality
  bool operator == (Lelem const & lel) const;
  private:
  std::string _name;
  Lex _lex;
  Conn _conn;
  bool _expandedbra;
};

std::ostream & operator << (std::ostream & o, Lelem const & lel);

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
};
/*!
    Equation
 */
class Equation {
public:
  // add string
  Equation & operator += (const Lelem & lel) { _eqn *= lel; return *this; };
  // clear _eqn
  void reseteq() { _eqn = Product<Lelem>();};
  // extract expression (remove parentheses and sums)
  bool extractit();
  // transform Product<Lelem> to Sum<Term>, if excopsonly: do only pure excitation operators (and bra/ket)
  bool do_sumterms(bool excopsonly=false);
  // add new operator
  void addnewop(const std::string& name, const Product<Lelem>& oper){ _newops[name] = oper; };
  Product<Lelem> eqn() const { return _eqn;};
  // get sum of terms
  Sum<Term,TFactor> sumterms() const { return _sumterms;};
private:
  // get position of the closing bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  lui closbrack(Product<Lelem> const & eqn, lui ipos);
  // get position of the opening bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  lui openbrack(Product<Lelem> const & eqn, lui ipos);
  // add connections: beg and end are positions of opening and closing parantheses in aterm
  Product<long int> addconnections(const Product< Lelem >& aterm, lui beg, lui end);
  // expand all newops
  Product<Lelem> expandnewops(Product<Lelem> const & eqn);
  // find the end position of current element in aterm, if bk==true: bra and ket are treated as brackets
  lui elem(Product<Lelem> const & aterm, lui beg, bool bk=false);
  // find the end position of current term
  lui term(Product<Lelem> const & eqn, lui beg);
  // expand equation
  Product<Lelem> expandeqn(Product<Lelem> const & eqn, std::vector< Product<long int> > & connections);
  // expand a term
  Product<Lelem> expandterm(Product<Lelem> const & aterm, std::vector< Product<long int> > & connections);
  // expand parantheses pair
  Product<Lelem>  expandpar(Product<Lelem> const & aterm, lui beg, std::vector< Product<long int> > & connections);
  // test if eqn is completely expanded
  bool expanded(Product<Lelem> const & eqn);
  // add connections to term, and term to _sumterms
  void addterm(Term& term, bool plus, lui beg, lui end, 
               Product< long int > const & indxoperterm, lui & nterm,bool excopsonly);
  // handle bra/ket
  Oper handle_braket(Lelem const & lel, Term & term);
  // handle explicit excitation index (like ^{ab}_{ij})
  Oper handle_explexcitation(std::string const & name, bool dg, Term & term);
  // handle excitation index
  Oper handle_excitation(std::string const & name, bool dg, Term & term);
  // handle factor
  TFactor handle_factor(Lelem const & lel);
  // handle operator
  Oper handle_operator(Lelem const & lel, Term & term, bool excopsonly=false);
  // handle sum
  void handle_sum(Lelem const & lel, Term & term);
  // handle parameter
  void handle_parameters(Term& term, bool excopsonly = false);
  // add nameadd to name (as superscript)
  void add2name(std::string & name, std::string const & nameadd);
  Product<Lelem> _eqn;
  Sum<Term, TFactor> _sumterms;
  // save names of pure excitation and deexciation operators
  Product<std::string> _excops;
  // save indices which the excitation (and deexcitation) operators have got.
  Product<Orbital> _occexcops, _virexcops;
  // save excclass
  Product<short> _exccls;
  // save spinsymmetry
  Product<Matrices::Spinsym> _spinsymexcs;
  // save position of a Matrix of the _excops in term, if -1: the _excops is not present in this term
  Product<int> _posexcopsterm;
  // save parameters in term
  Product<Lelem> _paramterm;
  // connections "map" in _eqn (starts from 1)
  // positive: connected; negative: disconnected
  // the earlier the connection comes the more important it is:
  // e.g. ((-2,-3),(2,3,4)) --> elements 2 and 3 are disconnected, and element 4 is connected to 2 and/or 3
  std::vector< Product<long int> > _connections;
  //custom operators from \newop
  std::map< std::string, Product<Lelem> > _newops;
};

/*!
    Input analyzer 
*/
class Finput {
public:
  // constructor
  Finput ();
  // constructor + init input-parameters
  Finput( std::string paramspath );
  // add string
  bool addline( const std::string& line );
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

