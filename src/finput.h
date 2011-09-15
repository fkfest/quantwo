#ifndef Finput_H
#define Finput_H
#include <string>
#include <vector>
#include <iostream>
#include "utilities.h"
#include "product.h"
#include "term.h"

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

/*!
    Input analyzer 
*/
class Finput {
  public:
  // word separators
  static const std::string separator;
  // neglect word separator after the gluer
  static const std::string gluer;
  // brackets
  static const std::string brackets;
  // define possible operator names (the order is important!)
  static const std::string commands[5]; // general commands
  static const std::string beqs[2]; // beq
  static const std::string eeqs[2]; // beq
  static const std::string comments[2]; // comment
  static const std::string skipops[5]; // operators to skip
  static const std::string refs[3]; // reference functions
  static const std::string dgs[2]; // daggers
  static const std::string bexcops[2]; // bare excitation operators
  static const char hms[4]; // hamiltonian parts
  // constructor
  Finput ();
  // add string
  Finput & operator += (std::string const & line);
  // get input
  std::string input() const;
  Product<Lelem> eqn() const;
  // get sum of terms
  Sum<Term,double> sumterms() const;
  // analyze input
  bool analyzeit();
  // analyze bra and ket
  
  private:
  // get positions of next word (begin and end)
  unsigned long int nextwordpos(std::string const & str, unsigned long int & ipos);
  // skip all characters in str beginning from ipos, which are present in what
  unsigned long int skip(std::string const & str, unsigned long int const & ipos, std::string const & what);
  // analyze command coming after backslash
  unsigned long int analyzecommand(std::string const & str, unsigned long int ipos);
  // get position of the closing bracket (which corresponds to the bracket on ipos)
  unsigned long int closbrack(std::string const & str, unsigned long int ipos);
  // get position of the closing bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  unsigned long int closbrack(Product<Lelem> const & eqn, unsigned long int ipos);
  // get position of the opening bracket in vector<Lelem> (which corresponds to the bracket on ipos)
  unsigned long int openbrack(Product<Lelem> const & eqn, unsigned long int ipos);
  // add connections: beg and end are positions of opening and closing parantheses in aterm
  Product<long int> addconnections(const Product< Lelem >& aterm, long unsigned int beg, long unsigned int end);
  // extract expression (remove parentheses and sums)
  bool extractit();
  // expand all Hamiltonians (\op H) to (\op F + \op W)
  Product<Lelem> expandH(Product<Lelem> const & eqn);
  // find the end position of current element in aterm, if bk==true: bra and ket are treated as brackets
  unsigned long int elem(Product<Lelem> const & aterm, unsigned long int beg, bool bk=false);
  // find the end position of current term
  unsigned long int term(Product<Lelem> const & eqn, unsigned long int beg);
  // expand equation
  Product<Lelem> expandeqn(Product<Lelem> const & eqn, std::vector< Product<long int> > & connections);
  // expand a term
  Product<Lelem> expandterm(Product<Lelem> const & aterm, std::vector< Product<long int> > & connections);
  // expand parantheses pair
  Product<Lelem>  expandpar(Product<Lelem> const & aterm, unsigned long int beg, std::vector< Product<long int> > & connections);
  // test if eqn is completely expanded
  bool expanded(Product<Lelem> const & eqn);
  // add connections to term, and term to _sumterms
  void addterm(Term& term, bool plus, long unsigned int beg, long unsigned int end, 
               Product< long int > const & indxoperterm, unsigned long int & nterm,bool excopsonly);
  // transform Product<Lelem> to Sum<Term>, if excopsonly: do only pure excitation operators (and bra/ket)
  bool do_sumterms(bool excopsonly=false);
  // handle bra/ket
  Oper handle_braket(Lelem const & lel, Term & term);
  // handle excitation index
  Oper handle_excitation(std::string const & name, bool dg, Term & term);
  // handle factor
  double handle_factor(Lelem const & lel);
  // handle operator
  Oper handle_operator(Lelem const & lel, Term & term, bool excopsonly=false);
  // handle sum
  void handle_sum(Lelem const & lel, Term & term);
  // handle parameter
  void handle_parameters(Term& term, bool excopsonly = false);
  // add nameadd to name (as superscript)
  void add2name(std::string & name, std::string const & nameadd);
  // variables
  std::string _input;
  Product<Lelem> _eqn;
  bool _eq;
  Sum<Term, double> _sumterms;
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
  // e.g. ((-2,-3),(2,3,4)) --> element 2 and 3 are disconnected, and element 4 is connected to 2 and/or 3
  std::vector< Product<long int> > _connections;
  
};

std::ostream & operator << (std::ostream & o, Finput const & inp);

#endif

