#ifndef Matrices_H
#define Matrices_H

#include <string>
#include <set>
#include <iostream>
#include <assert.h>
#include "product.h"
#include "orbital.h"
#include "globals.h"

namespace Ops {
  // enumerate operator types 
  enum Type {None, Exc, Exc0, Fock, FluctP, XPert, Deexc, Deexc0, Interm, Number};
  // generate Product<Orbital> from occupied and virtual orbitals and excitation class 
  Product<Orbital> genprodorb(short exccl,Orbital const & occ, Orbital const & virt);
};

typedef std::set<long int> TCon2;
// represents connecting line, i.e. tells to which index in which matrix it is connected
struct ConLine {
  ConLine() : imat(0),idx(0){};
  ConLine(lui imat_, lui idx_) : imat(imat_),idx(idx_){};
  lui imat;
  lui idx;
};
/*! 
    Implements class matrices (e.g. amplitudes, integrals etc)
*/
class Matrices {
  public:
  // enumerate vertices types
  enum Spinsym {Singlet, Triplet};
  Matrices();
  // construct from type and orbitals and name
  Matrices(Ops::Type t, Product<Orbital> p, std::string name="T", Spinsym matspinsym=Singlet);
  // return Type
  Ops::Type type() const;
  // return const reference to Product of orbitals
  const Product<Orbital>& orbitals() const;
  // return name
  std::string name() const;
  // return true if antisymmetrized form
  bool antisymform() const;
  // replace orbital orb1 with orb2
  void replace(Orbital orb1, Orbital orb2);
  // expand antisymmetrized matrix ( from antisymmetrized form < AB || CD > to the normal form < AB | CD > - < AB | DC > )
  // if firstpart=true : < AB | CD >, if firstpart=false : < AB | DC >
  // if return is true: expanded, if false: don't need to expand
  bool expandantisym(bool firstpart);
  // artificial ordering
  bool operator < (Matrices const & t) const;
  // equality of two Matrices (including symmetry properties)
  bool operator == (Matrices const & t) const;
  // set "kind" of matrix (_exccl, _intlines, _intvirt)
  void setkind(short exccl, short intlines, short intvirt);
  // return excitation class (has to be set previously!)
  short exccl() const {return _exccl;};
  // get the position of second orbital for the same electron (from position) 
  lui iorbel(lui ipos){ assert( ipos < _orbs.size() );
                       return (ipos%2==0?ipos+1:ipos-1); };
  // get the second orbital for the same electron (if not found return blank orbital)
  Orbital orbel(Orbital const & orb);
  // get the second orbital for the same electron (from position) 
  Orbital orbel(long int const & ipos){ return _orbs[iorbel(ipos)]; };
  // get the spin symmetry of the electron, which corresponds to the orbital on position ipos
  Spinsym spinsym(long int ipos);
  // set no spin for all orbitals
  void set_no_spin();
  // reset vertices
  void reset_vertices();
  // compare vertices
  bool vertices(long int ipos, Matrices & mat, long int ipos1, unsigned int indx);
  // set connections
  void set_connect(TCon2 connected2);
  // add connection
  void add_connect(long int con);
  // return connections
  TCon2 connected2() const;
  // set connections
  void set_conlines(Product<ConLine> conlines);
  // set connection
  void set_conline(lui iorb, lui imat, lui idx);
  // add connection
  void add_conline(lui imat, lui idx);
  // return connections
  const Product<ConLine>& conlines() const;
  // return ConLine for an orbital
  const ConLine& conline(lui iorb) const;
  
  private:
  Ops::Type _type;
  Product<Orbital> _orbs;
  std::string _name;
  Spinsym _matspinsym;
  bool _antisymform; // W is constructed in antisymmetrized form, but can be expanded later.
  // needed for comparison of terms:
  long int _indx;
  // connected to (index of matrix in term (start from 1!!!)):
  TCon2 _connected2;
  // connection lines (in the order of indices)
  Product<ConLine> _conlines;
  // following variables will be set by Term::matrixkind (according to Kallay and Surjan, JCP, 115 (2001), 2945)
  short _exccl, _intlines, _intvirt;
    
};

std::ostream & operator << (std::ostream & o, Matrices const & mat);

/*! 
    Implements class permutators (\Perm{ia,jb}(ia|jb)=(jb|ia))
*/
class Permut {
  public:
    Permut();
    // construct from product of indices
    Permut(Product<Orbital> p1,Product<Orbital> p2);
    // construct from orbitals
    Permut(Orbital o1,Orbital o2);
    // append permutator
    Permut & operator += (const Permut& p);
    // multiply permutation from RIGHT (e.g. P_this(i,j) P_p(j,k) => P_this(ij,jk); P_this(j,i) P_p(k,j) => P_this(k,i))
    Permut & operator *= (const Permut& p);
    // return orbitals "from"
    Product<Orbital> orbsfrom() const;
    // return orbitals "to"
    Product<Orbital> orbsto() const;
    // artificial ordering
    bool operator < (Permut const & p) const;
    // equality of permutators
    bool operator == (Permut const & p) const;
    bool is1() const { return _orbs.size() == 0; };
  private: 
    typedef std::map<Orbital,Orbital> TPerMap;
    // map _orbs[from] = to
    TPerMap _orbs;
};

std::ostream & operator << (std::ostream & o, Permut const & p);

#endif

