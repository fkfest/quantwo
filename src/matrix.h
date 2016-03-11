#ifndef Matrix_H
#define Matrix_H

#include <string>
#include <set>
#include <iostream>
#include <assert.h>
#include "globals.h"
#include "types.h"
#include "product.h"
#include "orbital.h"
#include "inpline.h" // for name-handling
#include "sum.h"
#include "kronecker.h"
#include "evertices.h"
// #include "operators.h" // for Creator/Annihilator types

namespace Ops {
  // generate Product<Orbital> from occupied and virtual orbitals and excitation class 
  Product<Orbital> genprodorb(short exccl,Orbital const & occ, Orbital const & virt);
}

typedef std::set<long int> TCon2;
// represents connecting line, i.e. tells to which index in which matrix it is connected
struct ConLine {
  ConLine() : imat(0),idx(0){};
  ConLine(lui imat_, lui idx_) : imat(imat_),idx(idx_){};
  lui imat;
  lui idx;
};
class Permut;
/*! 
    Implements class matrices (e.g. amplitudes, integrals etc)
*/
class Matrix {
  public:
  // enumerate vertices types
  enum Spinsym {Singlet, Triplet};
  Matrix();
  // construct from type and orbitals and name
  Matrix(Ops::Type t, const Product<Orbital>& p, uint npairs, short lmel=0, std::string name="T", 
           Spinsym matspinsym=Singlet, bool antisymW=true);
  // construct from type and orbitals and name
  Matrix(Ops::Type t, const Product<Orbital>& pcrea, const Product<Orbital>& panni, uint npairs, 
         short lmel=0, std::string name="T", Spinsym matspinsym=Singlet, bool antisymW=true);
  // construct from kronecker
  Matrix(const Kronecker& d);
  // return Type
  Ops::Type type() const;
  // return const reference to Product of orbitals
  const Product<Orbital>& orbitals() const;
  // return name
  std::string name() const;
  // return true if antisymmetrized form
  bool antisymform() const;
  // number of conserved electrons
  uint npairs() const { return _npairs;};
  // non-conserved electrons
  short lmel() const { return _lmel;};
  Spinsym matspinsym() const { return _matspinsym;};
  // combine with mat, in dontorbs: orbital indices in mat to let out
  void combine(const Matrix& mat, const Set<uint>& dontorbs = Set<uint>() );
  // replace orbital orb1 with orb2
  Return replace(Orbital orb1, Orbital orb2, bool smart);
  // replace spin spin1 with spin2
  Return replace(Spin spin1, Spin spin2, bool smart);
  // expand antisymmetrized matrix ( from antisymmetrized form < AB || CD > to the normal form < AB | CD > - < AB | DC > )
  // if firstpart=true : < AB | CD >, if firstpart=false : < AB | DC >
  // if return is true: expanded, if false: don't need to expand
  bool expandantisym(bool firstpart);
  // returns true if this is a non-singlet density matrix
  bool nonsingldm() const;
  // artificial ordering
  bool operator < (Matrix const & t) const;
  // equality of two Matrix (including symmetry properties)
  bool operator == (Matrix const & t) const;
  // equivalence of two matrices (i.e., without orbital names)
  bool equivalent( const Matrix& mat) const; 
  // number of vertices ("electrons") in the matrix
  uint nvertices() const { return _npairs+(_orbs.size()-2*_npairs);};
  // equivalent vertices (starting from 0+offs to nvertices-1+offs) (indistinguishability of electrons...)
  Equivalents equivertices(uint offs = 0) const;
  // creators orbitals (if anni=true - annihilators orbitals) (in the same order of electrons!)
  Product<Orbital> crobs(bool anni=false) const;
  // return true if the matrix is zero
  bool is0() const;
  // set cranorder for density matrix
  void set_cran( const Product<SQOpT::Gender>& cran );
  const Product<SQOpT::Gender>& get_cran() const {return _cranorder;};
  void calc_orbtypeshash();
  // set "kind" of matrix (_exccl, _intlines, _intvirt) and _orbtypeshash
  void setkind(short exccl, short intlines, short intvirt);
  // return excitation class (has to be set previously!)
  short exccl() const {return _exccl;};
  // get the position of second orbital for the same electron (from position) 
  long iorbel(lui ipos) const;
  // get the second orbital for the same electron (if not found return blank orbital)
  Orbital orbel(Orbital const & orb) const;
  // get the second orbital for the same electron (from position) 
  Orbital orbel(long int const & ipos) const;
  // get the spin symmetry of the electron, which corresponds to the orbital on position ipos
  Spinsym spinsym(long int ipos) const;
  // returns gender of the creation/annihilation operator that corresponds to the orbital on position ipos
  // information is either taken from _cranorder or guessed 
  SQOpT::Gender genderguess(uint ipos) const;
  // set no spin for all orbitals
  void set_no_spin();
  // reset vertices
  void reset_vertices();
  // compare vertices
  bool vertices(long int ipos, Matrix & mat, long int ipos1, unsigned int indx);
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
  // permute with p
  void permute(const Permut& p);
  // returns a plain name that can be used in algo's
  std::string plainname() const;
  private:
  // generate name from type or use the given name
  void gen_name(const std::string& name);
  // sets internal variables. it's used by constructors
  void create_Matrix(Ops::Type t, uint npairs, short lmel, 
                     std::string name, Spinsym matspinsym, bool antisymW);
  Ops::Type _type;
  // orbitals in the electron-order
  Product<Orbital> _orbs;
  std::string _name;
  Spinsym _matspinsym;
  Product<SQOpT::Gender> _cranorder; // sets creator or annihilator type for every orbital (for density matrices only!)
  bool _antisymform; // W is constructed in antisymmetrized form, but can be expanded later.
  // number of orbital pairs (== electrons for number-conserving operators and otherwise == conserved particles)
  uint _npairs;
  // non-conserved electrons
  short _lmel;
  // needed for comparison of terms:
  long int _indx;
  // connected to (index of matrix in term (start from 1!!!)):
  TCon2 _connected2;
  // connection lines (in the order of indices)
  Product<ConLine> _conlines;
  // following variables will be set by Term::matrixkind (according to Kallay and Surjan, JCP, 115 (2001), 2945)
  short _exccl, _intlines, _intvirt;
  // hash of types of all orbitals
  lui _orbtypeshash;
    
};

std::ostream & operator << (std::ostream & o, Matrix const & mat);

/*! 
    Implements class permutators (\Perm{ia,jb}(ia|jb)=(jb|ia))
*/
class Permut {
  public:
    Permut();
    // construct from product of indices
    Permut(Product<Orbital> p1,Product<Orbital> p2);
    // construct from orbitals (permutation from o1 to o2)
    Permut(Orbital o1,Orbital o2);
    // append permutator
    Permut & operator += (const Permut& p);
    // multiply permutation from RIGHT (e.g. P_this(i,j) P_p(j,k) => P_this(ij,jk); P_this(j,i) P_p(k,j) => P_this(k,i))
    Permut & operator *= (const Permut& p);
    // divide permutation from LEFT (e.g. P_this(ij,ji) / P_p(jk,kj) => P_this(ijk,kij))
    Permut & operator /= (const Permut& p);
    // return orbitals "from"
    Product<Orbital> orbsfrom() const;
    // return orbitals "to"
    Product<Orbital> orbsto() const;
    // permute orbital
    Orbital permutorb(const Orbital& orb) const;
    // artificial ordering
    bool operator < (Permut const & p) const;
    // equality of permutators
    bool operator == (Permut const & p) const;
    bool is1() const { return _orbs.size() == 0; };
  private: 
    typedef std::map<Orbital,Orbital> TPerMap;
    // map _orbs[from] = to
    TPerMap _orbs;
    uint dummy;
};

std::ostream & operator << (std::ostream & o, Permut const & p);

#endif

