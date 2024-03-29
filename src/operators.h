#ifndef Operators_H
#define Operators_H

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include "utilities.h"
#include "globals.h"
#include "types.h"
#include "orbital.h"
#include "matrix.h"
#include "product.h"
#include "term.h"
#include "assert.h"

/*!
    Implements a second quantized creation (\op a_i^\dagger)
    and annihilation ( \op a_i) operators
*/
class SQOp {
  public:
/*  // constructor from character (upper case: creator, lower case> annihilator)
  SQOp (std::string orb); */
  // constructor from gender and orbital
  SQOp (SQOpT::Gender gender, Orbital orb);
  // return gender
  SQOpT::Gender gender() const
    { return _gender;};
  // return gender in Particle/Hole formalism
  SQOpT::Gender genderPH() const
    { return _genderPH;};
  // return orbital
  Orbital orb() const;
  // check equality
  bool operator == (SQOp const & o) const
    { return _gender==o._gender && _orb==o._orb; };
  // check ordering relation (for sorting)
  bool operator < (SQOp const & o) const
    {
      if ( _gender<o._gender )
        return true;
      if (o._gender<_gender )
        return false;
      return _orb<o._orb;
    };
  // replace orbital orb1 with orb2
  Return replace(Orbital orb1, Orbital orb2, bool smart);
  // replace spin spin1 with spin2
  Return replace(Spin spin1, Spin spin2, bool smart);

  private:
  SQOpT::Gender _gender, _genderPH;
  Orbital _orb;
};

std::ostream & operator << (std::ostream & o, SQOp const & op);

class Term;
typedef Sum<Term,TFactor> TermSum;
/*!
    Implements Operators in second quantized form
*/
class Oper {
  public:
  //default constructor
  Oper ();
  // constructor from Type (for Hamiltonian)
  Oper (Ops::Type type, bool antisym = true, Term* pTerm = 0,
        const std::vector<OrbitalTypes>& orbtypes = std::vector<OrbitalTypes>());
  // constructor from excitation class and Type, lm: difference in number of electrons (#virts-#occs), exccl == #occs
  Oper (Ops::Type type, short exccl, std::string name="T", int lm=0, int pmsym=0, Term* pTerm = 0);
  // constructor from excitation class, Type and orbital-types (for occ. orbs: [0] and virts.: [1] -- for multireference).
  // lm: difference in number of electrons (#virts-#occs), exccl == #occs
  Oper (Ops::Type type, short exccl, const std::vector<OrbitalTypes>& orbtypes,
        std::string name="T", int lm=0, int pmsym=0, Term* pTerm = 0);
  // constructor from excitation class, Type and orbitals
  Oper (Ops::Type type, short exccl,Orbital occ, Orbital virt ,std::string name="T",
        int lm=0, int pmsym=0, Term* pTerm = 0);
  Oper (Ops::Type type, short exccl, const std::map<Orbital::Type,Orbital>& orbnames,
        const std::vector<OrbitalTypes>& orbtypes, std::string name="T", int lm=0, int pmsym=0, Term* pTerm = 0);
  // constructor from excitation class, Type and product of orbitals (virts.size - occs.size  == lm; occs.size = exccl)
  Oper (Ops::Type type, short exccl,const Product<Orbital>& occs, const Product<Orbital>& virts ,
        std::string name="T", int lm=0, int pmsym=0, Term* pTerm = 0);
  // constructor from orbitals (order and spins won't change)
  Oper (Ops::Type type, const Product<Orbital>& orbs, std::string name="T", int lm=0, int pmsym=0 );
  //return matrix (integral or amplitude)
  Matrix mat() const;
  //return operator
  Product<SQOp> SQprod() const;
  //return sum of operators/matrices - for the generalized normal ordering
  TermSum sumops() const { return _sumops; }
  //return prefactor
  TFactor prefac() const;
  // return orbitals
  const TOrbSet & orbs() const;
  // return summation indices
  TOrbSet sumorbs() const;

  private:
  // for hamiltonian-parts
  void create_Oper(const std::string& name, bool antisym, const std::vector<OrbitalTypes>& orbtypes);
  // for excitation operators
  void create_Oper(const short int& exccl, const Orbital& occ, const Orbital& virt,
                   const std::string& name, int lm, int pmsym);
  void create_Oper(const short int& exccl, const std::map<Orbital::Type,Orbital>& orbnames,
                   const std::vector<OrbitalTypes>& orbtypes, const std::string& name, int lm, int pmsym);
  void create_Oper(const Product< Orbital >& occs, const Product< Orbital >& virts,
                   const std::string& name, int pmsym);
  // create operator using orbitals in orbs (in the given order and spins!)
  void create_Oper(const Product<Orbital>& orbs, const std::string& name, int lm, int pmsym);
  // move _SQprod to _sumops and clear it
  void move_SQprod();
  // returns _SQprod in general normal ordering (_SQprod is cleared at the end)
  TermSum gennormord();
  typedef Product<uint> Inds;
  // term from the generalized normal ordering
  // in creas and annis: the indices (offsets) of creator and annihilator (from end) operators in _SQprod
  // to create density matrices
  TermSum gennormordterm(const Product<uint>& creas, const Product<uint>& annis);
  Ops::Type _type;
  Product<SQOp> _SQprod;
  Matrix _mat;
  TFactor _prefac;
  TOrbSet _orbs, _sumorbs;
  TermSum _sumops;
  Term * p_Term;
};

std::ostream & operator << (std::ostream & o, Oper const & op);

#endif

