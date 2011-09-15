#ifndef Operators_H
#define Operators_H

#include <iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include "utilities.h"
#include "globals.h"
#include "orbital.h"
#include "matrices.h"
#include "product.h"


/*!
    Implements a second quantized creation (\op a_i^\dagger)
    and annihilation ( \op a_i) operators (M. Hanrath)
*/
class SQOp {
  public:
  // enumerate operator types (Gen: for Particle/Hole formalism and general orbitals)
  enum Gender { Creator, Annihilator, Gen };
/*  // constructor from character (upper case: creator, lower case> annihilator)
  SQOp (std::string orb); */
  // constructor from gender and orbital
  SQOp (Gender gender, Orbital orb);
  // return gender
  Gender gender() const;
  // return gender in Particle/Hole formalism
  Gender genderPH() const;
  // return orbital
  Orbital orb() const;
  // check equality
  bool operator == (SQOp const & o) const;
  // check ordering relation (for sorting)
  bool operator < (SQOp const & o) const;
  // replace orbital orb1 with orb2
  void replace(Orbital orb1, Orbital orb2);

  private:
  Gender _gender;
  Orbital _orb;
};

std::ostream & operator << (std::ostream & o, SQOp const & op);

/*! 
    Implements Operators in second quantized form
*/
class Oper {
  public:
  //default constructor
  Oper ();
  // constructor from Type (for Hamiltonian)
  Oper (Ops::Type type);
  // constructor from excitation class and Type
  Oper (Ops::Type type, short exccl, std::string name="T");
  Oper (Ops::Type type, short exccl, void * term, Orbital (*freeorb)(void * , Orbital::Type), std::string name="T");
  // constructor from excitation class, Type and orbitals
  Oper (Ops::Type type, short exccl,Orbital occ, Orbital virt ,std::string name="T");
  //return matrix (integral or amplitude)
  Matrices mat() const;
  //return operator
  Product<SQOp> SQprod() const;
  //return prefactor
  double prefac() const;
  // return summation indices
  Product<Orbital> sumindx() const;
  // return real summation indices (without summations over bare excitations)
  Product<Orbital> realsumindx() const;

  private:
  void create_Oper(const short int& exccl, const Orbital& occ, const Orbital& virt, const std::string& name);
  Ops::Type _type;
  Product<SQOp> _SQprod;
  Matrices _mat;
  double _prefac;
  Product <Orbital> _sumindx, _fakesumindx;
};

std::ostream & operator << (std::ostream & o, Oper const & op);
  
#endif

