#ifndef Types_H
#define Types_H
#include <vector>
#include "product.h"
// various types

namespace SQOpT{
  // enumerate operator types (Gen: for Particle/Hole formalism and general orbitals)
  enum Gender { Creator, Annihilator, Gen };
}

namespace Ops {
  // enumerate operator types 
  enum Type 
  { None, 
    Exc, // excitation operators \op T_i
    Exc0, // bare excitation operators \op \tau_{\mu_i}
    Fock, // Fock 
    OneEl, // one-electron operator \op h
    FluctP, // fluctuation potential
    XPert, // external perturbation
    Deexc, // deexcitation operators \op T_i^\dg
    Deexc0, // bare deexcitation operators \op \tau_{\mu_i}^\dg
    Interm, // some intermediates
    DensM, // density matrix (for active orbitals)
    Delta, // kronecker
    Number};
}

// connections
typedef std::vector< Product<long int> > ConnectionsMap;

struct Order : public std::vector< unsigned int > {
  typedef std::vector< unsigned int > Base;
  using Base::Base;
};

#endif
