#ifndef Types_H
#define Types_H
#include <vector>
#include "product.h"
#include "arrays.h"
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
    Deexc, // deexcitation operators \op T_i^\dg
    Deexc0, // bare deexcitation operators \op \tau_{\mu_i}^\dg
    Fock, // Fock
    OneEl, // one-electron operator \op h
    FluctP, // fluctuation potential
    XPert, // external perturbation
    Interm, // some intermediates
    DensM, // density matrix (for active orbitals)
    Delta, // kronecker
    Overlap,
    Number};
}

// connections
typedef std::vector< Product<long int> > ConnectionsMap;

typedef Array<unsigned int> Order;
#endif
