#ifndef Orbital_H
#define Orbital_H

#include <string>
#include <iostream>
#include <assert.h>
#include "utilities.h"
#include "globals.h"

/*
    Orbitals (occ, virt, general) with spin (nospin, alpha, beta, general) 
*/

class Orbital {
  public:
  // enumerate orbital types
  enum Type { Occ, Virt, GenT };
  // enumerate Spin
  enum Spin { No, Up, Down, GenS };
  Orbital ();
  // constructor from name (a-h: virt, i-o: occ, p-z: general; lower case: no Spin, upper case: general)
  Orbital (std::string name);
  // constructor from name and Type 
  Orbital (std::string name, Type type);
  // constructor from name and Spin
  Orbital (std::string name, Spin spin);
  // full constructor
  Orbital (std::string name, Type type, Spin spin);
  // return orbital
  std::string name() const;
  Type type() const;
  Spin spin() const;
  // set spin
  void setspin(Spin spin);
  // check equality
  bool operator == (Orbital const & orb) const;
  // check inequality
  bool operator != (Orbital const & orb) const;
  // check ordering relation (for sorting)
  bool operator < (Orbital const & orb) const;
  // return letter-name of orbital
  std::string letname() const;
  // check main (i.e. letter) names (e.g. i<j; ii>j, i23==i42 )
  // -1: <; 0: ==; 1: >
  int comp_letname( const Orbital& orb ) const;
  
  private:
  Type _type;
  Spin _spin;
  std::string _name;
};

std::ostream & operator << (std::ostream & o, Orbital const & orb);
#endif

