#ifndef Orbital_H
#define Orbital_H

#include <string>
#include <iostream>
#include <assert.h>
#include "utilities.h"
#include "globals.h"
#include "product.h"

typedef uint Electron;

class Spin {
public:
  // enumerate Spin
  enum Type { No = 0, // no spin
    Up = 1, Down = 2, // alpha, beta
    Gen = 3 , // general spin (i.e. without any restrictions)
    GenS = 4, // general spin-sum, i.e. sum of spins (alpha + beta)
    GenD = 5}; // general spin-difference, i.e. (alpha - beta)
  // for hash and loops
  static const uint MaxType = 6;
  Spin (Type type = No) : _type(type), _el(0) {};
  // NOTE: we ignore electron labels for spin = Gen 
  Spin (const Electron& el, Type type = GenS) : _type(type), _el(el) { if (type == Gen) _el = 0;}; 
  
  Type type() const {return _type;}; 
  Electron el() const {return _el;}; 
  // returns hash for spin (type + electrons*MaxType)
  uint spinhash(bool distinguish_el = true) const {return (uint(_type) + (distinguish_el?(MaxType*uint(_el)):0));};
  void settype(Spin::Type type){ _type = type; };
  void setel(Electron el) { _el = el; };
  // check equality
  bool operator == (const Spin& spin) const;
  // check inequality
  bool operator != (const Spin& spin) const {return !(*this == spin);};
  // check ordering relation (for sorting)
  bool operator < (const Spin& spin) const;
  // replace s1 with s2. 
  // correct handling of Ups and Downs. If the current value is Up or Down and the new one is general - return Return::Repeat
  // By GenS and GenD intesection --> Return::Delete
  Return replace( const Spin& s1, const Spin& s2); 
private:
  Type _type;
  Electron _el;
  
};

std::ostream & operator << (std::ostream & o, const Spin& spin);

/*
    Orbitals (occ, virt, general) with spin (nospin, alpha, beta, general) 
*/

class Orbital {
  public:
  // enumerate orbital types. The numbers are important, since we iterate over them
  enum Type { 
    Occ = 1,  // occupied (i)
    Virt = 2, // virtual (a)
    GenT = 3, // general (p)
    Act = 4  // active (v)
  };
  // for hash and loops
  static const uint MaxType = 5;
  // enumerate Spin
//  enum Spin { No, Up, Down, GenS };
  Orbital ();
  // constructor from name (a-h: virt, i-o: occ, p-z: general; lower case: no Spin, upper case: general)
  Orbital (const std::string& name, Electron el = 0);
  // constructor from name and Type 
  Orbital (const std::string& name, Type type, Electron el = 0);
  // constructor from name and Spin
  Orbital (const std::string& name, Spin spin);
  Orbital (const std::string& name, Spin::Type spint, Electron el = 0);
  // full constructor
  Orbital (const std::string& name, Type type, Spin spin);
  Orbital (const std::string& name, Type type, Spin::Type spint, Electron el = 0);
  // return orbital
  std::string name() const;
  Type type() const;
  Spin spin() const;
  // set spin
  void setspin(Spin spin);
  // set electron (in _spin)
  void setel(Electron el){_spin.setel(el);};
  // check equality
  bool operator == (Orbital const & orb) const;
  // check inequality
  bool operator != (Orbital const & orb) const;
  // check ordering relation (for sorting)
  bool operator < (Orbital const & orb) const;
  // return letter-name of orbital
  std::string letname() const;
  // compare main (i.e. letter) names (e.g. i<j; ii>j, i23==i42 )
  // -1: <; 0: ==; 1: >
  int comp_letname( const Orbital& orb ) const;
  // replace orb1 with orb2
  Return replace( const Orbital& orb1, const Orbital& orb2, bool smart);
  // replace spin1 with spin2
  Return replace( const Spin& spin1, const Spin& spin2, bool smart);
  
  private:
  // generate orbital type from _name
  void gentype();
  Type _type;
  Spin _spin;
  std::string _name;
};

std::ostream & operator << (std::ostream & o, Orbital const & orb);

/*
    Electrons
*/
class Electrons {
public:
  Electrons() : _name(""){};
  Electrons(std::string name) : _name(name){};
  std::string name() const {return _name;};
  bool operator == (const Electrons& el) const { return _name == el._name; };
  bool operator != (const Electrons& el) const { return _name != el._name; };
  bool operator < (const Electrons& el) const { return _name < el._name; };
  
private:
  std::string _name;
};
std::ostream & operator << (std::ostream & o, const Electrons& el);

typedef Set<Orbital> TOrbSet;
typedef Set<Electrons> TElSet;

#endif

