#ifndef Orbital_H
#define Orbital_H

#include <string>
#include <iostream>
#include <assert.h>
#include "utilities.h"
#include "globals.h"

class Spin {
public:
  // enumerate Spin
  enum Type { No, // no spin
    Up, Down, // alpha, beta
    GenS, // general spin-sum, i.e. sum of spins (alpha + beta), also can represent a general spin!
    GenD}; // general spin-difference, i.e. (alpha - beta)
  Spin (Type type = No) : _type(type), _name("") {};
  // NOTE: we ignore names for spin != GenS or GenD 
  Spin (const std::string& name, Type type = GenS) : _type(type), _name(name) { 
    cleanname();};
  
  Type type() const {return _type;}; 
  std::string name() const {return _name;}; 
  void settype(Spin::Type type){ _type = type; cleanname();};
  // we ignore names for spin != GenS or GenD 
  void cleanname(){ if (_type != GenS && _type != GenD) _name = "";};
  // check equality
  bool operator == (const Spin& spin) const;
  // check inequality
  bool operator != (const Spin& spin) const {return !(*this == spin);};
  // check ordering relation (for sorting)
  bool operator < (const Spin& spin) const;
private:
  Type _type;
  std::string _name;
  
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
  Orbital (const std::string& name);
  // constructor from name and Type 
  Orbital (const std::string& name, Type type);
  // constructor from name and Spin
  Orbital (const std::string& name, Spin spin);
  Orbital (const std::string& name, Spin::Type spint);
  // full constructor
  Orbital (const std::string& name, Type type, Spin spin);
  Orbital (const std::string& name, Type type, Spin::Type spint);
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
  // compare main (i.e. letter) names (e.g. i<j; ii>j, i23==i42 )
  // -1: <; 0: ==; 1: >
  int comp_letname( const Orbital& orb ) const;
  
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
class Electron {
public:
  Electron() : _name(""){};
  Electron(std::string name) : _name(name){};
  std::string name() const {return _name;};
  bool operator == (const Electron& el) const { return _name == el._name; };
  bool operator != (const Electron& el) const { return _name != el._name; };
  bool operator < (const Electron& el) const { return _name < el._name; };
  
private:
  std::string _name;
};
std::ostream & operator << (std::ostream & o, const Electron& el);
#endif

