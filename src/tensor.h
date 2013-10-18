#ifndef Tensor_H
#define Tensor_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <stdint.h>
#include "globals.h"
#include "utilities.h"
#include "types.h"
//#include "product.h"
//#include "orbital.h"
//#include "inpline.h"
//#include "sum.h"

typedef uint64_t Length;

class SlotType {
public:
  enum Type {
    NoType = 0,
    Occ = 1,
    Virt = 2,
    GenT = 3,
    Act = 4,
    AO = 5,
    DF = 6,
    RI = 7
  };
  SlotType(Length n = 0, Type type = NoType) : _nIndices(n), _type(type) {};
  bool operator < ( const SlotType& st) const; 
  Length length() const { return _nIndices; };
  SlotType::Type type() const { return _type; };
  std::string typeLetter() const;
//private:
  Length _nIndices;
  Type _type;
};

//! output operator for slot types
std::ostream & operator << (std::ostream& o, const SlotType& st);

typedef std::set<SlotType> SlotTypes;
typedef std::vector<uint> Slots;

std::ostream & operator << (std::ostream& o, const Slots& ss);


class Cut {
public:
  enum Type {
    NoType = 0,
    // cut according to a list of indices (pairlist, tripleslist)
    List = 1,
    // cut according to a set of other indices (pairdomains etc)
    Domain = 2,
    // forced triangular mode (neglect all non-triangular blocks)
    Triang = 10
  };
  // cut strength
  enum Strength{
    DefSth = 0,
    // distant (e.g. strong+close+weak+distant pairs)
    Dist = 1,
    // weak (e.g. strong+close+weak pairs)
    Weak = 2,
    // close (e.g. strong+close pairs)
    Close = 3,
    // strong (e.g. strong pairs only)
    Strong = 4,
    // no strength, needed for search
    NoSth = 100
  };
  Cut::Type type() const { return _cutType; };
  Cut::Strength cutStrength() const { return _cutStrength; };
  std::string StrengthLetter() const;
  uint nSlotCutType() const { return _nSlotCutType; };
  const Slots& cutSlots() const { return _cutSlots; };
  const Slots& defSlots() const { return _defSlots; };
  bool operator < ( const Cut& cut ) const;
//private:
  Cut::Type _cutType;
  Cut::Strength 
      _cutStrength,
      // strength of the parent defining cut
      _PDCStrength;
  // type of the cut ( 1 = single; 2 = pair; 3 = triple ... )
  // eq. or larger (for united domains/lists) than number of cut (defining) slots in cutSlots and defSlots.
  uint _nSlotCutType;
  Slots _cutSlots;
  Slots _defSlots;
};

//! output operator for cuts
std::ostream & operator << (std::ostream& o, const Cut& cut);

class Action;
//slot types
typedef std::vector<const SlotType*> SlotTs;
typedef std::vector<Cut> Cuts;
typedef std::vector<const Action*> Actions;

class Tensor {
public:
  Tensor( const SlotTs& slots, std::string name = "T" ) : _slots(slots), _name(name) {};
  Tensor( const SlotTs& slots, const Cuts& cuts, std::string name = "T" ) : _slots(slots), _cuts(cuts), _name(name) {};
  const SlotTs& slots() const { return _slots; };
  const Cuts& cuts() const { return _cuts; };
  const Actions& parents() const { return _parents; };
  const std::string& name() const { return _name; };
  std::string SlotTypeLetters() const;
  bool operator < ( const Tensor& ten ) const;
//private:
  SlotTs _slots;
  Cuts _cuts;
  Actions _parents;
  std::string _name;
};

//! output operator for tensors
std::ostream & operator << (std::ostream& o, const Tensor& t);

typedef double Factor;
typedef double Cost;




#endif