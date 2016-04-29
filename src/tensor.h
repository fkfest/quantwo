#ifndef Tensor_H
#define Tensor_H

#include <string>
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <stdint.h>
#include <limits>
#include "globals.h"
#include "utilities.h"
#include "types.h"
#include "arrays.h"
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
  SlotType(const std::string& lettertype);
  bool operator < ( const SlotType& st) const; 
  Length length() const { return _nIndices; };
  SlotType::Type type() const { return _type; };
  std::string name(const std::string & oldname = "") const;
//private:
  Length _nIndices;
  Type _type;
  std::string _internalName;
};

//! output operator for slot types
std::ostream & operator << (std::ostream& o, const SlotType& st);

typedef std::set<SlotType> SlotTypes;
typedef Array<uint> Slots;
typedef std::set<uint> SlotUniqueSet;


std::ostream & operator << (std::ostream& o, const Slots& ss);
std::ostream & operator << (std::ostream& o, const SlotUniqueSet& ss);

class Symmetry {
public:
  void canonicalize();
  bool operator < ( const Symmetry& sym ) const;
  bool operator == ( const Symmetry& sym ) const;
  bool operator != ( const Symmetry& sym ) const { return !(*this == sym); };
  int _sign;
  // symmetric slots
  Slots _symSlots;
  // slots that has to be permuted simultaniously to the symSlots
  Slots _simSlots;
};

//! output operator for symmetry
std::ostream & operator << (std::ostream& o, const Symmetry& sym);

class Cut {
public:
  enum Type {
    NoType = 0,
    // cut according to a set of other indices (pairdomains etc)
    Domain = 1,
    // cut according to a list of indices (pairlist, tripleslist)
    List = 2,
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
  std::string strengthLetter() const;
  uint nSlotCutType() const { return _nSlotCutType; };
  const Slots& cutSlots() const { return _cutSlots; };
  const Slots& defSlots() const { return _defSlots; };
  bool operator < ( const Cut& cut ) const;
  bool operator == ( const Cut& cut ) const;
  bool operator != ( const Cut& cut ) const { return !(*this == cut); };
  void canonicalize();
//private:
  Cut::Type _cutType;
  Cut::Strength 
      _cutStrength,
      // strength of the parent defining cut
      _PDCStrength;
  // type of the cut ( 1 = single; 2 = pair; 3 = triple ... )
  // eq. or larger (for united domains/lists) than number of cut (defining) slots in cutSlots and defSlots.
  uint _nSlotCutType;
  // main slots
  Slots _cutSlots;
  // cut-defining slots. can be empty
  Slots _defSlots;
};

//! output operator for cuts
std::ostream & operator << (std::ostream& o, const Cut& cut);

class Action;
//slot types
typedef Array<const SlotType*> SlotTs;

// sort in decreasing order (e.g. [Fai])
void Canonicalize( SlotTs& sts, Slots& ref );

typedef Array<Symmetry> Symmetries;
typedef Array<Cut> Cuts;
typedef Array<const Action*> Actions;

class TensorBase {
public:
  TensorBase(std::string name = "T") : _name(name){};
  TensorBase( const Symmetries& syms, std::string name = "T" ) 
    : _syms(syms), _name(name) {};
  const Symmetries& syms() const { return _syms; };
  const SlotUniqueSet& phantom() const { return _phantomSlots; };
  const std::string& name() const { return _name; };
  bool phantomSlot( uint iSlot ) const { return _phantomSlots.count(iSlot); };
//  virtual bool operator < ( const TensorBase& ten ) const = 0;
  
//private:
  Symmetries _syms;
  SlotUniqueSet _phantomSlots;
  std::string _name;
};

const uint MAXNINDICES = 32;

typedef double Factor;
typedef double Cost;
const Cost MAXCOST = std::numeric_limits<double>::max()/10.0; 

// diagrammatic connections
// slotref 
struct Connections {
  // positions in bitmask correspond to a global reference to a connection line in this diagram 
  // (i.e., the orbital-index name in the term)
  std::bitset<MAXNINDICES> bitmask;
  // slotref[(bitmask>>ipos).count()] == iSlot in the tensor
  // usually is simply nSlots-(bitmask>>ipos).count()-1, but we will keep it for the moment...
  Array<unsigned short> slotref;
};

class DiagramTensor : public TensorBase {
public:
  DiagramTensor( std::string name = "" ) : TensorBase(name) {};
  DiagramTensor( const Connections& conns, std::string name = "" ) : TensorBase(name), _connect(conns) {};
//  DiagramTensor( const SlotTs& slots, const Symmetries& syms, const Cuts& cuts, std::string name = "T" ) 
//    : TensorBase(slots,syms,cuts,name) {};
//  bool operator < ( const DiagramTensor& ten ) const;
  std::string slotTypeLetters( const SlotTs& slottypes ) const;
  const Connections& connects() const { return _connect; };
//private:
  Connections _connect;
};

//! output operator for DiagramTensor
std::ostream & operator << (std::ostream& o, const DiagramTensor& t);

class Tensor : public TensorBase {
public:
  Tensor( const SlotTs& slots, std::string name = "T" ) : TensorBase(name), _slots(slots), _dummy(false) {};
  Tensor( const SlotTs& slots, const Symmetries& syms, const Cuts& cuts, std::string name = "T" ) 
    : TensorBase(syms,name), _slots(slots), _cuts(cuts), _dummy(false) {};
  const SlotTs& slots() const { return _slots; };
  const Actions& parents() const { return _parents; };
  const Cuts& cuts() const { return _cuts; };
  // add a parent action (if it's not there already)
  void add( const Action * pAct );
  std::string slotTypeLetters() const;
  bool equal( const Tensor& ten, bool considerprops = true ) const;
  bool operator < ( const Tensor& ten ) const;
  bool operator == ( const Tensor& ten ) const { return equal(ten); };
  /// Desc: A comma-separated string of cut specifications for local tensors:
  ///      "012/456": slots 012 depend on slots 456, e.g., triples-domains
  ///      "456": cut according to a list of orbitals, e.g., triples-list
  ///      "*": the orbital is already summed up, e.g. "4**": slot 4 comes from the triples-list
  ///      "t23": forced triangular order, i.e. all blocks not in the triangular order will be deleted
  ///             (useful for "dummy"-defining slots in intermediates)
  ///      "f23": phantom slots needed to define local approximations (block size eq 1)
  /// multiple declarations can be joined by ",". A total declaration
  /// might look like "12/34*,34*" for B[Faaii]: B^ij_abP = T^ijk_abc(kc|P), for example.
  /// one can use {} to define slots with more than one digit, e.g., 0123/45{10}{11}
  /// strength of cuts can be given by adding a character to the begin of the corresponding set of slots
  /// (otherwise the default strength is used):
  /// s/c/w/d: strong/s+close/s+c+weak/s+c+w+distant:
  ///      "d23": pair list for strong+close+weak+distant pairs
  ///      "0/1*,s1*": united domain for strong pairs
  ///      "0/s2*,1/d2*,s2*,d2*": united domains for slot 0 from strong pairs and for slot 1 from distant pairs.
  void CreateCutFromDesc( std::string const &desc );
  
//private:
  SlotTs _slots;
  Actions _parents;
  Cuts _cuts;
  bool _dummy;
};

//! output operator for tensors
std::ostream & operator << (std::ostream& o, const Tensor& t);

#endif
