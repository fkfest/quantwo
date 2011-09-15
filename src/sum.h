#ifndef Sum_H
#define Sum_H

#include <map>
#include <iostream>
#include <sstream>
#include <cmath>
#include "globals.h"


/*
    Implements a sum (M.Hanrath)
*/
template <class Object, class Field>
class Sum : public std::map<Object, Field> {
  public:
    // inplace addition of an Object
    Sum<Object, Field> &  operator += (Object const & o);
    // inplace subtraction of an Object
    Sum<Object, Field> &  operator -= (Object const & o);
    // inplace addition of a Sum of Objects
    Sum<Object, Field> &  operator += (Sum<Object,Field> const & s);
    // inplace subtraction of a Sum of Objects
    Sum<Object, Field> &  operator -= (Sum<Object,Field> const & s);
    // inplace addition of an Object-Field pair
    Sum<Object, Field> &  operator += (std::pair<Object,Field> const & p);
    // multiplication by a factor
    Sum<Object, Field> &  operator *= (Field const & f);
    // erase
  //  void erase();
    // artificial ordering
    bool operator < (Sum<Object, Field> const & s) const;
};

// output operator for Sums of Ts
template <class Object, class Field>
std::ostream & operator << (std::ostream & o, Sum<Object, Field> const & p);

#include "sum.cpp"

#endif

