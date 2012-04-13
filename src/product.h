#ifndef Product_H
#define Product_H
#include <vector>
#include <algorithm>
#include <iostream>
#include "utilities.h"


/*
    Implements a non-commutative product of Ts  (M. Hanrath)
*/
template <class T>
class Product : public std::vector<T> {
  public:
    Product<T> () : std::vector<T>(){};
    Product<T> ( typename Product<T>::const_iterator beg, typename Product<T>::const_iterator end)
     : std::vector<T>(beg,end){};
    // append t to product
    Product<T> & operator *= (T const & t);
    // append product to product
    Product<T> & operator *= (Product<T> const & p);
    // get sub product
    Product<T> subprod(unsigned long int beg, unsigned long int end) const;
    // search, if not found -> -1
    int find(T const & t) const;
    // bubble sort in increasing order
    void resort();
};
template <class T>
std::ostream & operator << (std::ostream & o, Product<T> const & p);

#include "product.cpp"

#endif

