#ifndef Array_H
#define Array_H
#include <vector>
#include <deque>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "utilities.h"


/*
some useful functions for vector arrays
*/
template <class T>
class Array : public std::vector<T> {
    typedef std::vector<T> Base;
  public:
    using Base::Base;
    void push_front(const T& t);
    // set from 0 to n
    void identity( unsigned long int n );
    bool is_identity() const;
    // append t to array
    Array<T> & operator += (T const & t);
    // append array to array
    Array<T> & operator += (Array<T> const & p);
    // get sub array
    Array<T> subarray(unsigned long int beg, unsigned long int end) const;
    // search (starting from position ipos), if not found -> -1
    int find(T const & t, uint ipos = 0) const;
    // quick sort in increasing order
    void resort();
    // referenced array ((*this)[ref[i]])
    Array<T> refarr( const Array<uint>& ref ) const;
    // canonicalize each element and sort in increasing order
    void canonicalize();
};
template <class T>
std::ostream & operator << (std::ostream & o, Array<T> const & p);

/*
some useful functions for deque arrays (useful for very large arrays and many insertions)
*/
template <class T>
class BigArray : public std::deque<T> {
    typedef std::deque<T> Base;
  public:
    using Base::Base;
    // append t to array
    BigArray<T> & operator += (T const & t);
    // append array to array
    BigArray<T> & operator += (BigArray<T> const & p);
    // get sub array
    BigArray<T> subarray(unsigned long int beg, unsigned long int end) const;
    // search (starting from position ipos), if not found -> -1
    int find(T const & t, uint ipos = 0) const;
    // quick sort in increasing order
    void resort();
};
template <class T>
std::ostream & operator << (std::ostream & o, BigArray<T> const & p);
#include "arrays.cpp"

#endif

