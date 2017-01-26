#ifndef Product_H
#define Product_H
#include <vector>
#include <list>
#include <set>
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
    Product<T> ( const T& a1 ) : std::vector<T>(){ 
      this->push_back(a1);
    }
    Product<T> ( const T& a1, const T& a2 ) : std::vector<T>(){ 
      this->push_back(a1); 
      this->push_back(a2);
    }
    Product<T> ( const T& a1, const T& a2, const T& a3 ) : std::vector<T>(){ 
      this->push_back(a1); 
      this->push_back(a2);
      this->push_back(a3);
    }
    // append t to product
    Product<T> & operator *= (T const & t);
    // append product to product
    Product<T> & operator *= (Product<T> const & p);
    // get sub product
    Product<T> subprod(unsigned long int beg, unsigned long int end) const;
    // search (starting from position ipos), if not found -> -1
    int find(T const & t, uint ipos = 0) const;
    // quick sort in increasing order
    void resort();
    // erase element at position ipos
    void eraseelem( uint ipos ) {this->erase(this->begin()+ipos);}
};
template <class T>
std::ostream & operator << (std::ostream & o, Product<T> const & p);

template <class T>
class Set;
/*
    Implements a non-commutative list of Ts  
*/
template <class T>
class List : public std::list<T> {
  public:
    List<T> () : std::list<T>(){};
    List<T> ( typename List<T>::const_iterator beg, typename List<T>::const_iterator end)
     : std::list<T>(beg,end){};
    List<T> ( const Product<T>& p)
     : std::list<T>(p.begin(),p.end()){};
    List<T> ( const Set<T>& p)
     : std::list<T>(p.begin(),p.end()){};
    // append t to list
    List<T> & operator *= (T const & t);
    // append product to list
    List<T> & operator *= (Product<T> const & p);
    // append list to list
    List<T> & operator *= (List<T> const & p);
    List<T> & operator *= (Set<T> const & p);
    // get sub product
    List<T> subprod(unsigned long int beg, unsigned long int end) const;
    // search, if not found -> -1
    int find(T const & t) const;
    // quick sort in increasing order
    void resort();
};
template <class T>
std::ostream & operator << (std::ostream & o, List<T> const & p);

/*
    Implements a set of Ts (commutative and unique)
*/
template <class T>
class Set : public std::set<T> {
  public:
    Set<T> () : std::set<T>(){};
    Set<T> ( typename Set<T>::const_iterator beg, typename Set<T>::const_iterator end)
     : std::set<T>(beg,end){};
    Set<T> ( const Product<T>& p)
     : std::set<T>(p.begin(),p.end()){};
    // add t to set
    Set<T> & operator *= (T const & t);
    // add product to set
    Set<T> & operator *= (Product<T> const & p);
    // add set to set
    Set<T> & operator *= (Set<T> const & p);
};
template <class T>
std::ostream & operator << (std::ostream & o, Set<T> const & p);
#include "product.cpp"

#endif

