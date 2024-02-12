#ifndef Util_H
#define Util_H

#include <string>
#include <map>
#include <list>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>
#include <bitset>
#ifdef __MACH__
#include <mach-o/dyld.h>    /* _NSGetExecutablePath */
#endif

#include "globals.h"
/*!
    Implements utilities (e.g. error function)
*/

// error function
void error(std::string what, std::string where="");

// log function
void say(std::string what, std::string where="");

// warning function
#define warning(x) { std::cerr << "WARNING: " << x << std::endl; }

// path of executable
std::string exepath();

// copied from IL namespace...
// find position of a substring what on the current level of the string str (i.e. don't search inside of {})
std::size_t curlyfind(const std::string& str, const std::string& what, std::size_t ipos = 0);

//Transform container of strings to a comma separated string
template <template<typename> class T>
inline 
std::string container2csstring(const T<std::string>& arr)
{
  std::string res;
  for( auto it = arr.begin(); it != arr.end()-1; ++it ){
    res += *it + ",";
  }
  res += arr.back();
  return res;
}

// string to number transformation
// call: str2num<double>(x,"3.14",std::dec), number will be in x
template <class T>
inline
bool str2num(T& t, const std::string& s,
             std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

// number to string transformation
// call: num2str(3.14,std::dec), string will be returned
template <class T>
inline
std::string num2str(const T& t,
             std::ios_base& (*f)(std::ios_base&))
{
  std::ostringstream oss;
  oss << f << t;
  return oss.str();
}
// anything to string transformation
template <class T>
inline
std::string any2str(const T& t)
{
  std::ostringstream oss;
  oss << t;
  return oss.str();
}

//return sign of a number as a char '-' or '+'
template <class T>
inline 
char sgnchar(const T& val) {
  if ( std::signbit(val) ) return '-' ;
  else return '+';
}

#define Error(what) error(what,"file "+std::string(__FILE__)+" line "+any2str<int>(__LINE__))

template<typename string_t>
string_t DirName(string_t source)
{
    source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
    return source;
}

template<typename string_t>
string_t FileName(string_t source, bool remove_ext = false)
{
    source.erase(source.begin(), std::find(source.rbegin(), source.rend(), '/').base());
    if ( remove_ext ) {
      size_t ix = source.rfind('.');
      if ( ix != std::string::npos )
        source.erase(source.begin() + ix, source.end());
    }
    return source;
}

// implementation of InSet commands (comparing one value with up to 10 other values)
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2)
{ return item==i1 || item==i2; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3)
{ return item==i1 || item==i2 || item==i3; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4 )
{ return item==i1 || item==i2 || item==i3 || item==i4; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5 || item==i6; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6,
           const T & i7 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5 || item==i6 || item==i7; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6,
           const T & i7, const T & i8 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5 || item==i6 || item==i7 || item==i8; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6,
           const T & i7, const T & i8, const T & i9 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5 || item==i6 || item==i7 || item==i8 || item==i9; }
template <typename T>
bool InSet(const T & item, const T & i1, const T & i2, const T & i3, const T & i4, const T & i5, const T & i6,
           const T & i7, const T & i8, const T & i9, const T & i10 )
{ return item==i1 || item==i2 || item==i3 || item==i4 || item==i5 || item==i6 || item==i7 || item==i8 || item==i9 || item==i10; }
// implementation of InSet command for search in an array
template <typename ValueType, size_t arraySize>
bool InSet(const ValueType& val, const ValueType (&arr)[arraySize])
{
  for (size_t i=0; i<arraySize; i++)
    if (val==arr[i]) return true;
  return false;
  /*return std::find(&arr[0], &arr[arraySize], val)!=&arr[arraySize]; }*/
}
// implementation of InSet command for search in an map::second
template <typename ValueType>
bool InSet(const ValueType& val, const std::map<std::string,ValueType>& m)
{
  for ( typename std::map<std::string,ValueType>::const_iterator it = m.begin(); it != m.end(); ++it )
    if ( val == it->second ) return true;
  return false;
}
// implementation of InSet command for search in an list
template < typename ValueType, class Cont >
bool InSet(const ValueType& val, const Cont& con)
{
  for ( typename Cont::const_iterator it = con.begin(); it != con.end(); ++it )
    if ( val == *it ) return true;
  return false;
}

// placeholder for __restrict__
#define RESTRICT
// insertion sort, the new order will be in pSel, returns the number of swaps
template<class T, class P>
long unsigned int InsertionSort( const T *RESTRICT pIn, P *RESTRICT pSel, uint N )
{
  long unsigned int nswaps = 0;
  for ( P* p = pSel+1; p != pSel+N; ++p ){
    P   s = *p,
      * q = p;
    const T& curr = pIn[s];
    for ( ; q != pSel && curr < pIn[*(q-1)]; --q, ++nswaps ) *q = *(q-1);
    *q = s;
  }
  return nswaps;
}

// insertion pointer sort, the new order will be in pSel, returns the number of swaps
template<class T, class P>
long unsigned int InsertionPSort( const T **RESTRICT pIn, P *RESTRICT pSel, uint N )
{
  long unsigned int nswaps = 0;
  for ( P* p = pSel+1; p != pSel+N; ++p ){
    P   s = *p,
      * q = p;
    const T& curr = *pIn[s];
    for ( ; q != pSel && curr < *pIn[*(q-1)]; --q, ++nswaps ) *q = *(q-1);
    *q = s;
  }
  return nswaps;
}

// insertion pointer sort, decreasing order, the new order will be in pSel, returns the number of swaps
// TODO: unite with outer InsertionSorts
template<class T, class P>
long unsigned int InsertionPSortD( const T **RESTRICT pIn, P *RESTRICT pSel, uint N )
{
  long unsigned int nswaps = 0;
  for ( P* p = pSel+1; p != pSel+N; ++p ){
    P   s = *p,
      * q = p;
    const T& curr = *pIn[s];
    for ( ; q != pSel && *pIn[*(q-1)] < curr; --q, ++nswaps ) *q = *(q-1);
    *q = s;
  }
  return nswaps;
}

// Creates new combination (analog to std::new_permutation) of m out of n (m <= n), sorted in triangular manner
// e.g. (0 1) 2, (0 2) 1, (1 2) 0
// first and last point begin and end of the array, and k corresponds to m
template <typename Iter>
bool next_combination(Iter first, Iter k, Iter last)
{
  if ((first == last) || (first == k) || (last == k)) return false;
  Iter i = k, last1 = last, first1 = first;
  --first1;
  --last1;
  while ( --i != first1 && !(*i < *last1) ){}
  if ( i == first1 ){
    std::rotate(first,k,last);
    return false;
  }
  Iter j = k;
  for ( ; !(*i < *j); ++j) {}
  std::iter_swap(i,j);
  std::rotate(++i,++j,last);
  for ( i = k; j != last; ++j, ++i) {}
  std::rotate(k,i,last);
  return true;
}

// Creates new bitset combination m out of n (m <= n)
// e.g. 1 1 0, 1 0 1, 0 1 1
//  std::bitset<8> x;
//  x[1]=true;
//  x[2]=true;
//  do {
//    xout << x << std::endl;
//  } while ( next_combination(x) );
template< std::size_t N >
bool next_combination( std::bitset<N>& bset)
{
  uint nbits = 0;
  for ( uint i = 0; i < N; ++i ) {
    if ( bset[i] ) {
      bset[i] = false;
      if ( i != nbits ) {
        // move the bit to the right and fill rest bits
        for ( uint j = i-nbits-1; j < i; ++j ) bset[j] = true;
        return true;
      } else {
        // all previous bits were full
        ++nbits;
      }
    }
  }
  return false;
}

#endif

