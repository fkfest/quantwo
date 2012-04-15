#ifndef Util_H
#define Util_H

#include <string>
#include <map>
#include <list>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#ifdef TARGET_OS_MAC
#include <mach-o/dyld.h>	/* _NSGetExecutablePath */
#endif

/*! 
    Implements utilities (e.g. error function)
*/

// error function
void error(std::string what, std::string where="");

// log function
void say(std::string what, std::string where="");

// path of executable
std::string exepath();

// string to number transformation
// call: str2num<double>(x,"3.14",std::dec), number will be in x
template <class T>
bool str2num(T& t, const std::string& s,std::ios_base& (*f)(std::ios_base&));

// number to string transformation
// call: num2str(3.14,std::dec), string will be returned
template <class T>
std::string num2str(const T& t,std::ios_base& (*f)(std::ios_base&));

template<typename string_t>
string_t DirName(string_t source)
{
    source.erase(std::find(source.rbegin(), source.rend(), '/').base(), source.end());
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

#include "utilities.cpp"

#endif

