#include "arrays.h"
template <class T>
inline
void Array<T>::identity(unsigned long int n)
{
  this->reserve(n);
  for ( uint i = 0; i < n; ++i ) this->push_back(i);
}
template <class T>
inline
bool Array<T>::is_identity() const
{
  for ( uint i = 0; i < this->size(); ++i ) 
    if ( (*this)[i] != i ) return false;
  return true;
}
template <class T>
inline
void Array<T>::push_front(const T& t)
{
  Array<T> temp;
  temp.reserve(this->size()+1);
  temp.push_back(t);
  for ( typename Array<T>::const_iterator it = this->begin(); it != this->end(); ++it )
    temp.push_back(*it);
  *this = temp;
}
template <class T>
inline
Array< T > & Array<T>::operator+= (T const & t) // append t
{
  this->push_back(t);
  return *this;
}
template <class T>
inline
Array< T >& Array<T>::operator+=(const Array<T>& p)
{
  for ( typename Array<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    this->push_back(*i);
  return *this;
}
template <class T>
inline
Array<T> Array<T>::subarray(unsigned long int beg, unsigned long int end) const
{
  unsigned long int end1=(end < this->size() ? end + 1 : this->size());
  return Array<T>(this->begin()+beg,this->begin()+end1);
}
template <class T>
inline
int Array<T>::find(const T& t, uint ipos) const
{
  long unsigned int i = ipos;
  for ( typename Array<T>::const_iterator it = this->begin()+ipos; it != this->end(); ++it, ++i )
    if ( t== *it ) return i;
  return -1;
}
template <class T>
inline
void Array<T>::resort()
{
  std::sort(this->begin(), this->end());
}
template <class T>
inline
Array<T> Array<T>::refarr( const Array<uint>& ref ) const
{
  Array<T> res;
  _foreach_cauto(Array<uint>,ir,ref){
    assert( *ir < this->size() );
    res.push_back((*this)[*ir]);
  }
  return res;
}
template <class T>
inline
void Array<T>::canonicalize()
{
  _foreach_auto(typename Array<T>,it,*this){
    it->canonicalize();
  }
  std::sort(this->begin(), this->end());
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, Array<T> const & p)
{
  for ( typename Array<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    o << *i;
  return o;
}

template <class T>
inline
BigArray< T > & BigArray<T>::operator+= (T const & t) // append t
{
  this->push_back(t);
  return *this;
}
template <class T>
inline
BigArray< T >& BigArray<T>::operator+=(const BigArray<T>& p)
{
  for ( typename BigArray<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    this->push_back(*i);
  return *this;
}
template <class T>
inline
BigArray<T> BigArray<T>::subarray(unsigned long int beg, unsigned long int end) const
{
  unsigned long int end1=(end < this->size() ? end + 1 : this->size());
  return BigArray<T>(this->begin()+beg,this->begin()+end1);
}
template <class T>
inline
int BigArray<T>::find(const T& t, uint ipos) const
{
  long unsigned int i = ipos;
  for ( typename BigArray<T>::const_iterator it = this->begin()+ipos; it != this->end(); ++it, ++i )
    if ( t== *it ) return i;
  return -1;
}
template <class T>
inline
void BigArray<T>::resort()
{
  std::sort(this->begin(), this->end());
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, BigArray<T> const & p)
{
  for ( typename BigArray<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    o << *i;
  return o;
}
