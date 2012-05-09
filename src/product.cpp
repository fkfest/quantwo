#include "product.h"
template <class T>
inline
Product< T > & Product<T>::operator *= (T const & t) // append t to product
{
  push_back(t);
  return *this;
}
template <class T>
inline
Product< T >& Product<T>::operator*=(const Product<T>& p)
{
  for ( typename Product<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    push_back(*i);
  return *this;
}
template <class T>
inline
Product<T> Product<T>::subprod(unsigned long int beg, unsigned long int end) const
{
  unsigned long int end1=(end < this->size() ? end + 1 : this->size());
  return Product<T>(this->begin()+beg,this->begin()+end1);
}
template <class T>
inline
int Product<T>::find(const T& t) const
{
  long unsigned int i = 0;
  for ( typename Product<T>::const_iterator it = this->begin(); it != this->end(); ++it, ++i )
    if ( t== *it ) return i;
  return -1;
}
template <class T>
inline
void Product<T>::resort()
{
  std::sort(this->begin(), this->end());
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, Product<T> const & p)
{
  for ( typename Product<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    o << *i;
  return o;
}

template <class T>
inline
List< T > & List<T>::operator *= (T const & t) // append t to product
{
  push_back(t);
  return *this;
}
template <class T>
inline
List< T >& List<T>::operator*=(const Product<T>& p)
{
  for ( typename Product<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    push_back(*i);
  return *this;
}
template <class T>
inline
List< T >& List<T>::operator*=(const List<T>& p)
{
  for ( typename List<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    push_back(*i);
  return *this;
}
template <class T>
inline
List< T >& List<T>::operator*=(const Set<T>& p)
{
  for ( typename Set<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    push_back(*i);
  return *this;
}
template <class T>
inline
List<T> List<T>::subprod(unsigned long int beg, unsigned long int end) const
{
  unsigned long int end1=(end < this->size() ? end + 1 : this->size());
  return List<T>(this->begin()+beg,this->begin()+end1);
}
template <class T>
inline
int List<T>::find(const T& t) const
{
  long unsigned int i = 0;
  for ( typename List<T>::const_iterator it = this->begin(); it != this->end(); ++it, ++i )
    if ( t== *it ) return i;
  return -1;
}
template <class T>
inline
void List<T>::resort()
{
  std::sort(this->begin(), this->end());
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, List<T> const & p)
{
  for ( typename List<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    o << *i;
  return o;
}

template <class T>
inline
Set< T > & Set<T>::operator *= (T const & t) 
{
  insert(t);
  return *this;
}
template <class T>
inline
Set< T >& Set<T>::operator*=(const Product<T>& p)
{
  insert(p.begin(),p.end());
  return *this;
}
template <class T>
inline
Set< T >& Set<T>::operator*=(const Set<T>& p)
{
  insert(p.begin(),p.end());
  return *this;
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, Set<T> const & p)
{
  for ( typename Set<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
    o << *i;
  return o;
}
