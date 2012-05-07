#include "sum.h"
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator += (Object const & o)
{
  (*this)[o] += 1;
  return *this;
}
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator -= (Object const & o)
{
  (*this)[o] -= 1;
  return *this;
}
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator += (Sum<Object, Field> const & s)
{
  for ( typename Sum<Object,Field>::const_iterator i=s.begin(); i!=s.end(); ++i )
    (*this)[i->first] += i->second;
  return *this;
}
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator -= (Sum<Object, Field> const & s)
{
  for ( typename Sum<Object,Field>::const_iterator i=s.begin(); i!=s.end(); ++i )
    (*this)[i->first] -= i->second;
  return *this;
}
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator += (std::pair<Object, Field> const & p)
{
  (*this)[p.first] += p.second;
  return *this;
}
template <class Object, class Field>
inline
Sum<Object, Field> &  Sum<Object, Field>::operator *= (Field const & f)
{
  for ( typename Sum<Object,Field>::const_iterator i=this->begin();i!=this->end(); ++i) {
    (*this)[i->first] *= f;
  }
  return *this;
}
template <class Object, class Field>
inline
bool Sum<Object, Field>::operator<(const Sum< Object, Field >& s) const
{
  if (this->size() < s.size())return true;
  if (s.size() < this->size())return false;
  typename Sum<Object,Field>::const_iterator j=s.begin();
  for ( typename Sum<Object,Field>::const_iterator i=this->begin();i!=this->end(); ++i) {
    if (i->first < j->first) return true;
    if (j->first < i->first) return false;
    if (i->second < j->second) return true;
    if (j->second < i->second) return false;
    ++j;
  }
  return false; // the same
}
template <class Object, class Field>
std::ostream & operator << (std::ostream & o, Sum<Object,Field> const & p)
{
  std::streampos ipos0;
  typename Sum<Object,Field>::const_iterator last=p.end(); 
  if (p.begin()!=last)last--;
  for ( typename Sum<Object,Field>::const_iterator i=p.begin();i!=p.end(); ++i) {
    ipos0=o.tellp();
    if ( i->second<0 ) {
      if ( i!=p.begin() )
        o << " - ";
      else
        o << "-";
    } else {
      if ( i!=p.begin() )
        o << " + ";
    }
    MyOut::pcurout->lenbuf += o.tellp()-ipos0;
    ipos0=o.tellp();
    if ( _todouble(_abs( _abs(i->second) - 1 )) > MyOut::pcurout->small ){
      o << _abs(i->second) << "*";
      MyOut::pcurout->lenbuf += 1; //times sign
      if ( typeid(Field) == typeid(double) )
        MyOut::pcurout->lenbuf += o.tellp()-ipos0-1;
    }
    o << i->first ; // lenline should be handled in other routines
    if ( ! MyOut::pcurout->breaklongline() && i!=last ) 
      MyOut::pcurout->flushbuf();
  }
  return o;
}

