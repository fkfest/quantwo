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
Sum<Object, Field> &  Sum<Object, Field>::operator /= (const Sum<Object, Field> & s)
{
  if ( s.empty() ){
    error("Can not divide by empty sum!");
  }
  Sum<Object, Field> rest(*this),result;
  typename Sum<Object,Field>::iterator it1;
  typename Sum<Object,Field>::const_iterator its1;
  Object o1,os;
  Field f1,fs;
  rest.clean();
  int n = 0;
  while( !rest.empty() && n < Numbers::big ){
    it1 = rest.begin();
    o1 = it1->first;
    f1 = it1->second;
    its1 = s.begin();
    o1 /= its1->first;
    f1 /= its1->second;
    rest.erase(it1);
    result[o1] += f1;
    for ( ++its1 ; its1 != s.end(); ++its1 ){
      os = its1->first;
      fs = its1->second;
      os *= o1;
      fs *= f1;
      rest += std::pair<Object,Field>(os,-fs);
    }
    rest.clean();
    ++n;
  }
  if ( n == Numbers::big ){
    error("Could not divide in 1000 iterations!");
  }
  result.clean();
  *this = result;
  return *this;
}
template <class Object, class Field>
inline
void  Sum<Object, Field>::clean()
{
  Sum<Object, Field> cleansum;
  for ( typename Sum<Object,Field>::const_iterator i=this->begin();i!=this->end(); ++i) {
    if ( _todouble(_abs(i->second)) >= Numbers::verysmall )
      cleansum[i->first] = i->second;
  }
  *this = cleansum;
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

