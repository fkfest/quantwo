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
  {
    push_back(*i);
  }
  return *this;
}
template <class T>
inline
Product<T> Product<T>::subprod(unsigned long int beg, unsigned long int end) const
{
  Product<T> result;
  unsigned long int end1=(end<this->size() ? end : this->size()-1);
  for ( typename Product<T>::const_iterator i=this->begin()+beg; i<=this->begin()+end1; ++i )
  {
    result *= *i;
  }
  return result;
}
template <class T>
inline
int Product<T>::find(const T& t) const
{
  for ( unsigned int i=0; i<this->size(); ++i)
  {
    if ( t==this->at(i) ) return i;
  }
  return -1;
}
template <class T>
inline
void Product<T>::resort()
{
  for (unsigned long int i=0; i<this->size(); i++)
    for ( unsigned int j=0; j<this->size()-i-1; ++j)
      if (this->at(j+1)<this->at(j)) std::swap(this->at(j+1),this->at(j));
}

template <class T>
inline
std::ostream & operator << (std::ostream & o, Product<T> const & p)
{
  if (p.size()>0)
  {
    for ( typename Product<T>::const_iterator i=p.begin(); i!=p.end(); ++i )
      o << *i;
  }
  return o;
}

