#include "unigraph.h"

UniGraph::UniGraph(const Term& term)
{
  pTerm = &term;
  const Product<Matrices> & mats = term.mat();
  for ( uint i = 0; i < mats.size(); ++i )
    _matsord.push_back(i);
  InsertionSort(&*mats.begin(),&*_matsord.begin(),mats.size());
 
  EquiVertices equimat;
  uint prev = _matsord.size();
  uint currvert = 0;
  _foreach_cauto(Order,im,_matsord){
    const Matrices mat = mats[*im];
    uint nextvert = currvert+mat.nvertices();
    // equivalent vertices?
    Equivalents equivs( mat.equivertices(currvert) );
    _equivs.insert(_equivs.end(),equivs.begin(),equivs.end());
    // equivalent matrices?
    if ( prev < _matsord.size() && mat.equivalent(mats[prev]) ) {
      JointVertices verts;
      for ( uint vert = currvert; vert < nextvert; ++vert) {
        verts.push_back(vert);
      }
      equimat.push_back(verts);
    } else if ( equimat.size() > 0 ){
      _equivs.push_back(equimat);
    }
    prev = *im;
    currvert = nextvert;
  }
  if ( equimat.size() > 0 ){
    _equivs.push_back(equimat);
  }
  
}
