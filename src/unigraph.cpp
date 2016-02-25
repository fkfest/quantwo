#include "unigraph.h"

UniGraph::UniGraph(const Term& term)
{
  pTerm = &term;
  const Product<Matrix> & mats = term.mat();
  for ( uint i = 0; i < mats.size(); ++i )
    _matsord.push_back(i);
  InsertionSort(&*mats.begin(),&*_matsord.begin(),mats.size());
 
  Product<Orbital> creators, annihilators;
  EquiVertices equimat;
  uint prev = 0;//just to give some value
  uint currvert = 0;
  // vertices of the previous matrix
  JointVertices verts;
  
  _foreach_cauto(Order,im,_matsord){
    const Matrix mat = mats[*im];
    uint nextvert = currvert+mat.nvertices();
    // creators and annihilators orbitals
    creators *= mat.crobs();
    annihilators *= mat.crobs(true);
    if ( InSet(mat.type(),Ops::Deexc0,Ops::Exc0) ){
      // sort external orbitals in order to have the same indices always on the same places
      std::sort(creators.begin()+currvert,creators.end());
      std::sort(annihilators.begin()+currvert,annihilators.end());
    }
    if ( creators.size() < nextvert ) creators.resize(nextvert);
    if ( annihilators.size() < nextvert ) annihilators.resize(nextvert);
    // equivalent matrices?
    if ( verts.size() > 0 && mat.equivalent(mats[prev]) ) {
      equimat.push_back(verts);
    } else if ( equimat.size() > 0 ){
      equimat.push_back(verts);
      _equivs.push_back(equimat);
      equimat.clear();
    }
    verts.clear();
    for ( uint vert = currvert; vert < nextvert; ++vert) {
      verts.push_back(vert);
    }
    prev = *im;
    // equivalent vertices?
    Equivalents equivs( mat.equivertices(currvert) );
    _equivs.insert(_equivs.end(),equivs.begin(),equivs.end());
    currvert = nextvert;
  }
  if ( equimat.size() > 0 ){
    equimat.push_back(verts);
    _equivs.push_back(equimat);
  }
  
  // connections
  assert(creators.size() == annihilators.size()); // even for non-conserving operators
  Orbital dummy;
  uint nverts = creators.size();
  _foreach_cauto( Product<Orbital>,itorb,creators ){
    if ( *itorb == dummy ) {
      // an invalid connection
      _vertconn.push_back(nverts);
    } else {
      // find in the annihilators
      int ipos = annihilators.find(*itorb);
      if ( ipos < 0 ) Error("Implementation error: no annihilator for a creator!");
      _vertconn.push_back(ipos);
    }
    _orbtypes.push_back(itorb->type());
  }
  
}

Product< Matrix > UniGraph::ordmats() const
{
  Product<Matrix> mats;
  assert( pTerm );
  const Product<Matrix>& origmats = pTerm->mat();
  assert( _matsord.size() == origmats.size());
  _foreach_cauto(Order,im,_matsord)
    mats.push_back(origmats[*im]);
  return mats;
}

void UniGraph::minimize()
{
  // order of vertices
  Order vertorder;
  vertorder.init(_vertconn.size());
  // last value for not connected vertices
  vertorder.push_back(_vertconn.size());
  Order connections(_vertconn), minconn(_vertconn);
  bool nextperm;
  uint minorder = 0, iord = 0;
  do {
    // next permutation of ieqv'th equivalent vertices
    nextperm = false;
    for ( uint ieqv = 0; ieqv < _equivs.size() && !nextperm; ++ieqv ) {
      nextperm = _equivs[ieqv].next_permutation(vertorder);
    }
    if ( nextperm ) {
      ++iord;
      // create new connection vector and compare to the old one
      for ( uint i = 0; i < _vertconn.size(); ++i )
        connections[vertorder[i]] = vertorder[_vertconn[i]];
//       xout << vertorder << " --> " << connections << std::endl;
      if ( connections < minconn ) {
        minconn = connections;
        minorder = iord;
      }
    }
  } while (nextperm); 
  xout << "Smallest connection vector (" << minorder<< "): " << minconn << std::endl;
}


std::ostream& operator<<(std::ostream& o, const UniGraph& ug)
{
  const Order& conns = ug.connections();
  Product<Matrix> mats = ug.ordmats();
  const Equivalents& equivs = ug.equivals();
  o << mats;
  o << equivs;
  o << "{";
  _foreach_cauto(Order,ic,conns){
    o << *ic << " ";
  }
  o << "}";
  return o;
}
