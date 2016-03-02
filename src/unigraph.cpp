#include "unigraph.h"

UniGraph::UniGraph(const Term& term)
{
  int permuteq = Input::iPars["prog"]["permuteq"];
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
    // equivalent matrices?
    if ( verts.size() > 0 && mat.equivalent(mats[prev]) ) {
      equimat.add(verts);
    } else if ( equimat.size() > 0 ){
      equimat.add(verts);
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
    // creators and annihilators orbitals
    creators *= mat.crobs();
    annihilators *= mat.crobs(true);
    if ( InSet(mat.type(),Ops::Deexc0,Ops::Exc0) ){
      // sort external orbitals in order to have the same indices always on the same places
      std::sort(creators.begin()+currvert,creators.end());
      std::sort(annihilators.begin()+currvert,annihilators.end());
      if ( permuteq > 0 ) {
        // allow permutations
        _eqperms.push_back(verts);
      }
    }
    // for non-conserved vertices
    if ( creators.size() < nextvert ) creators.resize(nextvert);
    if ( annihilators.size() < nextvert ) annihilators.resize(nextvert);
    currvert = nextvert;
  }
  if ( equimat.size() > 0 ){
    equimat.add(verts);
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
  // "from-verices" for allowed permutations
  for(JointVertices& pvs: _eqperms) {
    JointVertices fromvert;
    bool modify_pvs = false;
    for(uint& pv: pvs) {
      int ipos = _vertconn.find(pv);
      if ( ipos >= 0 ) 
        fromvert.push_back(ipos);
      if ( _vertconn[pv] == nverts ) {
        // no creator on this one, remove it from the list
        pv = nverts;
        modify_pvs = true;
      }
    }
    _eqperm_from.push_back(fromvert);
    if (modify_pvs) {
      pvs.erase(std::remove(pvs.begin(),pvs.end(),nverts),pvs.end());
    }
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

void UniGraph::apply_eqperms( Order& connections, const Order& vertorder,
                              PermVertices& ep, PermVertices& epf, PermVertices& epfo )
{
    ep = _eqperms;
    epf.set_with_order(_eqperm_from,vertorder);
    epfo = epf;
    ep.minpermute(connections,_eqperms);
    epf.minpermute(connections,epfo);
}

void UniGraph::minimize()
{
  int permuteq = Input::iPars["prog"]["permuteq"];
  bool permute4each_vertorder = (permuteq > 1);
  // order of vertices
  Order vertorder;
  vertorder.identity(_vertconn.size());
  // last value for not connected vertices
  vertorder.push_back(_vertconn.size());
  Order
    minvertorder(vertorder),
    connections(_vertconn), minconn(_vertconn);
  PermVertices eq_perms(_eqperms), eq_perm_from(_eqperm_from), eq_perm_from_orig(eq_perm_from);
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
      if (permute4each_vertorder) {
        // allowed permutations
        apply_eqperms(connections,vertorder,eq_perms,eq_perm_from,eq_perm_from_orig);
      }
//       xout << vertorder << " --> " << connections << std::endl;
      if ( connections < minconn ) {
        minconn = connections;
        minvertorder = vertorder;
        minorder = iord;
      }
    }
  } while (nextperm); 
  if ( !permute4each_vertorder ) {
    // try to minimize further by using allowed permutations
    apply_eqperms(minconn,minvertorder,eq_perms,eq_perm_from,eq_perm_from_orig);
  }
  xout << "Smallest connection vector (" << minorder<< "): " << minvertorder << " --> " << minconn << std::endl;
}


std::ostream& operator<<(std::ostream& o, const UniGraph& ug)
{
  const Order& conns = ug.connections();
  Product<Matrix> mats = ug.ordmats();
  const Equivalents& equivs = ug.equivals();
  const PermVertices& eqperms = ug.eqperms();
  const PermVertices& eqperm_from = ug.eqperm_from();
  o << mats;
  o << equivs;
  o << "{";
  for(auto& ic: conns){
    o << ic << " ";
  }
  o << "}";
  o << "Perm/";
  for(auto& jv: eqperms){
    o << jv << " ";
  }
  o <<"//";
  for(auto& jv: eqperm_from){
    o << jv << " ";
  }
  o <<"/";
  
  
  return o;
}
