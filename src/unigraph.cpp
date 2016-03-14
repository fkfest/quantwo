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
    const Matrix& mat = mats[*im];
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
      // store external orbitals for later 
      _extorbs_crea *= Product<Orbital>(creators.begin()+currvert,creators.end());
      _extorbs_anni *= Product<Orbital>(annihilators.begin()+currvert,annihilators.end());
      _extvertices += verts;
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
  _eqperms.sort();
  _eqperm_from.sort();
}

bool UniGraph::is_equal(const UniGraph& ug) const
{
  if ( _matsord.size() != ug._matsord.size() ) return false;
  const Product<Matrix>& mats = pTerm->mat();
  const Product<Matrix>& ugmats = ug.pTerm->mat();
  for ( uint i = 0; i < _matsord.size(); ++i ) {
    if ( ! mats[_matsord[i]].equivalent(ugmats[ug._matsord[i]]) ) return false;
  }
  return _orbtypes == ug._orbtypes && _vertconn == ug._vertconn;
  
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
  epf.sort();
  epfo = epf;
  ep.minpermute(connections,_eqperms);
  epf.minpermute(connections,epfo);
}
void UniGraph::gen_perms(const PermVertices& from, const PermVertices& to)
{
  assert( from.size() == to.size() );
  PermVertices::const_iterator itpv = to.begin();
  for( const JointVertices& fpv: from){
    assert( itpv->size() == fpv.size() );
    JointVertices::const_iterator it = itpv->begin();
    for ( const uint& fp: fpv ) {
      if ( *it != fp ) {
        assert( _perms.count(*it) == 0 );
        _perms[*it] = fp;
      }
      ++it;
    }
    ++itpv;
  }
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
  PermVertices eq_perms(_eqperms), eq_perm_from(_eqperm_from), eq_perm_from_orig(eq_perm_from),
               min_eq_perms, min_eq_perm_from;
  if ( permute4each_vertorder ) {
    min_eq_perms = eq_perms;
    min_eq_perm_from = eq_perm_from;
  }
  Order
    minvertorder(vertorder),
    connections(_vertconn);
  if (permute4each_vertorder) {
    // allowed permutations
    apply_eqperms(connections,vertorder,eq_perms,eq_perm_from,eq_perm_from_orig);
  }
  Order minconn(connections);
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
        if (permute4each_vertorder) {
          min_eq_perms = eq_perms;
          min_eq_perm_from = eq_perm_from;
        }
      }
    }
  } while (nextperm); 
  if ( !permute4each_vertorder ) {
    // try to minimize further by using allowed permutations
    apply_eqperms(minconn,minvertorder,eq_perms,eq_perm_from,eq_perm_from_orig);
    gen_perms(_eqperms,eq_perms);
    gen_perms(eq_perm_from_orig,eq_perm_from);
  } else {
    eq_perm_from_orig.set_with_order(_eqperm_from,minvertorder);
    gen_perms(_eqperms,min_eq_perms);
    gen_perms(eq_perm_from_orig,min_eq_perm_from);
  }
  _vertconn = minconn;
  xout << "Smallest connection vector (" << minorder<< "): " << minvertorder << " --> " << minconn << std::endl;
  if ( _perms.size() > 0 ) {
    xout << "Permutations: ";
    for ( const auto& perm: _perms ){
      xout << perm.first;
    }
    xout << " --> ";
    for ( const auto& perm: _perms ){
      xout << perm.second;
    }
    xout << std::endl;
  }
}

Term UniGraph::gen_term() const
{
  Term term;
  // it works for spinfree stuff only for now
  bool spinfree = true;
  // create annihilator-connection vector from the connection vector
  uint nverts = _vertconn.size();
  Order annicon(_vertconn.size(),nverts);
  for ( uint icon = 0; icon < _vertconn.size(); ++icon )
    if ( _vertconn[icon] < nverts )
      annicon[_vertconn[icon]] = icon;
  
  // create a list of orbitals from the type-list (replace external orbitals by the saved ones)
  for ( const auto& orb: _extorbs_crea )
    term.set_lastorb(orb,true);
  for ( const auto& orb: _extorbs_anni )
    term.set_lastorb(orb,true);
  Array<Orbital> orbs(nverts);
  // set external orbitals
  uint iexorbc = 0, iexorba = 0;
  for ( uint exvert: _extvertices ) {
    assert( exvert < nverts );
    if ( _vertconn[exvert] < nverts ) {
      orbs[exvert] = _extorbs_crea[iexorbc];
      term.addorb(orbs[exvert]);
      ++iexorbc;
    }
    if ( annicon[exvert] < nverts ) {
      orbs[annicon[exvert]] = _extorbs_anni[iexorba];
      term.addorb(orbs[annicon[exvert]]);
      ++iexorba;
    }
  }
  assert( _orbtypes.size() == nverts );
  for ( uint ivert = 0; ivert < nverts; ++ivert ) {
    if ( orbs[ivert].type() == Orbital::NoType && _orbtypes[ivert] != Orbital::NoType ) {
      // not set before
      orbs[ivert] = term.freeorbname(_orbtypes[ivert],spinfree);
      term.addorb(orbs[ivert]);
      term.addsummation(orbs[ivert]);
    }
  }
  // set permutations
  Permut permuts;
  for ( const auto& perm: _perms ){
    assert( perm.first < orbs.size() && perm.second < orbs.size() );
    permuts += Permut(orbs[perm.first],orbs[perm.second]); 
  }
  Sum<Permut,TFactor> sumpermut;
  sumpermut += std::make_pair(permuts,pTerm->prefac());
  term.setperm(sumpermut);
  // add matrices to the term
  const Product<Matrix>& mats = pTerm->mat();
  const Product<Matrix> newmats;
  uint currvert = 0;
  _foreach_cauto(Order,im,_matsord){
    const Matrix& mat = mats[*im];
    uint nextvert = currvert + mat.nvertices();
    assert( nextvert <= nverts );
    // use connection vectors and list of orbitals to create the product of orbitals of mat
    Product<Orbital> orbcre, orbani;
    for ( uint ivert = currvert; ivert < nextvert; ++ivert ){
      if ( _vertconn[ivert] < nverts ) {
        orbcre.push_back(orbs[ivert]);
      }
      if ( annicon[ivert] < nverts ) {
        orbani.push_back(orbs[annicon[ivert]]);
      }
    }
    // add the matrix to the term
    term *= Matrix(mat.type(),orbcre,orbani,mat.npairs(),mat.lmel(),mat.name(),mat.matspinsym(),
                   mat.antisymform());
    currvert = nextvert;
  }
  
  return term;
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
