#ifndef Term_H
#define Term_H

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "utilities.h"
#include "globals.h"
#include "types.h"
#include "product.h"
#include "operators.h"
#include "matrix.h"
#include "kronecker.h"
#include "sum.h"

#include <iostream>

class SQOp;
class Oper;
/*!
    A term consists of a Product of SQOperators and a Product of Kroneckers
*/
class Term {
  public:
    //! default constructor 
    Term();
    //! construct from Product<SQOp>, Product<Kronecker> will be empty
    Term(Product<SQOp> const & opProd);
    //! construct from Product<SQOp> and Product<Kronecker>
    Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd);
    //! construct from Product<SQOp>, Product<Kronecker>, Product<Matrix>, summation indices and prefactor
    Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd, 
         Product<Matrix> const & mat, const TOrbSet & orbs, const TOrbSet & sumorbs, 
         const TFactor& prefac, const ConnectionsMap& connections);
    //! validate term
    bool term_is_valid();
    //! append Operator
    Term & operator *= (Oper const & t);
    //! multiply by a factor
    Term & operator *= (const TFactor& fac);
    //! multiply by a permutation
    Term & operator *= (Permut const & perm);
    Term & operator *= (Sum<Permut,TFactor> const & perm);
    //! multiply by another term 
    // TODO: replace with a template
    Term & operator *= (Term const & t);
    //! multiply by a matrix
    Term & operator *= (const Matrix & mat);
    //! multiply a term by a sum
    Sum<Term,TFactor> times(const Sum<Term,TFactor>& s) const;
    //! multiply a term by a sum of matrices
    Sum<Term,TFactor> times(const Sum<Matrix,TFactor>& s) const;
    //! add permutator
    Term & operator += (Permut const & perm);
    //! add permutator with a factor
    Term & operator += (std::pair<Permut,TFactor> const & p);
    //! add connections
    void addconnection (Product<long int> const & connections);
    //! add summation indices
    void addsummation (Orbital const & orb, short excl);
    void addsummation (const Product<Orbital> & orbs);
    //! replace matrix on position ipos
    void replacematrix (Matrix const & mat, unsigned long int ipos);
    //! return contained Product<SQOp>
    Product<SQOp> opProd() const;
    //! return contained Product<Kronecker>
    Product<Kronecker>  kProd() const;
    //! return prefactor
    TFactor prefac() const;
    //! return matrices
    const Product<Matrix>& mat() const;
    //! return orbitals
    const TOrbSet& orbs() const;
    //! return summation indices
    const TOrbSet& sumorbs() const;
    //! generate set of external-lines orbitals
    TOrbSet extindx() const;
    //! generate set of orbitals that correspond to external-creation lines
    TOrbSet extcreaindx() const;
    //! return Sum of Permutators
    Sum<Permut,TFactor> perm() const;
    //! return connections
    ConnectionsMap connections() const;
    //! return true if term is zero
    bool term_is_0(double minfac) const;
    //! return true if term has to be removed 
    bool removeit() const;
    //! artificial ordering
    bool operator < (Term const & t) const;
    //! equal terms
    // terms will be not changed! (but const can't be applied) 
    // perm: permutation which brings t-term to this term (if true at return)
    bool equal(Term & t, Permut & perm);
    //! calculate normal ordering
    Sum<Term, TFactor>  normalOrder() const;
    //! calculate normal ordering, fully contracted terms only
    Sum<Term, TFactor>  normalOrder_fullyContractedOnly() const;
    //! calculate normal ordering in Particle/Hole formalism
    Sum<Term, TFactor>  normalOrderPH() const;
    //! calculate normal ordering in Particle/Hole formalism, fully contracted terms only
    Sum<Term, TFactor>  normalOrderPH_fullyContractedOnly() const;
    typedef std::list<int> TWMats;
    typedef std::list<TWMats> TWOps;
    //! Wick's theorem: call recursive routine wick
    // if genw == true: use the generalized Wick's theorem
    Sum<Term, TFactor>  wickstheorem(bool genw = false, int noord = 0) const;
    //! Wick's theorem, recursive: opers contains index of SQop in _opProd (divided into individual operators)
    Sum<Term, TFactor>  wick(TWOps& opers, TWMats& krons) const;
    //! generalized Wick's theorem, recursive: opers contains index of SQop in _opProd (divided into individual operators)
    Sum<Term, TFactor>  genwick(Term::TWOps& opers, const Term::TWMats& krons, Term::TWMats densmat) const;
    //! set connections for each matrix
    void setmatconnections();
    //! reduce equation (delete Kroneckers and summation indices)
    void reduceTerm();
    //! reduce electrons in equation according to Kroneckers that are left after reduceTerm()
    void reduceElectronsInTerm();
    //! transform Kroneckers to matrices
    void krons2mats();
    //! replace orbital orb1 with orb2, if smart == true : handle the spins in a smart way
    void replace(Orbital orb1, Orbital orb2, bool smart = true);
    //! replace spin spin1 with spin2
    void replace(Spin spin1, Spin spin2, bool smart = true);
    // delete "None" matrices (caution, the order of matrices can be important, so do it AFTER connection stuff!)
    // if unite_exc0 is true: remove all exc0 and dexc0 matrices and create one single exc0 (and/or dexc0) 
    void deleteNoneMats(bool unite_exc0 = true);
    // if pMat is zero - sets pMat to the it-matrix and increments it (it has to be from _mat)
    // if not zero - combines pMat and it (using skipping norbs-orbitals in *it) and deletes *it
    void combineMats(Matrix *& pMat, Product<Matrix>::iterator& it, const Set<uint>& norbs);
    // brilloin condition (return true for terms with occ-virt fock)
    bool brilloin() const;
    //! Determine connections (in reduced term!)
    void matrixkind();
    //! expand integral ( from antisymmetrized form < AB || CD > to the normal form < AB | CD > - < AB | DC > )
    //! if firstpart=true : < AB | CD >, if firstpart=false : < AB | DC >
    //! if return is true: expanded, if false: don't need to expand
    bool expandintegral(bool firstpart);
    //! check if we have any antisymmetrized matrices in term
    bool antisymmetrized();
    //! expand all antisymmetrical matrices in term 
    Sum<Term,TFactor> expand_antisym();
    //! one-electron matrices to fock
    Sum<Term, TFactor> oneel2fock(bool multiref = false);
    //! change matrix imat in _mat from oneel to fock
    Sum<Term, TFactor> change2fock(uint imat, bool multiref = false) const;
    //! true if has a non-singlet density matrix
    bool has_nonsingldm() const;
    //! reorder density matrices to singlet order
    Sum<Term, TFactor> dm2singlet();
    //! check electrons in DMs and make electron-deltas (returns true if replaced electrons in term)
    bool dmelectrons(uint imat);
    Sum<Term, TFactor> dmwickstheorem(const Matrix& dm) const;
    Sum<Term, TFactor> dmwick( Term::TWMats& opers, const Term::TWMats& krons, const Matrix& dm ) const;
    //! true if has general indices
    bool has_generalindices() const;
    //! replace remaining general indices by occupied (and active) orbitals
    Sum<Term, TFactor> removegeneralindices();
    //! Spin integration (if notfake false: calculate only _nloops, _nintloops, _nocc)
    void spinintegration(bool notfake);
    //! set prefactor of term to one
    void reset_prefac();
    // set permutations to p
    void setperm(const Sum<Permut,TFactor>& p){_perm = p;};
    // permute the term according to p
    void permute(const Permut& p);
    // resolve permutations
    Sum<Term, TFactor> resolve_permutations() const;
    //! compare actual connections with the needed (in _connections)
    //! return true if the term is ok
    bool properconnect() const;
    //! print diagram, which corresponds to this term 
    void printdiag(Output* pout) const;
    //! return free orbital name
    Orbital freeorbname(Orbital::Type type);
    //! static wrapper-function to be able to callback the member function freeorbname()
    static Orbital getfreeorbname(void * Obj, Orbital::Type type);
    //! return next free electron
    Electron nextelectron();
    //! static wrapper-function
    static Electron getnextelectron(void * Obj);
    //! set last orbital (onlylarger: only if it's larger than current one)
    void set_lastorb(Orbital orb, bool onlylarger = false);
    void set_lastel(Electron el, bool onlylarger = false);
    // set _lastorb using _orbs
    void set_lastorbs();
    // copy lastorbs and last electron from another term
    void copy_lastorbs( const Term & term ){ _lastorb = term._lastorb; _lastel = term._lastel;};
    // iterates over all orbitals in the term. Returns Orbital() if iorb > last one
    Orbital orb( uint iorb ) const;
    
  private:
    Sum<Term, TFactor>  normalOrder(bool fullyContractedOnly) const;
    Sum<Term, TFactor>  normalOrderPH(bool fullyContractedOnly) const;

    Product<SQOp> _opProd;
    Product<Kronecker>  _kProd;
    Product<Matrix> _mat;
    TOrbSet _orbs,_sumorbs;
    TFactor _prefac;
    Sum<Permut,TFactor> _perm;
    // connections of matrices in term
    // (abs(value)-1) gives the index of the corresponding matrix in _mat
    // positive sign: connected group of matrices Product<long int>
    // negative sign: disconnected --------------"-----------------
    ConnectionsMap _connections;

    std::map< Orbital::Type, Orbital > _lastorb;
    Electron _lastel;
    // for term comparison:
    // number of (all and internal only) loops and occupied orbitals (set in spinintegration)
    lui _nloops, _nintloops, _nocc;
};

//! output operator for Term
std::ostream & operator << (std::ostream & o, Term const & t);

class Termel : public Term {
public:
    Termel(){};
    Termel(const Term& term) : Term(term){};
private:
  TElSet _els, _sumels;
};

namespace Q2
{
  template <class T, class Q>
  Return replace(Product<T> &p, Q orb1, Q orb2, bool smart);
  template <class T, class Q>
  Return replace(Set<T> &p, Q orb1, Q orb2, bool smart);
  template <class T>
  Return replace(SQOp &op, T orb1, T orb2, bool smart);
  template <class T, class Q>
  Return replace(T &orb, Q orb1, Q orb2, bool smart);
  template <class T>
  Return replace(Matrix &mat, T orb1, T orb2, bool smart);
  template <class T>
  Return replace(Kronecker &kron, T orb1, T orb2, bool smart);
  template < class T >
  typename T::iterator findSpin(T orbs, const Spin& spin);
}
#endif


