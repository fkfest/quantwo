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
class Term;
typedef Sum<Term,TFactor> TermSum;
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
    //! validate term and finalize settings
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
    TermSum times(const TermSum& s) const;
    //! multiply a term by a sum of matrices
    TermSum times(const Sum<Matrix,TFactor>& s) const;
    //! add permutator
    Term & operator += (Permut const & perm);
    //! add permutator with a factor
    Term & operator += (std::pair<Permut,TFactor> const & p);
    //! add connections
    void addconnection (Product<long int> const & connections);
    //! add summation indices
    void addsummation (Orbital const & orb, short excl);
    void addsummation (const Product<Orbital> & orbs);
    void addsummation( const Orbital& orb );
    // add orbital index
    void addorb( const Orbital& orb );
    // add overlaps in virtual space
    void addoverlaps();
    //! replace matrix on position ipos
    void replacematrix (Matrix const & mat, unsigned long int ipos);
    //! return contained Product<SQOp>
    Product<SQOp> opProd() const;
    //! return contained Product<Kronecker>
    Product<Kronecker>  kProd() const;
    //! return prefactor
    TFactor prefac() const;
    // return termsum factors
    Product<TermSum> termsfacs() const { return _termsfacs; };
    void set_prefac(const TFactor& fac) { _prefac = fac; };
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
    // expand termsum-prefactors (from e.g. general normal ordered operators)
    TermSum expandtermsfacs();
    //! calculate normal ordering
    TermSum  normalOrder() const;
    //! calculate normal ordering, fully contracted terms only
    TermSum  normalOrder_fullyContractedOnly() const;
    //! calculate normal ordering in Particle/Hole formalism
    TermSum  normalOrderPH() const;
    //! calculate normal ordering in Particle/Hole formalism, fully contracted terms only
    TermSum  normalOrderPH_fullyContractedOnly() const;
    typedef std::list<int> TWMats;
    typedef std::list<TWMats> TWOps;
    //! Wick's theorem: call recursive routine wick
    // if genw == true: use the generalized Wick's theorem
    TermSum  wickstheorem(bool genw = false, int noord = 0) const;
    //! Wick's theorem, recursive: opers contains index of SQop in _opProd (divided into individual operators)
    TermSum  wick(TWOps& opers, TWMats& krons) const;
    //! generalized Wick's theorem, recursive: opers contains index of SQop in _opProd (divided into individual operators)
    TermSum  genwick(Term::TWOps& opers, const Term::TWMats& krons, Term::TWMats densmat) const;
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
    TermSum expand_antisym();
    //! one-electron matrices to fock
    TermSum oneel2fock(std::string decoration = "", bool multiref = false);
    //! change matrix imat in _mat from oneel to fock
    TermSum change2fock(uint imat, const std::string& decoration, bool multiref = false) const;
    //! replace E0 to <0|F|0>
    TermSum replaceE0(bool multiref = false);
    //! replace E0 at position imat in _mat to <0|F|0>
    TermSum replaceE0byfock(uint imat, bool multiref = false, bool replaceE0act = false) const;
    //! true if has a non-singlet density matrix
    bool has_nonsingldm() const;
    //! reorder density matrices to singlet order
    TermSum dm2singlet();
    //! check electrons in DMs and make electron-deltas (returns true if replaced electrons in term)
    bool dmelectrons(uint imat);
    TermSum dmwickstheorem(const Matrix& dm) const;
    TermSum dmwick( Term::TWMats& opers, const Term::TWMats& krons, const Matrix& dm ) const;
    //! true if has general indices
    bool has_generalindices() const;
    //! replace remaining general indices by occupied (and active) orbitals
    TermSum removegeneralindices();
    //! Spin integration (if notfake false: calculate only _nloops, _nintloops, _nocc)
    void spinintegration(bool notfake);
    //!Maximize overlap of two body integrals and amplitudes with <\Phi^{ij}_{ab}|
    void order();
    void maxloops();
     //!returns true if a loop with orbs1 and orbs2 is possible
    bool loop(Product<Orbital> orbs1, Product<Orbital> orbs2);
    bool loop(Array<Product<Orbital>>& elecorbs1, Array<Product<Orbital>>& elecorbs2);
    //!permute T using its antisymmetry and return the new terms as a TermSum
    TermSum addpermuteT(const TFactor fac);
    //! permute the first and third orbital in _orb of the j-th amplitude in _mat
    void permuteT(uint j);
    //! Spin expansion
    TermSum spinexpansion(const TFactor fac);
    //permute spins in TOrbSet in all possible (and partially unphysical ways)
    std::vector<TOrbSet> spinpermute(const TOrbSet& genorbs);
    //returns all possible permutations of spin in norbs orbitals
    std::vector<std::vector<Spin::Type>> spinpermutations(int norbs);
    //solves the combinatorical problem of spinpermutations() recursively
    void recspinperm(std::vector<Spin::Type> spintypes, std::vector<Spin::Type> comb, std::vector<std::vector<Spin::Type>>& spinperms, int nset, int ndraw);
    //check spin of term
    Return::Vals check_spin() const;
    //set spin of Term according to spin in TOrbSet
    void replace(TOrbSet& orbs);
    //select term if external orbitals have Spin::Type
    Return::Vals selectspin(const std::vector<Spin::Type>& spins) const;
    bool samespin() const;
    //! set prefactor of term to one
    void reset_prefac();
    // set permutations to p
    void setperm(const Sum<Permut,TFactor>& p){_perm = p;};
    // permute the term according to p
    void permute(const Permut& p);
    // resolve permutations
    TermSum resolve_permutations() const;
    //! compare actual connections with the needed (in _connections)
    //! return true if the term is ok
    bool properconnect() const;
    //! print diagram, which corresponds to this term
    void printdiag(Output* pout) const;
    //! return free orbital name
    Orbital freeorbname(Orbital::Type type, bool spinfree = false);
    //! static wrapper-function to be able to callback the member function freeorbname()
    static Orbital getfreeorbname(void * Obj, Orbital::Type type);
    //! return next free electron
    Electron nextelectron();
    //! static wrapper-function
    static Electron getnextelectron(void * Obj);
    //! set last orbital (onlylarger: only if it's larger than current one)
    void set_lastorb(Orbital orb, bool onlylarger = false);
    void set_lastel(Electron el, bool onlylarger = false);
    void set_isinput(bool isinputterm) {_isinputterm = isinputterm;};
    bool get_isinput() const {return _isinputterm;};
    // set _lastorb using _orbs
    void set_lastorbs();
    // copy lastorbs and last electron from another term
    void copy_lastorbs( const Term & term ){ _lastorb = term._lastorb; _lastel = term._lastel;};
    // iterates over all orbitals in the term. Returns Orbital() if iorb > last one
    Orbital orb( uint iorb ) const;
    void set_no_el();
    void clear_opProd() {_opProd.clear();}
    //! return matrices
    Product<Matrix>& get_mat(){return _mat;};
  private:
    TermSum  normalOrder(bool fullyContractedOnly) const;
    TermSum  normalOrderPH(bool fullyContractedOnly) const;

    Product<SQOp> _opProd;
    // termsum factors from e.g. generalized normal ordered SQ-operators
    // should be cleared before wickstheorem by calling expandtermsfacs
    Product<TermSum> _termsfacs;
    Product<Kronecker>  _kProd;
    Product<Matrix> _mat;
    // orbitals in term and summation orbitals
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
    // stores whether the matconnections have been set
    bool _matconnectionsset;
    bool _isinputterm = false;
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
  template <class T, class Q>
  Return replace(Sum<T,TFactor> &p, Q orb1, Q orb2);
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


