#ifndef Term_H
#define Term_H

#include <vector>
#include <stdlib.h>
#include <cmath>
#include "utilities.h"
#include "product.h"
#include "operators.h"
#include "kronecker.h"
#include "sum.h"
#include "globals.h"

#include <iostream>


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
    
    //! construct from Product<SQOp>, Product<Kronecker>, Product<Matrices>, summation indices and prefactor
    Term(Product<SQOp> const & opProd, Product<Kronecker> const & kProd, 
         Product<Matrices> const & mat, Product<Orbital> const & sumindx, const Product< Orbital >& realsumindx, 
         double const & prefac, const std::vector< Product<long int> >& connections);
    
    //! validate term
    bool term_is_valid();
    //! append Operator
    Term & operator *= (Oper const & t);
    
    //! multiply by a factor
    Term & operator *= (double const & fac);
    
    //! add permutator
    Term & operator += (Permut const & perm);
    
    //! add permutator with a factor
    Term & operator += (std::pair<Permut,double> const & p);
    
    //! add connections
    void addconnection (Product<long int> const & connections);
    
    //! add summation indices
    void addsummation (Orbital const & orb, short excl);
    
    //! add matrix
    void addmatrix (Matrices const & mat);
     
    //! replace matrix on position ipos
    void replacematrix (Matrices const & mat, unsigned long int ipos);
    
    //! return contained Product<SQOp>
    Product<SQOp> opProd() const;

    //! return contained Product<Kronecker>
    Product<Kronecker>  kProd() const;
    
    //! return prefactor
    double prefac() const;
    
    //! return matrices
    Product<Matrices> mat() const;
    
    //! return summation indices
    Product<Orbital> sumindx() const;
    
    //! return real summation indices
    Product<Orbital> realsumindx() const;
    
    //! generate Product of external-lines orbitals
    Product<Orbital> extindx() const;
    
    //! return Sum of Permutators
    Sum<Permut,double> perm() const;
    
    //! return connections
    std::vector< Product<long int> > connections() const;
    
    //! return true if term is zero
    bool term_is_0(double minfac) const;

    //! artificial ordering
    bool operator < (Term const & t) const;
    
    //! equal terms
    // terms will be not changed! (but const can't be applied) 
    bool equal(Term & t, Permut & perm);
    bool equal_old(Term & t, Permut & perm);
    
    //! calculate normal ordering
    Sum<Term, double>  normalOrder() const;

    //! calculate normal ordering, fully contracted terms only
    Sum<Term, double>  normalOrder_fullyContractedOnly() const;
    
    //! calculate normal ordering in Particle/Hole formalism
    Sum<Term, double>  normalOrderPH() const;

    //! calculate normal ordering in Particle/Hole formalism, fully contracted terms only
    Sum<Term, double>  normalOrderPH_fullyContractedOnly() const;
    
    //! Wick's theorem: call recursive routine wick
    Sum<Term, double>  wickstheorem() const;
    //! Wick's theorem, recursive: opers contains index of SQop in _opProd (divided into individual operators)
    Sum<Term, double>  wick(std::vector< std::vector<long int> > & opers, std::vector<long int> & krons) const;
    
    //! set connections for each matrix
    void setmatconnections();
    
    //! reduce equation (delete Kroneckers and summation indices)
    void reduceTerm();
    
    //! Determine connections (in reduced term!)
    void matrixkind();
    
    //! expand integral ( from antisymmetrized form < AB || CD > to the normal form < AB | CD > - < AB | DC > )
    //! if firstpart=true : < AB | CD >, if firstpart=false : < AB | DC >
    //! if return is true: expanded, if false: don't need to expand
    bool expandintegral(bool firstpart);
    
    //! check if we have any antisymmetrized matrices in term
    bool antisymmetrized();
    
    //! expand all antisymmetrical matrices in term, if spinintegr true: do spin integration
    Sum<Term,double> expand_antisym(bool spinintegr);
    
    //! Spin integration (if notfake false: calculate only _nloops, _nintloops, _nocc)
    void spinintegration(bool notfake);

    //! set prefactor of term to one
    void reset_prefac();
    
    //! compare actual connections with the needed (in _connections)
    //! return true if the term is ok
    bool properconnect() const;
    
    //! return free orbital name
    Orbital freeorbname(Orbital::Type type);
    //! static wrapper-function to be able to callback the member function freeorbname()
    static Orbital getfreeorbname(void * Obj, Orbital::Type type);
//    //! get last orbital
//    Orbital get_lastorb(Orbital::Type type) const { return _lastorb[type]; };
    //! set last orbital (onlylarger: only if it's larger than current one)
    void set_lastorb(Orbital orb, bool onlylarger = false);

  private:
    Sum<Term, double>  normalOrder(bool fullyContractedOnly) const;
    Sum<Term, double>  normalOrderPH(bool fullyContractedOnly) const;

    Product<SQOp> _opProd;
    Product<Kronecker>  _kProd;
    Product<Matrices> _mat;
    Product <Orbital> _sumindx,_realsumindx;
    double _prefac;
    // connections of matrices in term
    // (abs(value)-1) gives the index of the corresponding matrix in _mat
    // positive sign: connected group of matrices Product<long int>
    // negative sign: disconnected --------------"-----------------
    std::vector< Product<long int> > _connections;

    Sum<Permut,double> _perm;
    std::map< Orbital::Type, Orbital > _lastorb;
    // for term comparison:
    // number of (all and internal only) loops and occupied orbitals (set in spinintegration)
    unsigned long int _nloops, _nintloops, _nocc;
};

//! output operator for Term
std::ostream & operator << (std::ostream & o, Term const & t);

namespace Q2
{
  Sum<Term,double> reduceSum(Sum<Term,double> s);
  Sum<Term,double> normalOrderPH(Sum<Term,double> s);
  Sum<Term,double> wick(Sum<Term,double> s);
  template <class T>
  void replace(Product<T> &p, Orbital orb1, Orbital orb2);
  void replace(SQOp &op, Orbital orb1, Orbital orb2);
  void replace(Orbital &orb, Orbital orb1, Orbital orb2);
  void replace(Matrices &mat,Orbital orb1, Orbital orb2);
  void replace(Kronecker &kron, Orbital orb1, Orbital orb2);
}
#endif

