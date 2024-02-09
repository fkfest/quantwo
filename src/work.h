#ifndef Work_H
#define Work_H

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include "utilities.h"
#include "product.h"
#include "operators.h"
#include "kronecker.h"
#include "sum.h"
#include "arrays.h"
#include "globals.h"
#include "term.h"
#include "unigraph.h"
#include "finput.h"
#include "factorizer.h"
#include <iostream>


namespace Q2
{
  TermSum evalEq(Finput& finput);
  TermSum reduceSum(TermSum s);
  TermSum Kroneckers(const TermSum& s);
  TermSum SingletDM(TermSum s);
  bool has_nonsingldm(const TermSum& s);
  TermSum OneEl2Fock(const TermSum& s);
  TermSum ReplaceE0(const TermSum& s);
  TermSum GeneralIndices(TermSum s);
  bool has_generalindices(const TermSum& s);
  TermSum ZeroTerms(const TermSum& s);
  TermSum EqualTerms(const TermSum& s, double minfac);
  TermSum SmallTerms(const TermSum& s, double minfac);
  TermSum VirtSpace(const TermSum& s);
  TermSum ResolvePermutaions(const TermSum& s, bool inputterms = false);

  TermSum normalOrderPH(const TermSum& s);
  TermSum wick(const TermSum& s);
  //!helps EqualTerms to find equal terms. Assumes full anti-symmetry of amplitudes!
  TermSum PreConditioner(const TermSum& s);
  TermSum postaction(const TermSum& s);
  void printdiags(Output* pout, const TermSum& s);

  void printalgo(std::ofstream& out, const TermSum& s);
}
#endif

