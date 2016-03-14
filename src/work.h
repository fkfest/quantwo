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
#include "globals.h"
#include "term.h"
#include "unigraph.h"
#include "finput.h"
#include "factorizer.h"
#include <iostream>


namespace Q2
{
  TermSum reduceSum(TermSum s);
  TermSum Kroneckers(TermSum s);
  TermSum SingletDM(TermSum s);
  bool has_nonsingldm(const TermSum& s);
  TermSum OneEl2Fock(TermSum s);
  TermSum GeneralIndices(TermSum s);
  bool has_generalindices(const TermSum& s);
  TermSum ZeroTerms(TermSum s);
  TermSum EqualTerms(TermSum s, double minfac);
  TermSum SmallTerms(TermSum s, double minfac);
  TermSum ResolvePermutaions(TermSum s);
  
  TermSum normalOrderPH(TermSum s);
  TermSum wick(TermSum s);
  TermSum postaction(TermSum s);
  void printdiags(Output* pout, TermSum s);
  
  void printalgo(std::ofstream& out, TermSum s);
}
#endif

