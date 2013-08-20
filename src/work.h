#ifndef Work_H
#define Work_H

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include "utilities.h"
#include "product.h"
#include "operators.h"
#include "kronecker.h"
#include "sum.h"
#include "globals.h"
#include "term.h"
#include "finput.h"

#include <iostream>


namespace Q2
{
  Sum<Term,TFactor> reduceSum(Sum<Term,TFactor> s);
  Sum<Term,TFactor> Kroneckers(Sum<Term,TFactor> s);
  Sum<Term,TFactor> EqualTerms(Sum<Term,TFactor> s, double minfac);
  Sum<Term,TFactor> SmallTerms(Sum<Term,TFactor> s, double minfac);
  Sum<Term,TFactor> ResolvePermutaions(Sum<Term,TFactor> s);
  
  Sum<Term,TFactor> normalOrderPH(Sum<Term,TFactor> s);
  Sum<Term,TFactor> wick(Sum<Term,TFactor> s);
  Sum<Term,TFactor> postaction(Sum<Term,TFactor> s);
  void printdiags(Output* pout, Sum<Term,TFactor> s);
}
#endif