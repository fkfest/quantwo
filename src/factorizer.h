#ifndef Factorizer_H
#define Factorizer_H

#include <string>
#include <set>
#include <iostream>
#include <assert.h>
#include <stdint.h>
#include "globals.h"
#include "types.h"
#include "product.h"
#include "orbital.h"
#include "inpline.h" // for name-handling
#include "sum.h"
#include "term.h"
#include "tensor.h"
#include "action.h"
#include "expression.h"

namespace Translators{
  SlotType orb2slot(const Orbital& orb);
  Tensor mat2tensor(const Matrix& mat, const std::map<Orbital,const SlotType*>& slotorbs);
  Diagram term2diagram(const Term& term, Factor fact, const std::map<Orbital,const SlotType*>& slotorbs, const Expression& expr);
}

class Factorizer {
public:
  Factorizer(const TermSum& s);



  Expression _expression;
};



#endif

