#ifndef Equation_H
#define Equation_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>
#include "utilities.h"
#include "product.h"
#include "term.h"
#include "matrices.h"
#include "globals.h"

/*!
    Equation
 */
class Equation : public Sum<Term, TFactor> {
public:
  
};

#endif
