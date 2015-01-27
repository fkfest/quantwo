#ifndef KRONECKER_H
#define KRONECKER_H

#include <iostream>
#include "orbital.h"
#include "globals.h"


/*!
    Implements a Kronecker symbol \[ \delta_{ij} \]
*/

class Kronecker {
public:
    //! construct from given orbitals
    Kronecker(Orbital orb1, Orbital orb2);
    //! get first index
    Orbital orb1() const;
    //! get second index
    Orbital orb2() const;
    //! "artificial" ordering relation (for sorting purposes)
    bool operator < (Kronecker const & k) const;
    // replace orbital orb1 with orb2
    Return replace(Orbital orb1, Orbital orb2, bool smart);
    // replace spin spin1 with spin2
    Return replace(Spin spin1, Spin spin2, bool smart);

private:
Orbital    _orb1;
Orbital    _orb2;
};


//! output operator for Kronecker
std::ostream & operator << (std::ostream & o, Kronecker const & k);

#endif

