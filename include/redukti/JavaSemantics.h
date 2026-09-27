// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// Small helpers reproducing Java semantics that have no direct C++ equivalent.
#ifndef REDUKTI_JAVASEMANTICS_H
#define REDUKTI_JAVASEMANTICS_H

#include <cmath>

namespace redukti {

/**
 * Mirrors `Double.compare(a, b) == 0`, which -- unlike `==` -- treats NaN as
 * equal to itself and -0.0 as distinct from 0.0. Java's generated equals()
 * methods use it, so the ported ones do too.
 */
inline bool doubleCompareEquals(double a, double b) {
    if (std::isnan(a) || std::isnan(b))
        return std::isnan(a) && std::isnan(b);
    if (a == b)
        return std::signbit(a) == std::signbit(b);
    return false;
}

} // namespace redukti

#endif // REDUKTI_JAVASEMANTICS_H
