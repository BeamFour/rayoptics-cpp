// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Goptical: Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
// Licensed under the GNU General Public License, version 3 or later; see LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.mathlib.Quaternion
#include "redukti/mathlib/Quaternion.h"

#include "redukti/JavaSemantics.h"
#include "redukti/Text.h"

namespace redukti::mathlib {

std::string Quaternion::toString() const {
    return "[" + doubleToString(x) + "," + doubleToString(y) + "," + doubleToString(z) +
           "," + doubleToString(w) + "]";
}

bool Quaternion::equals(const Quaternion &other) const {
    return doubleCompareEquals(other.x, x) && doubleCompareEquals(other.y, y) &&
           doubleCompareEquals(other.z, z) && doubleCompareEquals(other.w, w);
}

} // namespace redukti::mathlib
