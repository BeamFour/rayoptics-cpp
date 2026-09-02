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
