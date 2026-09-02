// C++ port of org.redukti.mathlib.Vector3
#include "redukti/mathlib/Vector3.h"

#include "redukti/Text.h"

#include <cmath>

namespace {
/** Mirrors `Double.compare(a, b) == 0`, which distinguishes -0.0 from 0.0. */
bool javaDoubleEquals(double a, double b) {
    if (std::isnan(a) || std::isnan(b))
        return std::isnan(a) && std::isnan(b);
    if (a == b)
        return std::signbit(a) == std::signbit(b);
    return false;
}
} // namespace

namespace redukti::mathlib {

std::string Vector3::toString() const {
    return "[" + doubleToString(x) + "," + doubleToString(y) + "," +
           doubleToString(z) + "]";
}

bool Vector3::equals(const Vector3 &other) const {
    return javaDoubleEquals(other.x, x) && javaDoubleEquals(other.y, y) &&
           javaDoubleEquals(other.z, z);
}

} // namespace redukti::mathlib
