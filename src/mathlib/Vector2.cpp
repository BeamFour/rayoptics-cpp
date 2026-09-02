// C++ port of org.redukti.mathlib.Vector2
#include "redukti/mathlib/Vector2.h"

#include "redukti/Text.h"

namespace redukti::mathlib {

std::string Vector2::toString() const {
    return "[" + doubleToString(x) + "," + doubleToString(y) + "]";
}

} // namespace redukti::mathlib
