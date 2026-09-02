// C++ port of org.redukti.mathlib.RootResult
#include "redukti/mathlib/RootResult.h"

#include "redukti/Text.h"

namespace redukti::mathlib {

std::string RootResult::toString() const {
    return "RootResult{root=" + doubleToString(root) +
           ", converged=" + (converged ? "true" : "false") +
           ", iterations=" + intToString(iterations) + "}";
}

} // namespace redukti::mathlib
