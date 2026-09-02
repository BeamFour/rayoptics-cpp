// C++ port of org.redukti.mathlib.BrentSolver
#ifndef REDUKTI_MATHLIB_BRENTSOLVER_H
#define REDUKTI_MATHLIB_BRENTSOLVER_H

#include "redukti/mathlib/RootResult.h"
#include "redukti/mathlib/ScalarObjectiveFunction.h"

namespace redukti::mathlib {

class BrentSolver {
public:
    static const double TOL;

    static RootResult find_root(double a, double b, ScalarObjectiveFunction &fn);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_BRENTSOLVER_H
