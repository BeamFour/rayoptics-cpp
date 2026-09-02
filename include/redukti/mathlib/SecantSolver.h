// C++ port of org.redukti.mathlib.SecantSolver
#ifndef REDUKTI_MATHLIB_SECANTSOLVER_H
#define REDUKTI_MATHLIB_SECANTSOLVER_H

#include "redukti/mathlib/RootResult.h"
#include "redukti/mathlib/ScalarObjectiveFunction.h"

namespace redukti::mathlib {

class SecantSolver {
public:
    static RootResult find_root(ScalarObjectiveFunction &f, double x0, int maxiter,
                                double tol);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_SECANTSOLVER_H
