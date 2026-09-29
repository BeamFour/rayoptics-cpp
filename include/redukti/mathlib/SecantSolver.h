// C++ port of org.redukti.mathlib.SecantSolver
#ifndef REDUKTI_MATHLIB_SECANTSOLVER_H
#define REDUKTI_MATHLIB_SECANTSOLVER_H

#include "redukti/mathlib/RootResult.h"
#include "redukti/mathlib/ScalarObjectiveFunction.h"

#include <optional>

namespace redukti::mathlib {

class SecantSolver {
public:
    static RootResult find_root(ScalarObjectiveFunction &f, double x0, int maxiter,
                                double tol);

    /**
     * The second point the secant starts from. Without one it sits just outside x0, so
     * the search moves outward; a caller that must search inward passes it. The Java
     * takes a nullable Double here.
     */
    static RootResult find_root(ScalarObjectiveFunction &f, double x0,
                                std::optional<double> x1, int maxiter, double tol);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_SECANTSOLVER_H
