// C++ port of org.redukti.mathlib.SecantSolver
#include "redukti/mathlib/SecantSolver.h"

#include <cmath>

namespace redukti::mathlib {

RootResult SecantSolver::find_root(ScalarObjectiveFunction &f, double x0, int maxiter,
                                   double tol) {
    const double eps = 1e-4;
    double p0 = x0;
    double p1 = x0 * (1 + eps);
    p1 += (p1 >= 0 ? eps : -eps);
    double q0 = f.eval(p0).value();
    double q1 = f.eval(p1).value();
    if (std::abs(q1) < std::abs(q0)) {
        double tmp = p0;
        p0 = p1;
        p1 = tmp;
        tmp = q0;
        q0 = q1;
        q1 = tmp;
    }
    double p = 0.0;
    for (int i = 0; i < maxiter; i++) {
        if (q1 == q0) {
            return RootResult((p1 + p0) / 2.0, false, i);
        } else {
            if (std::abs(q1) > std::abs(q0)) {
                p = (-q0 / q1 * p1 + p0) / (1.0 - q0 / q1);
            } else {
                p = (-q1 / q0 * p0 + p1) / (1.0 - q1 / q0);
            }
            if (std::abs(p - p1) < tol) {
                return RootResult(p, true, i);
            }
        }
        p0 = p1;
        q0 = q1;
        p1 = p;
        q1 = f.eval(p1).value();
    }
    return RootResult(p, false, maxiter);
}

} // namespace redukti::mathlib
