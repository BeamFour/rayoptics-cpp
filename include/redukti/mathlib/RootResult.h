// C++ port of org.redukti.mathlib.RootResult
#ifndef REDUKTI_MATHLIB_ROOTRESULT_H
#define REDUKTI_MATHLIB_ROOTRESULT_H

#include <string>

namespace redukti::mathlib {

class RootResult {
public:
    double root;
    bool converged;
    int iterations;

    RootResult(double root_, bool converged_, int iterations_)
        : root(root_), converged(converged_), iterations(iterations_) {}

    std::string toString() const;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_ROOTRESULT_H
