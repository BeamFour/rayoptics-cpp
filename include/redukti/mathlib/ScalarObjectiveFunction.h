// C++ port of org.redukti.mathlib.ScalarObjectiveFunction
#ifndef REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H
#define REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H

#include <optional>

namespace redukti::mathlib {

/**
 * eval() returns Java's boxed Double, and the null really is used: VigCalc's
 * Fn_r_pupil_coordinate and Wideangle's wrappers return null when the ray fails
 * before the surface of interest, and Wideangle::find_edge branches on
 * `fc == null` to decide which half of the bracket to keep. So this is
 * std::optional<double>, not double.
 *
 * The root solvers (Brent, Secant) auto-unbox in the Java and would throw a
 * NullPointerException on a null; they call .value() here, which throws
 * std::bad_optional_access in the same situation.
 *
 * Modelled as an abstract base rather than std::function: every implementation
 * in the codebase is a named class (Trace.SecantFunction, VigCalc's pupil
 * coordinate functions, Wideangle.Eval_Z_Enp_Function), never a lambda.
 */
class ScalarObjectiveFunction {
public:
    virtual ~ScalarObjectiveFunction() = default;
    virtual std::optional<double> eval(double x) = 0;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H
