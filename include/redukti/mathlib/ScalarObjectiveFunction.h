// C++ port of org.redukti.mathlib.ScalarObjectiveFunction
#ifndef REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H
#define REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H

namespace redukti::mathlib {

/**
 * Java declares eval() as returning the boxed Double, but every call site
 * immediately unboxes it and no implementation returns null, so this returns a
 * plain double.
 *
 * Modelled as an abstract base rather than std::function: every implementation
 * in the codebase is a named class (Trace.SecantFunction, VigCalc's pupil
 * coordinate functions, Wideangle.Eval_Z_Enp_Function), never a lambda.
 */
class ScalarObjectiveFunction {
public:
    virtual ~ScalarObjectiveFunction() = default;
    virtual double eval(double x) = 0;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_SCALAROBJECTIVEFUNCTION_H
