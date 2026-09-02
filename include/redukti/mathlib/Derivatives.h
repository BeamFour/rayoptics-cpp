// C++ port of org.redukti.mathlib.Derivatives / DerivFunction / DerivResult
//
// Ported from deriv/deriv.c in GNU Scientific Library
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2004, 2007 Brian Gough
// Licensed under the GNU General Public License, version 3 or later.
#ifndef REDUKTI_MATHLIB_DERIVATIVES_H
#define REDUKTI_MATHLIB_DERIVATIVES_H

#include <functional>

namespace redukti::mathlib {

/**
 * Java's DerivFunction is a functional interface with no named implementations
 * anywhere in the codebase, so it maps to std::function rather than to an
 * abstract base the way ScalarObjectiveFunction does.
 */
using DerivFunction = std::function<double(double)>;

/** C++ port of org.redukti.mathlib.DerivResult */
class DerivResult {
public:
    double result;
    double abserr;

    DerivResult(double result_, double abserr_) : result(result_), abserr(abserr_) {}
};

/**
 * NOTE: nothing outside mathlib calls central_derivative. Ported for
 * completeness; the optimizer builds its Jacobians in optim instead.
 */
class Derivatives {
public:
    static const double GSL_DBL_EPSILON;

    struct EvalResult {
        double result;
        double abserr_round;
        double abserr_trunc;
    };

    static EvalResult central_deriv(const DerivFunction &f, double x, double h);

    static DerivResult central_derivative(const DerivFunction &f, double x, double h);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_DERIVATIVES_H
