// C++ port of org.redukti.optim.SolverTolerances
#ifndef REDUKTI_OPTIM_SOLVERTOLERANCES_H
#define REDUKTI_OPTIM_SOLVERTOLERANCES_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/MinPack.h"

#include <cmath>
#include <string>

namespace redukti::optim {

/**
 * The lmder stopping tolerances, and the defaults a trial starts from. A trial's
 * 'solver' settings replace them one at a time; everything else keeps the value
 * here, so these four fields are the only place a default is stated.
 */
class SolverTolerances {
public:
    /** Don't stop on how far the variables moved; ray-trace noise makes late steps tiny. */
    static constexpr double DEFAULT_XTOL = 0.0;
    /** maxEvaluations placeholder for 100 * (variables + 1), which needs the variable count. */
    static constexpr int FROM_VARIABLE_COUNT = 0;

    /** Relative reduction in the sum of squares: the square root of the machine epsilon. */
    static double defaultFtol() { return std::sqrt(mathlib::MinPack::dpmpar(1)); }
    /** Stop when the gradient is genuinely flat. */
    static double defaultGtol() { return std::sqrt(mathlib::MinPack::dpmpar(1)); }

    SolverTolerances()
        : ftol_(defaultFtol()), xtol_(DEFAULT_XTOL), gtol_(defaultGtol()),
          maxEvaluations_(FROM_VARIABLE_COUNT) {}

    SolverTolerances(double ftol, double xtol, double gtol, int maxEvaluations)
        : ftol_(ftol), xtol_(xtol), gtol_(gtol), maxEvaluations_(maxEvaluations) {
        check(ftol_, "ftol");
        check(xtol_, "xtol");
        check(gtol_, "gtol");
        if (maxEvaluations_ < 0)
            throw IllegalArgumentException("max-evaluations must be positive");
    }

    /** relative reduction in the sum of squares that stops the solve */
    double ftol() const { return ftol_; }
    /**
     * Relative change in the varied values themselves - the curvatures, thicknesses
     * and aspheric terms - between two iterations, below which the solve stops.
     * 0 disables the test.
     */
    double xtol() const { return xtol_; }
    /** gradient flatness that stops the solve */
    double gtol() const { return gtol_; }
    /** the lmder maxfev as set, which may be FROM_VARIABLE_COUNT */
    int maxEvaluations() const { return maxEvaluations_; }

    /** The trial step limit for a solve over `variables` variables. */
    int maxEvaluations(int variables) const {
        return maxEvaluations_ == FROM_VARIABLE_COUNT ? (variables + 1) * 100 : maxEvaluations_;
    }

private:
    static void check(double value, const char *what) {
        if (!std::isfinite(value) || value < 0.0)
            throw IllegalArgumentException(std::string(what) +
                                           " must be finite and non-negative");
    }

    double ftol_;
    double xtol_;
    double gtol_;
    int maxEvaluations_;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_SOLVERTOLERANCES_H
