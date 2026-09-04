// C++ port of org.redukti.optim.LMDerMeritFunction and LMDerSolver
#include "redukti/optim/LMDer.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/LMLSolver.h"
#include "redukti/mathlib/M.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

namespace redukti::optim {

const double LMDerMeritFunction::BIGVAL = mathlib::LMLSolver::BIGVAL;

LMDerMeritFunction::LMDerMeritFunction(Analysis *analysis,
                                       std::vector<std::shared_ptr<Var>> vars_,
                                       std::vector<std::shared_ptr<Goal>> functions_,
                                       bool use_native)
    : analysis(analysis), vars(std::move(vars_)), functions(std::move(functions_)),
      use_native(use_native) {
    weights.resize(functions.size());
    for (std::size_t i = 0; i < functions.size(); i++) {
        if (!std::isfinite(functions[i]->_weight) || functions[i]->_weight < 0.0)
            throw IllegalArgumentException("Goal weight must be finite and non-negative");
        weights[i] = std::sqrt(functions[i]->_weight);
    }
    for (std::size_t i = 0; i < vars.size(); i++)
        vars[i]->read_from_prescription();
}

LMDerMeritFunction::LMDerMeritFunction(std::shared_ptr<Analysis> analysis_,
                                       std::vector<std::shared_ptr<Var>> vars_,
                                       std::vector<std::shared_ptr<Goal>> functions_,
                                       bool use_native)
    : analysis_owner(std::move(analysis_)), analysis(analysis_owner.get()),
      vars(std::move(vars_)), functions(std::move(functions_)), use_native(use_native) {
    weights.resize(functions.size());
    for (std::size_t i = 0; i < functions.size(); i++) {
        if (!std::isfinite(functions[i]->_weight) || functions[i]->_weight < 0.0)
            throw IllegalArgumentException("Goal weight must be finite and non-negative");
        weights[i] = std::sqrt(functions[i]->_weight);
    }
    for (const auto &var : vars)
        var->read_from_prescription();
}

int LMDerMeritFunction::apply(int m, int n, std::vector<double> &x,
                              std::vector<double> &fvec, std::vector<double> &fjac,
                              int ldfjac, int iflag) {
    // m should be size of outs
    // n should be size of vars
    // x is current guess for vars
    // fvec is the result of outs
    // fjac is jacobian

    assert(m == static_cast<int>(functions.size()));
    assert(n == static_cast<int>(vars.size()));
    (void)m;
    (void)n;

    // if the nprint parameter to lmder is positive, the function is
    // called every nprint iterations with iflag=0, so that the
    // function may perform special operations, such as printing
    // residuals.
    if (iflag == 0)
        return 0;
    if (iflag != 2) {
        computeResiduals(x, fvec);
    } else {
        // compute jacobian
        if (!buildJacobian(x, fjac, ldfjac))
            return -99;
    }
    return 0;
}

int LMDerMeritFunction::apply(int m, int n, std::vector<double> &x,
                              std::vector<double> &fvec, int iflag) {
    assert(m == static_cast<int>(functions.size()));
    assert(n == static_cast<int>(vars.size()));
    (void)m;
    (void)n;

    if (iflag == 0)
        return 0;
    computeResiduals(x, fvec);
    return 0;
}

void LMDerMeritFunction::computeResiduals(std::vector<double> &x,
                                          std::vector<double> &fvec) {
    bool okay = true;
    try {
        for (std::size_t i = 0; i < x.size(); i++) {
            vars[i]->set_scaled_value(x[i]);
            vars[i]->write_to_prescription();
        }
        analysis->compute();
    } catch (const Exception &) {
        okay = false;
    }
    for (std::size_t i = 0; i < functions.size(); i++) {
        double value = okay ? functions[i]->value() : BIGVAL;
        double r;
        if (!std::isfinite(value) || value >= BIGVAL) {
            r = BIGVAL;
        } else {
            r = (value - functions[i]->_target) * weights[i];
        }
        fvec[i] = std::isfinite(r) ? r : BIGVAL;
    }
}

bool LMDerMeritFunction::buildJacobian(std::vector<double> &x, std::vector<double> &fjac,
                                       int ldfjac) {
    const int n = static_cast<int>(vars.size());
    const int m = static_cast<int>(functions.size());
    std::vector<double> base(static_cast<std::size_t>(m), 0.0);
    std::vector<double> forward(static_cast<std::size_t>(m), 0.0);
    std::vector<double> backward(static_cast<std::size_t>(m), 0.0);
    std::vector<double> restored(static_cast<std::size_t>(m), 0.0);
    std::vector<double> delta(static_cast<std::size_t>(n), 0.0);
    bool success = false;

    // The Java body sits in a try/finally; the restoration below is the finally.
    // Nothing here throws -- evaluate() swallows everything -- so a plain
    // sequence reproduces it without needing a scope guard.
    {
        // A derivative is meaningful only around a fully valid accepted point.
        // This also re-establishes x if the preceding LM trial was rejected.
        success = evaluate(x, delta, base);
        for (int j = 0; success && j < n; j++) {
            // Fixed absolute step per variable (scaled units). A relative step
            // degenerates at zero-valued parameters such as fresh aspherics.
            double step = vars[static_cast<std::size_t>(j)]->_d_delta;
            if (!std::isfinite(step) || step <= 0.0) {
                success = false;
                break;
            }

            bool usable = false;
            for (int attempt = 0; attempt <= MAX_JACOBIAN_STEP_REDUCTIONS; attempt++) {
                delta[static_cast<std::size_t>(j)] = step;
                bool forwardComplete = evaluate(x, delta, forward);
                delta[static_cast<std::size_t>(j)] = -step;
                bool backwardComplete = evaluate(x, delta, backward);
                delta[static_cast<std::size_t>(j)] = 0.0;

                usable = everyResidualHasPerturbedValue(forward, backward);
                if (forwardComplete && backwardComplete)
                    break;
                if (attempt < MAX_JACOBIAN_STEP_REDUCTIONS)
                    step *= 0.5;
            }
            if (!usable) {
                success = false;
                break;
            }

            for (int i = 0; i < m; i++) {
                auto ui = static_cast<std::size_t>(i);
                double derivative;
                if (isUsable(forward[ui]) && isUsable(backward[ui]))
                    derivative = (forward[ui] - backward[ui]) / (2.0 * step);
                else if (isUsable(forward[ui]))
                    derivative = (forward[ui] - base[ui]) / step;
                else
                    derivative = (base[ui] - backward[ui]) / step;
                derivative *= weights[ui];
                if (!std::isfinite(derivative)) {
                    success = false;
                    break;
                }
                fjac[static_cast<std::size_t>(i + j * ldfjac)] = derivative;
            }
        }
    }
    // A failed probe must never leak its prescription or partial analysis
    // state to the caller. Restoration failure invalidates the Jacobian.
    std::fill(delta.begin(), delta.end(), 0.0);
    // The Java writes `success &= evaluate(...)`, a non-short-circuiting
    // bitwise AND on booleans: the restore runs even when the Jacobian has
    // already failed, which is the whole point of doing it in a finally.
    // `success && evaluate(...)` would skip it and leave the prescription
    // and analysis describing the last probe point.
    const bool restoredOk = evaluate(x, delta, restored);
    success = success && restoredOk;
    return success;
}

bool LMDerMeritFunction::everyResidualHasPerturbedValue(
    const std::vector<double> &forward, const std::vector<double> &backward) {
    for (std::size_t i = 0; i < forward.size(); i++)
        if (!isUsable(forward[i]) && !isUsable(backward[i]))
            return false;
    return true;
}

bool LMDerMeritFunction::isUsable(double value) {
    return std::isfinite(value) && value < BIGVAL;
}

bool LMDerMeritFunction::evaluate(std::vector<double> &x, std::vector<double> &delta,
                                  std::vector<double> &values) {
    try {
        for (std::size_t i = 0; i < delta.size(); i++) {
            vars[i]->set_scaled_value(x[i] + delta[i]);
            vars[i]->write_to_prescription();
        }
        analysis->compute();
    } catch (const Exception &) {
        std::fill(values.begin(), values.end(),
                  std::numeric_limits<double>::quiet_NaN());
        return false;
    }
    bool complete = true;
    for (std::size_t i = 0; i < functions.size(); i++) {
        double value = functions[i]->value();
        if (isUsable(value)) {
            values[i] = value;
        } else {
            values[i] = std::numeric_limits<double>::quiet_NaN();
            complete = false;
        }
    }
    return complete;
}

void LMDerMeritFunction::validateInitialContrastSamples() {
    int invalidCount = 0;
    bool haveFirst = false;
    std::string first;
    for (const auto &contrast : analysis->_contrasts) {
        for (const auto &field : contrast.fields) {
            for (const auto &wavelength : field.wavelengths) {
                for (const auto &sample : wavelength.samples) {
                    if (sample.valid)
                        continue;
                    invalidCount++;
                    if (!haveFirst) {
                        haveFirst = true;
                        const auto &failure = *sample.failure;
                        first = "first failure: " + failure.ray + " ray encountered " +
                                failure.exceptionType + " at surface " +
                                intToString(failure.surface) +
                                ", field=" + doubleToString(field.field->y) +
                                ", wavelength=" + doubleToString(wavelength.wavelength) +
                                " nm" + ", frequency=" +
                                doubleToString(contrast.spatialFrequency) +
                                " cycles/mm" + ", pupil=" + sample.pupil.toString();
                    }
                }
            }
        }
    }
    if (invalidCount > 0)
        throw IllegalStateException("Cannot start optimization: " +
                                    intToString(invalidCount) +
                                    " contrast samples contain failed rays; " + first);
}

void LMDerMeritFunction::validateInputs() {
    analysis->compute();
    validateInitialContrastSamples();
}

std::unique_ptr<Solver> LMDerMeritFunction::getSolver() {
    validateInputs();
    if (analysis_owner)
        return std::make_unique<LMDerSolver>(analysis_owner, vars, functions, use_native);
    return std::make_unique<LMDerSolver>(analysis, vars, functions, use_native);
}

std::string LMDerMeritFunction::toString() {
    std::string sb;
    sb += "Vars:\n";
    for (std::size_t i = 0; i < vars.size(); i++)
        sb += vars[i]->toString() + "\n";
    sb += "Values:\n";
    for (std::size_t i = 0; i < functions.size(); i++)
        sb += functions[i]->toString() + "\n";
    sb += "RMS: " + doubleToString(getRMS()) + "\n";
    return sb;
}

double LMDerMeritFunction::getRMS() {
    double sos = 0.0;
    std::vector<double> resid(functions.size(), 0.0);
    for (std::size_t i = 0; i < functions.size(); i++) {
        resid[i] = (functions[i]->_target - functions[i]->value()) * weights[i];
        sos += mathlib::M::square(resid[i]);
    }
    if (mathlib::M::isZero(sos))
        return sos;
    return std::sqrt(sos / static_cast<double>(functions.size()));
}

// ---------------------------------------------------------------------------
// LMDerSolver
// ---------------------------------------------------------------------------

int LMDerSolver::solve() {
    LMDerMeritFunction fcn = analysis_owner
                                 ? LMDerMeritFunction(analysis_owner, vars, functions, use_native)
                                 : LMDerMeritFunction(analysis, vars, functions, use_native);
    int m = static_cast<int>(functions.size()); // number of functions
    int n = static_cast<int>(vars.size());      // number of variables, must not exceed m
    if (m < n)
        throw IllegalArgumentException("Number of goals must be >= number of variables");

    // Setup initial solution vector - this is just read from the values in
    // the prescription
    std::vector<double> x(static_cast<std::size_t>(n), 0.0);
    for (std::size_t i = 0; i < x.size(); i++) {
        // x[] will be scaled values
        x[i] = vars[i]->read_from_prescription();
    }
    std::vector<double> diag(static_cast<std::size_t>(n), 1.0);

    int info = 0;
    if (use_native) {
        // The Java's native MinpackFFM path is commented out there too.
    } else {
        std::vector<double> fvec(static_cast<std::size_t>(m), 0.0); // Results of goals
        std::vector<double> fjac(static_cast<std::size_t>(m) * static_cast<std::size_t>(n),
                                 0.0); // Space for jacobian
        int ldfjac = m;
        double ftol = std::sqrt(mathlib::MinPack::dpmpar(1));
        double xtol = 0.;      // don't stop on step size; ray-trace noise makes late steps tiny
        double gtol = 1.0e-12; // stop when the gradient is genuinely flat
        int maxfev = (n + 1) * 100;
        int mode = 1; // 1=scale internally 2=scale using diag
        double factor = 100;
        int nprint = 1;

        std::vector<int> nfev(1, 0), njev(1, 0);
        std::vector<int> ipvt(static_cast<std::size_t>(n), 0);
        double epsfcn = 0.0;

        std::vector<double> qtf(static_cast<std::size_t>(n), 0.0);
        std::vector<double> wa1(static_cast<std::size_t>(n), 0.0);
        std::vector<double> wa2(static_cast<std::size_t>(n), 0.0);
        std::vector<double> wa3(static_cast<std::size_t>(n), 0.0);
        std::vector<double> wa4(static_cast<std::size_t>(m), 0.0);
        info = mathlib::MinPack::lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol,
                                       maxfev, diag, mode, factor, nprint, nfev, njev,
                                       ipvt, qtf, wa1, wa2, wa3, wa4, epsfcn);
        std::cout << "lmder: info=" << info << " nfev=" << nfev[0] << " njev=" << njev[0]
                  << " (each njev costs 2n=" << (2 * n)
                  << " evaluations via central differences)" << std::endl;
    }

    // Set final solution vector
    for (std::size_t i = 0; i < x.size(); i++) {
        vars[i]->set_scaled_value(x[i]);
        // Update prescription
        vars[i]->write_to_prescription();
    }
    // A failed Jacobian or rejected trial can leave Analysis describing the
    // last perturbed point, including partially populated contrast results.
    // Recompute at lmder's accepted x so reporting is coherent even for info < 0.
    analysis->compute();

    return info;
}

} // namespace redukti::optim
