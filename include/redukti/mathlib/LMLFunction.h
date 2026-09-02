// C++ port of org.redukti.mathlib.LMLFunction
//
// @author: M.Lampton (c) 2005 Stellar Software
// Original License: GPL v2
#ifndef REDUKTI_MATHLIB_LMLFUNCTION_H
#define REDUKTI_MATHLIB_LMLFUNCTION_H

#include <vector>

namespace redukti::mathlib {

class LMLFunction {
public:
    virtual ~LMLFunction() = default;

    /** Returns sos, or BIGVAL if parms failed. */
    virtual double computeResiduals() = 0;

    /** Allows LM to request a new Jacobian. false if parms failed. */
    virtual bool buildJacobian() = 0;

    /** Get residual at i */
    virtual double getResidual(int i) = 0;

    /** Get Jacobian at i,j */
    virtual double getJacobian(int i, int j) = 0;

    /**
     * Allows LM to modify parms[] and reevaluate its fit.
     * Returns sum-of-squares for nudged params.
     * This is the only place that parms[] are modified.
     * Moves parms, builds resid[], returns sos or BIGVAL
     */
    virtual double nudge(const std::vector<double> &delta) = 0;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_LMLFUNCTION_H
