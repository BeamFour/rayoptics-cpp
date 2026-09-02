// C++ port of org.redukti.mathlib.LMLSolver
//
// Levenberg Marquardt Lampton solver
// From BeamFour project - AutoRay.
// M.Lampton, 1997 Computers In Physics v.11 #10 110-115.
//
// @author: M.Lampton (c) 2005 Stellar Software
// Original License: GPL v2
#ifndef REDUKTI_MATHLIB_LMLSOLVER_H
#define REDUKTI_MATHLIB_LMLSOLVER_H

#include "redukti/mathlib/LMLFunction.h"

#include <vector>

namespace redukti::mathlib {

/**
 * Constructor is used to set up all parms including host for callback.
 * Sole public method is iLMiter() which performs one iteration.
 * Arrays parms[], resid[], jac[][] are unknown here; instead the host must
 * fake these and provide results through the LMLFunction callbacks.
 *
 * Exit leaves host with parms[] optimized through its sequence of nudges.
 */
class LMLSolver {
public:
    static const double BIGVAL;
    static const int DOWNITER = 0;  // iteration ok, want more
    static const int LEVELITER = 1; // iteration ok, all done.
    static const int MAXITER = 2;   // did enough iterations.
    static const int BADITER = 3;   // ray killed bail out

    /**
     * Builds an instance of the solver.
     *
     * @param gH User defined function
     * @param gtol Tolerance
     * @param gnparms # of parameters to adjust
     * @param gnpts # of points to fit to
     */
    LMLSolver(LMLFunction &gH, double gtol, int gnparms, int gnpts);

    /**
     * Called repeatedly by LMhost to perform each LM iteration.
     * Returns BADITER to shut down ray failed;
     * Returns DOWNITER if iteration went OK, more needed;
     * Returns LEVELITER if iteration went OK, all done.
     * Ref: M.Lampton, Computers in Physics v.11 pp.110-115 1997.
     */
    int iLMiter();

private:
    /** inverts the array a[N][N] by Gauss-Jordan method */
    double gaussj(std::vector<std::vector<double>> &a, int N);

    const double LMBOOST = 2.0;     // damping increase per bad step
    const double LMSHRINK = 0.10;   // damping decrease per good step
    const double LAMBDAZERO = 0.001;// initial damping
    const double LAMBDAMAX = 1E3;   // max damping

    int niter = 0;                  // local diagnostic only
    double sos = 0, sosinit = 0, lambda = 0; // local diagnostic only

    LMLFunction *myH = nullptr;     // overwritten by constructor
    double lmtol = 1E-6;            // overwritten by constructor
    int lmiter = 200;               // overwritten by constructor
    int nparms = 0;                 // overwritten by constructor
    int npts = 0;                   // overwritten by constructor

    std::vector<double> delta;                  // local
    std::vector<double> beta;                   // local
    std::vector<std::vector<double>> alpha;     // local
    std::vector<std::vector<double>> amatrix;   // local
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_LMLSOLVER_H
