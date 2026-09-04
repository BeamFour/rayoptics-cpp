// C++ port of org.redukti.optim.Solver, LMDerMeritFunction and LMDerSolver.
#ifndef REDUKTI_OPTIM_LMDER_H
#define REDUKTI_OPTIM_LMDER_H

#include "redukti/mathlib/MinPack.h"
#include "redukti/optim/Goal.h"
#include "redukti/optim/Var.h"

#include <memory>
#include <string>
#include <vector>

namespace redukti::optim {

class Solver {
public:
    virtual ~Solver() = default;
    virtual int solve() = 0;
};

class LMDerMeritFunction : public mathlib::Lmder_Function {
public:
    /**
     * The vars and goals are shared with the OptimizationSetup that built them
     * and with the solver this hands out, which is what the Java's shared array
     * references amount to.
     */
    LMDerMeritFunction(Analysis *analysis, std::vector<std::shared_ptr<Var>> vars,
                       std::vector<std::shared_ptr<Goal>> functions, bool use_native);
    LMDerMeritFunction(std::shared_ptr<Analysis> analysis,
                       std::vector<std::shared_ptr<Var>> vars,
                       std::vector<std::shared_ptr<Goal>> functions, bool use_native);

    bool hasJacobian() override { return true; }

    using mathlib::Lmder_Function::apply;

    int apply(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
              std::vector<double> &fjac, int ldfjac, int iflag) override;

    int apply(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
              int iflag) override;

    bool buildJacobian(std::vector<double> &x, std::vector<double> &fjac, int ldfjac);

    /** Retained for callers that need to evaluate a complete perturbed goal vector. */
    bool nudge(std::vector<double> &x, std::vector<double> &delta,
               std::vector<double> &resid) {
        return evaluate(x, delta, resid);
    }

    std::unique_ptr<Solver> getSolver();

    std::string toString();

    const std::vector<std::shared_ptr<Goal>> &goals() const { return functions; }
    const std::vector<std::shared_ptr<Var>> &variables() const { return vars; }

    double getRMS();

private:
    static const double BIGVAL;
    static constexpr int MAX_JACOBIAN_STEP_REDUCTIONS = 8;

    std::vector<double> weights;
    std::shared_ptr<Analysis> analysis_owner;
    Analysis *analysis;
    std::vector<std::shared_ptr<Var>> vars;
    std::vector<std::shared_ptr<Goal>> functions;
    bool use_native;

    /**
     * Evaluates the weighted residuals (value - target)*weight at x.
     * On any failure (killed ray, NaN in x, exception in the analysis)
     * all residuals are set to BIGVAL so that lmder rejects the trial step.
     * lmder minimizes ||fvec||^2, so fvec must be the deviation from target:
     * using the raw goal value would drive EFL/Fno/MTF towards zero.
     */
    void computeResiduals(std::vector<double> &x, std::vector<double> &fvec);

    /** Evaluate raw goal values at x + delta; invalid goals are stored as NaN. */
    bool evaluate(std::vector<double> &x, std::vector<double> &delta,
                  std::vector<double> &values);

    static bool everyResidualHasPerturbedValue(const std::vector<double> &forward,
                                               const std::vector<double> &backward);
    static bool isUsable(double value);

    void validateInitialContrastSamples();
    void validateInputs();
};

class LMDerSolver : public Solver {
public:
    LMDerSolver(Analysis *analysis, std::vector<std::shared_ptr<Var>> vars,
                std::vector<std::shared_ptr<Goal>> functions, bool use_native)
        : analysis(analysis), vars(std::move(vars)), functions(std::move(functions)),
          use_native(use_native) {}
    LMDerSolver(std::shared_ptr<Analysis> analysis_, std::vector<std::shared_ptr<Var>> vars,
                std::vector<std::shared_ptr<Goal>> functions, bool use_native)
        : analysis_owner(std::move(analysis_)), analysis(analysis_owner.get()),
          vars(std::move(vars)), functions(std::move(functions)), use_native(use_native) {}

    int solve() override;

private:
    std::shared_ptr<Analysis> analysis_owner;
    Analysis *analysis;
    /** number of vars in lmder parlance */
    std::vector<std::shared_ptr<Var>> vars;
    /** number of functions in lmder parlance */
    std::vector<std::shared_ptr<Goal>> functions;
    bool use_native = false;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_LMDER_H
