// Expectations captured by running the Java classes in rayoptics/target/classes
// on JDK 25. As in MatrixTest, values are asserted exactly -- see the note
// there. The LMLSolver fixture is the LMDemo from TestLMLSolver.java.
#include "TestHarness.h"

#include "redukti/mathlib/BrentSolver.h"
#include "redukti/mathlib/Derivatives.h"
#include "redukti/mathlib/LMLSolver.h"
#include "redukti/mathlib/SecantSolver.h"

#include <cmath>
#include <vector>

using namespace redukti::mathlib;

namespace {

// f(x) = x^3 - 2x - 5, the classic Brent test; root near 2.0945514815423265.
struct Cubic : ScalarObjectiveFunction {
    int calls = 0;
    double eval(double x) override {
        calls++;
        return x * x * x - 2.0 * x - 5.0;
    }
};

struct Cosine : ScalarObjectiveFunction {
    double eval(double x) override { return std::cos(x) - x; }
};

struct Quadratic : ScalarObjectiveFunction {
    double eval(double x) override { return x * x - 4.0; }
};

struct NoRoot : ScalarObjectiveFunction {
    double eval(double x) override { return x * x + 1.0; }
};

} // namespace

TEST(brent_finds_roots) {
    Cubic c;
    RootResult r = BrentSolver::find_root(1.0, 3.0, c);
    CHECK_STR_EQ(r.toString(),
                 "RootResult{root=2.0945514815423265, converged=true, iterations=10}");
    CHECK_EQ(c.calls, 11);

    Cosine cosf;
    CHECK_STR_EQ(BrentSolver::find_root(0.0, 2.0, cosf).toString(),
                 "RootResult{root=0.7390851332151607, converged=true, iterations=7}");

    Quadratic q;
    CHECK_STR_EQ(BrentSolver::find_root(0.0, 5.0, q).toString(),
                 "RootResult{root=1.999999999999977, converged=true, iterations=10}");
}

TEST(brent_exact_root_at_bracket_endpoint) {
    Quadratic q;
    // fa == 0 and fb == 0 both short-circuit before any iteration.
    CHECK_STR_EQ(BrentSolver::find_root(2.0, 5.0, q).toString(),
                 "RootResult{root=2.0, converged=true, iterations=0}");
    CHECK_STR_EQ(BrentSolver::find_root(0.0, 2.0, q).toString(),
                 "RootResult{root=2.0, converged=true, iterations=0}");
}

TEST(brent_rejects_bad_bracket) {
    NoRoot n;
    // Same sign at both ends: reported as root 0.0, not converged.
    CHECK_STR_EQ(BrentSolver::find_root(0.0, 1.0, n).toString(),
                 "RootResult{root=0.0, converged=false, iterations=0}");
}

TEST(secant_finds_roots) {
    Cubic c;
    CHECK_STR_EQ(SecantSolver::find_root(c, 2.0, 100, 1e-12).toString(),
                 "RootResult{root=2.0945514815423265, converged=true, iterations=5}");
    Cosine cosf;
    CHECK_STR_EQ(SecantSolver::find_root(cosf, 0.5, 100, 1e-12).toString(),
                 "RootResult{root=0.7390851332151606, converged=true, iterations=5}");
    Quadratic q;
    CHECK_STR_EQ(SecantSolver::find_root(q, 1.0, 100, 1e-12).toString(),
                 "RootResult{root=2.0, converged=true, iterations=7}");
    // x0 == 0 makes the initial bracket depend on the eps nudge alone.
    CHECK_STR_EQ(SecantSolver::find_root(q, 0.0, 100, 1e-12).toString(),
                 "RootResult{root=2.0, converged=true, iterations=28}");
}

TEST(secant_gives_up_at_maxiter) {
    NoRoot n;
    CHECK_STR_EQ(SecantSolver::find_root(n, 1.0, 5, 1e-12).toString(),
                 "RootResult{root=1.0016011211104183, converged=false, iterations=5}");
}

TEST(central_derivative) {
    DerivResult d1 = Derivatives::central_derivative([](double x) { return std::sin(x); },
                                                     1.0, 1e-4);
    CHECK_CLOSE(d1.result, 0.5403023058670838, 0.0);
    CHECK_CLOSE(d1.abserr, 9.636876518932678e-11, 0.0);

    DerivResult d2 = Derivatives::central_derivative([](double x) { return x * x * x; },
                                                     2.0, 1e-3);
    CHECK_CLOSE(d2.result, 11.999999999971859, 0.0);
    CHECK_CLOSE(d2.abserr, 1.2437794525465432e-9, 0.0);

    DerivResult d3 = Derivatives::central_derivative([](double x) { return std::exp(x); },
                                                     0.5, 1e-2);
    CHECK_CLOSE(d3.result, 1.6487212707011694, 0.0);
    CHECK_CLOSE(d3.abserr, 2.186528896904694e-10, 0.0);

    DerivResult d4 = Derivatives::central_derivative(
        [](double x) { return 3.0 * x + 1.0; }, 7.0, 1e-3);
    CHECK_CLOSE(d4.result, 3.0000000000001132, 0.0);
    CHECK_CLOSE(d4.abserr, 3.397325823439896e-11, 0.0);
}

// ---------------------------------------------------------------------------
// LMDemo, ported from TestLMLSolver.java: fits an inverse even polynomial
// p0 / (1 + p1*x^2 + p2*x^4 + p3*x^6) to 67 sampled points.
// ---------------------------------------------------------------------------
namespace {

class LMDemo : public LMLFunction {
public:
    static constexpr double DELTAP = 1e-6; // parm step
    static constexpr double BIGVAL = 9e99; // fault flag

    // {x, y} data pairs [row][col]
    std::vector<std::vector<double>> data = {
        {0.00, 0.6793}, {0.03, 0.6787}, {0.06, 0.6768}, {0.09, 0.6736}, {0.12, 0.6691},
        {0.15, 0.6634}, {0.18, 0.6565}, {0.21, 0.6482}, {0.24, 0.6388}, {0.27, 0.6280},
        {0.30, 0.6161}, {0.33, 0.6030}, {0.36, 0.5887}, {0.39, 0.5733}, {0.42, 0.5568},
        {0.45, 0.5394}, {0.48, 0.5210}, {0.51, 0.5019}, {0.54, 0.4820}, {0.57, 0.4614},
        {0.60, 0.4404}, {0.63, 0.4191}, {0.66, 0.3975}, {0.69, 0.3758}, {0.72, 0.3542},
        {0.75, 0.3328}, {0.78, 0.3117}, {0.81, 0.2910}, {0.84, 0.2710}, {0.87, 0.2515},
        {0.90, 0.2328}, {0.93, 0.2149}, {0.96, 0.1979}, {0.99, 0.1818}, {1.02, 0.1667},
        {1.05, 0.1524}, {1.08, 0.1392}, {1.11, 0.1268}, {1.14, 0.1154}, {1.17, 0.1048},
        {1.20, 0.0951}, {1.23, 0.0862}, {1.26, 0.0781}, {1.29, 0.0706}, {1.32, 0.0639},
        {1.35, 0.0577}, {1.38, 0.0521}, {1.41, 0.0471}, {1.44, 0.0425}, {1.47, 0.0384},
        {1.50, 0.0347}, {1.53, 0.0313}, {1.56, 0.0283}, {1.59, 0.0255}, {1.62, 0.0231},
        {1.65, 0.0209}, {1.68, 0.0189}, {1.71, 0.0171}, {1.74, 0.0155}, {1.77, 0.0140},
        {1.80, 0.0127}, {1.83, 0.0115}, {1.86, 0.0105}, {1.89, 0.0095}, {1.92, 0.0087},
        {1.95, 0.0079}, {1.98, 0.0072}};

    int NPTS = static_cast<int>(data.size());
    double WEIGHT = 1.0;
    std::vector<double> parms = {1.0, 1.0, 1.0, 1.0};
    int NPARMS = static_cast<int>(parms.size());
    std::vector<double> resid = std::vector<double>(static_cast<std::size_t>(NPTS), 0.0);
    std::vector<std::vector<double>> jac = std::vector<std::vector<double>>(
        static_cast<std::size_t>(NPTS),
        std::vector<double>(static_cast<std::size_t>(NPARMS), 0.0));

    // inverse even polynomial. Called only by computeResiduals().
    double func(int i, const std::vector<std::vector<double>> &d,
                const std::vector<double> &p) {
        double x = d[i][0];
        double x2 = x * x;
        double x4 = x2 * x2;
        double x6 = x4 * x2;
        double denom = 1 + p[1] * x2 + p[2] * x4 + p[3] * x6;
        return p[0] / denom;
    }

    double computeResiduals() override {
        double sumsq = 0.0;
        for (int i = 0; i < NPTS; i++) {
            double y = data[i][1]; // row i, col 1
            resid[i] = (func(i, data, parms) - y) * WEIGHT;
            sumsq += resid[i] * resid[i];
        }
        return sumsq;
    }

    bool buildJacobian() override {
        std::vector<double> delta(static_cast<std::size_t>(NPARMS), 0.0);
        double FACTOR = 0.5 / DELTAP;
        double d = 0;

        for (int j = 0; j < NPARMS; j++) {
            for (int k = 0; k < NPARMS; k++)
                delta[k] = (k == j) ? DELTAP : 0.0;

            d = nudge(delta); // resid at pplus
            if (d == BIGVAL)
                return false;
            for (int i = 0; i < NPTS; i++)
                jac[i][j] = getResidual(i);

            for (int k = 0; k < NPARMS; k++)
                delta[k] = (k == j) ? -2 * DELTAP : 0.0;

            d = nudge(delta); // resid at pminus
            if (d == BIGVAL)
                return false;

            for (int i = 0; i < NPTS; i++)
                jac[i][j] -= getResidual(i); // fetches resid[]

            for (int i = 0; i < NPTS; i++)
                jac[i][j] *= FACTOR;

            for (int k = 0; k < NPARMS; k++)
                delta[k] = (k == j) ? DELTAP : 0.0;

            d = nudge(delta);
            if (d == BIGVAL)
                return false;
        }
        return true;
    }

    double getResidual(int i) override { return resid[i]; }

    double getJacobian(int i, int j) override { return jac[i][j]; }

    double nudge(const std::vector<double> &delta) override {
        for (int j = 0; j < NPARMS; j++)
            parms[j] += delta[j];
        return computeResiduals();
    }
};

int solve(LMLSolver &solver) {
    int status;
    do {
        status = solver.iLMiter();
    } while (status == LMLSolver::DOWNITER);
    return status;
}

} // namespace

TEST(lmlsolver_baseline) {
    LMDemo function;
    LMLSolver solver(function, 1E-12, function.NPARMS, function.NPTS);
    int status = solve(solver);
    CHECK_EQ(status, LMLSolver::LEVELITER);
    // TestLMLSolver.testBaseline only asserts to 1e-6; these are the full
    // values the Java solver actually lands on, so the whole iteration
    // sequence has to match, not just the converged neighbourhood.
    CHECK_CLOSE(function.parms[0], 0.6804239765015776, 0.0);
    CHECK_CLOSE(function.parms[1], 1.109957043602709, 0.0);
    CHECK_CLOSE(function.parms[2], 0.7185315615803544, 0.0);
    CHECK_CLOSE(function.parms[3], 1.0486691656282272, 0.0);
}
