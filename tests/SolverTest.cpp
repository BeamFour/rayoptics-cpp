// Expectations captured by running the Java classes in rayoptics/target/classes
// on JDK 25. As in MatrixTest, values are asserted exactly -- see the note
// there.
#include "TestHarness.h"

#include "redukti/mathlib/BrentSolver.h"
#include "redukti/mathlib/SecantSolver.h"

#include <cmath>
#include <optional>
#include <vector>

using namespace redukti::mathlib;

namespace {

// f(x) = x^3 - 2x - 5, the classic Brent test; root near 2.0945514815423265.
struct Cubic : ScalarObjectiveFunction {
    int calls = 0;
    std::optional<double> eval(double x) override {
        calls++;
        return x * x * x - 2.0 * x - 5.0;
    }
};

struct Cosine : ScalarObjectiveFunction {
    std::optional<double> eval(double x) override { return std::cos(x) - x; }
};

struct Quadratic : ScalarObjectiveFunction {
    std::optional<double> eval(double x) override { return x * x - 4.0; }
};

struct NoRoot : ScalarObjectiveFunction {
    std::optional<double> eval(double x) override { return x * x + 1.0; }
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
