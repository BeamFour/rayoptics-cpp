// C++ port of org.redukti.mathlib.BrentSolver
#include "redukti/mathlib/BrentSolver.h"

#include <cmath>

namespace redukti::mathlib {

const double BrentSolver::TOL = 5e-14;

RootResult BrentSolver::find_root(double a, double b, ScalarObjectiveFunction &fn)
//  Given a bracket (t[0], t[1]), Brent() sets t[0] to root,
//  and returns the number of calls to zDiff() taken.
//  Press et al NUMERICAL RECIPES IN C 2nd edition 1992 p.361
//  R.P.Brent ALGORITHMS... Prentice-Hall, NJ 1973.
{
    int MAXIT = 50;
    double c, d = 0, e = 0, min1, min2;
    double fc, p, q, r, s, toler, xm, fa, fb;
    c = b;
    fa = fn.eval(a).value();
    if (fa == 0.)
        return RootResult(a, true, 0);
    fb = fn.eval(b).value();
    if (fb == 0.)
        return RootResult(b, true, 0);
    if (fa * fb > 0)
        // bad starting bracket
        return RootResult(0., false, 0);
    fc = fb;
    for (int iter = 1; iter < MAXIT; iter++) {
        if (fb * fc > 0) {
            c = a;
            fc = fa;
            e = d = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        toler = 2.0 * TOL * std::abs(b) + TOL;
        xm = 0.5 * (c - b);
        if ((std::abs(xm) <= toler) || (fb == 0.0))
            return RootResult(b, true, iter);
        if ((std::abs(e) >= toler) && (std::abs(fa) > std::abs(fb))) {
            s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0)
                q = -q;
            p = std::abs(p);
            min1 = 3.0 * xm * q - std::abs(toler * q);
            min2 = std::abs(e * q);
            if (2.0 * p < (min1 < min2 ? min1 : min2)) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        if (std::abs(d) > toler)
            b += d;
        else
            b += (xm > 0.0) ? std::abs(toler) : -std::abs(toler);
        fb = fn.eval(b).value();
    }
    return RootResult(b, false, MAXIT); // SNH
}

} // namespace redukti::mathlib
