// Ported from TestHybrd.java and TestLmdir.java -- the standard MINPACK driver
// examples, with the same expected values and tolerances.
#include "TestHarness.h"

#include "redukti/mathlib/MinPack.h"

#include <cmath>
#include <vector>

using namespace redukti::mathlib;

namespace {

/** subroutine fcn for hybrd example. */
struct HybrdExample : Hybrd_Function {
    void apply(int n, std::vector<double> &x, std::vector<double> &fvec,
               std::vector<int> &iflag) override {
        int k;
        double temp, temp1, temp2;

        if (iflag[0] == 0) {
            // insert print statements here when nprint is positive.
            return;
        }

        /* compute residuals */
        for (k = 1; k <= n; k++) {
            temp = (3 - 2 * x[k - 1]) * x[k - 1];
            temp1 = 0;
            if (k != 1) {
                temp1 = x[k - 1 - 1];
            }
            temp2 = 0;
            if (k != n) {
                temp2 = x[k + 1 - 1];
            }
            fvec[k - 1] = temp - temp1 - 2 * temp2 + 1;
        }
    }
};

} // namespace

TEST(minpack_hybrd_example) {
    int j, n, maxfev, ml, mu, mode, nprint, info, ldfjac, lr;
    std::vector<int> nfev(1, 0);
    double xtol, epsfcn, factor, fnorm;
    std::vector<double> x(9, 0.0), fvec(9, 0.0), diag(9, 0.0), fjac(9 * 9, 0.0),
        r(45, 0.0), qtf(9, 0.0), wa1(9, 0.0), wa2(9, 0.0), wa3(9, 0.0), wa4(9, 0.0);

    n = 9;

    /* the following starting values provide a rough solution. */
    for (j = 1; j <= 9; j++) {
        x[j - 1] = -1.;
    }

    ldfjac = 9;
    lr = 45;

    /* set xtol to the square root of the machine precision. */
    xtol = std::sqrt(MinPack::dpmpar(1));

    maxfev = 2000;
    ml = 1;
    mu = 1;
    epsfcn = 0.;
    mode = 2;
    for (j = 1; j <= 9; j++) {
        diag[j - 1] = 1.;
    }

    factor = 1.e2;
    nprint = 0;

    HybrdExample fcn;
    info = MinPack::hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, diag, mode,
                          factor, nprint, nfev, fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3,
                          wa4);
    fnorm = MinPack::enorm(n, 0, fvec);

    CHECK_CLOSE(fnorm, 0.1192636E-07, 1e-7);
    CHECK_EQ(nfev[0], 14);
    CHECK_EQ(info, 1);

    const double expected[9] = {-0.5706545, -0.6816283, -0.7017325,
                                -0.7042129, -0.7013690, -0.6918656,
                                -0.6657920, -0.5960342, -0.4164121};
    for (std::size_t i = 0; i < x.size(); i++) {
        CHECK_CLOSE(x[i], expected[i], 1e-7);
    }
}

namespace {

/** subroutine fcn for lmder example. */
struct LmderExample : Lmder_Function {
    using Lmder_Function::apply; // the five-argument overload stays visible

    std::vector<double> y = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1,
                             3.2e-1, 3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1,
                             7.3e-1, 9.6e-1, 1.34,   2.1,    4.39};

    int apply(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
              std::vector<double> &fjac, int ldfjac, int iflag) override {
        int i;
        double tmp1, tmp2, tmp3, tmp4;
        (void)m;
        (void)n;

        if (iflag == 0) {
            // insert print statements here when nprint is positive.
            return 0;
        }

        if (iflag != 2) {
            /* compute residuals */
            for (i = 0; i < 15; ++i) {
                tmp1 = i + 1;
                tmp2 = 15 - i;
                tmp3 = (i > 7) ? tmp2 : tmp1;
                fvec[i] = y[i] - (x[0] + tmp1 / (x[1] * tmp2 + x[2] * tmp3));
            }
        } else {
            /* compute Jacobian */
            for (i = 0; i < 15; ++i) {
                tmp1 = i + 1;
                tmp2 = 15 - i;
                tmp3 = (i > 7) ? tmp2 : tmp1;
                tmp4 = (x[1] * tmp2 + x[2] * tmp3);
                tmp4 = tmp4 * tmp4;
                fjac[i + ldfjac * 0] = -1.;
                fjac[i + ldfjac * 1] = tmp1 * tmp2 / tmp4;
                fjac[i + ldfjac * 2] = tmp1 * tmp3 / tmp4;
            }
        }
        return 0;
    }

    bool hasJacobian() override { return true; }
};

} // namespace

TEST(minpack_lmder_example) {
    int ldfjac, maxfev, mode, nprint, info;
    std::vector<int> nfev(1, 0), njev(1, 0);
    std::vector<int> ipvt(3, 0);
    double ftol, xtol, gtol, factor, fnorm;
    std::vector<double> x(3, 0.0);
    std::vector<double> fvec(15, 0.0);
    std::vector<double> fjac(15 * 3, 0.0);
    std::vector<double> diag(3, 0.0);
    std::vector<double> qtf(3, 0.0);
    std::vector<double> wa1(3, 0.0);
    std::vector<double> wa2(3, 0.0);
    std::vector<double> wa3(3, 0.0);
    std::vector<double> wa4(15, 0.0);
    const int m = 15;
    const int n = 3;

    /* the following starting values provide a rough fit. */
    x[0] = 1.;
    x[1] = 1.;
    x[2] = 1.;

    ldfjac = 15;

    ftol = std::sqrt(MinPack::dpmpar(1));
    xtol = std::sqrt(MinPack::dpmpar(1));
    gtol = 0.;

    maxfev = 400;
    mode = 1;
    factor = 1.e2;
    nprint = 0;

    LmderExample fcn;
    info = MinPack::lmder(fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev,
                          diag, mode, factor, nprint, nfev, njev, ipvt, qtf, wa1, wa2,
                          wa3, wa4);

    fnorm = MinPack::enorm(m, 0, fvec);
    CHECK_CLOSE(fnorm, 0.09063596, 1e-8);
    CHECK_EQ(nfev[0], 6);
    CHECK_EQ(njev[0], 5);
    CHECK_EQ(info, 1);
    CHECK_CLOSE(x[0], 0.08241058, 1e-6);
    CHECK_CLOSE(x[1], 1.133037, 1e-6);
    CHECK_CLOSE(x[2], 2.343695, 1e-6);
}

TEST(minpack_dpmpar) {
    CHECK_CLOSE(MinPack::dpmpar(1), 2.22044604926e-16, 0.0);
    CHECK_CLOSE(MinPack::dpmpar(2), 2.22507385852e-308, 0.0);
    CHECK_CLOSE(MinPack::dpmpar(3), 1.79769313485e+308, 0.0);
}

TEST(minpack_enorm) {
    std::vector<double> v = {3.0, 4.0, 0.0, 12.0};
    CHECK_CLOSE(MinPack::enorm(2, 0, v), 5.0, 0.0);
    CHECK_CLOSE(MinPack::enorm(2, 2, v), 12.0, 0.0);
    CHECK_CLOSE(MinPack::enorm(0, 0, v), 0.0, 0.0);
}
