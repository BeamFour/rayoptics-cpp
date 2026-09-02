// Ported from TestHybrd2.java: the standard MINPACK hybrd1 test battery
// (22 problems, 55 cases) with the exact expected fnorm/nfev/info values.
// The Java only asserts fnorm; nfev and info are asserted here as well, since
// the expected table carries them and they pin the whole iteration sequence.
#include "TestHarness.h"

#include "redukti/mathlib/MinPack.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace redukti::mathlib;

namespace {

struct TestHybrd2 : Hybrd_Function {

    static constexpr double epsmch = 2.22044604926e-16;
    static constexpr double zero =    0.0;
    static constexpr double half =     .5;
    static constexpr double one =     1.0;
    static constexpr double two =     2.0;
    static constexpr double three =   3.0;
    static constexpr double four =    4.0;
    static constexpr double five =    5.0;
    static constexpr double seven =   7.0;
    static constexpr double eight =   8.0;
    static constexpr double ten =    10.0;
    static constexpr double twenty = 20.0;
    static constexpr double twntf =  25.0;

    int nfev = 0;
    int nprob =0;

    struct TestProblem {
        int nprob, n, ntries;
    };

    struct ExpectedResult {
        int nprob;
        int n;
        int nfev;
        int info;
        double fnorm;
    };
    void run () {

        const std::vector<TestProblem> problems = {
                {1,    2,    3},
                {2,    4,    3},
                {3,    2,    2},
                {4,    4,    3},
                {5,    3,    3},
                {6,    6,    2},
                {6,    9,    2},
                {7,    5,    3},
                {7,    6,    3},
                {7,    7,    3},
                {7,    8,    1},
                {7,    9,    1},
                {8,   10,    3},
                {8,   30,    1},
                {8,   40,    1},
                {9,   10,    3},
                {10,    1,    3},
                {10,   10,    3},
                {11,   10,    3},
                {12,   10,    3},
                {13,   10,    3},
                {14,   10,    3}
        };

        const std::vector<ExpectedResult> expectedResults = {
                {1,2,22,1,0.0},
                {1,2,9,1,0.0},
                {1,2,9,1,5.373479439185758E-13},
                {2,4,106,4,1.630462542097203E-33},
                {2,4,108,4,6.076260751859731E-33},
                {2,4,132,4,1.6826101886919138E-33},
                {3,2,181,1,1.7124935913770545E-9},
                {3,2,11,1,3.7444845557566E-8},
                {4,4,94,1,4.005066887716473E-11},
                {4,4,234,1,5.154548231991567E-10},
                {4,4,514,1,2.4911329353710294E-11},
                {5,3,27,1,2.753457427134646E-13},
                {5,3,31,1,4.4279070194841094E-14},
                {5,3,40,1,1.610826463283445E-13},
                {6,6,96,1,3.1916461155034835E-13},
                {6,6,310,1,2.3361280345344312E-11},
                {6,9,167,1,1.817744267637861E-14},
                {6,9,167,1,1.1388692433566192E-12},
                {7,5,17,1,3.540291609809256E-12},
                {7,5,256,1,1.9023836575828113E-10},
                {7,5,522,1,2.2691301265379657E-11},
                {7,6,25,1,6.841873686727006E-10},
                {7,6,172,1,2.281995876562388E-11},
                {7,6,280,1,7.373849407048581E-11},
                {7,7,20,1,2.398700313876028E-9},
                {7,7,525,1,1.6243564328339643E-10},
                {7,7,426,4,1.1214382813465733},
                {7,8,120,4,0.06440508384820999},
                {7,9,41,1,1.8898077704177782E-9},
                {8,10,31,1,1.312131443737577E-14},
                {8,10,31,1,1.6011864169946884E-14},
                {8,10,38,1,4.213000162292041E-15},
                {8,30,113,1,2.0734139515353978E-13},
                {8,40,196,1,8.885114241808755E-14},
                {9,10,16,1,2.533550868211982E-15},
                {9,10,19,1,1.7259388604461542E-13},
                {9,10,52,1,4.177635246055282E-10},
                {10,1,7,1,5.551115123125783E-17},
                {10,1,9,1,5.551115123125783E-17},
                {10,1,16,1,2.7755575615628914E-17},
                {10,10,16,1,5.009131711331394E-15},
                {10,10,19,1,2.1883100011823385E-13},
                {10,10,39,1,3.0348043981948507E-15},
                {11,10,130,4,0.005296401898597665},
                {11,10,84,1,5.915975632815838E-11},
                {11,10,85,1,1.8562840902695322E-9},
                {12,10,31,1,5.1194673452253975E-12},
                {12,10,35,1,7.399639659366749E-11},
                {12,10,66,1,0.0},
                {13,10,21,1,1.4938793750003843E-8},
                {13,10,59,1,5.3374415147502625E-9},
                {13,10,42,1,9.878879503809141E-11},
                {14,10,30,1,2.0579005008925252E-9},
                {14,10,45,1,7.953614269648959E-9},
                {14,10,58,1,4.526425424093767E-10},
        };

        int IC,INFO,K,LWA;
        std::vector<int> NA(60, 0), NF(60, 0), NP(60, 0), NX(60, 0);
        double FACTOR,FNORM1,FNORM2,TOL;
        std::vector<double> FNM(60, 0.0), FVEC(40, 0.0), WA(2660, 0.0), X(40, 0.0);

        TOL = std::sqrt(MinPack::dpmpar(1));
        LWA = 2660;
        IC = 0;

        for (const TestProblem &p : problems) {
            FACTOR = one;
            nprob = p.nprob;
            for (K = 1; K <= p.ntries; K++) {
                IC = IC + 1;
                initpt(p.n,X,p.nprob,FACTOR);
                vecfcn(p.n,X,FVEC,p.nprob);
                FNORM1 = MinPack::enorm(p.n, 0, FVEC);
                nfev = 0;
                INFO = MinPack::hybrd1(*this,p.n,X,FVEC,TOL,WA,LWA);
                FNORM2 = MinPack::enorm(p.n,0,FVEC);
                NP[IC-1] = p.nprob;
                NA[IC-1] = p.n;
                NF[IC-1] = nfev;
                NX[IC-1] = INFO;
                FNM[IC-1] = FNORM2;
                const ExpectedResult &want = expectedResults[IC-1];

                // The iteration path must match the JVM exactly on every
                // problem: same number of function evaluations, same exit code.
                CHECK_EQ(nfev, want.nfev);
                CHECK_EQ(INFO, want.info);

                // Problem 11 is the only one in the battery whose residual
                // function calls sin/cos. Those are not correctly rounded and
                // the JVM's differ from libm's by an ulp on some arguments,
                // which the residual norm amplifies. nfev and info still match
                // exactly there, so the algorithm is taking the identical path;
                // only the last few digits of the residual move. Every other
                // problem is polynomial or rational and matches to 1e-15.
                // An absolute bound, not a relative one: two of the three
                // problem-11 cases converge to a residual of ~1e-9, where a
                // relative comparison is measuring noise against noise. The
                // largest residual in the battery is 5.3e-3, so 1e-6 absolute
                // is still a real constraint. Observed deltas: 2.5e-8, 2e-15,
                // 4.2e-11.
                if (p.nprob == 11) {
                    CHECK(std::abs(FNORM2 - want.fnorm) < 1e-6);
                } else {
                    CHECK_CLOSE(FNORM2, want.fnorm, 1e-15);
                }
                FACTOR = ten*FACTOR;
            }
        }

    }
    void apply(int n, std::vector<double> &x, std::vector<double> &fvec, std::vector<int> &iflag) override {
        vecfcn(n,x,fvec,nprob);
        nfev = nfev + 1;
    }
    static void vecfcn(int n, std::vector<double> &x, std::vector<double> &fvec,
                              int nprob) {

        int i,iev,j,k,k1,k2,kp1,ml,mu;
        double c1=1e4,c2=1.0001,c3=2.0e2,c4=2.02e1,c5=1.98e1,c6=1.8e2,c7=2.5e-1,c8=5.0e-1,c9=2.9e1;

        switch (nprob) {

            case 1:

                fvec[0] = one - x[0];
                fvec[1] = ten*(x[1] - x[0]*x[0]);

                return;

            case 2:

                fvec[0] = x[0] + ten*x[1];
                fvec[1] = std::sqrt(five)*(x[2] - x[3]);
                fvec[2] = (x[1] - two*x[2])*(x[1] - two*x[2]);
                fvec[3] = std::sqrt(ten)*(x[0] - x[3])*(x[0] - x[3]);

                return;

            case 3:

                fvec[0] = c1*x[0]*x[1] - one;
                fvec[1] = std::exp(-x[0]) + std::exp(-x[1]) - c2;

                return;

            case 4: {

                double temp1 = x[1] - x[0] * x[0];
                double temp2 = x[3] - x[2] * x[2];
                fvec[0] = -c3 * x[0] * temp1 - (one - x[0]);
                fvec[1] = c3 * temp1 + c4 * (x[1] - one) + c5 * (x[3] - one);
                fvec[2] = -c6 * x[2] * temp2 - (one - x[2]);
                fvec[3] = c6 * temp2 + c4 * (x[3] - one) + c5 * (x[1] - one);

                return;
            }

            case 5: {

                double tpi = eight * std::atan(one);
                double temp1;
                if (x[1] < 0.0) {
                    temp1 = -c7;
                } else {
                    temp1 = c7;
                }
                if (x[0] > zero) temp1 = std::atan(x[1] / x[0]) / tpi;
                if (x[0] < zero) temp1 = std::atan(x[1] / x[0]) / tpi + c8;
                double temp2 = std::sqrt(x[0] * x[0] + x[1] * x[1]);
                fvec[0] = ten * (x[2] - ten * temp1);
                fvec[1] = ten * (temp2 - one);
                fvec[2] = x[2];

                return;

            }

            case 6: {

                for (k = 1; k <= n; k++) {
                    fvec[k-1] = zero;
                }

                for (i = 1; i <= 29; i++) {

                    double ti = i/c9;
                    double sum1 = zero;
                    double temp = one;

                    for (j = 2; j <= n; j++) {
                        sum1 += (double)(j-1)*temp*x[j-1];
                        temp *= ti;
                    }

                    double sum2 = zero;
                    temp = one;

                    for (j = 1; j <= n; j++) {
                        sum2 += temp*x[j-1];
                        temp *= ti;
                    }

                    double temp1 = sum1 - sum2*sum2 - one;
                    double temp2 = two*ti*sum2;
                    temp = one/ti;

                    for (k = 1; k <=n; k++) {
                        fvec[k-1] = fvec[k-1] + temp*((double)(k-1) - temp2)*temp1;
                        temp = ti*temp;
                    }
                }

                double temp = x[1] - x[0]*x[0] - one;
                fvec[0] = fvec[0] + x[0]*(one - two*temp);
                fvec[1] = fvec[1] + temp;

                return;

            }

            case 7: {

                for (k = 1; k <= n; k++) {
                    fvec[k-1] = zero;
                }
                for (j = 1; j <= n; j++) {
                    double tmp1 = one;
                    double tmp2 = two*x[j-1] - one;
                    double temp = two*tmp2;

                    for (i = 1; i <= n; i++) {
                        fvec[i-1] += tmp2;
                        double ti = temp*tmp2 - tmp1;
                        tmp1 = tmp2;
                        tmp2 = ti;
                    }
                }

                double dx = one/(double)n;
                iev = -1;

                for (k = 1; k <= n; k++) {
                    fvec[k-1] *= dx;
                    if (iev > 0) fvec[k-1] += one/((double)(k*k) - one);
                    iev = -iev;
                }

                return;

            }

            case 8: {

                double sum = -((double)n+1.0);
                double prod = one;

                for (j = 1; j <= n; j++) {
                    sum += x[j-1];
                    prod *= x[j-1];
                }

                for (k = 1; k <= n; k++) {
                    fvec[k-1] = x[k-1] + sum;
                }

                fvec[n-1] = prod - one;

                return;

            }

            case 9: {

                double h = one/(double)(n+1);
                for (k =1; k <= n; k++) {
                    double t = x[k-1] + (double)(k) * h + one;
                    double temp = t*t*t;
                    double temp1 = zero;
                    if (k != 1)
                        temp1 = x[k-2];
                    double temp2 = zero;
                    if (k != n)
                        temp2 = x[k];
                    fvec[k-1] = two * x[k-1] - temp1 - temp2 + temp * h*h / two;
                }

                return;

            }

            case 10: {

                double h = one/(double)(n+1);
                for (k = 1; k <= n; k++) {
                    double tk = (double)(k)*h;
                    double sum1 = zero;
                    for (j = 1; j<= k; j++) {
                        double tj = (double)(j)*h;
                        double t = x[j-1] + tj + one;
                        double temp = t*t*t;
                        sum1 = sum1 + tj*temp;
                    }
                    double sum2 = zero;
                    kp1 = k + 1;
                    if (n >= kp1) {
                        for (j = kp1; j <= n; j++) {
                            double tj = (double)(j)*h;
                            double t = x[j-1] + tj + one;
                            double temp = t*t*t;
                            sum2 = sum2 + (one - tj)*temp;
                        }
                    }
                    fvec[k-1] = x[k-1] + h*((one - tk)*sum1 + tk*sum2)/two;
                }

                return;

            }

            case 11: {

                double sum = zero;
                for (j = 1; j <= n; j++) {
                    fvec[j-1] = std::cos(x[j-1]);
                    sum = sum + fvec[j-1];
                }
                for (k = 1; k <= n; k++) {
                    fvec[k-1] = (double)(n+k) - std::sin(x[k-1]) - sum - (double)(k)*fvec[k-1];
                }

                return;
            }

            case 12: {

                double sum = zero;
                for (j = 1; j <= n; j++) {
                    sum = sum + (double)(j)*(x[j-1] - one);
                }
                double temp = sum*(one + two*sum*sum);
                for (k = 1; k <= n; k++) {
                    fvec[k-1] = x[k-1] - one + (double)(k)*temp;
                }

                return;
            }

            case 13: {

                for (k = 1; k <= n; k++) {
                    double temp = (three - two*x[k-1])*x[k-1];
                    double temp1 = zero;
                    if (k != 1)
                        temp1 = x[k-2];
                    double temp2 = zero;
                    if (k != n)
                        temp2 = x[k];
                    fvec[k-1] = temp - temp1 - two*temp2 + one;
                }

                return;

            }

            case 14: {

                ml = 5;
                mu = 1;
                for (k = 1; k<= n; k++) {
                    k1 = std::max(1,k-ml); // MAX0?
                    k2 = std::min(k+mu,n);     // MIN0
                    double temp = zero;
                    for (j = k1; j <= k2; j++) {
                        if (j != k)
                            temp = temp + x[j-1]*(one + x[j-1]);
                    }
                    fvec[k-1] = x[k-1]*(two + five*x[k-1]*x[k-1]) + one - temp;
                }

                return;

            }

        }

    }
    static void initpt(int n, std::vector<double> &x, int nprob, double factor) {

        int j;

        double c1 = 1.2,h;

        switch (nprob) {

            case 1:

                x[0] = -c1;
                x[1] = one;

                break;

            case 2:

                x[0] = three;
                x[1] = -one;
                x[2] = zero;
                x[3] = one;

                break;

            case 3:

                x[0] = zero;
                x[1] = one;

                break;

            case 4:

                x[0] = -three;
                x[1] = -one;
                x[2] = -three;
                x[3] = -one;

                break;

            case 5:

                x[0] = -one;
                x[1] = zero;
                x[2] = zero;

                break;

            case 6:

                for (j = 1; j <= n; j++) {
                    x[j-1] = zero;
                }

                break;

            case 7:

                h = one/(double)(n+1);
                for (j = 1; j <= n; j++) {
                    x[j-1] = j*h;
                }

                break;

            case 8:

                for (j = 1; j <= n; j++) {
                    x[j-1] = half;
                }

                break;

            case 9:
            case 10:

                h = one/(double)(n+1);
                for (j = 1; j <= n; j++) {
                    double tj = (double) j * h;
                    x[j-1] = tj * (tj - one);
                }

                break;

            case 11:

                h = one/(double)(n);
                for (j = 1; j <= n; j++) {
                    x[j-1] = h;
                }

                break;

            case 12:

                h = one/(double)(n);
                for (j = 1; j <= n; j++) {
                    x[j-1] = one - (double)(j) * h;
                }

                break;

            case 13:
            case 14:

                for (j = 1; j <= n; j++) {
                    x[j-1] = -one;
                }

                break;

        }

        if (factor == one) return;

        if (nprob != 6) {

            for (j = 1; j <= n; j++) {
                x[j-1] *= factor;
            }

        } else {

            for (j = 1; j <= n; j++) {
                x[j-1] = factor;
            }

        }

        return;

    }

};

} // namespace

TEST(minpack_hybrd1_battery) {
    TestHybrd2 instance;
    instance.run();
}
