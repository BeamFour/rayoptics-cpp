// C++ port of org.redukti.mathlib.MinPack
//
// argonne national laboratory. minpack project. march 1980.
// burton s. garbow, kenneth e. hillstrom, jorge j. more
//
// Transliterated from the Java, which preserves the original Fortran
// structure: flat arrays with explicit leading dimensions, explicit start
// offsets in place of pointer arithmetic, and single-element arrays used as
// in/out parameters. Java's labelled breaks become gotos to a label placed
// immediately after the loop, and lmpar's try/catch -- which the Java used to
// emulate Fortran's "go to 220" -- becomes a plain goto.
#include "redukti/mathlib/MinPack.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace redukti::mathlib {



    double MinPack::dpmpar(int i) {
        switch (i) {
            case 1:
                return 2.22044604926e-16;
            case 2:
                return 2.22507385852e-308;
            default:
                return 1.79769313485e+308;
        }
    }

    
    int MinPack::lmder1(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                             std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac, double tol,
                             std::vector<int> &ipvt, std::vector<double> &wa, int lwa) {
        // Setting epsfcn to 0.0 causes machine precision to be used
        return lmder1(fcn, m, n, x, fvec, fjac, ldfjac, tol, ipvt, wa, lwa, 0.0);
    }
    
    int MinPack::lmder1(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                             std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac, double tol,
                             std::vector<int> &ipvt, std::vector<double> &wa, int lwa, double epsfcn) {
        
        const double factor = 100.;

        int mode;
        std::vector<int> nfev(1, 0), njev(1, 0);
        double ftol, gtol, xtol;
        int maxfev, nprint;
        int info;

        if (n <= 0 || m < n || ldfjac < m || tol < 0. || lwa < n * 5 + m) {
            return 0;
        }

        maxfev = (n + 1) * 100;
        ftol = tol;
        xtol = tol;
        gtol = 0.;
        mode = 1;
        nprint = 0;

        std::vector<double> wa1(n, 0.0);
        std::vector<double> wa2(n, 0.0);
        std::vector<double> wa3(n, 0.0);
        std::vector<double> wa4(n, 0.0);
        std::vector<double> wa5(m, 0.0);

        info = lmder(fcn, m, n, x, fvec, fjac, ldfjac,
                ftol, xtol, gtol, maxfev, wa, mode, factor, nprint,
                nfev, njev, ipvt, wa1, wa2, wa3, wa4, wa5, epsfcn);
        if (info == 8) {
            info = 4;
        }
        return info;

    }

    int MinPack::lmder(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                            std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac, double ftol,
                            double xtol, double gtol, int maxfev, std::vector<double> &diag, int mode, double factor, int nprint,
                            std::vector<int> &nfev, std::vector<int> &njev, std::vector<int> &ipvt, std::vector<double> &qtf,
                            std::vector<double> &wa1, std::vector<double> &wa2, std::vector<double> &wa3, std::vector<double> &wa4) {
        // Setting epsfcn to 0.0 causes machine precision to be used
        return lmder(fcn, m, n, x, fvec, fjac, ldfjac,
                ftol, xtol, gtol, maxfev, diag, mode, factor, nprint,
                nfev, njev, ipvt, qtf, wa1, wa2, wa3, wa4, 0.0);
    }

    int MinPack::lmder(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                            std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac, double ftol,
                            double xtol, double gtol, int maxfev, std::vector<double> &diag, int mode, double factor, int nprint,
                            std::vector<int> &nfev, std::vector<int> &njev, std::vector<int> &ipvt, std::vector<double> &qtf,
                            std::vector<double> &wa1, std::vector<double> &wa2, std::vector<double> &wa3, std::vector<double> &wa4, double epsfcn)
    {
        
        const double p1 = .1;
        const double p5 = .5;
        const double p25 = .25;
        const double p75 = .75;
        const double p0001 = 1e-4;

        double d1, d2;

        int i, j, l;
        std::vector<double> par(1, 0.0);
        double sum;
        int iter;
        double temp, temp1, temp2;
        int iflag;
        double delta = 0.;
        double ratio;
        double fnorm, gnorm, pnorm, xnorm = 0., fnorm1, actred, dirder,
                epsmch, prered;
        int info;

        epsmch = dpmpar(1);

        info = 0;
        iflag = 0;
        nfev[0] = 0;
        njev[0] = 0;

        do {
            if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. ||
                    gtol < 0. || maxfev <= 0 || factor <= 0.) {
                goto processing_end;
            }
            if (mode == 2) {
                for (j = 0; j < n; ++j) {
                    if (diag[j] <= 0.) {
                        goto processing_end;
                    }
                }
            }

            iflag = fcn.hasJacobian() ?
                    fcn.apply(m, n, x, fvec, fjac, ldfjac, 1) :
                    fcn.apply(m, n, x, fvec, 1);
            nfev[0] = 1;
            if (iflag < 0) {
                goto processing_end;
            }
            fnorm = enorm(m, 0, fvec);

            par[0] = 0.;
            iter = 1;

            for (; ; ) {

                iflag = fcn.hasJacobian() ?
                        fcn.apply(m, n, x, fvec, fjac, ldfjac, 2) :
                        fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, epsfcn, wa4);
                njev[0] = njev[0]+1;
                if (iflag < 0) {
                    goto processing_end;
                }

                if (nprint > 0) {
                    iflag = 0;
                    if ((iter - 1) % nprint == 0) {
                        iflag = fcn.hasJacobian() ?
                                fcn.apply(m, n, x, fvec, fjac, ldfjac, 0) :
                                fcn.apply(m, n, x, fvec, 0);
                    }
                    if (iflag < 0) {
                        goto processing_end;
                    }
                }

                qrfac(m, n, fjac, ldfjac, 1, ipvt, n,
                        wa1, wa2, wa3);

                if (iter == 1) {
                    if (mode != 2) {
                        for (j = 0; j < n; ++j) {
                            diag[j] = wa2[j];
                            if (wa2[j] == 0.) {
                                diag[j] = 1.;
                            }
                        }
                    }

                    for (j = 0; j < n; ++j) {
                        wa3[j] = diag[j] * x[j];
                    }
                    xnorm = enorm(n, 0, wa3);
                    delta = factor * xnorm;
                    if (delta == 0.) {
                        delta = factor;
                    }
                }

                for (i = 0; i < m; ++i) {
                    wa4[i] = fvec[i];
                }
                for (j = 0; j < n; ++j) {
                    if (fjac[j + j * ldfjac] != 0.) {
                        sum = 0.;
                        for (i = j; i < m; ++i) {
                            sum += fjac[i + j * ldfjac] * wa4[i];
                        }
                        temp = -sum / fjac[j + j * ldfjac];
                        for (i = j; i < m; ++i) {
                            wa4[i] += fjac[i + j * ldfjac] * temp;
                        }
                    }
                    fjac[j + j * ldfjac] = wa1[j];
                    qtf[j] = wa4[j];
                }

                gnorm = 0.;
                if (fnorm != 0.) {
                    for (j = 0; j < n; ++j) {
                        l = ipvt[j]-1;
                        if (wa2[l] != 0.) {
                            sum = 0.;
                            for (i = 0; i <= j; ++i) {
                                sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
                            }
                            
                            d1 = std::abs(sum / wa2[l]);
                            gnorm = std::max(gnorm, d1);
                        }
                    }
                }

                if (gnorm <= gtol) {
                    info = 4;
                }
                if (info != 0) {
                    goto processing_end;
                }

                if (mode != 2) {
                    for (j = 0; j < n; ++j) {
                        
                        d1 = diag[j];
                        d2 = wa2[j];
                        diag[j] = std::max(d1, d2);
                    }
                }

                do {

                    lmpar (n, fjac, ldfjac, ipvt, diag, qtf, delta,
                            par, wa1, wa2, wa3, wa4);

                    for (j = 0; j < n; ++j) {
                        wa1[j] = -wa1[j];
                        wa2[j] = x[j] + wa1[j];
                        wa3[j] = diag[j] * wa1[j];
                    }
                    pnorm = enorm(n, 0, wa3);

                    if (iter == 1) {
                        delta = std::min(delta, pnorm);
                    }

                    iflag = fcn.hasJacobian() ?
                            fcn.apply(m, n, wa2, wa4, fjac, ldfjac, 1) :
                            fcn.apply(m, n, wa2, wa4, 1);
                    nfev[0] = nfev[0] + 1;
                    if (iflag < 0) {
                        goto processing_end;
                    }
                    fnorm1 = enorm(m, 0, wa4);

                    actred = -1.;
                    if (p1 * fnorm1 < fnorm) {
                        
                        d1 = fnorm1 / fnorm;
                        actred = 1 - d1 * d1;
                    }

                    for (j = 0; j < n; ++j) {
                        wa3[j] = 0.;
                        l = ipvt[j] - 1;
                        temp = wa1[l];
                        for (i = 0; i <= j; ++i) {
                            wa3[i] += fjac[i + j * ldfjac] * temp;
                        }
                    }
                    temp1 = enorm(n, 0, wa3) / fnorm;
                    temp2 = (std::sqrt(par[0]) * pnorm) / fnorm;
                    prered = temp1 * temp1 + temp2 * temp2 / p5;
                    dirder = -(temp1 * temp1 + temp2 * temp2);

                    ratio = 0.;
                    if (prered != 0.) {
                        ratio = actred / prered;
                    }

                    if (ratio <= p25) {
                        if (actred >= 0.) {
                            temp = p5;
                        } else {
                            temp = p5 * dirder / (dirder + p5 * actred);
                        }
                        if (p1 * fnorm1 >= fnorm || temp < p1) {
                            temp = p1;
                        }
                        
                        d1 = pnorm / p1;
                        delta = temp * std::min(delta, d1);
                        par[0] = par[0] / temp;
                    } else {
                        if (par[0] == 0. || ratio >= p75) {
                            delta = pnorm / p5;
                            par[0] = p5 * par[0];
                        }
                    }

                    if (ratio >= p0001) {

                        for (j = 0; j < n; ++j) {
                            x[j] = wa2[j];
                            wa2[j] = diag[j] * x[j];
                        }
                        for (i = 0; i < m; ++i) {
                            fvec[i] = wa4[i];
                        }
                        xnorm = enorm(n, 0, wa2);
                        fnorm = fnorm1;
                        ++iter;
                    }

                    if (std::abs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
                        info = 1;
                    }
                    if (delta <= xtol * xnorm) {
                        info = 2;
                    }
                    if (std::abs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info == 2) {
                        info = 3;
                    }
                    if (info != 0) {
                        goto processing_end;
                    }

                    if (nfev[0] >= maxfev){
                        info = 5;
                    }
                    if (std::abs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.) {
                        info = 6;
                    }
                    if (delta <= epsmch * xnorm) {
                        info = 7;
                    }
                    if (gnorm <= epsmch) {
                        info = 8;
                    }
                    if (info != 0) {
                        goto processing_end;
                    }

                } while (ratio < p0001);

            }
        } while (false);
        processing_end:;

        if (iflag < 0) {
            info = iflag;
        }
        if (nprint > 0) {
            int t = fcn.hasJacobian() ?
                    fcn.apply(m, n, x, fvec, fjac, ldfjac, 0) :
                    fcn.apply(m, n, x, fvec, 0);
        }
        return info;

    }

    double MinPack::enorm(int n, int start, std::vector<double> &x) {
        int i;
        double agiant, floatn, s1, s2, s3, xabs,
                x1max, x3max;
        double enorm;

        const double one = 1.0;
        const double zero = 0.0;
        const double rdwarf = 3.834e-20;
        const double rgiant = 1.304e+19;

        s1 = 0.0;
        s2 = 0.0;
        s3 = 0.0;
        x1max = 0.0;
        x3max = 0.0;
        floatn = n;
        agiant = rgiant / floatn;

        for (i = start; i < n + start; i++) {
            xabs = std::abs(x[i]);
            if (xabs <= rdwarf || xabs >= agiant) {
                if (xabs > rdwarf) {
                    // Sum for large components.
                    if (xabs > x1max) {
                        s1 = one + s1 * (x1max / xabs) * (x1max / xabs);
                        x1max = xabs;
                    } else {
                        s1 += (xabs / x1max) * (xabs / x1max);
                    }
                } else {
                    // Sum for small components.
                    if (xabs > x3max) {
                        s3 = one + s3 * (x3max / xabs) * (x3max / xabs);
                        x3max = xabs;
                    } else {
                        if (xabs != zero) s3 += (xabs / x3max) * (xabs / x3max);
                    }
                }
            } else {
                // Sum for intermediate components.
                s2 += xabs * xabs;
            }
        }

        // Calculation of norm.
        if (s1 != zero) {
            enorm = x1max * std::sqrt(s1 + (s2 / x1max) / x1max);
        } else {
            if (s2 != zero) {
                if (s2 >= x3max) {
                    enorm = std::sqrt(s2 * (one + (x3max / s2) * (x3max * s3)));
                } else {
                    enorm = std::sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
                }
            } else {
                enorm = x3max * std::sqrt(s3);
            }
        }
        return enorm;
    }

    void MinPack::qrfac(int m, int n, std::vector<double> &a, int lda,
                             int pivot, std::vector<int> &ipvt, int lipvt, std::vector<double> &rdiag,
                             std::vector<double> &acnorm, std::vector<double> &wa) {
        
        const double p05 = .05;

        double d1;

        int i, j, k, jp1;
        double sum;
        double temp;
        int minmn;
        double epsmch;
        double ajnorm;

        epsmch = dpmpar(1);

        for (j = 0; j < n; ++j) {
            acnorm[j] = enorm(m, j * lda, a);
            rdiag[j] = acnorm[j];
            wa[j] = rdiag[j];
            if (pivot != 0) {
                ipvt[j] = j + 1;
            }
        }

        minmn = std::min(m, n);
        for (j = 0; j < minmn; j++) {
            if (pivot != 0) {

                int kmax = j;
                for (k = j; k < n; ++k) {
                    if (rdiag[k] > rdiag[kmax]) {
                        kmax = k;
                    }
                }
                if (kmax != j) {
                    for (i = 0; i < m; ++i) {
                        temp = a[i + j * lda];
                        a[i + j * lda] = a[i + kmax * lda];
                        a[i + kmax * lda] = temp;
                    }
                    rdiag[kmax] = rdiag[j];
                    wa[kmax] = wa[j];
                    k = ipvt[j];
                    ipvt[j] = ipvt[kmax];
                    ipvt[kmax] = k;
                }
            }

            ajnorm = enorm(m - j, j + j * lda, a);
            if (ajnorm != 0.) {
                if (a[j + j * lda] < 0.) {
                    ajnorm = -ajnorm;
                }
                for (i = j; i < m; ++i) {
                    a[i + j * lda] /= ajnorm;
                }
                a[j + j * lda] += 1;

                jp1 = j + 1;
                if (n > jp1) {
                    for (k = jp1; k < n; ++k) {
                        sum = 0.;
                        for (i = j; i < m; ++i) {
                            sum += a[i + j * lda] * a[i + k * lda];
                        }
                        temp = sum / a[j + j * lda];
                        for (i = j; i < m; ++i) {
                            a[i + k * lda] -= temp * a[i + j * lda];
                        }
                        if (pivot != 0 && rdiag[k] != 0.) {
                            temp = a[j + k * lda] / rdiag[k];
                            
                            d1 = 1 - temp * temp;
                            rdiag[k] *= std::sqrt((std::max(0., d1)));
                            
                            d1 = rdiag[k] / wa[k];
                            if (p05 * (d1 * d1) <= epsmch) {
                                rdiag[k] = enorm(m - (j + 1), jp1 + k * lda, a);
                                wa[k] = rdiag[k];
                            }
                        }
                    }
                }
            }
            rdiag[j] = -ajnorm;
        }
    }

    void MinPack::qrsolv(int n, std::vector<double> &r, int ldr,
                              std::vector<int> &ipvt, std::vector<double> &diag, std::vector<double> &qtb, std::vector<double> &x,
                              std::vector<double> &sdiag, std::vector<double> &wa) {
        
        const double p5 = .5;
        const double p25 = .25;

        int i, j, k, l;
        double cos, sin, sum, temp;
        int nsing;
        double qtbpj;

        for (j = 0; j < n; ++j) {
            for (i = j; i < n; ++i) {
                r[i + j * ldr] = r[j + i * ldr];
            }
            x[j] = r[j + j * ldr];
            wa[j] = qtb[j];
        }

        for (j = 0; j < n; ++j) {

            l = ipvt[j] - 1;
            if (diag[l] != 0.) {
                for (k = j; k < n; ++k) {
                    sdiag[k] = 0.;
                }
                sdiag[j] = diag[l];

                qtbpj = 0.;
                for (k = j; k < n; ++k) {

                    if (sdiag[k] != 0.) {
                        if (std::abs(r[k + k * ldr]) < std::abs(sdiag[k])) {
                            double cotan;
                            cotan = r[k + k * ldr] / sdiag[k];
                            sin = p5 / std::sqrt(p25 + p25 * (cotan * cotan));
                            cos = sin * cotan;
                        } else {
                            double tan;
                            tan = sdiag[k] / r[k + k * ldr];
                            cos = p5 / std::sqrt(p25 + p25 * (tan * tan));
                            sin = cos * tan;
                        }

                        temp = cos * wa[k] + sin * qtbpj;
                        qtbpj = -sin * wa[k] + cos * qtbpj;
                        wa[k] = temp;

                        r[k + k * ldr] = cos * r[k + k * ldr] + sin * sdiag[k];
                        if (n > k + 1) {
                            for (i = k + 1; i < n; ++i) {
                                temp = cos * r[i + k * ldr] + sin * sdiag[i];
                                sdiag[i] = -sin * r[i + k * ldr] + cos * sdiag[i];
                                r[i + k * ldr] = temp;
                            }
                        }
                    }
                }
            }

            sdiag[j] = r[j + j * ldr];
            r[j + j * ldr] = x[j];
        }

        nsing = n;
        for (j = 0; j < n; ++j) {
            if (sdiag[j] == 0. && nsing == n) {
                nsing = j;
            }
            if (nsing < n) {
                wa[j] = 0.;
            }
        }
        if (nsing >= 1) {
            for (k = 1; k <= nsing; ++k) {
                j = nsing - k;
                sum = 0.;
                if (nsing > j + 1) {
                    for (i = j + 1; i < nsing; ++i) {
                        sum += r[i + j * ldr] * wa[i];
                    }
                }
                wa[j] = (wa[j] - sum) / sdiag[j];
            }
        }

        for (j = 0; j < n; ++j) {
            l = ipvt[j] - 1;
            x[l] = wa[j];
        }
    }

    void MinPack::lmpar(int n, std::vector<double> &r, int ldr,
                             std::vector<int> &ipvt, std::vector<double> &diag, std::vector<double> &qtb, double delta,
                             std::vector<double> &par, std::vector<double> &x, std::vector<double> &sdiag, std::vector<double> &wa1,
                             std::vector<double> &wa2) {
        
        const double p1 = .1;
        const double p001 = .001;

        double d1, d2;

        int j, l;
        double fp;
        double parc, parl;
        int iter;
        double temp, paru, dwarf;
        int nsing;
        double gnorm;
        double dxnorm;

        dwarf = dpmpar(2);

        nsing = n;
        for (j = 0; j < n; ++j) {
            wa1[j] = qtb[j];
            if (r[j + j * ldr] == 0. && nsing == n) {
                nsing = j;
            }
            if (nsing < n) {
                wa1[j] = 0.;
            }
        }
        if (nsing >= 1) {
            int k;
            for (k = 1; k <= nsing; ++k) {
                j = nsing - k;
                wa1[j] /= r[j + j * ldr];
                temp = wa1[j];
                if (j >= 1) {
                    int i;
                    for (i = 0; i < j; ++i) {
                        wa1[i] -= r[i + j * ldr] * temp;
                    }
                }
            }
        }
        for (j = 0; j < n; ++j) {
            l = ipvt[j] - 1;
            x[l] = wa1[j];
        }

        iter = 0;
        for (j = 0; j < n; ++j) {
            wa2[j] = diag[j] * x[j];
        }
        dxnorm = enorm(n, 0, wa2);
        fp = dxnorm - delta;

        {

            if (fp <= p1 * delta) {
                goto terminate;
            }

            parl = 0.;
            if (nsing >= n) {
                for (j = 0; j < n; ++j) {
                    l = ipvt[j] - 1;
                    wa1[j] = diag[l] * (wa2[l] / dxnorm);
                }
                for (j = 0; j < n; ++j) {
                    double sum = 0.;
                    if (j >= 1) {
                        int i;
                        for (i = 0; i < j; ++i) {
                            sum += r[i + j * ldr] * wa1[i];
                        }
                    }
                    wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
                }
                temp = enorm(n, 0, wa1);
                parl = fp / delta / temp / temp;
            }

            for (j = 0; j < n; ++j) {
                double sum;
                int i;
                sum = 0.;
                for (i = 0; i <= j; ++i) {
                    sum += r[i + j * ldr] * qtb[i];
                }
                l = ipvt[j] - 1;
                wa1[j] = sum / diag[l];
            }
            gnorm = enorm(n, 0, wa1);
            paru = gnorm / delta;
            if (paru == 0.) {
                paru = dwarf / std::min(delta, p1) ;
            }

            par[0] = std::max(par[0], parl);
            par[0] = std::min(par[0], paru);
            if (par[0] == 0.) {
                par[0] = gnorm / dxnorm;
            }

            for (; ; ) {
                ++iter;

                if (par[0] == 0.) {
                    
                    d1 = dwarf;
                    d2 = p001 * paru;
                    par[0] = std::max(d1, d2);
                }
                temp = std::sqrt(par[0]);
                for (j = 0; j < n; ++j) {
                    wa1[j] = temp * diag[j];
                }
                qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
                for (j = 0; j < n; ++j) {
                    wa2[j] = diag[j] * x[j];
                }
                dxnorm = enorm(n, 0, wa2);
                temp = fp;
                fp = dxnorm - delta;

                if (std::abs(fp) <= p1 * delta || (parl == 0. && fp <= temp && temp < 0.) || iter == 10) {
                    goto terminate;
                }

                for (j = 0; j < n; ++j) {
                    l = ipvt[j] - 1;
                    wa1[j] = diag[l] * (wa2[l] / dxnorm);
                }
                for (j = 0; j < n; ++j) {
                    wa1[j] /= sdiag[j];
                    temp = wa1[j];
                    if (n > j + 1) {
                        int i;
                        for (i = j + 1; i < n; ++i) {
                            wa1[i] -= r[i + j * ldr] * temp;
                        }
                    }
                }
                temp = enorm(n, 0, wa1);
                parc = fp / delta / temp / temp;

                if (fp > 0.) {
                    parl = std::max(parl, par[0]);
                }
                if (fp < 0.) {
                    paru = std::min(paru, par[0]);
                }

                d1 = parl;
                d2 = par[0] + parc;
                par[0] = std::max(d1, d2);

            }
        }
        terminate:;


        if (iter == 0) {
            par[0] = 0.;
        }
    }

    int MinPack::fdjac2(Lmder_Function &fcn_mn, int m, int n, std::vector<double> &x,
                             std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                             double epsfcn, std::vector<double> &wa)
    {
        
        double h;
        int i, j;
        double eps, temp, epsmch;
        int iflag;

        epsmch = dpmpar(1);

        eps = std::sqrt((std::max(epsfcn,epsmch)));
        for (j = 0; j < n; ++j) {
            temp = x[j];
            h = eps * std::abs(temp);
            if (h == 0.) {
                h = eps;
            }
            x[j] = temp + h;
        
            iflag = fcn_mn.apply(m, n, x, wa, 2);
            if (iflag < 0) {
                return iflag;
            }
            x[j] = temp;
            for (i = 0; i < m; ++i) {
                fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
            }
        }
        return 0;
    }

    void MinPack::r1mpyq(int m, int n, std::vector<double> &a, int lda, std::vector<double> &v, std::vector<double> &w) {
        double cos,sin,temp;
        int i, j, nmj, nm1;

        nm1 = n - 1;
        if (nm1 < 1)
            return;
        for (nmj = 1; nmj <= nm1; nmj++) {
            j = n - nmj;
            if (std::abs(v[j-1]) > 1.0) {
                cos = 1.0 / v[j-1];
                sin = std::sqrt(1.0 - cos * cos);
            } else {
                sin = v[j-1];
                cos = std::sqrt(1.0 - sin * sin);
            }
            for (i = 1; i <= m; i++) {
                temp = cos * a[(i-1) + (j-1) * lda] - sin * a[(i-1) + (n - 1) * lda];
                a[(i-1) + (n - 1) * lda] = sin * a[(i-1) + (j-1) * lda] + cos * a[(i-1) + (n - 1) * lda];
                a[(i-1) + (j-1) * lda] = temp;
            }
        }

        for (j = 1; j <= nm1; j++) {
            if (std::abs(w[j-1]) > 1.0) {
                cos = 1.0 / w[j-1];
                sin = std::sqrt(1.0 - cos * cos);
            } else {
                sin = w[j-1];
                cos = std::sqrt(1.0 - sin * sin);
            }
            for (i = 1; i <= m; i++) {
                temp = cos * a[(i-1) + (j-1) * lda] + sin * a[(i-1) + (n - 1) * lda];
                a[(i-1) + (n - 1) * lda] = -sin * a[(i-1) + (j-1) * lda] + cos * a[(i-1) + (n - 1) * lda];
                a[(i-1) + (j-1) * lda] = temp;
            }
        }
    }

    bool MinPack::r1updt(int m, int n, std::vector<double> &s, int ls, std::vector<double> &u, std::vector<double> &v,
                                 std::vector<double> &w) {
        const double p25 = 0.25;
        const double p5 = 0.5;

        int i,j,jj,l,nm1;
        double cos,cotan,giant;
        double sin,tan,tau,temp;

        bool sing;

        giant = dpmpar(3);

        jj = (n * (2 * m - n + 1)) / 2 - (m - n);

        l = jj;
        for (i = n; i <= m; i++) {
            w[i - 1] = s[l - 1];
            l = l + 1;
        }

        nm1 = n - 1;
        if (nm1 >= 1) {
            for (int nmj = 1; nmj <= nm1; nmj++) {
                j = n - nmj;
                jj = jj - (m - j + 1);
                w[j - 1] = 0.0;

                if (v[j - 1] != 0.0) {

                    if (!(std::abs(v[n - 1]) >= std::abs(v[j - 1]))) {
                        cotan = v[n - 1] / v[j - 1];
                        sin = p5 / std::sqrt(p25 + p25 * cotan * cotan);
                        cos = sin * cotan;
                        tau = 1.0;
                        if (std::abs(cos) * giant > 1.0) {
                            tau = 1.0 / cos;
                        }
                    } else {
                        tan = v[j - 1] / v[n - 1];
                        cos = p5 / std::sqrt(p25 + p25 * tan * tan);
                        sin = cos * tan;
                        tau = sin;
                    }

                    v[n - 1] = sin * v[j - 1] + cos * v[n - 1];
                    v[j - 1] = tau;

                    l = jj;
                    for (i = j; i <= m; i++) {
                        temp = cos * s[l - 1] - sin * w[i - 1];
                        w[i - 1] = sin * s[l - 1] + cos * w[i - 1];
                        s[l - 1] = temp;
                        l = l + 1;
                    }
                }
            }
        }
        
        for (i = 1; i <= m; i++) {
            w[i - 1] = w[i - 1] + v[n - 1] * u[i - 1];
        }
        
        sing = false;
        if (nm1 >= 1) {
            for (j = 1; j <= nm1; j++) {

                if (w[j - 1] != 0.0) {

                    if (!(std::abs(s[jj - 1]) >= std::abs(w[j - 1]))) {
                        cotan = s[jj - 1] / w[j - 1];
                        sin = p5 / std::sqrt(p25 + p25 * cotan * cotan);
                        cos = sin * cotan;
                        tau = 1.0;
                        if (std::abs(cos) * giant > 1.0) {
                            tau = 1.0 / cos;
                        }
                    } else {
                        tan = w[j - 1] / s[jj - 1];
                        cos = p5 / std::sqrt(p25 + p25 * tan * tan);
                        sin = cos * tan;
                        tau = sin;
                    }

                    l = jj;
                    for (i = j; i <= m; i++) {
                        temp = cos * s[l - 1] + sin * w[i - 1];
                        w[i - 1] = -sin * s[l - 1] + cos * w[i - 1];
                        s[l - 1] = temp;
                        l = l + 1;
                    }
                    
                    w[j - 1] = tau;
                }
                
                if (s[jj - 1] == 0.0) {
                    sing = true;
                }
                jj = jj + (m - j + 1);
            }
        }
        
        l = jj;
        for (i = n; i <= m; i++) {
            s[l - 1] = w[i - 1];
            l = l + 1;
        }
        if (s[jj - 1] == 0.0) {
            sing = true;
        }
        return sing;
    }

    void MinPack::dogleg(int n, std::vector<double> &r, int lr, std::vector<double> &diag, std::vector<double> &qtb,
                              double delta, std::vector<double> &x, std::vector<double> &wa1, std::vector<double> &wa2) {

        int i,j,jj,jp1,k,l;
        double alpha,bnorm,epsmch,gnorm,qnorm,sgnorm,sum,temp;
        
        epsmch = dpmpar(1);
        
        jj = (n * (n + 1)) / 2 + 1;

        for (k = 1; k <= n; k++) {
            j = n - k + 1;
            jp1 = j + 1;
            jj = jj - k;
            l = jj + 1;
            sum = 0.0;
            if (n >= jp1) {
                for (i = jp1; i <= n; i++) {
                    sum = sum + r[l - 1] * x[i - 1];
                    l = l + 1;
                }
            }
            temp = r[jj - 1];
            if (temp == 0.0) {
                l = j;
                for (i = 1; i <= j; i++) {
                    temp = std::max(temp, std::abs(r[l - 1]));
                    l = l + n - i;
                }
                temp = epsmch * temp;
                if (temp == 0.0) {
                    temp = epsmch;
                }
            }
            x[j - 1] = (qtb[j - 1] - sum) / temp;
        }
        
        for (j = 1; j <= n; j++) {
            wa1[j-1] = 0.0;
            wa2[j-1] = diag[j-1] * x[j-1];
        }
        qnorm = enorm(n, 0, wa2);

        if (qnorm <= delta) {
            return;
        }
        
        l = 1;
        for (j = 1; j <= n; j++) {
            temp = qtb[j-1];
            for (i = j; i <= n; i++) {
                wa1[i-1] = wa1[i-1] + r[l-1] * temp;
                l = l + 1;
            }
            wa1[j-1] = wa1[j-1] / diag[j-1];
        }

        gnorm = enorm(n, 0, wa1);
        sgnorm = 0.0;
        alpha = delta / qnorm;

        if (gnorm != 0.0) {

            for (j = 1; j <= n; j++) {
                wa1[j-1] = (wa1[j-1] / gnorm) / diag[j-1];
            }
            l = 1;
            for (j = 1; j <= n; j++) {
                sum = 0.0;
                for (i = j; i <= n; i++) {
                    sum = sum + r[l-1] * wa1[i-1];
                    l = l + 1;
                }
                wa2[j-1] = sum;
            }
            temp = enorm(n, 0, wa2);
            
            if (temp == 0.)
                sgnorm = delta;
            else
                sgnorm = (gnorm / temp) / temp;

            alpha = 0.0;

            if (sgnorm < delta) {

                bnorm = enorm(n, 0, qtb);
                temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
                temp = temp - (delta / qnorm) * (sgnorm / delta) * (sgnorm / delta)
                        + std::sqrt((temp - (delta / qnorm))*(temp - (delta / qnorm))
                        + (1.0 - (delta / qnorm) * (delta / qnorm))
                        * (1.0 - (sgnorm / delta) * (sgnorm / delta)));
                alpha = ((delta / qnorm)
                        * (1.0 - (sgnorm / delta) * (sgnorm / delta))) / temp;
            }
        }
        
        temp = (1.0 - alpha) * std::min(sgnorm, delta);
        for (j = 1; j <= n; j++) {
            x[j-1] = temp * wa1[j-1] + alpha * x[j-1];
        }
    }

    void MinPack::qform(int m, int n, std::vector<double> &q, int ldq, std::vector<double> &wa) {
        int i,j,k,minmn;
        double sum,temp;
        
        minmn = std::min(m, n);
        if (minmn >= 2) {
            for (j = 2; j <= minmn; j++) {
                int jm1 = j - 1;
                for (i = 1; i <= jm1; i++) {
                    q[(i-1) + (j-1) * ldq] = 0.0;
                }
            }
        }
        
        int np1 = n + 1;
        if (m >= np1) {
            for (j = np1; j <= m; j++) {
                for (i = 1; i <= m; i++) {
                    q[(i-1) + (j-1) * ldq] = 0.0;
                }
                q[(j-1) + (j-1) * ldq] = 1.0;
            }
        }

        for (int l = 1; l <= minmn; l++) {
            k = minmn - l + 1;
            for (i = k; i <= m; i++) {
                wa[i-1] = q[(i-1) + (k-1) * ldq];
                q[(i-1) + (k-1) * ldq] = 0.0;
            }
            q[(k-1) + (k-1) * ldq] = 1.0;

            if (wa[k-1] != 0.0) {
                for (j = k; j <= m; j++) {
                    sum = 0.0;
                    for (i = k; i <= m; i++) {
                        sum = sum + q[(i-1) + (j-1) * ldq] * wa[i-1];
                    }
                    temp = sum / wa[k-1];
                    for (i = k; i <= m; i++) {
                        q[(i-1) + (j-1) * ldq] = q[(i-1) + (j-1) * ldq] - temp * wa[i-1];
                    }
                }
            }
        }
    }

    
    void MinPack::fdjac1(Hybrd_Function &fcn,
                              int n, std::vector<double> &x, std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac, std::vector<int> &iflag,
                              int ml, int mu, double epsfcn, std::vector<double> &wa1, std::vector<double> &wa2) {
        int i,j,k,msum;
        double eps,epsmch,h,temp;

        epsmch = dpmpar(1);

        eps = std::sqrt(std::max(epsfcn, epsmch));
        msum = ml + mu + 1;

        if (msum >= n) {
            for (j = 1; j <= n; j++) {
                temp = x[j-1];
                h = eps * std::abs(temp);
                if (h == 0.0) {
                    h = eps;
                }
                x[j-1] = temp + h;
                fcn.apply(n, x, wa1, iflag);
                if (iflag[0] < 0) {
                    break;
                }
                x[j-1] = temp;
                for (i = 1; i <= n; i++) {
                    fjac[(i-1) + (j-1) * ldfjac] = (wa1[i-1] - fvec[i-1]) / h;
                }
            }
            return;
        }

        for (k = 1; k <= msum; k++) {
            for (j = k; j <= n; j = j + msum) {
                wa2[j-1] = x[j-1];
                h = eps * std::abs(wa2[j-1]);
                if (h == 0.0) {
                    h = eps;
                }
                x[j-1] = wa2[j-1] + h;
            }
            fcn.apply(n, x, wa1, iflag);
            if (iflag[0] < 0) {
                break;
            }
            for (j = k; j <= n; j = j + msum) {
                x[j-1] = wa2[j-1];
                h = eps * std::abs(wa2[j-1]);
                if (h == 0.0) {
                    h = eps;
                }
                for (i = 1; i <= n; i++) {
                    fjac[(i-1) + (j-1) * ldfjac] = 0.0;
                    if (i >= j - mu && i <= j + ml) {
                        fjac[(i-1) + (j-1) * ldfjac] = (wa1[i-1] - fvec[i-1]) / h;
                    }
                }
            }
        }
    }

    int MinPack::hybrd(Hybrd_Function &fcn,
                            int n, std::vector<double> &x,
                            std::vector<double> &fvec, double xtol, int maxfev, int ml, int mu, double epsfcn,
                            std::vector<double> &diag, int mode, double factor, int nprint, std::vector<int> &nfev,
                            std::vector<double> &fjac, int ldfjac, std::vector<double> &r, int lr, std::vector<double> &qtf, std::vector<double> &wa1,
                            std::vector<double> &wa2, std::vector<double> &wa3, std::vector<double> &wa4) {
        int i,iflag,info,iter,j,l,msum,ncfail,ncsuc,nslow1,nslow2;
        std::vector<int> iwa(1, 0);
        bool jeval;
        const double p001 = 0.001;
        const double p0001 = 0.0001;
        const double p1 = 0.1;
        const double p5 = 0.5;
        double actred,delta,epsmch,fnorm,fnorm1,pnorm,
                prered,ratio,sum,temp,xnorm;
        std::vector<int> iflag_(1, 0);

        epsmch = dpmpar(1);

        info = 0;
        iflag = 0;
        nfev[0] = 0;

        if (n <= 0
            ||xtol < 0.0
            ||maxfev <= 0
            ||ml < 0
            ||mu < 0
            ||factor <= 0.0
            ||ldfjac < n
            ||lr < (n * (n + 1)) / 2) {
            info = 0;
            return info;
        }
        if (mode == 2) {
            for (j = 1; j <= n; j++) {
                if (diag[j-1] <= 0.0) {
                    info = 0;
                    return info;
                }
            }
        }

        iflag = 1;
        iflag_[0] = iflag;
        fcn.apply(n, x, fvec, iflag_);
        iflag = iflag_[0];
        nfev[0] = 1;
        if (iflag < 0) {
            info = iflag;
            return info;
        }

        fnorm = enorm(n, 0, fvec);

        msum = std::min(ml + mu + 1, n);

        iter = 1;
        ncsuc = 0;
        ncfail = 0;
        nslow1 = 0;
        nslow2 = 0;
        delta = 0;
        xnorm = 0;

        for (;;) {
            jeval = true;

            iflag = 2;
            iflag_[0] = iflag;
            fdjac1(fcn, n, x, fvec, fjac, ldfjac, iflag_, ml, mu, epsfcn, wa1, wa2);
            iflag = iflag_[0];

            nfev[0] = nfev[0] + msum;
            if (iflag < 0) {
                goto outerloop_end;
            }

            qrfac(n, n, fjac, ldfjac, 0, iwa, 1, wa1, wa2, wa3);

            if (iter == 1) {
                if (mode != 2) {
                    for (j = 1; j <= n; j++) {
                        diag[j - 1] = wa2[j - 1];
                        if (wa2[j - 1] == 0.0) {
                            diag[j - 1] = 1.0;
                        }
                    }
                }

                for (j = 1; j <= n; j++) {
                    wa3[j - 1] = diag[j - 1] * x[j - 1];
                }
                xnorm = enorm(n, 0, wa3);
                delta = factor * xnorm;
                if (delta == 0.0) delta = factor;
            }

            for (i = 1; i <= n; i++) {
                qtf[i - 1] = fvec[i - 1];
            }
            for (j = 1; j <= n; j++) {
                if (fjac[(j - 1) + (j - 1) * ldfjac] != 0.0) {
                    sum = 0.0;
                    for (i = j; i <= n; i++) {
                        sum = sum + fjac[(i - 1) + (j - 1) * ldfjac] * qtf[i - 1];
                    }
                    temp = -sum / fjac[(j - 1) + (j - 1) * ldfjac];
                    for (i = j; i <= n; i++) {
                        qtf[i - 1] = qtf[i - 1] + fjac[(i - 1) + (j - 1) * ldfjac] * temp;
                    }
                }
            }
            
            for (j = 1; j <= n; j++) {
                l = j;
                int jm1 = j - 1;
                if (jm1 >= 1) {
                    for (i = 1; i <= jm1; i++) {
                        r[l - 1] = fjac[(i - 1) + (j - 1) * ldfjac];
                        l = l + n - i;
                    }
                }
                r[l - 1] = wa1[j - 1];
                if (wa1[j - 1] == 0.0) {
                    std::fprintf(stderr, "  Matrix is singular.\n");
                }
            }

            qform(n, n, fjac, ldfjac, wa1);

            if (mode != 2) {
                for (j = 1; j <= n; j++) {
                    diag[j - 1] = std::max(diag[j - 1], wa2[j - 1]);
                }
            }

            for (; ; ) {
                
                if (nprint > 0) {
                    if ((iter - 1) % nprint == 0) {
                        iflag = 0;
                        iflag_[0] = iflag;
                        fcn.apply(n, x, fvec, iflag_);
                        iflag = iflag_[0];
                        if (iflag < 0) {
                            goto outerloop_end;
                        }
                    }
                }

                dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);

                for (j = 1; j <= n; j++) {
                    wa1[j - 1] = -wa1[j - 1];
                    wa2[j - 1] = x[j - 1] + wa1[j - 1];
                    wa3[j - 1] = diag[j - 1] * wa1[j - 1];
                }
                pnorm = enorm(n, 0, wa3);

                if (iter == 1) {
                    delta = std::min(delta, pnorm);
                }

                iflag = 1;
                iflag_[0] = iflag;
                fcn.apply(n, wa2, wa4, iflag_);
                iflag = iflag_[0];
                nfev[0] = nfev[0] + 1;
                if (iflag < 0) {
                    goto outerloop_end;
                }
                fnorm1 = enorm(n, 0, wa4);

                actred = -1.0;
                if (fnorm1 < fnorm) {
                    actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);
                }

                l = 1;
                for (i = 1; i <= n; i++) {
                    sum = 0.0;
                    for (j = i; j <= n; j++) {
                        sum = sum + r[l - 1] * wa1[j - 1];
                        l = l + 1;
                    }
                    wa3[i - 1] = qtf[i - 1] + sum;
                }
                temp = enorm(n, 0, wa3);
                prered = 0.0;
                if (temp < fnorm) {
                    prered = 1.0 - (temp / fnorm) * (temp / fnorm);
                }

                ratio = 0.0;
                if (prered > 0.0) {
                    ratio = actred / prered;
                }

                if (!(ratio >= p1)) {
                    ncsuc = 0;
                    ncfail = ncfail + 1;
                    delta = p5 * delta;
                } else {
                    ncfail = 0;
                    ncsuc = ncsuc + 1;
                    if (ratio >= p5 || ncsuc > 1) {
                        delta = std::max(delta, pnorm / p5);
                    }
                    if (std::abs(ratio - 1.0) <= p1) {
                        delta = pnorm / p5;
                    }
                }
                
                if (!(ratio < p0001)) {

                    for (j = 1; j <= n; j++) {
                        x[j - 1] = wa2[j - 1];
                        wa2[j - 1] = diag[j - 1] * x[j - 1];
                        fvec[j - 1] = wa4[j - 1];
                    }
                    xnorm = enorm(n, 0, wa2);
                    fnorm = fnorm1;
                    iter = iter + 1;
                }

                nslow1 = nslow1 + 1;
                if (actred >= p001) nslow1 = 0;
                if (jeval) nslow2 = nslow2 + 1;
                if (actred >= p1) nslow2 = 0;
                
                if (delta <= xtol * xnorm || fnorm == 0.0) info = 1;
                if (info != 0) goto outerloop_end;

                if (nfev[0] >= maxfev) info = 2;
                if (p1 * std::max(p1 * delta, pnorm) <= epsmch * xnorm) info = 3;
                if (nslow2 == 5) info = 4;
                if (nslow1 == 10) info = 5;
                if (info != 0) goto outerloop_end;

                if (ncfail == 2) goto innerloop_end;

                for (j = 1; j <= n; j++) {
                    sum = 0.0;
                    for (i = 1; i <= n; i++) {
                        sum = sum + fjac[(i - 1) + (j - 1) * ldfjac] * wa4[i - 1];
                    }
                    wa2[j - 1] = (sum - wa3[j - 1]) / pnorm;
                    wa1[j - 1] = diag[j - 1] * ((diag[j - 1] * wa1[j - 1]) / pnorm);
                    if (ratio >= p0001) {
                        qtf[j - 1] = sum;
                    }
                }
                
                r1updt(n, n, r, lr, wa1, wa2, wa3);
                r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
                r1mpyq(1, n, qtf, 1, wa2, wa3);

                jeval = false;
            }
        innerloop_end:;
            
        }
        outerloop_end:;

        if (iflag < 0) {
            info = iflag;
        }
        if (nprint > 0) {
            iflag_[0] = 0;
            fcn.apply(n, x, fvec, iflag_);
        }
        return info;

    }

    int MinPack::hybrd1(Hybrd_Function &fcn, int n,
                              std::vector<double> &x, std::vector<double> &fvec, double tol, std::vector<double> &wa, int lwa) {
        // Setting epsfcn to 0 means default machine precision will be used
        return hybrd1(fcn, n, x, fvec, tol, wa, lwa, 0.0);
    }

    int MinPack::hybrd1(Hybrd_Function &fcn, int n,
                              std::vector<double> &x, std::vector<double> &fvec, double tol, std::vector<double> &wa, int lwa, double epsfcn )
    {
        int info,j,lr,maxfev,ml,mode,mu,nprint;
        std::vector<int> nfev(1, 0);
        double factor,xtol;

        info = 0;

        if ( n <= 0 || tol <= 0.0 || lwa < ( n * ( 3 * n + 13 ) ) / 2 )
        {
            return info;
        }

        maxfev = 200 * ( n + 1 );
        xtol = tol;
        ml = n - 1;
        mu = n - 1;
        mode = 2;
        for ( j = 0; j < n; j++ )
        {
            wa[j] = 1.0;
        }
        factor = 100.0;
        nprint = 0;
        lr = ( n * ( n + 1 ) ) / 2;

        std::vector<double> fjac(n*n, 0.0);
        std::vector<double> r(lr, 0.0);
        std::vector<double> qtf(n, 0.0);
        std::vector<double> wa1(n, 0.0);
        std::vector<double> wa2(n, 0.0);
        std::vector<double> wa3(n, 0.0);
        std::vector<double> wa4(n, 0.0);
        info = hybrd ( fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn, wa, mode,
                factor, nprint, nfev, fjac, n, r, lr,
                qtf, wa1, wa2, wa3, wa4 );

        if ( info == 5 )
        {
            info = 4;
        }
        return info;
    }

    void MinPack::chkder(int m, int n, std::vector<double> &x, std::vector<double> &fvec, std::vector<double> &fjac,
                int ldfjac, std::vector<double> &xp, std::vector<double> &fvecp, int mode, std::vector<double> &err) {
        
        int i, j;
        double eps, epsf, temp, epsmch;
        double epslog;

        const double factor = 100.0;

        epsmch = dpmpar(1);

        eps = std::sqrt(epsmch);

        if (mode == 1) {
            for (j = 0; j < n; j++) {
                if (x[j] == 0.0) {
                    temp = eps;
                } else {
                    temp = eps * std::abs(x[j]);
                }
                xp[j] = x[j] + temp;
            }
            return;
        }
        
        epsf = factor * epsmch;
        epslog = std::log10(eps);
        for (i = 0; i < m; i++) {
            err[i] = 0.;
        }

        for (j = 0; j < n; j++) {
            if (x[j] == 0.0) {
                temp = 1.;
            } else {
                temp = std::abs(x[j]);
            }
            for (i = 0; i < m; i++) {
                err[i] = err[i] + temp * fjac[i + j * ldfjac];
            }
        }

        for (i = 0; i < m; i++) {
            temp = 1.;
            if (fvec[i] != 0.0 &&
                    fvecp[i] != 0.0 &&
                    epsf * std::abs(fvec[i]) <= std::abs(fvecp[i] - fvec[i])) {
                temp = eps * std::abs((fvecp[i] - fvec[i]) / eps - err[i])
                        / (std::abs(fvec[i]) + std::abs(fvecp[i]));

                if (temp <= epsmch) {
                    err[i] = 1.;
                } else if (temp < eps) {
                    err[i] = (std::log10(temp) - epslog) / epslog;
                } else {
                    err[i] = 0.0;
                }
            }
        }
    }
} // namespace redukti::mathlib
