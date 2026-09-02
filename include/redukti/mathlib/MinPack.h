// C++ port of org.redukti.mathlib.MinPack
//
// argonne national laboratory. minpack project. march 1980.
// burton s. garbow, kenneth e. hillstrom, jorge j. more
#ifndef REDUKTI_MATHLIB_MINPACK_H
#define REDUKTI_MATHLIB_MINPACK_H

#include <vector>

namespace redukti::mathlib {

/**
 * Callback for lmder1 / lmder.
 *
 * Two overloads, exactly as in the Java interface: implementations that can
 * supply their own jacobian override the seven-argument form and return true
 * from hasJacobian(); implementations that cannot override the five-argument
 * form and leave hasJacobian() false, in which case lmder computes the
 * jacobian by finite differences via fdjac2.
 *
 * NOTE: C++ name hiding means a subclass that overrides only one overload
 * hides the other. Derived classes should carry `using Lmder_Function::apply;`.
 */
class Lmder_Function {
public:
    virtual ~Lmder_Function() = default;

    /**
     * @param iflag if 1, calculate the functions at x and return them in fvec,
     *              do not alter fjac; if 2, calculate the jacobian at x and
     *              return it in fjac, do not alter fvec.
     * @return a negative value to terminate lmder1/lmder
     */
    virtual int apply(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
                      std::vector<double> &fjac, int ldfjac, int iflag) {
        (void)m; (void)n; (void)x; (void)fvec; (void)fjac; (void)ldfjac; (void)iflag;
        return -1;
    }

    /**
     * @param iflag if 1, calculate the functions at x and return them in fvec.
     * @return a negative value to terminate lmder1/lmder
     */
    virtual int apply(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
                      int iflag) {
        (void)m; (void)n; (void)x; (void)fvec; (void)iflag;
        return -1;
    }

    virtual bool hasJacobian() { return false; }
};

/** Callback for hybrd / hybrd1. */
class Hybrd_Function {
public:
    virtual ~Hybrd_Function() = default;
    virtual void apply(int n, std::vector<double> &x, std::vector<double> &fvec,
                       std::vector<int> &iflag) = 0;
};

class MinPack {
public:
    /**
     * Double precision machine parameters.
     * dpmpar(1) = machine precision, dpmpar(2) = smallest magnitude,
     * dpmpar(3) = largest magnitude.
     */
    static double dpmpar(int i);

    static int lmder1(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                      std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                      double tol, std::vector<int> &ipvt, std::vector<double> &wa,
                      int lwa);

    static int lmder1(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                      std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                      double tol, std::vector<int> &ipvt, std::vector<double> &wa,
                      int lwa, double epsfcn);

    static int lmder(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                     std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                     double ftol, double xtol, double gtol, int maxfev,
                     std::vector<double> &diag, int mode, double factor, int nprint,
                     std::vector<int> &nfev, std::vector<int> &njev,
                     std::vector<int> &ipvt, std::vector<double> &qtf,
                     std::vector<double> &wa1, std::vector<double> &wa2,
                     std::vector<double> &wa3, std::vector<double> &wa4);

    static int lmder(Lmder_Function &fcn, int m, int n, std::vector<double> &x,
                     std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                     double ftol, double xtol, double gtol, int maxfev,
                     std::vector<double> &diag, int mode, double factor, int nprint,
                     std::vector<int> &nfev, std::vector<int> &njev,
                     std::vector<int> &ipvt, std::vector<double> &qtf,
                     std::vector<double> &wa1, std::vector<double> &wa2,
                     std::vector<double> &wa3, std::vector<double> &wa4, double epsfcn);

    static double enorm(int n, int start, std::vector<double> &x);

    static void qrfac(int m, int n, std::vector<double> &a, int lda, int pivot,
                      std::vector<int> &ipvt, int lipvt, std::vector<double> &rdiag,
                      std::vector<double> &acnorm, std::vector<double> &wa);

    static void qrsolv(int n, std::vector<double> &r, int ldr, std::vector<int> &ipvt,
                       std::vector<double> &diag, std::vector<double> &qtb,
                       std::vector<double> &x, std::vector<double> &sdiag,
                       std::vector<double> &wa);

    static void lmpar(int n, std::vector<double> &r, int ldr, std::vector<int> &ipvt,
                      std::vector<double> &diag, std::vector<double> &qtb, double delta,
                      std::vector<double> &par, std::vector<double> &x,
                      std::vector<double> &sdiag, std::vector<double> &wa1,
                      std::vector<double> &wa2);

    static int fdjac2(Lmder_Function &fcn_mn, int m, int n, std::vector<double> &x,
                      std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                      double epsfcn, std::vector<double> &wa);

    static void r1mpyq(int m, int n, std::vector<double> &a, int lda,
                       std::vector<double> &v, std::vector<double> &w);

    static bool r1updt(int m, int n, std::vector<double> &s, int ls,
                       std::vector<double> &u, std::vector<double> &v,
                       std::vector<double> &w);

    static void dogleg(int n, std::vector<double> &r, int lr, std::vector<double> &diag,
                       std::vector<double> &qtb, double delta, std::vector<double> &x,
                       std::vector<double> &wa1, std::vector<double> &wa2);

    static void qform(int m, int n, std::vector<double> &q, int ldq,
                      std::vector<double> &wa);

    static void fdjac1(Hybrd_Function &fcn, int n, std::vector<double> &x,
                       std::vector<double> &fvec, std::vector<double> &fjac, int ldfjac,
                       std::vector<int> &iflag, int ml, int mu, double epsfcn,
                       std::vector<double> &wa1, std::vector<double> &wa2);

    static int hybrd(Hybrd_Function &fcn, int n, std::vector<double> &x,
                     std::vector<double> &fvec, double xtol, int maxfev, int ml, int mu,
                     double epsfcn, std::vector<double> &diag, int mode, double factor,
                     int nprint, std::vector<int> &nfev, std::vector<double> &fjac,
                     int ldfjac, std::vector<double> &r, int lr,
                     std::vector<double> &qtf, std::vector<double> &wa1,
                     std::vector<double> &wa2, std::vector<double> &wa3,
                     std::vector<double> &wa4);

    static int hybrd1(Hybrd_Function &fcn, int n, std::vector<double> &x,
                      std::vector<double> &fvec, double tol, std::vector<double> &wa,
                      int lwa);

    static int hybrd1(Hybrd_Function &fcn, int n, std::vector<double> &x,
                      std::vector<double> &fvec, double tol, std::vector<double> &wa,
                      int lwa, double epsfcn);

    static void chkder(int m, int n, std::vector<double> &x, std::vector<double> &fvec,
                       std::vector<double> &fjac, int ldfjac, std::vector<double> &xp,
                       std::vector<double> &fvecp, int mode, std::vector<double> &err);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_MINPACK_H
