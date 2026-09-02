// C++ port of org.redukti.mathlib.LMLSolver
#include "redukti/mathlib/LMLSolver.h"

#include <cmath>
#include <cstdio>

namespace redukti::mathlib {

const double LMLSolver::BIGVAL = 9.876543e+99;

LMLSolver::LMLSolver(LMLFunction &gH, double gtol, int gnparms, int gnpts)
// Constructor sets up private fields, including host for callbacks.
{
    myH = &gH;
    lmtol = gtol;
    nparms = gnparms;
    npts = gnpts;
    niter = 0;
    delta.assign(static_cast<std::size_t>(nparms), 0.0);
    beta.assign(static_cast<std::size_t>(nparms), 0.0);
    alpha.assign(static_cast<std::size_t>(nparms),
                 std::vector<double>(static_cast<std::size_t>(nparms), 0.0));
    amatrix.assign(static_cast<std::size_t>(nparms),
                   std::vector<double>(static_cast<std::size_t>(nparms), 0.0));
    lambda = LAMBDAZERO;
}

int LMLSolver::iLMiter() {
    sosinit = myH->computeResiduals();
    if (sosinit == BIGVAL)  // failed ray?
        return BADITER;     // cannot proceed, request host OUTER LOOP exit.

    if (!myH->buildJacobian()) // ask host for new Jacobian.
        return BADITER;        // cannot proceed, request host OUTER LOOP exit.

    for (int k = 0; k < nparms; k++) // get downhill gradient beta
    {
        beta[k] = 0.0;
        for (int i = 0; i < npts; i++)
            beta[k] -= myH->getResidual(i) * myH->getJacobian(i, k);
    }
    for (int k = 0; k < nparms; k++) // get undamped curvature matrix alpha
        for (int j = 0; j < nparms; j++) {
            alpha[j][k] = 0.0;
            for (int i = 0; i < npts; i++)
                alpha[j][k] += myH->getJacobian(i, j) * myH->getJacobian(i, k);
        }

    double rise = 0;
    do /// LMinner damping loop searches for one downhill step
    {
        niter++;                         // local diagnostic only
        for (int k = 0; k < nparms; k++) // copy and damp it
            for (int j = 0; j < nparms; j++)
                amatrix[j][k] = alpha[j][k] + ((j == k) ? lambda : 0.0);
        gaussj(amatrix, nparms); // invert

        for (int k = 0; k < nparms; k++) // compute delta[]
        {
            delta[k] = 0.0;
            for (int j = 0; j < nparms; j++)
                delta[k] += amatrix[j][k] * beta[j];
        }

        sos = myH->nudge(delta); // try it out.
        rise = (sos - sosinit) / (1 + sosinit);

        //---four possibilities and three exits---------

        if (rise <= -lmtol)     // good downhill step!
        {
            lambda *= LMSHRINK; // shrink lambda
            return DOWNITER;    // request another OUTER LOOP iteration.
        }

        if (rise <= 0.0)        // good step but level; all done.
        {
            lambda *= LMSHRINK; // no need to shrink lambda?
            return LEVELITER;   // return to host: OUTER LOOP exit.
        }

        for (int k = 0; k < nparms; k++) // reverse course!
            delta[k] *= -1.0;
        myH->nudge(delta); // sosprev is still OK

        if (rise < lmtol)     // finished but keep prev parms
        {
            return LEVELITER; // return to host: OUTER LOOP exit.
        }
        // Diagnostic print carried over verbatim from the Java.
        std::printf("niter = %d\n", niter);
        if (niter >= lmiter)
            return MAXITER;

        lambda *= LMBOOST;         // UPITER:  apply more damping.
    } while (lambda < LAMBDAMAX);  // and stay in this INNER LOOP.

    return BADITER; // exceeded LAMBDAMAX, so request host OUTER LOOP exit.
                    // usual cause is ray damage during minimization.
}

double LMLSolver::gaussj(std::vector<std::vector<double>> &a, int N) {
    double det = 1.0, big, save;
    int i, j, k, L;
    // Fixed 100-element pivot arrays, as in the Java. Overflows silently for
    // N > 100; no caller comes close (nparms is the number of optimizer
    // variables).
    int ik[100];
    int jk[100];
    for (k = 0; k < N; k++) {
        big = 0.0;
        for (i = k; i < N; i++)
            for (j = k; j < N; j++) // find biggest element
                if (std::abs(big) <= std::abs(a[i][j])) {
                    big = a[i][j];
                    ik[k] = i;
                    jk[k] = j;
                }
        if (big == 0.0)
            return 0.0;
        i = ik[k];
        if (i > k)
            for (j = 0; j < N; j++) // exchange rows
            {
                save = a[k][j];
                a[k][j] = a[i][j];
                a[i][j] = -save;
            }
        j = jk[k];
        if (j > k)
            for (i = 0; i < N; i++) {
                save = a[i][k];
                a[i][k] = a[i][j];
                a[i][j] = -save;
            }
        for (i = 0; i < N; i++) // build the inverse
            if (i != k)
                a[i][k] = -a[i][k] / big;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                if ((i != k) && (j != k))
                    a[i][j] += a[i][k] * a[k][j];
        for (j = 0; j < N; j++)
            if (j != k)
                a[k][j] /= big;
        a[k][k] = 1.0 / big;
        det *= big; // bomb point
    }               // end k loop
    for (L = 0; L < N; L++) {
        k = N - L - 1;
        j = ik[k];
        if (j > k)
            for (i = 0; i < N; i++) {
                save = a[i][k];
                a[i][k] = -a[i][j];
                a[i][j] = save;
            }
        i = jk[k];
        if (i > k)
            for (j = 0; j < N; j++) {
                save = a[k][j];
                a[k][j] = -a[i][j];
                a[i][j] = save;
            }
    }
    return det;
}

} // namespace redukti::mathlib
