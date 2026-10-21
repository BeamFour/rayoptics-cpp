/*
 * Gaussian quadrature utilities.
 *
 * C++ port of org.redukti.mathlib.GaussQ.
 *
 * The direct Gauss-Legendre implementation follows W. H. Press et al.,
 * Numerical Recipes: The Art of Scientific Computing, 3rd ed., section 4.6
 * (Cambridge University Press, 2007).
 *
 * The alternative implementation is a port of the Legendre, no-fixed-
 * endpoint path through gaussq.f from Netlib's "go" (Golden Oldies)
 * collection: https://netlib.org/go/gaussq.f
 * The Jacobi-matrix construction and the calculation of weights from the first
 * eigenvector components are the Golub-Welsch algorithm in reference 1. The
 * Netlib routine replaces that paper's shifted QR procedure with GAUSQ2, a
 * modified EISPACK IMTQL2 implicit QL eigensolver.
 *
 * Original gaussq.f notices:
 *
 *     original version 20 jan 1975 from stanford
 *     modified 21 dec 1983 by eric grosse
 *       imtql2 => gausq2
 *       hex constant => d1mach (from core library)
 *       compute pi using datan
 *       removed accuracy claims, description of method
 *       added single precision version
 *
 * References retained from gaussq.f:
 *
 * 1. G. H. Golub and J. H. Welsch, "Calculation of Gaussian Quadrature
 *    Rules," Mathematics of Computation 23 (1969), 221-230.
 * 2. G. H. Golub, "Some Modified Matrix Eigenvalue Problems," SIAM Review
 *    15 (1973), 318-334, section 7.
 * 3. A. H. Stroud and D. Secrest, Gaussian Quadrature Formulas,
 *    Prentice-Hall, Englewood Cliffs, N.J., 1966.
 */
#ifndef REDUKTI_MATHLIB_GAUSSQ_H
#define REDUKTI_MATHLIB_GAUSSQ_H

#include <vector>

namespace redukti::mathlib {

/** Generates Gauss-Legendre nodes and weights on `[-1, 1]`. */
class GaussQ {
public:
    GaussQ() = delete;

    /**
     * Computes an `order`-point Gauss-Legendre rule by refining the roots of the
     * Legendre polynomial with Newton's method.
     *
     * @return {nodes, weights}, with nodes in ascending order
     */
    static std::vector<std::vector<double>> gauss_legendre_nodes_and_weights(int order);

    /**
     * Computes the same rule using the Golub-Welsch Jacobi-matrix method. This is the
     * `kind=1`, `kpts=0` path through Netlib `gaussq.f`.
     *
     * @return {nodes, weights}, with nodes in ascending order
     */
    static std::vector<std::vector<double>>
    gauss_legendre_nodes_and_weights_golub_welsch(int order);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_GAUSSQ_H
