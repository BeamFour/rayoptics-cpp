// C++ port of org.redukti.mathlib.GaussQ. The attributions and references are in the
// header.
#include "redukti/mathlib/GaussQ.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"

#include <cmath>

namespace redukti::mathlib {

namespace {

void requirePositiveOrder(int order) {
    if (order < 1)
        throw IllegalArgumentException("Quadrature order must be at least 1");
}

/**
 * This routine is a translation of an ALGOL procedure, Numerische Mathematik 12,
 * 377-383 (1968), by Martin and Wilkinson, as modified in Numerische Mathematik 15,
 * 450 (1970), by Dubrulle; Handbook for Automatic Computation, vol. II, Linear
 * Algebra, 241-248 (1971). This is a modified version of the EISPACK routine IMTQL2.
 *
 * It finds the eigenvalues and first components of the eigenvectors of a symmetric
 * tridiagonal matrix by the implicit QL method.
 *
 * @param diagonal on input, the diagonal elements of the symmetric tridiagonal
 *                 matrix; on output, its eigenvalues in ascending order
 * @param offDiagonal on input, the subdiagonal elements in positions 0 through
 *                    n - 2; position n - 1 is unused. The array is destroyed during
 *                    the calculation.
 * @param firstComponents on input, the first row of the identity matrix,
 *                        [1, 0, ..., 0]; on output, the first component of each
 *                        normalized eigenvector, in the same order as the returned
 *                        eigenvalues
 * @return zero for normal completion, or the one-based index of an eigenvalue not
 *         determined after 30 iterations
 */
int gausq2(std::vector<double> &diagonal, std::vector<double> &offDiagonal,
           std::vector<double> &firstComponents) {
    int order = static_cast<int>(diagonal.size());
    if (order == 1)
        return 0;

    offDiagonal[static_cast<std::size_t>(order - 1)] = 0.0;
    // Java's Math.ulp(1.0).
    double machinePrecision = std::nextafter(1.0, 2.0) - 1.0;

    for (int l = 0; l < order; l++) {
        int iterations = 0;
        while (true) {
            /* Look for a small sub-diagonal element. */
            int m;
            for (m = l; m < order; m++) {
                auto k = static_cast<std::size_t>(m);
                if (m == order - 1 ||
                    std::abs(offDiagonal[k]) <=
                        machinePrecision * (std::abs(diagonal[k]) + std::abs(diagonal[k + 1])))
                    break;
            }

            double p = diagonal[static_cast<std::size_t>(l)];
            if (m == l)
                break;
            if (iterations == 30)
                return l + 1;
            iterations++;

            /* Form shift. */
            auto li = static_cast<std::size_t>(l);
            double g = (diagonal[li + 1] - p) / (2.0 * offDiagonal[li]);
            double r = std::hypot(g, 1.0);
            g = diagonal[static_cast<std::size_t>(m)] - p +
                offDiagonal[li] / (g + std::copysign(r, g));
            double sine = 1.0;
            double cosine = 1.0;
            p = 0.0;

            for (int i = m - 1; i >= l; i--) {
                auto k = static_cast<std::size_t>(i);
                double f = sine * offDiagonal[k];
                double b = cosine * offDiagonal[k];
                if (std::abs(f) >= std::abs(g)) {
                    cosine = g / f;
                    r = std::hypot(cosine, 1.0);
                    offDiagonal[k + 1] = f * r;
                    sine = 1.0 / r;
                    cosine *= sine;
                } else {
                    sine = f / g;
                    r = std::hypot(sine, 1.0);
                    offDiagonal[k + 1] = g * r;
                    cosine = 1.0 / r;
                    sine *= cosine;
                }
                g = diagonal[k + 1] - p;
                r = (diagonal[k] - g) * sine + 2.0 * cosine * b;
                p = sine * r;
                diagonal[k + 1] = g + p;
                g = cosine * r - b;

                /* Form first component of vector. */
                f = firstComponents[k + 1];
                firstComponents[k + 1] = sine * firstComponents[k] + cosine * f;
                firstComponents[k] = cosine * firstComponents[k] - sine * f;
            }

            diagonal[li] -= p;
            offDiagonal[li] = g;
            offDiagonal[static_cast<std::size_t>(m)] = 0.0;
        }
    }

    /* Order eigenvalues and eigenvectors. */
    for (int i = 0; i < order - 1; i++) {
        int smallest = i;
        double value = diagonal[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < order; j++) {
            if (diagonal[static_cast<std::size_t>(j)] < value) {
                smallest = j;
                value = diagonal[static_cast<std::size_t>(j)];
            }
        }
        if (smallest != i) {
            auto a = static_cast<std::size_t>(i);
            auto b = static_cast<std::size_t>(smallest);
            diagonal[b] = diagonal[a];
            diagonal[a] = value;
            double component = firstComponents[a];
            firstComponents[a] = firstComponents[b];
            firstComponents[b] = component;
        }
    }
    return 0;
}

} // namespace

std::vector<std::vector<double>> GaussQ::gauss_legendre_nodes_and_weights(int order) {
    requirePositiveOrder(order);

    auto n = static_cast<std::size_t>(order);
    std::vector<double> nodes(n, 0.0);
    std::vector<double> weights(n, 0.0);
    int rootsToFind = (order + 1) / 2;

    for (int i = 0; i < rootsToFind; i++) {
        double x = std::cos(M::PI * (i + 0.75) / (order + 0.5));
        double derivative;
        double delta;
        do {
            double p0 = 1.0;
            double p1 = x;
            for (int degree = 2; degree <= order; degree++) {
                double p2 = ((2.0 * degree - 1.0) * x * p1 - (degree - 1.0) * p0) / degree;
                p0 = p1;
                p1 = p2;
            }
            double polynomial = p1;
            double previousPolynomial = p0;
            derivative = order * (x * polynomial - previousPolynomial) / (x * x - 1.0);
            delta = polynomial / derivative;
            x -= delta;
        } while (std::abs(delta) > 1.0e-15);

        double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
        auto k = static_cast<std::size_t>(i);
        nodes[k] = -x;
        nodes[n - 1 - k] = x;
        weights[k] = weight;
        weights[n - 1 - k] = weight;
    }
    return {nodes, weights};
}

std::vector<std::vector<double>>
GaussQ::gauss_legendre_nodes_and_weights_golub_welsch(int order) {
    requirePositiveOrder(order);

    auto n = static_cast<std::size_t>(order);
    std::vector<double> nodes(n, 0.0);
    std::vector<double> offDiagonal(n, 0.0);
    std::vector<double> firstComponents(n, 0.0);

    /*
     * kind = 1: Legendre polynomials p(x) on (-1, +1), w(x) = 1.
     * The matrix is symmetric.  nodes contains its diagonal elements and
     * offDiagonal its off-diagonal elements.
     */
    for (int i = 1; i < order; i++) {
        double index = i;
        offDiagonal[static_cast<std::size_t>(i - 1)] =
            index / std::sqrt(4.0 * index * index - 1.0);
    }

    /* The first row of the identity matrix. */
    firstComponents[0] = 1.0;

    int errorIndex = gausq2(nodes, offDiagonal, firstComponents);
    if (errorIndex != 0)
        throw ArithmeticException("Gauss-Legendre eigensolver did not converge at index " +
                                  intToString(errorIndex));

    /* muzero = integral from -1 to 1 of w(x) dx = 2. */
    std::vector<double> weights(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
        weights[i] = 2.0 * firstComponents[i] * firstComponents[i];
    return {nodes, weights};
}

} // namespace redukti::mathlib
