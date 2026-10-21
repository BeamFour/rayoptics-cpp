// C++ port of org.redukti.mathlib.GaussQTest.
//
// The two Gauss-Legendre implementations are checked against each other, against
// published order-16 data, and against the polynomial moments a rule of each order
// must integrate exactly.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/GaussQ.h"

#include <cmath>
#include <string>
#include <vector>

namespace {

using redukti::mathlib::GaussQ;

/*
 * Order-16 reference rule from John Burkardt's QUADRATURE_RULES_LEGENDRE
 * dataset, distributed under the GNU LGPL:
 * https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_legendre/
 */
const double ORDER_16_NODES[] = {
    -0.9894009349916499,  -0.9445750230732326,  -0.8656312023878318,
    -0.7554044083550030,  -0.6178762444026438,  -0.4580167776572274,
    -0.2816035507792589,  -0.09501250983763744, 0.09501250983763744,
    0.2816035507792589,   0.4580167776572274,   0.6178762444026438,
    0.7554044083550030,   0.8656312023878318,   0.9445750230732326,
    0.9894009349916499};

const double ORDER_16_WEIGHTS[] = {
    0.02715245941175407, 0.06225352393864787, 0.09515851168249277, 0.1246289712555339,
    0.1495959888165767,  0.1691565193950026,  0.1826034150449236,  0.1894506104550686,
    0.1894506104550686,  0.1826034150449236,  0.1691565193950026,  0.1495959888165767,
    0.1246289712555339,  0.09515851168249277, 0.06225352393864787, 0.02715245941175407};

double sum(const std::vector<double> &values) {
    double result = 0.0;
    for (double value : values)
        result += value;
    return result;
}

void assertMoments(const std::vector<std::vector<double>> &rule, int order) {
    CHECK_CLOSE(sum(rule[1]), 2.0, 2.0e-14);
    for (int degree = 0; degree < 2 * order; degree++) {
        double actual = 0.0;
        for (int i = 0; i < order; i++) {
            auto k = static_cast<std::size_t>(i);
            actual += rule[1][k] * std::pow(rule[0][k], degree);
        }
        double expected = degree % 2 == 0 ? 2.0 / (degree + 1.0) : 0.0;
        CHECK_CLOSE(actual, expected, 2.0e-14);
    }
    for (int i = 1; i < order; i++)
        CHECK(rule[0][static_cast<std::size_t>(i - 1)] <
              rule[0][static_cast<std::size_t>(i)]);
}

void assertReferenceRule(const std::vector<std::vector<double>> &rule) {
    for (std::size_t i = 0; i < 16; i++) {
        CHECK_CLOSE(rule[0][i], ORDER_16_NODES[i], 5.0e-15);
        CHECK_CLOSE(rule[1][i], ORDER_16_WEIGHTS[i], 5.0e-15);
    }
}

} // namespace

TEST(gaussq_direct_rule_has_known_three_point_values) {
    auto rule = GaussQ::gauss_legendre_nodes_and_weights(3);

    CHECK_CLOSE(rule[0][0], -std::sqrt(3.0 / 5.0), 1.0e-15);
    CHECK_CLOSE(rule[0][1], 0.0, 1.0e-15);
    CHECK_CLOSE(rule[0][2], std::sqrt(3.0 / 5.0), 1.0e-15);
    CHECK_CLOSE(rule[1][0], 5.0 / 9.0, 5.0e-15);
    CHECK_CLOSE(rule[1][1], 8.0 / 9.0, 5.0e-15);
    CHECK_CLOSE(rule[1][2], 5.0 / 9.0, 5.0e-15);
}

TEST(gaussq_golub_welsch_agrees_with_direct_root_refinement) {
    for (int order = 1; order <= 64; order++) {
        auto direct = GaussQ::gauss_legendre_nodes_and_weights(order);
        auto golubWelsch = GaussQ::gauss_legendre_nodes_and_weights_golub_welsch(order);

        for (int i = 0; i < order; i++) {
            auto k = static_cast<std::size_t>(i);
            CHECK_CLOSE(golubWelsch[0][k], direct[0][k], 2.0e-14);
            CHECK_CLOSE(golubWelsch[1][k], direct[1][k], 2.0e-14);
        }
    }
}

TEST(gaussq_both_rules_agree_with_published_order_sixteen_data) {
    assertReferenceRule(GaussQ::gauss_legendre_nodes_and_weights(16));
    assertReferenceRule(GaussQ::gauss_legendre_nodes_and_weights_golub_welsch(16));
}

TEST(gaussq_both_rules_integrate_polynomial_moments) {
    for (int order = 1; order <= 8; order++) {
        assertMoments(GaussQ::gauss_legendre_nodes_and_weights(order), order);
        assertMoments(GaussQ::gauss_legendre_nodes_and_weights_golub_welsch(order), order);
    }
}

TEST(gaussq_rules_reject_non_positive_order) {
    CHECK_THROWS(GaussQ::gauss_legendre_nodes_and_weights(0),
                 redukti::IllegalArgumentException);
    CHECK_THROWS(GaussQ::gauss_legendre_nodes_and_weights_golub_welsch(-1),
                 redukti::IllegalArgumentException);
}
