// C++ port of org.redukti.mathlib.Sphere3
#ifndef REDUKTI_MATHLIB_SPHERE3_H
#define REDUKTI_MATHLIB_SPHERE3_H

#include "redukti/mathlib/Line3.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector3.h"

#include <array>
#include <cmath>
#include <optional>

namespace redukti::mathlib {

class Sphere3 {
public:
    Vector3 center;
    double radius;

    Sphere3(const Vector3 &center_, double radius_) : center(center_), radius(radius_) {}

    /**
     * Based on Geometric Tools for Computer Graphics.
     * Returns up to two solutions; Java returns a Double[] whose entries may be
     * null, which maps to std::optional here.
     */
    std::array<std::optional<double>, 2> intersect(const Line3 &line) const {
        auto diff = line.origin.minus(center);
        auto a0 = diff.dot(diff) - radius * radius;
        auto a1 = line.direction.dot(diff);
        auto discr = a1 * a1 - a0;
        if (M::isZero(discr)) {
            return {-a1, std::nullopt};
        } else if (discr > 0.0) {
            auto root = std::sqrt(discr);
            return {-a1 - root, -a1 + root};
        }
        return {std::nullopt, std::nullopt};
    }
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_SPHERE3_H
