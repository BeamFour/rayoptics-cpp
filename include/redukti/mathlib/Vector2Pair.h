// C++ port of org.redukti.mathlib.Vector2Pair
//
// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
#ifndef REDUKTI_MATHLIB_VECTOR2PAIR_H
#define REDUKTI_MATHLIB_VECTOR2PAIR_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/mathlib/Vector3Pair.h"

#include <cmath>
#include <string>

namespace redukti::mathlib {

class Vector2Pair {
public:
    Vector2 v0;
    Vector2 v1;

    Vector2Pair(const Vector2 &v0_, const Vector2 &b) : v0(v0_), v1(b) {}

    bool isEquals(const Vector2Pair &other, double tolerance) const {
        return v0.isEqual(other.v0, tolerance) && v1.isEqual(other.v1, tolerance);
    }

    double ln_intersect_ln_scale(const Vector2Pair &line) const {
        // based on
        // http://geometryalgorithms.com/Archive/algorithm_0104/algorithm_0104B.htm

        Vector2 w = v0.minus(line.v0);

        double d = v1.x * line.v1.y - v1.y * line.v1.x;

        if (std::abs(d) < 1e-10)
            throw IllegalArgumentException("ln_intersect_ln_scale: lines are parallel");

        double s = (line.v1.x * w.y - line.v1.y * w.x) / d;

        return s;
    }

    Vector2 ln_intersect_ln(const Vector2Pair &line) const {
        return v0.plus(v1.times(ln_intersect_ln_scale(line)));
    }

    /**
     * Create a 2d vector pair and initialize vectors from
     * specified components of vectors from an other pair.
     */
    static Vector2Pair from(const Vector3Pair &v, int c0, int c1) {
        return Vector2Pair(Vector2::from(v.v0, c0, c1), Vector2::from(v.v1, c0, c1));
    }

    std::string toString() const { return "[" + v0.toString() + "," + v1.toString() + "]"; }

    static const Vector2Pair vector2_pair_00;
};

inline const Vector2Pair Vector2Pair::vector2_pair_00{Vector2::vector2_0,
                                                      Vector2::vector2_0};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_VECTOR2PAIR_H
