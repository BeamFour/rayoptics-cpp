// C++ port of org.redukti.mathlib.Vector3Pair
//
// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
#ifndef REDUKTI_MATHLIB_VECTOR3PAIR_H
#define REDUKTI_MATHLIB_VECTOR3PAIR_H

#include "redukti/mathlib/Vector3.h"

#include <string>

namespace redukti::mathlib {

/*
Notes:
Plane can be represented as two vectors: point and normal
Ray can be represented as two vectors: origin and direction
*/
class Vector3Pair {
public:
    Vector3 v0;
    Vector3 v1;

    /**
     * @param v0 First vector, origin / point
     * @param v1 Second vector, direction / normal
     */
    Vector3Pair(const Vector3 &v0_, const Vector3 &v1_) : v0(v0_), v1(v1_) {}

    const Vector3 &point() const { return v0; }
    const Vector3 &origin() const { return v0; }
    const Vector3 &direction() const { return v1; }
    const Vector3 &normal() const { return v1; }

    double z0() const { return v0.z; }
    double z1() const { return v1.z; }

    bool isEquals(const Vector3Pair &other, double tolerance) const {
        return v0.isEqual(other.v0, tolerance) && v1.isEqual(other.v1, tolerance);
    }

    double pl_ln_intersect_scale(const Vector3Pair &line) const {
        // See https://en.wikipedia.org/wiki/Line-plane_intersection
        return (origin().dot(normal()) - normal().dot(line.origin())) /
               (line.normal().dot(normal()));
    }

    Vector3 pl_ln_intersect(const Vector3Pair &line) const {
        return line.v0.plus(line.v1.times(pl_ln_intersect_scale(line)));
    }

    /**
     * Swap the given element between the member vectors and return a new pair.
     *
     * NOTE: the Java builds n0 and n1 but then constructs both vectors of the
     * returned pair from n0, so the swap is not actually observable. Kept
     * verbatim -- the method has no callers anywhere in the codebase.
     */
    static Vector3Pair swapElement(const Vector3Pair &p, int j) {
        double n0[3];
        double n1[3];

        for (int i = 0; i < 3; i++) {
            if (i == j) {
                // swap
                n0[i] = p.v1.v(i);
                n1[i] = p.v0.v(i);
            } else {
                // retain original
                n0[i] = p.v0.v(i);
                n1[i] = p.v1.v(i);
            }
        }
        (void)n1;
        return Vector3Pair(Vector3(n0[0], n0[1], n0[2]), Vector3(n0[0], n0[1], n0[2]));
    }

    std::string toString() const { return "[" + v0.toString() + "," + v1.toString() + "]"; }

    double x1() const { return v1.x; }

    static const Vector3Pair position_000_001;
};

inline const Vector3Pair Vector3Pair::position_000_001{Vector3::vector3_0,
                                                       Vector3::vector3_001};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_VECTOR3PAIR_H
