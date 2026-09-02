// C++ port of org.redukti.mathlib.Quaternion
//
// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
#ifndef REDUKTI_MATHLIB_QUATERNION_H
#define REDUKTI_MATHLIB_QUATERNION_H

#include "redukti/mathlib/Vector3.h"

#include <cmath>
#include <string>

namespace redukti::mathlib {

class Quaternion {
public:
    double x, y, z, w;

    Quaternion(double x_, double y_, double z_, double w_) : x(x_), y(y_), z(z_), w(w_) {}

    /**
     * Get shortest arc rotation between for a to be in the same direction as b.
     * Both vectors must be unit vectors.
     */
    static Quaternion get_rotation_between(const Vector3 &a, const Vector3 &b) {
        // Do not know the source of following equation
        // Believe it generates a Quaternion representing the rotation
        // of vector a to vector b
        // Closest match of the algo:
        // https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
        // FIXME It seems this implementation is not safe
        // See QuaternionBase<Derived>::setFromTwoVectors in eigen library
        // Also stackoverflow discussion

        Vector3 cp = a.cross(b);
        double _x = cp.x;
        double _y = cp.y;
        double _z = cp.z;
        double _w = a.dot(b) + 1.0;
        double n = norm(_x, _y, _z, _w);
        _x = _x / n;
        _y = _y / n;
        _z = _z / n;
        _w = _w / n;
        return Quaternion(_x, _y, _z, _w);
    }

    static double norm(double x, double y, double z, double w) {
        return std::sqrt(x * x + y * y + z * z + w * w);
    }

    std::string toString() const;

    /** Mirrors Java's equals(), which compares with Double.compare(). */
    bool equals(const Quaternion &other) const;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_QUATERNION_H
