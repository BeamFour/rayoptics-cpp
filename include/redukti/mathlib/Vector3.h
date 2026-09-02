// C++ port of org.redukti.mathlib.Vector3
#ifndef REDUKTI_MATHLIB_VECTOR3_H
#define REDUKTI_MATHLIB_VECTOR3_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2.h"

#include <cmath>
#include <string>

namespace redukti::mathlib {

/** See the note in Vector2.h on final fields and the dropped `x()` accessors. */
class Vector3 {
public:
    double x;
    double y;
    double z;

    Vector3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {
        if (std::isnan(x_) || std::isnan(y_) || std::isnan(z_)) {
            throw IllegalArgumentException("NaN");
        }
    }

    explicit Vector3(double v) : Vector3(v, v, v) {}

    Vector3 deg2rad() const {
        return Vector3(M::toRadians(x), M::toRadians(y), M::toRadians(z));
    }

    Vector3 tan() const { return Vector3(std::tan(x), std::tan(y), std::tan(z)); }

    Vector3 times(double scale) const { return Vector3(x * scale, y * scale, z * scale); }

    Vector3 divide(double scale) const { return Vector3(x / scale, y / scale, z / scale); }

    Vector3 negate() const { return Vector3(-x, -y, -z); }

    double length() const { return std::sqrt(x * x + y * y + z * z); }

    Vector3 normalize() const;

    double dot(const Vector3 &vector) const {
        return x * vector.x + y * vector.y + z * vector.z;
    }

    /**
     * The cross product a x b is defined as a vector c that is perpendicular
     * (orthogonal) to both a and b, with a direction given by the right-hand
     * rule and a magnitude equal to the area of the parallelogram that the
     * vectors span.
     *
     * https://en.wikipedia.org/wiki/Cross_product
     */
    Vector3 cross(const Vector3 &b) const {
        return Vector3(y * b.z - z * b.y,
                       z * b.x - x * b.z,
                       x * b.y - y * b.x);
    }

    Vector3 add(const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }

    bool isZero() const { return M::isZero(x * x + y * y + z * z); }

    bool any() const { return !isZero(); }

    Vector3 plus(const Vector3 &v) const { return Vector3(x + v.x, y + v.y, z + v.z); }

    Vector3 minus(const Vector3 &v) const { return Vector3(x - v.x, y - v.y, z - v.z); }

    std::string toString() const;

    Vector3 withX(double v) const { return Vector3(v, y, z); }
    Vector3 withY(double v) const { return Vector3(x, v, z); }
    Vector3 withZ(double v) const { return Vector3(x, y, v); }

    Vector3 sin() const { return Vector3(std::sin(x), std::sin(y), std::sin(z)); }

    Vector2 project_xy() const { return Vector2(x, y); }

    Vector2 project_zy() const { return Vector2(z, y); }

    double v(int i) const {
        switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        case 2:
            return z;
        default:
            throw IllegalArgumentException("Invalid offset " + std::to_string(i));
        }
    }

    Vector3 v(int i, double value) const {
        double x1 = this->x;
        double y1 = this->y;
        double z1 = this->z;
        switch (i) {
        case 0:
            x1 = value;
            break;
        case 1:
            y1 = value;
            break;
        case 2:
            z1 = value;
            break;
        default:
            throw IllegalArgumentException("Invalid offset " + std::to_string(i));
        }
        return Vector3(x1, y1, z1);
    }

    bool isEqual(const Vector3 &other, double tolerance) const {
        return std::abs(this->x - other.x) < tolerance &&
               std::abs(this->y - other.y) < tolerance &&
               std::abs(this->z - other.z) < tolerance;
    }

    /** Mirrors Java's equals(), which compares with Double.compare(). */
    bool equals(const Vector3 &other) const;

    bool effectivelyEqual(const Vector3 &o) const { return isEqual(o, 1e-13); }

    static const Vector3 ZERO;
    static const Vector3 vector3_0;
    static const Vector3 vector3_1;
    static const Vector3 vector3_001;
    static const Vector3 vector3_010;
    static const Vector3 vector3_100;
};

inline const Vector3 Vector3::ZERO{0, 0, 0};
inline const Vector3 Vector3::vector3_0{0.0, 0.0, 0.0};
inline const Vector3 Vector3::vector3_1{1.0, 1.0, 1.0};
inline const Vector3 Vector3::vector3_001{0.0, 0.0, 1.0};
inline const Vector3 Vector3::vector3_010{0.0, 1.0, 0.0};
inline const Vector3 Vector3::vector3_100{1.0, 0.0, 0.0};

inline Vector3 Vector3::normalize() const {
    double lengthsq = x * x + y * y + z * z;
    if (M::isZero(lengthsq)) {
        return ZERO;
    } else {
        double denom = std::sqrt(lengthsq);
        return Vector3(x / denom, y / denom, z / denom);
    }
}

// Deferred from Vector2.h, which only had a forward declaration of Vector3.
inline Vector2 Vector2::from(const Vector3 &v3, int a, int b) {
    double x = v3.v(a);
    double y = v3.v(b);
    return Vector2(x, y);
}

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_VECTOR3_H
