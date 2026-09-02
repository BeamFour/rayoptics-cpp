// C++ port of org.redukti.mathlib.Vector2
#ifndef REDUKTI_MATHLIB_VECTOR2_H
#define REDUKTI_MATHLIB_VECTOR2_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"

#include <array>
#include <cmath>
#include <string>

namespace redukti::mathlib {

class Vector3;

/**
 * Java's fields are `final`; here they are plain members and every method is
 * const and returns a new instance, so the type is immutable by convention.
 * They are not declared const because Java *variables* are rebindable and the
 * ported code relies on that (`v = v.plus(w)`).
 *
 * Java's `x()` accessor is dropped -- the public field `x` is identical and a
 * C++ member cannot share its name with a method. Java's `x(double)` "wither"
 * becomes `withX(double)` for the same reason.
 */
class Vector2 {
public:
    double x;
    double y;

    Vector2(double x_, double y_) : x(x_), y(y_) {
        if (std::isnan(x_) || std::isnan(y_)) {
            throw IllegalArgumentException("NaN");
        }
    }

    explicit Vector2(double v) : Vector2(v, v) {}

    Vector2 plus(const Vector2 &v) const { return Vector2(x + v.x, y + v.y); }

    Vector2 minus(const Vector2 &v) const { return Vector2(x - v.x, y - v.y); }

    Vector2 divide(double scalar) const { return Vector2(x / scalar, y / scalar); }

    Vector2 times(double scalar) const { return Vector2(x * scalar, y * scalar); }

    std::array<double, 2> as_array() const { return {x, y}; }

    Vector2 withX(double value) const { return Vector2(value, y); }
    Vector2 withY(double value) const { return Vector2(x, value); }

    double v(int i) const {
        switch (i) {
        case 0:
            return x;
        case 1:
            return y;
        default:
            throw IllegalArgumentException("Invalid offset " + std::to_string(i));
        }
    }

    Vector2 set(int i, double value) const {
        double x1 = this->x;
        double y1 = this->y;
        switch (i) {
        case 0:
            x1 = value;
            break;
        case 1:
            y1 = value;
            break;
        default:
            throw IllegalArgumentException("Invalid offset " + std::to_string(i));
        }
        return Vector2(x1, y1);
    }

    /** element by element divide */
    Vector2 ebeDivide(const Vector2 &v) const { return Vector2(x / v.x, y / v.y); }

    /** element by element multiply */
    Vector2 ebeTimes(const Vector2 &v) const { return Vector2(x * v.x, y * v.y); }

    Vector2 negate() const { return Vector2(-x, -y); }

    double len() const {
        double r = x * x + y * y;
        return std::sqrt(r);
    }

    static Vector2 from(const Vector3 &v3, int a, int b);

    Vector2 normalize() const {
        double lengthsq = x * x + y * y;
        if (M::isZero(lengthsq)) {
            return *this;
        } else {
            double denom = std::sqrt(lengthsq);
            return Vector2(x / denom, y / denom);
        }
    }

    std::string toString() const;

    bool isEqual(const Vector2 &other, double tolerance) const {
        return std::abs(this->x - other.x) < tolerance &&
               std::abs(this->y - other.y) < tolerance;
    }

    static const Vector2 vector2_0;
    static const Vector2 vector2_1;
    static const Vector2 vector2_10;
    static const Vector2 vector2_01;
};

inline const Vector2 Vector2::vector2_0{0.0, 0.0};
inline const Vector2 Vector2::vector2_1{1.0, 1.0};
inline const Vector2 Vector2::vector2_10{1.0, 0.0};
inline const Vector2 Vector2::vector2_01{0.0, 1.0};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_VECTOR2_H
