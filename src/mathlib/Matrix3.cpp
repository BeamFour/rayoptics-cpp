// C++ port of org.redukti.mathlib.Matrix3
#include "redukti/mathlib/Matrix3.h"

#include "redukti/Text.h"

namespace redukti::mathlib {

Matrix3 Matrix3::get_rotation_matrix(int axis, double angleInRadians) {
    assert(axis < 3 && axis >= 0);

    /*
     * Note on convention used below.
     *
     * See https://mathworld.wolfram.com/RotationMatrix.html
     * coordinate system rotations of the x-, y-, and z-axes in a
     * counterclockwise direction when looking towards the origin give the
     * matrices.
     *
     * This appears to correspond to xyz convention described in appendix A,
     * Classical Mechanics, Goldstein, 3rd Ed.
     */
    double n00, n01, n02, n10, n11, n12, n20, n21, n22;
    switch (axis) {
    case 0:
        // rotation counter clockwise around the X axis
        n00 = 1;
        n01 = 0;
        n02 = 0;
        n10 = 0;
        n11 = std::cos(angleInRadians);
        n12 = std::sin(angleInRadians);
        n20 = 0;
        n21 = -std::sin(angleInRadians);
        n22 = std::cos(angleInRadians);
        break;

    case 1:
        // rotation counter clockwise around the Y axis
        n00 = std::cos(angleInRadians);
        n01 = 0;
        n02 = -std::sin(angleInRadians);
        n10 = 0;
        n11 = 1;
        n12 = 0;
        n20 = std::sin(angleInRadians);
        n21 = 0;
        n22 = std::cos(angleInRadians);
        break;

    case 2:
        // rotation counter clockwise around the Z axis
        n00 = std::cos(angleInRadians);
        n01 = std::sin(angleInRadians);
        n02 = 0;
        n10 = -std::sin(angleInRadians);
        n11 = std::cos(angleInRadians);
        n12 = 0;
        n20 = 0;
        n21 = 0;
        n22 = 1;
        break;
    default:
        throw IllegalArgumentException("Invalid rotation axis, must be 0=x, 1=y or 2=z");
    }
    return Matrix3(n00, n01, n02,
                   n10, n11, n12,
                   n20, n21, n22);
}

std::string Matrix3::toString() const {
    return "[[" + doubleToString(m00) + "," + doubleToString(m01) + "," +
           doubleToString(m02) + "],\n [" + doubleToString(m10) + "," +
           doubleToString(m11) + "," + doubleToString(m12) + "],\n [" +
           doubleToString(m20) + "," + doubleToString(m21) + "," +
           doubleToString(m22) + "]]";
}

} // namespace redukti::mathlib
