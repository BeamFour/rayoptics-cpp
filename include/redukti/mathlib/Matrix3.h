// C++ port of org.redukti.mathlib.Matrix3
// Code derived from https://github.com/jvanverth/essentialmath
// Portions Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_MATHLIB_MATRIX3_H
#define REDUKTI_MATHLIB_MATRIX3_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Quaternion.h"
#include "redukti/mathlib/Vector3.h"

#include <cassert>
#include <cmath>
#include <string>

namespace redukti::mathlib {

/**
 * Column major 3d matrix where
 *
 *    0=m00 3=m01 6=m02
 *    1=m10 4=m11 7=m12
 *    2=m20 5=m21 8=m22
 */
class Matrix3 {
public:
    // Package-private and final in Java; public here because the tests read them.
    double m00, m01, m02;
    double m10, m11, m12;
    double m20, m21, m22;

    Matrix3(double m00_, double m01_, double m02_, double m10_, double m11_, double m12_,
            double m20_, double m21_, double m22_)
        : m00(m00_), m01(m01_), m02(m02_), m10(m10_), m11(m11_), m12(m12_), m20(m20_),
          m21(m21_), m22(m22_) {}

    /* Create a diagonal matrix with given values */
    static Matrix3 diag(double x, double y, double z) {
        return Matrix3(x, 0.0, 0.0,
                       0.0, y, 0.0,
                       0.0, 0.0, z);
    }

    static Matrix3 identity() { return diag(1., 1., 1.); }

    Matrix3 inverse() const {
        // compute determinant
        double cofactor0 = m11 * m22 - m21 * m12;
        double cofactor3 = m20 * m12 - m10 * m22;
        double cofactor6 = m10 * m21 - m20 * m11;
        double det = m00 * cofactor0 + m01 * cofactor3 + m02 * cofactor6;

        if (M::isZero(det)) {
            throw RuntimeException("Determinant is 0; singular matrix");
        }

        // create adjoint matrix and multiply by 1/det to get inverse
        double invDet = 1.0f / det;
        double n00 = invDet * cofactor0;
        double n10 = invDet * cofactor3;
        double n20 = invDet * cofactor6;

        double n01 = invDet * (m21 * m02 - m01 * m22);
        double n11 = invDet * (m00 * m22 - m20 * m02);
        double n21 = invDet * (m20 * m01 - m00 * m21);

        double n02 = invDet * (m01 * m12 - m11 * m02);
        double n12 = invDet * (m10 * m02 - m00 * m12);
        double n22 = invDet * (m00 * m11 - m10 * m01);

        return Matrix3(n00, n01, n02,
                       n10, n11, n12,
                       n20, n21, n22);
    }

    Matrix3 transpose() const {
        double n10 = m01;
        double n01 = m10;
        double n20 = m02;
        double n02 = m20;
        double n21 = m12;
        double n12 = m21;
        return Matrix3(m00, n01, n02,
                       n10, m11, n12,
                       n20, n21, m22);
    }

    /*
    // Quaternion to Rotation Matrix
    // Q = (x, y, z, w)
    // 1-2y^2-2z^2      2xy-2wz         2xz+2wy
    // 2xy+2zw          1-2x^2-2z^2     2yz-2xw
    // 2xz-2yw          2yz+2xw         1-2x^2-2y^2
    //
    // see https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
    */
    static Matrix3 to_rotation_matrix(const Quaternion &q) {
        double n00 = 1.0 - 2.0 * (q.y * q.y + q.z * q.z);
        double n10 = 2.0 * (q.x * q.y + q.z * q.w);
        double n20 = 2.0 * (q.x * q.z - q.y * q.w);

        double n01 = 2.0 * (q.x * q.y - q.z * q.w);
        double n11 = 1.0 - 2.0 * (q.x * q.x + q.z * q.z);
        double n21 = 2.0 * (q.z * q.y + q.x * q.w);

        double n02 = 2.0 * (q.x * q.z + q.y * q.w);
        double n12 = 2.0 * (q.y * q.z - q.x * q.w);
        double n22 = 1.0 - 2.0 * (q.x * q.x + q.y * q.y);
        return Matrix3(n00, n01, n02,
                       n10, n11, n12,
                       n20, n21, n22);
    }

    /**
     * Get a rotation matrix that describes the rotation of vector a
     * to obtain vector 3.
     * @param from Unit vector
     * @param to Unit vector
     */
    static Matrix3 get_rotation_between(const Vector3 &from, const Vector3 &to) {
        // Do not know the source of following equation
        // Believe it generates a Quaternion representing the rotation
        // of vector a to vector b
        Quaternion q = Quaternion::get_rotation_between(from, to);
        return to_rotation_matrix(q);
    }

    static Matrix3 euler2mat(double roll_angle, double pitch_angle, double yaw_angle) {
        double si = std::sin(roll_angle), sj = std::sin(pitch_angle),
               sk = std::sin(yaw_angle);
        double ci = std::cos(roll_angle), cj = std::cos(pitch_angle),
               ck = std::cos(yaw_angle);
        double cc = ci * ck, cs = ci * sk;
        double sc = si * ck, ss = si * sk;

        // gamma (roll), beta (pitch), alpha (yaw)
        // see https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        // More formally, it is an intrinsic rotation whose Tait-Bryan angles are
        // alpha, beta, gamma, about axes z, y, x, respectively.
        // The formula below corresponds to yaw.multiply(pitch.multiply(roll))
        // which means roll followed by pitch followed by yaw
        double n00 = cj * ck;      // Cos(y) * Cos(z)
        double n01 = sj * sc - cs; // Sin(y) * Sin(x) * Cos(z) - Cos(x) * Sin(z)
        double n02 = sj * cc + ss; // Sin(y) * Cos(x) * Cos(z) + Sin(x) * Sin(z)
        double n10 = cj * sk;      // Cos(y) * Sin(z)
        double n11 = sj * ss + cc; // Sin(y) * Sin(x) * Sin(z) + Cos(x) * Cos(z)
        double n12 = sj * cs - sc; // Sin(y) * Cos(x) * Sin(z) - Sin(x) * Cos(z)
        double n20 = -sj;          // -Sin(y)
        double n21 = cj * si;      // Cos(y) * Sin(x)
        double n22 = cj * ci;      // Cos(y) * Cos(x)

        return Matrix3(n00, n01, n02,
                       n10, n11, n12,
                       n20, n21, n22);
    }

    static Matrix3 euler2mat(const Vector3 &euler) {
        return euler2mat(euler.x, euler.y, euler.z);
    }

    /**
     * Rotating (intrinsic) frame x-y-z euler angles to a rotation matrix,
     * i.e. the transforms3d euler2mat(ai, aj, ak, axes=rxyz), which is
     * Rx(ai) * Ry(aj) * Rz(ak).
     *
     * euler2mat(Vector3) is the static (extrinsic) frame form, the transforms3d
     * default axes=sxyz = Rz(ak) * Ry(aj) * Rx(ai). The two agree when only one
     * angle is non-zero and differ for compound rotations. They are related by
     * rxyz(e) == sxyz(-e).transpose().
     */
    static Matrix3 euler2mat_rxyz(const Vector3 &euler) {
        return euler2mat(euler.negate()).transpose();
    }

    /** rotate around vert axis */
    static Matrix3 yaw(double angle) {
        double sine_theta = std::sin(angle), cos_theta = std::cos(angle);
        return Matrix3(cos_theta, -sine_theta, 0.0f,
                       sine_theta, cos_theta, 0.0f,
                       0.0f, 0.0f, 1.0f);
    }

    /** rotate around sideways axis - i.e. tilt */
    static Matrix3 pitch(double angle) {
        double sine_theta = std::sin(angle), cos_theta = std::cos(angle);
        return Matrix3(cos_theta, 0.0f, sine_theta,
                       0.0f, 1.0f, 0.0f,
                       -sine_theta, 0.0f, cos_theta);
    }

    /** Rotate around forward axis, i.e. turn */
    static Matrix3 roll(double angle) {
        double sine_theta = std::sin(angle), cos_theta = std::cos(angle);
        return Matrix3(1.0f, 0.0f, 0.0f,
                       0.0f, cos_theta, -sine_theta,
                       0.0f, sine_theta, cos_theta);
    }

    /**
     * Get rotation matrix for rotation about axis.
     * @param axis the axis of rotation, x=0, y=1, z=2
     * @param angleInRadians the angle to rotate in radians
     */
    static Matrix3 get_rotation_matrix(int axis, double angleInRadians);

    /**
     * rotate v1 into v2 using equivalent angle rotation.
     *
     * Compute a rotation matrix from v1 to v2. Take the cross product of the
     * input vectors to get the rotation axis. The equivalent angle rotation is
     * equation 2.80 from Introduction to Robotics, 2nd ed, by John J Craig.
     *
     * NOTE: n20 and n21 below repeat the expressions used for n02 and n12
     * rather than negating their s terms as the textbook equivalent-angle
     * matrix does. Carried over verbatim from the Java: Wideangle, FieldSpec
     * and OpticalSpecs are calibrated against this output, so changing it here
     * would silently move ray aiming results.
     */
    static Matrix3 rot_v1_into_v2(const Vector3 &v1, const Vector3 &v2) {
        auto rot_axis = v1.cross(v2).negate();
        auto s = rot_axis.length();
        auto cosine_ang = v1.dot(v2);
        auto c = cosine_ang;
        auto v = 1.0 - cosine_ang;
        auto ax = rot_axis.normalize();
        auto n00 = ax.x * ax.x * v + c;
        auto n01 = ax.x * ax.y * v - ax.z * s;
        auto n02 = ax.x * ax.z * v + ax.y * s;
        auto n10 = ax.x * ax.y * v + ax.z * s;
        auto n11 = ax.y * ax.y * v + c;
        auto n12 = ax.y * ax.z * v + ax.x * s;
        auto n20 = ax.x * ax.z * v + ax.y * s;
        auto n21 = ax.y * ax.z * v + ax.x * s;
        auto n22 = ax.z * ax.z * v + c;
        return Matrix3(n00, n01, n02,
                       n10, n11, n12,
                       n20, n21, n22);
    }

    Matrix3 multiply(const Matrix3 &other) const {
        double n00 = m00 * other.m00 + m01 * other.m10 + m02 * other.m20;
        double n10 = m10 * other.m00 + m11 * other.m10 + m12 * other.m20;
        double n20 = m20 * other.m00 + m21 * other.m10 + m22 * other.m20;
        double n01 = m00 * other.m01 + m01 * other.m11 + m02 * other.m21;
        double n11 = m10 * other.m01 + m11 * other.m11 + m12 * other.m21;
        double n21 = m20 * other.m01 + m21 * other.m11 + m22 * other.m21;
        double n02 = m00 * other.m02 + m01 * other.m12 + m02 * other.m22;
        double n12 = m10 * other.m02 + m11 * other.m12 + m12 * other.m22;
        double n22 = m20 * other.m02 + m21 * other.m12 + m22 * other.m22;
        return Matrix3(n00, n01, n02,
                       n10, n11, n12,
                       n20, n21, n22);
    }

    Vector3 multiply(const Vector3 &other) const {
        double x = m00 * other.x + m01 * other.y + m02 * other.z;
        double y = m10 * other.x + m11 * other.y + m12 * other.z;
        double z = m20 * other.x + m21 * other.y + m22 * other.z;
        return Vector3(x, y, z);
    }

    std::string toString() const;

    static const Matrix3 IDENTITY;
};

inline const Matrix3 Matrix3::IDENTITY = Matrix3::identity();

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_MATRIX3_H
