// C++ port of org.redukti.mathlib.Transform3
//
// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
#ifndef REDUKTI_MATHLIB_TRANSFORM3_H
#define REDUKTI_MATHLIB_TRANSFORM3_H

#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/mathlib/Vector3Pair.h"

#include <string>

namespace redukti::mathlib {

class Transform3 {
public:
    Vector3 translation;
    /* Rotation matrix to rotate a unit vector toward z to the required direction */
    Matrix3 rotation_matrix;
    /** Whether to use rotation matrix */
    bool use_rotation_matrix;

    Transform3()
        : translation(Vector3::vector3_0), rotation_matrix(Matrix3::diag(1.0, 1.0, 1.0)),
          use_rotation_matrix(false) {}

    explicit Transform3(const Vector3Pair &position)
        : translation(position.point()), rotation_matrix(Matrix3::diag(1.0, 1.0, 1.0)),
          use_rotation_matrix(false) {
        if (position.direction().x == 0 && position.direction().y == 0) {
            if (position.direction().z < 0.0) {
                this->rotation_matrix = Matrix3::diag(1.0, 1.0, -1.0);
                this->use_rotation_matrix = true;
            } else {
                this->rotation_matrix = Matrix3::diag(1.0, 1.0, 1.0);
                this->use_rotation_matrix = false;
            }
        } else {
            // Get a rotation matrix representing the rotation of unit vector in z
            // to the direction vector.
            this->rotation_matrix =
                Matrix3::get_rotation_between(Vector3::vector3_001, position.direction());
            this->use_rotation_matrix = true;
        }
    }

    Transform3(const Vector3 &translation_, const Matrix3 &rotation_matrix_,
               bool use_rotation_matrix_)
        : translation(translation_), rotation_matrix(rotation_matrix_),
          use_rotation_matrix(use_rotation_matrix_) {}

    /** Apply this transforms rotation to given vector */
    Vector3 apply_rotation(const Vector3 &v) const {
        if (use_rotation_matrix)
            return this->rotation_matrix.multiply(v);
        else
            return v;
    }

    /**
     * Composition. New translation is set to: apply parent's linear
     * transformation on child translation and add parent translation. New
     * linear matrix is the product of the parent and child matrices.
     *
     * @param p Parent component
     * @param c Child component
     */
    static Transform3 compose(const Transform3 &p, const Transform3 &c) {
        Vector3 translation = p.apply_rotation(c.translation).plus(p.translation);
        bool use_rotation_matrix = p.use_rotation_matrix || c.use_rotation_matrix;
        Matrix3 rotation_matrix = p.rotation_matrix.multiply(c.rotation_matrix);
        return Transform3(translation, rotation_matrix, use_rotation_matrix);
    }

    /** Transform given vector - vector is rotated and then translated. */
    Vector3 transform(const Vector3 &v) const {
        return apply_rotation(v).plus(translation);
    }

    /** Create an inverse of the transform */
    Transform3 inverse() const {
        Matrix3 rotation_matrix = this->rotation_matrix.inverse();
        Vector3 translation = rotation_matrix.multiply(this->translation.negate());
        return Transform3(translation, rotation_matrix, true);
    }

    /**
     * Rotate by x, y, and z axis.
     * @param v Vector with angles per axis
     */
    Transform3 rotate_axis_by_angles(const Vector3 &v) const {
        Transform3 t = *this;
        for (int i = 0; i < 3; i++) { // i stands for x,y,z axis
            if (v.v(i) != 0.0) {
                t = t.rotate_axis_by_angle(i, v.v(i));
            }
        }
        return t;
    }

    /**
     * Rotate around specified axis
     * @param axis   0=x, 1=y, 2=z
     * @param dangle Angle of rotation in degrees
     */
    Transform3 rotate_axis_by_angle(int axis, double dangle) const {
        return rotate_axis_by_radian(axis, M::toRadians(dangle));
    }

    /**
     * Rotate around specified axis
     * @param axis   0=x, 1=y, 2=z
     * @param rangle Angle of rotation in radians
     */
    Transform3 rotate_axis_by_radian(int axis, double rangle) const {
        Matrix3 r = Matrix3::get_rotation_matrix(axis, rangle);
        r = r.multiply(this->rotation_matrix);
        return Transform3(this->translation, r, true);
    }

    Vector3Pair transform_pair(const Vector3Pair &p) const {
        return Vector3Pair(transform(p.v0), transform(p.v1));
    }

    Transform3 set_translation(const Vector3 &translation_) const {
        return Transform3(translation_, rotation_matrix, use_rotation_matrix);
    }

    Transform3 set_direction(const Vector3 &direction) const {
        if (direction.x == 0.0 && direction.y == 0.0) {
            if (direction.z < 0.0) {
                return Transform3(translation, Matrix3::diag(1.0, 1.0, -1.0), true);
            } else {
                return Transform3(translation, Matrix3::diag(1.0, 1.0, 1.0), false);
            }
        } else {
            return Transform3(
                translation, Matrix3::get_rotation_between(Vector3::vector3_001, direction),
                true);
        }
    }

    Vector3Pair transform_line(const Vector3Pair &v) const {
        return Vector3Pair(transform(v.origin()), apply_rotation(v.direction()));
    }

    std::string toString() const;
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_TRANSFORM3_H
