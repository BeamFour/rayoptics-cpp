// C++ port of org.redukti.rayoptics.math.Tfm3d
#ifndef REDUKTI_RAYOPTICS_MATH_TFM3D_H
#define REDUKTI_RAYOPTICS_MATH_TFM3D_H

#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/Vector3.h"

#include <optional>

namespace redukti::rayoptics::math {

class Tfm3d {
public:
    /**
     * Nullable in the Java: DecenterData::tform_before_surf and
     * tform_after_surf both construct a Tfm3d with a null rotation matrix to
     * mean "no rotation", and Transform tests for it.
     */
    std::optional<mathlib::Matrix3> rt;
    mathlib::Vector3 t;

    Tfm3d(const std::optional<mathlib::Matrix3> &rt_, const mathlib::Vector3 &t_)
        : rt(rt_), t(t_) {}
};

} // namespace redukti::rayoptics::math

#endif // REDUKTI_RAYOPTICS_MATH_TFM3D_H
