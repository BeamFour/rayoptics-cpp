// C++ port of org.redukti.rayoptics.elem.surface.DecenterData
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_SURFACE_DECENTERDATA_H
#define REDUKTI_RAYOPTICS_ELEM_SURFACE_DECENTERDATA_H

#include "redukti/mathlib/Matrix3.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/math/Tfm3d.h"

#include <optional>
#include <string>

namespace redukti::rayoptics::elem::surface {

/** Fields are package-private in the Java; public here. */
class DecenterData {
public:
    std::string dtype;
    mathlib::Vector3 dec;
    mathlib::Vector3 euler;
    mathlib::Vector3 rot_pt;
    /** Null in the Java when there is no rotation. */
    std::optional<mathlib::Matrix3> rot_mat;

    DecenterData(std::string dtype_, double x, double y, double alpha, double beta,
                 double gamma)
        : dtype(std::move(dtype_)), dec(x, y, 0.0), euler(alpha, beta, gamma),
          rot_pt(mathlib::Vector3::ZERO), rot_mat(std::nullopt) {}

    explicit DecenterData(std::string dtype_)
        : DecenterData(std::move(dtype_), 0.0, 0.0, 0.0, 0.0, 0.0) {}

    void update() {
        if (euler.any())
            rot_mat = mathlib::Matrix3::euler2mat_rxyz(convertl2r().deg2rad());
        else
            rot_mat = std::nullopt;
    }

    void apply_scale_factor(double scale_factor) {
        dec = dec.times(scale_factor);
        rot_pt = rot_pt.times(scale_factor);
    }

    math::Tfm3d tform_before_surf() const {
        if (dtype != "reverse")
            return math::Tfm3d(rot_mat, dec);
        else
            return math::Tfm3d(std::nullopt, mathlib::Vector3::ZERO);
    }

    math::Tfm3d tform_after_surf() const {
        if (dtype == "reverse" || dtype == "dec and return") {
            std::optional<mathlib::Matrix3> rt = rot_mat;
            if (rot_mat.has_value())
                rt = rot_mat->transpose();
            return math::Tfm3d(rt, dec.negate());
        } else if (dtype == "bend")
            return math::Tfm3d(rot_mat, mathlib::Vector3::ZERO);
        else
            return math::Tfm3d(std::nullopt, mathlib::Vector3::ZERO);
    }

private:
    mathlib::Vector3 convertl2r() const {
        return mathlib::Vector3(-euler.x, -euler.y, euler.z);
    }
};

} // namespace redukti::rayoptics::elem::surface

#endif // REDUKTI_RAYOPTICS_ELEM_SURFACE_DECENTERDATA_H
