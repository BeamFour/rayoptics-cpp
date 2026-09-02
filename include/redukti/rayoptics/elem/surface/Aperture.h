// C++ port of org.redukti.rayoptics.elem.surface.{Aperture,Circular,InteractionMode,IntersectionResult}
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ELEM_SURFACE_APERTURE_H
#define REDUKTI_RAYOPTICS_ELEM_SURFACE_APERTURE_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <cmath>
#include <optional>

namespace redukti::rayoptics::elem::surface {

enum class InteractionMode {
    Transmit,
    Reflect,
};

class IntersectionResult {
public:
    /** distance to intersection point */
    double distance;
    /** intersection point p */
    mathlib::Vector3 intersection_point;

    IntersectionResult(double x, const mathlib::Vector3 &v)
        : distance(x), intersection_point(v) {}
};

class Aperture {
public:
    double x_offset = 0.0;
    double y_offset = 0.0;
    double rotation = 0.0;
    bool is_obscuration = false;

    Aperture(double x_offset_, double y_offset_, double rotation_, bool is_obscuration_)
        : x_offset(x_offset_), y_offset(y_offset_), rotation(rotation_),
          is_obscuration(is_obscuration_) {}

    virtual ~Aperture() = default;

    virtual mathlib::Vector2 dimension() const = 0;
    virtual void set_dimension(double x, double y) = 0;

    virtual double max_dimension() const {
        mathlib::Vector2 d = dimension();
        return std::sqrt(d.x * d.x + d.y * d.y);
    }

    /** `fuzz` is a nullable Double in the Java; nullopt selects the default. */
    virtual bool point_inside(double x, double y, std::optional<double> fuzz) const = 0;

    virtual mathlib::Vector2 edge_pt_target(const mathlib::Vector2 &rel_dir) const = 0;

    util::Pair<mathlib::Vector2, mathlib::Vector2> bounding_box() const {
        mathlib::Vector2 center(x_offset, y_offset);
        mathlib::Vector2 extent = dimension();
        return util::Pair<mathlib::Vector2, mathlib::Vector2>(center.minus(extent),
                                                              center.plus(extent));
    }

    virtual void apply_scale_factor(double scale_factor) {
        x_offset *= scale_factor;
        y_offset *= scale_factor;
    }

    mathlib::Vector2 tform(double x, double y) const {
        x -= x_offset;
        y -= y_offset;
        return mathlib::Vector2(x, y);
    }
};

class Circular : public Aperture {
public:
    double radius = 1.0;

    Circular(double x_offset_, double y_offset_, double rotation_, double radius_,
             bool is_obscuration_)
        : Aperture(x_offset_, y_offset_, rotation_, is_obscuration_), radius(radius_) {}

    mathlib::Vector2 dimension() const override { return mathlib::Vector2(radius, radius); }

    void set_dimension(double x, double y) override {
        (void)y;
        radius = x;
    }

    double max_dimension() const override { return radius; }

    bool point_inside(double x, double y, std::optional<double> fuzz) const override {
        double f = fuzz.has_value() ? *fuzz : 1e-5;
        mathlib::Vector2 v = tform(x, y);
        auto ans = std::sqrt(v.x * v.x + v.y * v.y) <= radius + f;
        return !is_obscuration ? ans : !ans;
    }

    mathlib::Vector2 edge_pt_target(const mathlib::Vector2 &rel_dir) const override {
        return rel_dir.normalize().times(radius);
    }

    void apply_scale_factor(double scale_factor) override {
        Aperture::apply_scale_factor(scale_factor);
        radius *= scale_factor;
    }
};

} // namespace redukti::rayoptics::elem::surface

#endif // REDUKTI_RAYOPTICS_ELEM_SURFACE_APERTURE_H
