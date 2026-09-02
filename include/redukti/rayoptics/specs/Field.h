// C++ port of org.redukti.rayoptics.specs.{Field,ReadOnlyField}
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SPECS_FIELD_H
#define REDUKTI_RAYOPTICS_SPECS_FIELD_H

#include "redukti/mathlib/Vector2.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::raytr {
class ChiefRayPkg;
class ReferenceSphere;
class RayPkg;
} // namespace redukti::rayoptics::raytr

namespace redukti::rayoptics::specs {

class FieldSpec;

/**
 * One field point.
 *
 * Field sits at the centre of the raytr/specs dependency cycle: it refers to
 * ChiefRayPkg and ReferenceSphere in raytr, and RayPkg refers back to a
 * ReadOnlyField snapshot of it. The cyclic members are shared_ptr to types
 * that are only forward-declared here, so the destructor and update() are
 * defined out of line in Field.cpp where those types are complete.
 */
class Field {
public:
    double x = 0.0;   // x field component
    double y = 0.0;   // y field component
    double vux = 0.0; // +x vignetting factor
    double vuy = 0.0; // +y vignetting factor
    double vlx = 0.0; // -x vignetting factor
    double vly = 0.0; // -y vignetting factor
    double wt = 0.0;  // field weight

    /** Nullable in the Java; two elements when present. */
    std::optional<std::vector<double>> aim_info;
    std::optional<double> z_enp;
    std::shared_ptr<raytr::ChiefRayPkg> chief_ray;
    std::shared_ptr<raytr::ReferenceSphere> ref_sphere;
    std::map<std::string, std::shared_ptr<const raytr::RayPkg>> pupil_rays;
    /** Borrowed: the FieldSpec that owns this Field's list. Nullable. */
    FieldSpec *fov = nullptr;

    explicit Field(FieldSpec *fov_);
    ~Field();

    void update();

    void apply_scale_factor(double scale_factor) {
        x *= scale_factor;
        y *= scale_factor;
    }

    /** Applies the vignetting factors to a pupil coordinate pair. */
    std::vector<double> apply_vignetting(const std::vector<double> &pupil) const {
        std::vector<double> vig_pupil = pupil;
        vig_pupil[0] *= vignetting_scale_x(pupil[0]);
        vig_pupil[1] *= vignetting_scale_y(pupil[1]);
        return vig_pupil;
    }

    double vignetting_scale_x(double x_) const { return vignetting_scale(x_, vlx, vux); }
    double vignetting_scale_y(double y_) const { return vignetting_scale(y_, vly, vuy); }

    void clear_vignetting() { vux = vuy = vlx = vly = 0.; }

    std::string toString() const;

    void list_str(std::string &sb, const std::string &fmtstr) const;

    bool is_relative() const;
    double max_field() const;

    // NOTE: the x and y accessors are not symmetric in the Java. xv() returns
    // the field value and xf() the fraction, but yf() mirrors xv() and yv()
    // mirrors xf() -- so for y the two are swapped relative to x. Carried over
    // verbatim. The only caller outside list_str is a diagnostic message in
    // Wideangle, so the asymmetry shows up in report text and nowhere else.
    double xv() const {
        if (is_relative())
            return _get_x_by_fref();
        return x;
    }
    double xf() const {
        if (is_relative())
            return x;
        return _get_x_by_vref();
    }
    double yf() const {
        if (is_relative())
            return _get_y_by_fref();
        return y;
    }
    double yv() const {
        if (is_relative())
            return y;
        return _get_y_by_vref();
    }

    double _get_x_by_fref() const { return x * max_field(); }
    double _get_y_by_fref() const { return y * max_field(); }

private:
    double _get_x_by_vref() const { return x / max_field(); }
    double _get_y_by_vref() const { return y / max_field(); }

    static double vignetting_scale(double coordinate, double lower, double upper) {
        double factor = coordinate < 0.0 ? lower : upper;
        return factor == 0.0 ? 1.0 : 1.0 - factor;
    }
};

/** An immutable snapshot of a Field, taken when a RayPkg records its field. */
class ReadOnlyField {
public:
    double x;
    double y;
    double vux;
    double vuy;
    double vlx;
    double vly;
    double wt;
    std::optional<mathlib::Vector2> aim_info;
    std::optional<double> z_enp;
    std::shared_ptr<raytr::ChiefRayPkg> chief_ray;
    std::shared_ptr<raytr::ReferenceSphere> ref_sphere;
    FieldSpec *fov;

    explicit ReadOnlyField(const Field &fld);
    ~ReadOnlyField();
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_FIELD_H
