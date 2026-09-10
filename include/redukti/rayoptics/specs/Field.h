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
 * a single field point, chief ray pkg and pupil limits
 *
 *     The Field class manages several types of data:
 *
 *     - the field coordinates, unscaled and fractional
 *     - aim info for tracing through the stop surface
 *     - the vignetting factors for the pupil definition
 *     - pkgs for the chief ray and reference sphere
 *
 *     The Field can have a reference to a fov/FieldSpec (recommended!) which is used to support the fractional and value interfaces simultaneously. If
 *     no fov is given, a max_field may be specified, with the default being unit field size.
 *
 *     Attributes:
 *         vux: +x vignetting factor
 *         vuy: +y vignetting factor
 *         vlx: -x vignetting factor
 *         vly: -y vignetting factor
 *         wt: field weight
 *         aim_info: x, y chief ray coords on the paraxial entrance pupil plane,
 *                   or z_enp for wide angle fovs
 *         chief_ray: ray package for the ray from the field point throught the
 *                    center of the aperture stop, traced in the central
 *                    wavelength
 *         ref_sphere: a tuple containing (image_pt, ref_dir, ref_sphere_radius)
 *         fov: :class:`~.FieldSpec` to be used as reference or None
 */
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

    /**
     * Populated for non wide-angle system
     * x, y chief ray coords on the paraxial entrance pupil plane
     * When this is populated z_enp should be null.
     */
    /** Nullable in the Java; two elements when present. */
    std::optional<std::vector<double>> aim_info;
    /**
     * The z center of the real pupil for `fld`, wrt 1st ifc
     * Populated for wide-angle system
     * When this is populated aim_info should be null;
     */
    std::optional<double> z_enp;
    /**
     * ray package for the ray from the field point through the
     * center of the aperture stop, traced in the central wavelength
     */
    std::shared_ptr<raytr::ChiefRayPkg> chief_ray;
    /**
     * a tuple containing (image_pt, ref_dir, ref_sphere_radius)
     */
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

    /**
     * Scale relative pupil coordinates by this field's vignetting factors,
     * returning a new array. The argument is not modified.
     *
     * This differs from upstream, deliberately. Upstream writes
     * vig_pupil = pupil[:], which for a numpy array is a view rather
     * than a copy, so scaling vig_pupil also scales the caller's array
     * in place. Callers that read their pupil array back after tracing see the
     * vignetted value there: analyses.trace_ray_fan records the pupil
     * after calling trace_safe and so reports vignetted fan
     * coordinates, where this implementation reports the nominal ones.
     *
     * The rays traced are the same either way - only what a caller observes in
     * its own array differs - but it is visible when comparing fan data against
     * upstream, where a first ray at nominal -1.0 shows up there as
     * -1 * (1 - vlx). See tools/src/main/python/README.md.
     *
     * @param pupil relative pupil coordinates, unmodified by this call
     * @return a new array with the vignetting factors applied
     */
    /** Applies the vignetting factors to a pupil coordinate pair. */
    std::vector<double> apply_vignetting(const std::vector<double> &pupil) const {
        std::vector<double> vig_pupil = pupil;
        vig_pupil[0] *= vignetting_scale_x(pupil[0]);
        vig_pupil[1] *= vignetting_scale_y(pupil[1]);
        return vig_pupil;
    }

    /**
     * Factor by which apply_vignetting scales an x pupil coordinate.
     * The upper and lower factors differ, so the scale depends on the sign of
     * the coordinate and the map has a kink at the axis.
     */
    double vignetting_scale_x(double x_) const { return vignetting_scale(x_, vlx, vux); }
    /** Factor by which apply_vignetting scales a y pupil coordinate. */
    double vignetting_scale_y(double y_) const { return vignetting_scale(y_, vly, vuy); }

    /**
     * Resets vignetting values to 0.
     */
    void clear_vignetting() { vux = vuy = vlx = vly = 0.; }

    std::string toString() const;

    void list_str(std::string &sb, const std::string &fmtstr) const;

    bool is_relative() const;
    /**
     * the maximum field value used for the fractional field calculation.
     */
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

/** Analysis metadata captured at result creation, with no model dependencies. */
/** Analysis metadata with no model dependencies; results own it as const. */
class FieldSnapshot {
public:
    double x, y, vux, vuy, vlx, vly, wt;
    explicit FieldSnapshot(const Field &field)
        : x(field.x), y(field.y), vux(field.vux), vuy(field.vuy),
          vlx(field.vlx), vly(field.vly), wt(field.wt), label(field.toString()) {}
    std::string toString() const { return label; }
private:
    std::string label;
};

/**
 * A readonly snapshot of a Field
 */
/** Snapshot taken when a RayPkg records its field. */
class ReadOnlyField {
public:
    double x;  // x field component
    double y;  // y field component
    double vux;  // +x vignetting factor
    double vuy;  // +y vignetting factor
    double vlx;  // -x vignetting factor
    double vly;  // -y vignetting factor
    double wt;  //  field weight
    std::optional<mathlib::Vector2> aim_info;
    /**
     * Populated for wide-angle system
     */
    std::optional<double> z_enp;
    std::shared_ptr<raytr::ChiefRayPkg> chief_ray;
    /**
     * a tuple containing (image_pt, ref_dir, ref_sphere_radius)
     */
    std::shared_ptr<raytr::ReferenceSphere> ref_sphere;
    FieldSpec *fov;

    explicit ReadOnlyField(const Field &fld);
    ~ReadOnlyField();
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_FIELD_H
