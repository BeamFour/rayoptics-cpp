// C++ port of the raytr data types:
//   RaySeg, RayPkg, RayData, RayDataWithZ_Enp, ReferenceSphere,
//   ChiefRayExitPupilSegment, ChiefRayPkg, RefSphereCR, RayResult and its
//   variants, AimInfo, VigResult, GridItem, the trace definitions and enums.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_RAYTYPES_H
#define REDUKTI_RAYOPTICS_RAYTR_RAYTYPES_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/seq/Interface.h"
#include "redukti/rayoptics/specs/Field.h"

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::raytr {

enum class RayFanType {
    TransverseRayFan,
    OpticalPathDifference,
};

enum class PupilType {
    REL_PUPIL, // relative pupil coordinates
    AIM_PT,    // aim point on pupil plane
    AIM_DIR,   // aim direction in object space
};

/** One segment of a traced ray. */
class RaySeg {
public:
    mathlib::Vector3 p;
    mathlib::Vector3 d;
    double dst;
    mathlib::Vector3 nrml;

    // TODO phase

    RaySeg(const mathlib::Vector3 &p_, const mathlib::Vector3 &d_, double dst_,
           const mathlib::Vector3 &nrml_)
        : p(p_), d(d_), dst(dst_), nrml(nrml_) {}

    RaySeg(const RaySeg &other, double dst_delta)
        : RaySeg(other.p, other.d, other.dst + dst_delta, other.nrml) {}

    std::string toString() const;
};

/**
 * A traced ray and its optical path length.
 *
 * Immutable once constructed, and shared widely -- it is stored in RayResult,
 * in TraceException's payload, and in Field's pupil_rays -- so it is passed as
 * shared_ptr<const RayPkg> throughout, matching Java's reference semantics
 * without copying the segment list.
 */
class RayPkg {
public:
    std::vector<RaySeg> ray;
    double op_delta;
    double wvl;
    /** Null unless the ray was traced for a specific field. */
    std::shared_ptr<const specs::ReadOnlyField> fld;
    std::optional<mathlib::Vector2> input_pupil;
    std::optional<mathlib::Vector2> vig_pupil;

    RayPkg(std::vector<RaySeg> ray_, double op_delta_, double wvl_)
        : ray(std::move(ray_)), op_delta(op_delta_), wvl(wvl_) {}

    RayPkg(std::vector<RaySeg> ray_, double op_delta_, double wvl_,
           const specs::Field *fld_, std::optional<mathlib::Vector2> input_pupil_,
           std::optional<mathlib::Vector2> vig_pupil_);

    /** Java's with(fld, input_pupil, vig_pupil). */
    std::shared_ptr<const RayPkg> with(const specs::Field *fld_,
                                       std::optional<mathlib::Vector2> input_pupil_,
                                       std::optional<mathlib::Vector2> vig_pupil_) const;

    std::string toString() const;
};

/** A point and a direction cosine. */
class RayData {
public:
    mathlib::Vector3 pt;
    mathlib::Vector3 dir;

    RayData(const mathlib::Vector3 &pt_, const mathlib::Vector3 &dir_)
        : pt(pt_), dir(dir_) {}
};

class RayDataWithZ_Enp {
public:
    RayData ray_data;
    double z_enp;

    RayDataWithZ_Enp(const RayData &ray_data_, double z_enp_)
        : ray_data(ray_data_), z_enp(z_enp_) {}
};

class ReferenceSphere {
public:
    mathlib::Vector3 image_pt;
    mathlib::Vector3 ref_dir;
    double ref_sphere_radius;
    math::Tfm3d lcl_tfrm_last;

    ReferenceSphere(const mathlib::Vector3 &image_pt_, const mathlib::Vector3 &ref_dir_,
                    double ref_sphere_radius_, const math::Tfm3d &lcl_tfrm_last_)
        : image_pt(image_pt_), ref_dir(ref_dir_),
          ref_sphere_radius(ref_sphere_radius_), lcl_tfrm_last(lcl_tfrm_last_) {}
};

class ChiefRayExitPupilSegment {
public:
    mathlib::Vector3 exp_pt;
    mathlib::Vector3 exp_dir;
    double exp_dst;
    /** Shared with the SequentialModel's interface list. */
    std::shared_ptr<seq::Interface> ifc;
    mathlib::Vector3 b4_pt;
    mathlib::Vector3 b4_dir;

    ChiefRayExitPupilSegment(const mathlib::Vector3 &exp_pt_,
                             const mathlib::Vector3 &exp_dir_, double exp_dst_,
                             std::shared_ptr<seq::Interface> ifc_,
                             const mathlib::Vector3 &b4_pt_,
                             const mathlib::Vector3 &b4_dir_)
        : exp_pt(exp_pt_), exp_dir(exp_dir_), exp_dst(exp_dst_), ifc(std::move(ifc_)),
          b4_pt(b4_pt_), b4_dir(b4_dir_) {}
};

class ChiefRayPkg {
public:
    std::shared_ptr<const RayPkg> chief_ray;
    std::shared_ptr<const ChiefRayExitPupilSegment> cr_exp_seg;

    ChiefRayPkg(std::shared_ptr<const RayPkg> chief_ray_,
                std::shared_ptr<const ChiefRayExitPupilSegment> cr_exp_seg_)
        : chief_ray(std::move(chief_ray_)), cr_exp_seg(std::move(cr_exp_seg_)) {}
};

class RefSphereCR {
public:
    std::shared_ptr<const ReferenceSphere> ref_sphere;
    std::shared_ptr<const ChiefRayPkg> chief_ray_pkg;

    RefSphereCR(std::shared_ptr<const ReferenceSphere> ref_sphere_,
                std::shared_ptr<const ChiefRayPkg> chief_ray_pkg_)
        : ref_sphere(std::move(ref_sphere_)), chief_ray_pkg(std::move(chief_ray_pkg_)) {}
};

/**
 * A traced ray or the reason it failed.
 *
 * This is the shape the Java already uses to carry trace failures as data
 * rather than as an in-flight exception: callers catch a TraceException, pull
 * the partial ray out of it, and store both here. Exactly one of pkg and err
 * is typically set.
 */
class RayResult {
public:
    std::shared_ptr<const RayPkg> pkg;
    std::shared_ptr<exceptions::TraceException> err;

    RayResult() = default;
    RayResult(std::shared_ptr<const RayPkg> pkg_,
              std::shared_ptr<exceptions::TraceException> err_)
        : pkg(std::move(pkg_)), err(std::move(err_)) {}
};

class RayResultWithStartCoord {
public:
    std::optional<std::vector<double>> start_coords;
    RayResult rr;
};

class RayResultWithStopCoord {
public:
    mathlib::Vector3 stop_coord;
    RayResult rr;
    int stop_idx;

    RayResultWithStopCoord(const mathlib::Vector3 &stop_coord_, RayResult rr_,
                           int stop_idx_)
        : stop_coord(stop_coord_), rr(std::move(rr_)), stop_idx(stop_idx_) {}
};

class RayResultWithZEnp {
public:
    std::optional<double> z_enp;
    RayResult rr;

    RayResultWithZEnp(std::optional<double> z_enp_, RayResult rr_)
        : z_enp(z_enp_), rr(std::move(rr_)) {}
};

class AimInfo {
public:
    /** aim_pt is used for paraxial aiming */
    std::vector<double> aim_pt;
    /** the actual entrance pupil distance wrt the 1st ifc for a field */
    std::optional<double> z_enp;

    AimInfo(std::vector<double> aim_pt_, std::optional<double> z_enp_)
        : aim_pt(std::move(aim_pt_)), z_enp(z_enp_) {}
};

class VigResult {
public:
    double vig;
    std::optional<int> clip_indx;
    std::shared_ptr<const RayPkg> ray_pkg;

    VigResult(double vig_, std::optional<int> clip_indx_,
              std::shared_ptr<const RayPkg> ray_pkg_)
        : vig(vig_), clip_indx(clip_indx_), ray_pkg(std::move(ray_pkg_)) {}
};

class GridItem {
public:
    mathlib::Vector2 pupil;
    std::optional<double> result;
    std::shared_ptr<const RayPkg> ray_pkg;
    double weight;
    bool valid;

    GridItem(const mathlib::Vector2 &pupil_, std::shared_ptr<const RayPkg> ray_pkg_)
        : GridItem(pupil_, std::move(ray_pkg_), std::nullopt, 1.0, true) {}

    GridItem(const mathlib::Vector2 &pupil_, std::shared_ptr<const RayPkg> ray_pkg_,
             double result_)
        : GridItem(pupil_, std::move(ray_pkg_), result_, 1.0, true) {}

    static GridItem failed(const mathlib::Vector2 &pupil_) {
        return GridItem(pupil_, nullptr, std::nullopt, 1.0, false);
    }

    GridItem withWeight(double weight_) const {
        return GridItem(pupil, ray_pkg, result, weight_, valid);
    }

    std::string toString() const { return pupil.toString(); }

private:
    GridItem(const mathlib::Vector2 &pupil_, std::shared_ptr<const RayPkg> ray_pkg_,
             std::optional<double> result_, double weight_, bool valid_)
        : pupil(pupil_), result(result_), ray_pkg(std::move(ray_pkg_)), weight(weight_),
          valid(valid_) {}
};

class TraceFanDef {
public:
    mathlib::Vector2 start;
    mathlib::Vector2 stop;
    int num_rays;

    TraceFanDef(const mathlib::Vector2 &fan_start, const mathlib::Vector2 &fan_stop,
                int num_rays_)
        : start(fan_start), stop(fan_stop), num_rays(num_rays_) {}
};

class TraceGridDef {
public:
    mathlib::Vector2 grid_start;
    mathlib::Vector2 grid_stop;
    int num_rays;

    TraceGridDef(const mathlib::Vector2 &grid_start_, const mathlib::Vector2 &grid_stop_,
                 int num_rays_)
        : grid_start(grid_start_), grid_stop(grid_stop_), num_rays(num_rays_) {}
};

class TraceRingsDef {
public:
    double cx = 0;
    double cy = 0;
    int num_rings = 21;
    /** Normalized inner radius; zero selects a filled circular pupil. */
    double min_radius = 0.0;
    double max_radius = 1.0;
    int num_points_in_ring_one = 6;
    bool hexapolar = true;
};

class TraceGridByWvl {
public:
    double wvl;
    std::vector<GridItem> grid;

    TraceGridByWvl(double wvl_, std::vector<GridItem> grid_)
        : wvl(wvl_), grid(std::move(grid_)) {}
};

/** Java's `record ContrastRayTriplet(...)`. */
class ContrastRayTriplet {
public:
    mathlib::Vector2 pupil;
    std::shared_ptr<const RayPkg> reference;
    std::shared_ptr<const RayPkg> sagittal;
    std::shared_ptr<const RayPkg> tangential;
    std::shared_ptr<exceptions::TraceException> referenceError;
    std::shared_ptr<exceptions::TraceException> sagittalError;
    std::shared_ptr<exceptions::TraceException> tangentialError;
    double weight;

    ContrastRayTriplet(const mathlib::Vector2 &pupil_,
                       std::shared_ptr<const RayPkg> reference_,
                       std::shared_ptr<const RayPkg> sagittal_,
                       std::shared_ptr<const RayPkg> tangential_,
                       std::shared_ptr<exceptions::TraceException> referenceError_,
                       std::shared_ptr<exceptions::TraceException> sagittalError_,
                       std::shared_ptr<exceptions::TraceException> tangentialError_,
                       double weight_)
        : pupil(pupil_), reference(std::move(reference_)), sagittal(std::move(sagittal_)),
          tangential(std::move(tangential_)), referenceError(std::move(referenceError_)),
          sagittalError(std::move(sagittalError_)),
          tangentialError(std::move(tangentialError_)), weight(weight_) {}
};

/** Java's `record ContrastTraceByWvl<T>(double wavelength, List<T> samples)`. */
template <typename T> class ContrastTraceByWvl {
public:
    double wavelength;
    std::vector<T> samples;

    ContrastTraceByWvl(double wavelength_, std::vector<T> samples_)
        : wavelength(wavelength_), samples(std::move(samples_)) {}
};

/**
 * Java's `interface ImageFilter`.
 *
 * The return is nullable on purpose. Trace::trace_grid, trace_rings and
 * trace_gaussian_quadrature call the filter with a null `pkg` for a ray that
 * failed or fell outside the pupil, and the analysis callbacks answer null
 * there rather than fabricating a grid point. The empty optional is what makes
 * `append_if_none` mean anything, so callers must test it.
 */
class ImageFilter {
public:
    virtual ~ImageFilter() = default;
    virtual std::optional<GridItem> apply(const mathlib::Vector2 &pupil,
                                          const std::shared_ptr<const RayPkg> &pkg) = 0;
};

/**
 * Java's `interface TraceGridCallback`. Null (an empty optional) means the ray
 * did not produce a grid point -- see the note on ImageFilter.
 */
using TraceGridCallback = std::function<std::optional<GridItem>(
    const mathlib::Vector2 &p, int wi, const std::shared_ptr<const RayPkg> &ray_pkg,
    specs::Field &fld, double wvl, double foc)>;

/** Java's `interface RayFanCallback`. Null means "no value at this pupil point". */
using RayFanCallback = std::function<std::optional<double>(
    optical::OpticalModel *opt_model, const mathlib::Vector2 &p, int wi,
    const std::shared_ptr<const RayPkg> &ray_pkg, specs::Field &fld, double wvl,
    double foc)>;

/** Java's `interface ContrastTraceCallback<T>`. */
template <typename T>
using ContrastTraceCallback = std::function<T(
    const ContrastRayTriplet &rays, specs::Field &field, double wavelength, double focus)>;

/** Java's `class TraceFanPoints`. */
class TraceFanPoints {
public:
    double wvl;
    std::vector<double> fan_x;
    /** Nullable: a traced ray whose callback returned null contributes a null. */
    std::vector<std::optional<double>> fan_y;
    std::vector<GridItem> fan;

    TraceFanPoints(double wvl_, std::vector<double> fan_x_,
                   std::vector<std::optional<double>> fan_y_, std::vector<GridItem> fan_)
        : wvl(wvl_), fan_x(std::move(fan_x_)), fan_y(std::move(fan_y_)),
          fan(std::move(fan_)) {}
};

/** Java's `class TraceFanResult`. */
class TraceFanResult {
public:
    /** Left unset until setFanType; Java leaves the field null. */
    std::optional<RayFanType> type;
    std::shared_ptr<const specs::FieldSnapshot> fld;
    int fi;
    int xy;
    std::vector<TraceFanPoints> fans;
    double max_rho_val;
    double max_y_val;

    TraceFanResult(specs::Field *fld_, int fi_, int xy_, std::vector<TraceFanPoints> fans_,
                   double max_rho_val_, double max_y_val_)
        : fld(fld_ ? std::make_shared<const specs::FieldSnapshot>(*fld_) : nullptr), fi(fi_), xy(xy_), fans(std::move(fans_)),
          max_rho_val(max_rho_val_), max_y_val(max_y_val_) {}

    TraceFanResult &setFanType(RayFanType type_) {
        this->type = type_;
        return *this;
    }
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_RAYTYPES_H
