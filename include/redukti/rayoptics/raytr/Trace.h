// C++ port of org.redukti.rayoptics.raytr.Trace
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_TRACE_H
#define REDUKTI_RAYOPTICS_RAYTR_TRACE_H

#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/seq/SurfaceData.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}
namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::raytr {

/** Java's `record GaussianQuadraturePoint(Vector2 pupil, double weight)`. */
class GaussianQuadraturePoint {
public:
    mathlib::Vector2 pupil;
    double weight;

    GaussianQuadraturePoint(const mathlib::Vector2 &pupil_, double weight_)
        : pupil(pupil_), weight(weight_) {}
};

/** Java's `record RayDataFrame(...)` neighbours, kept with their owner. */
class RayDataFrame {
public:
    std::vector<mathlib::Vector3> inc_pt;
    std::vector<mathlib::Vector3> after_dir;
    std::vector<double> after_dst;
    std::vector<mathlib::Vector3> normal;

    explicit RayDataFrame(const std::vector<RaySeg> &raySegList);
};

class RayDataFrameByField {
public:
    specs::Field *fld;
    std::vector<RayDataFrame> frames;

    RayDataFrameByField(specs::Field *fld_, std::vector<RayDataFrame> frames_)
        : fld(fld_), frames(std::move(frames_)) {}
};

class Trace {
public:
    static RayResult trace_ray(optical::OpticalModel *opt_model,
                               const mathlib::Vector2 &pupil, specs::Field &fld,
                               double wvl, TraceOptions &trace_options);

    /**
     * Trace one pupil coordinate, returning the failure as data rather than
     * letting the exception escape.
     */
    static RayResult trace_safe(optical::OpticalModel *opt_model,
                                const mathlib::Vector2 &pupil, specs::Field &fld,
                                double wvl, const TraceOptions &trace_options);

    static std::shared_ptr<const RayPkg> trace(seq::SequentialModel *seq_model,
                                               const mathlib::Vector3 &pt0,
                                               const mathlib::Vector3 &dir0, double wvl,
                                               const TraceOptions &trace_options);

    static std::shared_ptr<const RayPkg> trace_base(optical::OpticalModel *opt_model,
                                                    const std::vector<double> &pupil,
                                                    specs::Field &fld, double wvl,
                                                    const TraceOptions &trace_options);

    static RayResultWithStartCoord get_1d_solution(seq::SequentialModel *seq_model,
                                                   std::optional<int> ifcx,
                                                   const mathlib::Vector3 &pt0,
                                                   double dist, double wvl,
                                                   double y_target, bool not_wa);

    static RayResultWithStartCoord get_2d_solution(seq::SequentialModel *seq_model,
                                                   std::optional<int> ifcx,
                                                   const mathlib::Vector3 &pt0,
                                                   double dist, double wvl,
                                                   const std::vector<double> &xy_target,
                                                   bool not_wa);

    static RayResultWithStartCoord iterate_ray(optical::OpticalModel *opt_model,
                                               std::optional<int> ifcx,
                                               const std::vector<double> &xy_target,
                                               specs::Field &fld, double wvl);

    static std::vector<std::shared_ptr<const RayPkg>> trace_boundary_rays_at_field(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl,
        TraceOptions &trace_options);

    static std::map<std::string, std::shared_ptr<const RayPkg>> boundary_ray_dict(
        optical::OpticalModel *opt_model,
        const std::vector<std::shared_ptr<const RayPkg>> &rim_rays);

    static std::vector<std::vector<std::shared_ptr<const RayPkg>>> trace_boundary_rays(
        optical::OpticalModel *opt_model, TraceOptions &trace_options);

    static std::vector<RayDataFrame> trace_ray_list_at_field(
        optical::OpticalModel *opt_model, const std::vector<std::vector<double>> &ray_list,
        specs::Field &fld, double wvl, double foc, TraceOptions &trace_options);

    static RayDataFrameByField trace_field(optical::OpticalModel *opt_model,
                                           specs::Field &fld, double wvl, double foc);

    static std::vector<RayDataFrameByField> trace_all_fields(
        optical::OpticalModel *opt_model);

    static std::shared_ptr<const ChiefRayPkg> trace_chief_ray(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc);

    static void apply_paraxial_vignetting(optical::OpticalModel *opt_model);

    static std::shared_ptr<const ChiefRayPkg> get_chief_ray_pkg(
        optical::OpticalModel *opt_model, specs::Field &fld, double wvl, double foc);

    static RefSphereCR setup_pupil_coords(optical::OpticalModel *opt_model,
                                          specs::Field &fld, double wvl, double foc,
                                          std::optional<mathlib::Vector2> image_pt,
                                          std::optional<mathlib::Vector2> image_delta);

    /** Iterate the chief ray for a field to the stop surface. */
    static AimInfo aim_chief_ray(optical::OpticalModel *opt_model, specs::Field &fld,
                                 std::optional<double> wvl);

    static std::vector<GridItem> trace_fan(optical::OpticalModel *opt_model,
                                           const TraceFanDef &fan_rng, specs::Field &fld,
                                           double wvl, double foc, bool append_if_none,
                                           ImageFilter *img_filter,
                                           const TraceOptions &trace_options);

    static std::vector<GridItem> trace_grid(optical::OpticalModel *opt_model,
                                            const TraceGridDef &grid_rng,
                                            specs::Field &fld, double wvl, double foc,
                                            ImageFilter *img_filter, bool append_if_none,
                                            const TraceOptions &trace_options);

    static std::vector<GridItem> trace_rings(optical::OpticalModel *opt_model,
                                             const TraceRingsDef &grid_rng,
                                             specs::Field &fld, double wvl, double foc,
                                             ImageFilter *img_filter,
                                             bool append_if_none,
                                             const TraceOptions &trace_options);

    static std::vector<GridItem> trace_gaussian_quadrature(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, specs::Field &fld, double wvl, double foc,
        ImageFilter *img_filter, bool append_if_none,
        const TraceOptions &trace_options);

    static std::vector<ContrastRayTriplet> trace_contrast(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, specs::Field &fld, double wvl,
        const TraceOptions &trace_options);

    static std::vector<ContrastRayTriplet> trace_contrast(
        optical::OpticalModel *opt_model, const TraceRingsDef &grid_rng,
        std::optional<int> num_spokes, const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift,
        std::optional<mathlib::Vector2> sagittal_exit_shift,
        std::optional<mathlib::Vector2> tangential_exit_shift, specs::Field &fld,
        double wvl, const TraceOptions &trace_options, bool aim_exit_pupil);

    static std::vector<GaussianQuadraturePoint> generate_contrast_quadrature(
        const TraceRingsDef &grid_rng, std::optional<int> num_spokes,
        const mathlib::Vector2 &sagittal_shift,
        const mathlib::Vector2 &tangential_shift, specs::Field &fld);

    static bool inside_vignetted_pupil(const mathlib::Vector2 &pupil,
                                       const specs::Field &fld);

    static std::vector<GaussianQuadraturePoint> generate_gaussian_quadrature(
        const TraceRingsDef &grid_rng, int num_rings, std::optional<int> num_spokes);

    // ---- "raw" variants, which iterate over an explicit path -----------------

    static RayResultWithStartCoord get_1d_solution_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const mathlib::Vector3 &pt0, double dist, double wvl, double y_target,
        bool not_wa);

    static RayResultWithStartCoord get_2d_solution_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const mathlib::Vector3 &pt0, double dist, double wvl,
        const std::vector<double> &xy_target, bool not_wa);

    static RayResultWithStartCoord iterate_ray_raw(
        const std::vector<seq::PathSeg> &pthlist, std::optional<int> ifcx,
        const std::vector<double> &xy_target, const mathlib::Vector3 &pt0,
        const mathlib::Vector3 &d0, double obj2pup_dist, double eprad, double wvl,
        bool not_wa);

    static void list_ray(std::string &sb, const RayPkg &ray_pkg,
                         const std::optional<math::Tfm3d> &tfrms,
                         std::optional<int> start);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_TRACE_H
