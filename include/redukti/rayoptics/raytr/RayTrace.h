// C++ port of org.redukti.rayoptics.raytr.{RayTrace,RayTraceOptions,TraceOptions}
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_RAYTRACE_H
#define REDUKTI_RAYOPTICS_RAYTR_RAYTRACE_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/seq/SurfaceData.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::raytr {

class TraceOptions {
public:
    std::optional<double> pt_inside_fuzz;
    bool check_apertures = false;
    bool apply_vignetting = true;
    PupilType pupil_type = PupilType::REL_PUPIL;
    std::optional<std::string> output_filter;
    std::optional<std::string> rayerr_filter;
    std::optional<mathlib::Vector2> image_pt_2d;
    std::optional<mathlib::Vector2> image_delta;

    TraceOptions copy() const { return *this; }
};

class RayTraceOptions {
public:
    std::optional<int> first_surf;
    std::optional<int> last_surf;
    bool print_details = false;
    double eps = 1.0e-12;
    bool check_apertures = false;
    bool intersect_obj = true;
    bool filter_out_phantoms = false;
    std::optional<double> pt_inside_fuzz;

    RayTraceOptions() = default;
    explicit RayTraceOptions(const TraceOptions &trace_options) {
        this->check_apertures = trace_options.check_apertures;
        this->pt_inside_fuzz = trace_options.pt_inside_fuzz;
    }
};

class RayTrace {
public:
    /**
     * Trace a ray through the sequential model.
     *
     * Throws TraceMissedSurfaceException, TraceTIRException or
     * TraceRayBlockedException, each carrying the ray traced so far in its
     * ray_pkg -- callers routinely catch one and use that partial ray.
     */
    static std::shared_ptr<const RayPkg> trace(seq::SequentialModel *seq_model,
                                               const mathlib::Vector3 &pt0,
                                               const mathlib::Vector3 &dir0, double wvl,
                                               RayTraceOptions &options);

    static std::shared_ptr<const RayPkg> trace(seq::SequentialModel *seq_model,
                                               const mathlib::Vector3 &pt0,
                                               const mathlib::Vector3 &dir0, double wvl);

    static std::shared_ptr<const RayPkg> trace_raw(const std::vector<seq::PathSeg> &path,
                                                   const mathlib::Vector3 &pt0,
                                                   const mathlib::Vector3 &dir0,
                                                   double wvl,
                                                   const RayTraceOptions &options);

    /** Refract the ray at the interface. Throws TraceTIRException. */
    static mathlib::Vector3 bend(const mathlib::Vector3 &d_in,
                                 const mathlib::Vector3 &normal, double n_in,
                                 double n_out);

    /** Reflect the ray at the interface. */
    static mathlib::Vector3 reflect(const mathlib::Vector3 &d_in,
                                    const mathlib::Vector3 &normal);

    /** Distance from the axis, per eq 3.20/3.21. */
    static double eic_distance_from_axis(const mathlib::Vector3 &p,
                                         const mathlib::Vector3 &d, util::ZDir z_dir);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_RAYTRACE_H
