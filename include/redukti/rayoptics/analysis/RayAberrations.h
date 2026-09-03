// C++ port of org.redukti.rayoptics.analysis.RayAberrationResult,
// TransverseRayAberrationAnalysis and WavefrontAberrationAnalysis.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_RAYABERRATIONS_H
#define REDUKTI_RAYOPTICS_ANALYSIS_RAYABERRATIONS_H

#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::analysis {

/** The fans traced for every field and meridian, in trace order. */
class RayAberrationResult {
public:
    std::vector<raytr::TraceFanResult> results;

    void add(raytr::TraceFanResult fan_result) {
        results.push_back(std::move(fan_result));
    }

    /** Null when no fan matches; borrowed from `results`. */
    const raytr::TraceFanPoints *get_fans(int fi, int xy, double wvl) const;

    std::string list_ray_fans() const;
};

class TransverseRayAberrationAnalysis {
public:
    /** The fan callback. Null means the ray produced no aberration value. */
    static std::optional<double> ray_abr(optical::OpticalModel *opt_model,
                                         const mathlib::Vector2 &p, int xy,
                                         const std::shared_ptr<const raytr::RayPkg> &ray_pkg,
                                         specs::Field &fld, double wvl, double foc);

    static raytr::TraceFanResult eval_abr_fan(optical::OpticalModel *opt_model, int fi,
                                              int xy, int num_rays, bool append_if_none,
                                              const raytr::TraceOptions &trace_options);

    static RayAberrationResult eval(optical::OpticalModel *opt_model, int num_rays,
                                    bool append_if_none,
                                    const raytr::TraceOptions &trace_options);
};

class WavefrontAberrationAnalysis {
public:
    /** Optical path difference in waves. Null means no value at this point. */
    static std::optional<double> opd(optical::OpticalModel *opt_model,
                                     const mathlib::Vector2 &p, int xy,
                                     const std::shared_ptr<const raytr::RayPkg> &ray_pkg,
                                     specs::Field &fld, double wvl, double foc);

    static raytr::TraceFanResult eval_opd_fan(optical::OpticalModel *opt_model, int fi,
                                              int xy, int num_rays, bool append_if_none,
                                              const raytr::TraceOptions &trace_options);

    static RayAberrationResult eval(optical::OpticalModel *opt_model, int num_rays,
                                    bool append_if_none,
                                    const raytr::TraceOptions &trace_options);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_RAYABERRATIONS_H
