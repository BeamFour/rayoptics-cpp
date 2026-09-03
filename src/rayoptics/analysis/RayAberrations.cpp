// C++ port of RayAberrationResult, TransverseRayAberrationAnalysis and
// WavefrontAberrationAnalysis.
#include "redukti/rayoptics/analysis/RayAberrations.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/raytr/WaveAbr.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"
#include "redukti/rayoptics/util/Orientation.h"

namespace redukti::rayoptics::analysis {

namespace Orientation = util::Orientation;

using mathlib::Vector2;
using raytr::RayFanType;
using raytr::RayPkg;
using raytr::TraceFanPoints;
using raytr::TraceFanResult;

// ---------------------------------------------------------------------------
// RayAberrationResult
// ---------------------------------------------------------------------------

const TraceFanPoints *RayAberrationResult::get_fans(int fi, int xy, double wvl) const {
    for (const auto &result : results) {
        if (result.fi == fi && result.xy == xy) {
            for (const auto &fan : result.fans) {
                if (fan.wvl == wvl)
                    return &fan;
            }
        }
    }
    return nullptr;
}

std::string RayAberrationResult::list_ray_fans() const {
    std::string sb;
    for (const auto &result : results) {
        sb += result.fld->toString();
        sb += " xy=";
        sb += std::to_string(result.xy);
        sb += "\n";
        for (std::size_t i = 0; i < result.fans.size(); i++) {
            const auto &fan = result.fans[i];
            if (i > 0)
                sb += ",";
            sb += doubleToString(fan.wvl);
        }
        sb += "\n";
        std::size_t row_count = result.fans[0].fan_y.size();
        for (std::size_t i = 0; i < row_count; i++) {
            for (std::size_t j = 0; j < result.fans.size(); j++) {
                const auto &fan = result.fans[j];
                if (j > 0)
                    sb += ",";
                if (fan.fan_y.size() >= row_count)
                    // Java appends the Double, so a null prints as "null".
                    sb += fan.fan_y[i].has_value() ? doubleToString(*fan.fan_y[i])
                                                   : std::string("null");
                else
                    sb += "ERR";
            }
            sb += "\n";
        }
    }
    return sb;
}

// ---------------------------------------------------------------------------
// TransverseRayAberrationAnalysis
// ---------------------------------------------------------------------------

std::optional<double> TransverseRayAberrationAnalysis::ray_abr(
    optical::OpticalModel *opt_model, const Vector2 &p, int xy,
    const std::shared_ptr<const RayPkg> &ray_pkg, specs::Field &fld, double wvl,
    double foc) {
    (void)opt_model;
    (void)p;
    (void)wvl;
    // Java guards on `ray_pkg.ray != null`. A RayPkg is only ever built with a
    // segment list, and here it holds a vector, so the guard is always true.
    auto image_pt = fld.ref_sphere->image_pt;
    const auto &ray = ray_pkg->ray;
    auto dist = foc / util::Lists::get(ray, -1).d.z;
    auto defocused_pt =
        util::Lists::get(ray, -1).p.plus(util::Lists::get(ray, -1).d.times(dist));
    auto t_abr = defocused_pt.minus(image_pt);
    return t_abr.v(xy);
}

TraceFanResult TransverseRayAberrationAnalysis::eval_abr_fan(
    optical::OpticalModel *opt_model, int fi, int xy, int num_rays, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    auto seq_model = opt_model->seq_model.get();
    return seq_model
        ->trace_fan(TransverseRayAberrationAnalysis::ray_abr, fi, xy, num_rays,
                    append_if_none, trace_options)
        .setFanType(RayFanType::TransverseRayFan);
}

RayAberrationResult TransverseRayAberrationAnalysis::eval(
    optical::OpticalModel *opt_model, int num_rays, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    RayAberrationResult result;
    auto fov = opt_model->optical_spec->fov.get();
    for (std::size_t fi = 0; fi < fov->fields.size(); fi++) {
        for (int xy = 0; xy < Orientation::COUNT; xy++) {
            result.add(eval_abr_fan(opt_model, static_cast<int>(fi), xy, num_rays,
                                    append_if_none, trace_options));
        }
    }
    return result;
}

// ---------------------------------------------------------------------------
// WavefrontAberrationAnalysis
// ---------------------------------------------------------------------------

std::optional<double> WavefrontAberrationAnalysis::opd(
    optical::OpticalModel *opt_model, const Vector2 &p, int xy,
    const std::shared_ptr<const RayPkg> &ray_pkg, specs::Field &fld, double wvl,
    double foc) {
    (void)p;
    (void)xy;
    auto convert_to_waves = 1.0 / opt_model->nm_to_sys_units(wvl);
    // See the note in ray_abr: the Java null check on the segment list cannot
    // fail here.
    const auto &fod = opt_model->optical_spec->parax_data->fod;
    auto ops = raytr::WaveAbr::wave_abr_full_calc(fod, fld, wvl, foc, ray_pkg,
                                                  fld.chief_ray, fld.ref_sphere);
    return convert_to_waves * ops;
}

TraceFanResult WavefrontAberrationAnalysis::eval_opd_fan(
    optical::OpticalModel *opt_model, int fi, int xy, int num_rays, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    auto seq_model = opt_model->seq_model.get();
    return seq_model
        ->trace_fan(WavefrontAberrationAnalysis::opd, fi, xy, num_rays, append_if_none,
                    trace_options)
        .setFanType(RayFanType::OpticalPathDifference);
}

RayAberrationResult WavefrontAberrationAnalysis::eval(
    optical::OpticalModel *opt_model, int num_rays, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    RayAberrationResult result;
    auto fov = opt_model->optical_spec->fov.get();
    for (std::size_t fi = 0; fi < fov->fields.size(); fi++) {
        for (int xy = 0; xy < Orientation::COUNT; xy++) {
            result.add(eval_opd_fan(opt_model, static_cast<int>(fi), xy, num_rays,
                                    append_if_none, trace_options));
        }
    }
    return result;
}

} // namespace redukti::rayoptics::analysis
