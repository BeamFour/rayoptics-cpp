// C++ port of SpotOptions, SpotAnalysisResult and SpotAnalysis.
#include "redukti/rayoptics/analysis/SpotAnalysis.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"

#include <cmath>
#include <limits>

namespace redukti::rayoptics::analysis {

using mathlib::Vector2;
using mathlib::Vector3;
using raytr::GridItem;
using raytr::RayPkg;
using raytr::TraceGridByWvl;

// ---------------------------------------------------------------------------
// SpotOptions
// ---------------------------------------------------------------------------

SpotOptions::SpotOptions(bool useGaussGuadrature) {
    _trace_options.check_apertures = true;
    if (useGaussGuadrature)
        use_gaussian_quadrature();
    else
        use_hexapolar();
}

SpotOptions &SpotOptions::num_rays(int rays) {
    this->_num_rays_or_rings = rays;
    return *this;
}

SpotOptions &SpotOptions::num_rings(int rings) {
    this->_num_rays_or_rings = rings;
    return *this;
}

SpotOptions &SpotOptions::use_hexapolar() {
    _pattern = PATTERN_HEXAPOLAR;
    _num_rays_or_rings = 64;
    return *this;
}

SpotOptions &SpotOptions::use_grid() {
    _pattern = PATTERN_GRID;
    _num_rays_or_rings = 64;
    return *this;
}

SpotOptions &SpotOptions::use_gaussian_quadrature() {
    _pattern = PATTERN_GAUSS_QUADRATURE;
    _num_rays_or_rings = 14;
    _num_spokes = 20;
    return *this;
}

SpotOptions &SpotOptions::num_spokes(std::optional<int> spokes) {
    if (spokes.has_value() && *spokes < 3)
        throw IllegalArgumentException("Gaussian quadrature requires at least 3 spokes");
    this->_num_spokes = spokes;
    return *this;
}

SpotOptions &SpotOptions::inner_pupil_radius(double radius) {
    if (!std::isfinite(radius) || radius < 0.0 || radius >= 1.0)
        throw IllegalArgumentException(
            "Inner pupil radius must be finite and in [0, 1)");
    this->_inner_pupil_radius = radius;
    return *this;
}

SpotOptions &SpotOptions::use_centroid(bool value) {
    this->_use_centroid = value;
    return *this;
}

SpotOptions &SpotOptions::append_failed_rays(bool value) {
    this->_append_failed_rays = value;
    return *this;
}

SpotOptions &SpotOptions::check_apertures(bool value) {
    this->_trace_options.check_apertures = value;
    return *this;
}

// ---------------------------------------------------------------------------
// SpotAnalysisResult::SpotResultsForField
// ---------------------------------------------------------------------------

SpotAnalysisResult::SpotResultsForField::SpotResultsForField(
    specs::Field *fld_, std::vector<TraceGridByWvl> trace_results_, double ref_wvl,
    bool use_centroid)
    : fld(fld_), image_pt(fld_->ref_sphere->image_pt),
      trace_results(std::move(trace_results_)) {
    // The traced grids have to be in their final home before the intercepts are
    // built: each SpotIntercepts keeps a pointer into this vector.
    std::optional<Vector2> centroid;
    for (const auto &result : trace_results) {
        SpotIntercepts s(result);
        if (result.wvl == ref_wvl && use_centroid)
            centroid = s.compute_centroid();
        intercepts.push_back(std::move(s));
    }
    if (centroid.has_value()) {
        for (auto &intercept : intercepts)
            intercept.adjust_to_centroid(*centroid);
    }
    computeMeanMax();
}

void SpotAnalysisResult::SpotResultsForField::computeMeanMax() {
    max_radius = 0;
    mean_radius = 0;
    double totalWeight = 0.0;
    for (const auto &results : intercepts) {
        for (std::size_t i = 0; i < results.x.size(); i++) {
            if (!results.valid[i])
                continue;
            double r = results.x[i] * results.x[i] + results.y[i] * results.y[i];
            double l = std::sqrt(r);
            if (l > max_radius) {
                max_radius = l;
            }
            mean_radius += results.weights[i] * l * l;
            totalWeight += results.weights[i];
        }
    }
    mean_radius = totalWeight > 0.0 ? std::sqrt(mean_radius / totalWeight)
                                    : std::numeric_limits<double>::quiet_NaN();
}

std::string SpotAnalysisResult::SpotResultsForField::toString() const {
    return "Field " + fld->toString() + " mean radius " +
           doubleToString(get_mean_radius()) + " max radius " +
           doubleToString(get_max_radius());
}

// ---------------------------------------------------------------------------
// SpotAnalysisResult
// ---------------------------------------------------------------------------

SpotAnalysisResult &SpotAnalysisResult::add(specs::Field *fld,
                                            std::vector<TraceGridByWvl> trace_results,
                                            double ref_wvl) {
    spot_results.push_back(
        SpotResultsForField(fld, std::move(trace_results), ref_wvl, use_centroid));
    return *this;
}

std::vector<double> SpotAnalysisResult::fields() const {
    std::vector<double> fields(spot_results.size());
    for (std::size_t i = 0; i < spot_results.size(); i++)
        fields[i] = spot_results[i].fld->y;
    return fields;
}

std::vector<MTFResultByFreq> SpotAnalysisResult::computeMTFs(
    const std::vector<int> &freqs) const {
    std::vector<PolyMTF> mtfs;
    for (std::size_t i = 0; i < spot_results.size(); i++) {
        const auto &spotFld = spot_results[i];
        auto cfg = spotFld.mtfHistogramConfig();
        std::optional<PolyMTF> polyMtfForField;
        for (const auto &intercepts : spotFld.intercepts) {
            MonochromaticGeometricMTF mtf(intercepts, cfg);
            if (!polyMtfForField.has_value())
                polyMtfForField.emplace(mtf.mtf.fft_size, mtf.h2d.pixel_size);
            polyMtfForField->add(mtf.mtf, 1.0);
        }
        if (polyMtfForField.has_value()) {
            polyMtfForField->compute();
            mtfs.push_back(std::move(*polyMtfForField));
        }
    }
    std::vector<MTFResultByFreq> mtfResults;
    for (auto freq : freqs)
        mtfResults.push_back(MTFResultByFreq(mtfs, freq));
    return mtfResults;
}

std::string SpotAnalysisResult::toString() const {
    std::string sb;
    for (const auto &result : spot_results) {
        sb += result.toString();
        sb += "\n";
    }
    return sb;
}

// ---------------------------------------------------------------------------
// SpotAnalysis
// ---------------------------------------------------------------------------

std::optional<GridItem> SpotAnalysis::spot(const Vector2 &p, int wi,
                                           const std::shared_ptr<const RayPkg> &ray_pkg,
                                           specs::Field &fld, double wvl, double foc) {
    (void)p;
    (void)wi;
    (void)wvl;
    if (ray_pkg != nullptr) {
        auto image_pt = fld.ref_sphere->image_pt;
        const auto &ray = ray_pkg->ray;
        auto dist = foc / util::Lists::get(ray, -1).d.z;
        auto defocussed_pt =
            util::Lists::get(ray, -1).p.plus(util::Lists::get(ray, -1).d.times(dist));
        auto t_abr = defocussed_pt.minus(image_pt);
        return GridItem(t_abr.project_xy(), ray_pkg);
    } else
        return std::nullopt;
}

std::vector<TraceGridByWvl> SpotAnalysis::eval_grid(
    optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rays,
    const raytr::TraceOptions &trace_options) {
    auto seq_model = opt_model->seq_model.get();
    return seq_model->trace_grid(SpotAnalysis::spot, fi, wl, num_rays, false,
                                 trace_options);
}

std::vector<TraceGridByWvl> SpotAnalysis::eval_rings(
    optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rays,
    const raytr::TraceOptions &trace_options) {
    auto seq_model = opt_model->seq_model.get();
    return seq_model->trace_rings(SpotAnalysis::spot, fi, wl, num_rays, false,
                                  trace_options);
}

std::vector<TraceGridByWvl> SpotAnalysis::eval_gaussian_quadrature(
    optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rings,
    std::optional<int> num_spokes, double innerPupilRadius, bool append_if_none,
    const raytr::TraceOptions &trace_options) {
    return opt_model->seq_model->trace_gaussian_quadrature(
        SpotAnalysis::spot, fi, wl, num_rings, num_spokes, innerPupilRadius,
        append_if_none, trace_options);
}

SpotAnalysisResult SpotAnalysis::eval(optical::OpticalModel *opt_model,
                                      const SpotOptions &options) {
    auto num_rays = options._num_rays_or_rings;
    auto &trace_options = options._trace_options;
    SpotAnalysisResult result(options._use_centroid);
    auto fov = opt_model->optical_spec->fov.get();
    auto ref_wvl = fov->optical_spec->wvls->central_wvl();
    for (std::size_t fi = 0; fi < fov->fields.size(); fi++) {
        specs::Field *f = fov->fields[fi].get();
        if (options.is_grid())
            result.add(f,
                       eval_grid(opt_model, static_cast<int>(fi), std::nullopt, num_rays,
                                 trace_options),
                       ref_wvl);
        else if (options.is_gauss_quadrature())
            result.add(f,
                       eval_gaussian_quadrature(opt_model, static_cast<int>(fi),
                                                std::nullopt, num_rays,
                                                options._num_spokes,
                                                options._inner_pupil_radius,
                                                options._append_failed_rays,
                                                trace_options),
                       ref_wvl);
        else
            result.add(f,
                       eval_rings(opt_model, static_cast<int>(fi), std::nullopt, num_rays,
                                  trace_options),
                       ref_wvl);
    }
    return result;
}

} // namespace redukti::rayoptics::analysis
