// C++ port of org.redukti.rayoptics.analysis.SpotOptions, SpotAnalysisResult
// and SpotAnalysis.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_SPOTANALYSIS_H
#define REDUKTI_RAYOPTICS_ANALYSIS_SPOTANALYSIS_H

#include "redukti/rayoptics/analysis/Histogram.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/SpotIntercepts.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::analysis {

/** Fluent configuration for a spot diagram, mirroring the Java builder. */
class SpotOptions {
public:
    static constexpr int PATTERN_HEXAPOLAR = 1;
    static constexpr int PATTERN_GAUSS_QUADRATURE = 2;
    static constexpr int PATTERN_GRID = 3;

    raytr::TraceOptions _trace_options;
    int _pattern = 0;
    bool _use_centroid = true;
    bool _append_failed_rays = false;
    int _num_rays_or_rings = 0;
    std::optional<int> _num_spokes;
    double _inner_pupil_radius = 0.0;

    explicit SpotOptions(bool useGaussGuadrature);

    SpotOptions() : SpotOptions(false) {}

    SpotOptions &num_rays(int rays);
    SpotOptions &num_rings(int rings);
    SpotOptions &use_hexapolar();
    SpotOptions &use_grid();
    SpotOptions &use_gaussian_quadrature();
    SpotOptions &num_spokes(std::optional<int> spokes);
    SpotOptions &inner_pupil_radius(double radius);
    SpotOptions &use_centroid(bool value);
    SpotOptions &append_failed_rays(bool value);
    SpotOptions &check_apertures(bool value);

    bool is_gauss_quadrature() const { return _pattern == PATTERN_GAUSS_QUADRATURE; }
    bool is_hexapolar() const { return _pattern == PATTERN_HEXAPOLAR; }
    bool is_grid() const { return _pattern == PATTERN_GRID; }
};

class SpotAnalysisResult {
public:
    bool use_centroid;

    class SpotResultsForField {
    public:
        std::shared_ptr<const specs::FieldSnapshot> fld;
        mathlib::Vector3 image_pt;
        std::vector<raytr::TraceGridByWvl> trace_results;
        std::vector<SpotIntercepts> intercepts;
        double max_radius = 0;
        double mean_radius = 0;

        SpotResultsForField(specs::Field *fld_,
                            std::vector<raytr::TraceGridByWvl> trace_results_,
                            double ref_wvl, bool use_centroid);

        std::string toString() const;

        double get_max_radius() const { return max_radius * 1000; }
        double get_mean_radius() const { return mean_radius * 1000; }

        Histogram::Config mtfHistogramConfig() const {
            return Histogram::adaptiveConfig(max_radius);
        }

    private:
        void computeMeanMax();
    };

    std::vector<SpotResultsForField> spot_results;

    explicit SpotAnalysisResult(bool use_centroid_) : use_centroid(use_centroid_) {}

    SpotAnalysisResult &add(specs::Field *fld,
                            std::vector<raytr::TraceGridByWvl> trace_results,
                            double ref_wvl);

    std::vector<double> fields() const;

    std::vector<MTFResultByFreq> computeMTFs(const std::vector<int> &freqs) const;

    std::string toString() const;
};

class SpotAnalysis {
public:
    /**
     * The grid callback. Returns null when the ray failed, which is what stops
     * a failed trace from contributing a point -- see raytr::ImageFilter.
     */
    static std::optional<raytr::GridItem> spot(const mathlib::Vector2 &p, int wi,
                                               const std::shared_ptr<const raytr::RayPkg> &ray_pkg,
                                               specs::Field &fld, double wvl, double foc);

    static std::vector<raytr::TraceGridByWvl> eval_grid(
        optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rays,
        const raytr::TraceOptions &trace_options);

    static std::vector<raytr::TraceGridByWvl> eval_rings(
        optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rays,
        const raytr::TraceOptions &trace_options);

    static std::vector<raytr::TraceGridByWvl> eval_gaussian_quadrature(
        optical::OpticalModel *opt_model, int fi, std::optional<int> wl, int num_rings,
        std::optional<int> num_spokes, double innerPupilRadius, bool append_if_none,
        const raytr::TraceOptions &trace_options);

    static SpotAnalysisResult eval(optical::OpticalModel *opt_model,
                                   const SpotOptions &options);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_SPOTANALYSIS_H
