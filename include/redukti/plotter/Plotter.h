// C++ port of org.redukti.plotter: Colors, SpotDiagram, GeoMTFPlot,
// GeoMTFByFieldPlot and RayAberrationPlot.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_PLOTTER_PLOTTER_H
#define REDUKTI_PLOTTER_PLOTTER_H

#include "redukti/Text.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/render/Renderer.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::plotter {

class Colors {
public:
    /** Maps a wavelength in nm to the colour it is drawn in. */
    static render::Rgb get_wavelen_color(double wl);
};

class SpotDiagram {
public:
    /** Borrowed; the caller owns the analysis result. */
    const rayoptics::analysis::SpotAnalysisResult::SpotResultsForField *result;

    explicit SpotDiagram(
        const rayoptics::analysis::SpotAnalysisResult::SpotResultsForField &result_)
        : result(&result_) {}

    /** Java takes a nullable Double; absent means "use the spot's max radius". */
    std::string plot(std::optional<double> radius) const;
};

class GeoMTFPlot {
public:
    const rayoptics::specs::Field *fld;
    const rayoptics::analysis::MonochromaticGeometricMTF *geo_mtf;

    GeoMTFPlot(const rayoptics::specs::Field &fld_,
               const rayoptics::analysis::MonochromaticGeometricMTF &geo_mtf_)
        : fld(&fld_), geo_mtf(&geo_mtf_) {}

    std::string plot() const;
};

class GeoMTFByFieldPlot {
public:
    std::vector<rayoptics::analysis::MTFResultByFreq> mtfs_by_freq;
    std::vector<double> fields;

    GeoMTFByFieldPlot(std::vector<rayoptics::analysis::MTFResultByFreq> mtfs_by_freq_,
                      const std::vector<double> &fields_)
        : mtfs_by_freq(std::move(mtfs_by_freq_)), fields(fields_) {}

    static std::string freq_legend(const std::vector<int> &freqs);

    std::string plot() const;

    /** The CSV the tool writes alongside the SVG. */
    std::string toString() const;
};

class RayAberrationPlot {
public:
    const rayoptics::analysis::RayAberrationResult *ray_aberration_results;

    explicit RayAberrationPlot(
        const rayoptics::analysis::RayAberrationResult &ray_aberration_results_)
        : ray_aberration_results(&ray_aberration_results_) {}

    std::string plot(const rayoptics::raytr::TraceFanResult &fan_result,
                     double yscale) const;

private:
    double auto_y_scale() const;
};

} // namespace redukti::plotter

#endif // REDUKTI_PLOTTER_PLOTTER_H
