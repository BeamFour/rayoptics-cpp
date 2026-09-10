// C++ port of org.redukti.plotter: Colors, SpotDiagram, GeoMTFPlot,
// GeoMTFByFieldPlot and RayAberrationPlot.
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
    /** get rgb color associated with wavelen */
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

    /**
     * Supply a radius for the plt - if not supplied this is set to
     * 1000 * computed max_radius.
     * Example value is 600. that covers most lenses, past and present
     * but modern lenses tend to have much smaller spot diagrams.
     */
    /** Java takes a nullable Double; absent means "use the spot's max radius". */
    std::string plot(std::optional<double> radius) const;
};

class GeoMTFPlot {
public:
    std::shared_ptr<const rayoptics::specs::FieldSnapshot> fld;
    const rayoptics::analysis::MonochromaticGeometricMTF *geo_mtf;

    GeoMTFPlot(const rayoptics::specs::FieldSnapshot &fld_,
               const rayoptics::analysis::MonochromaticGeometricMTF &geo_mtf_)
        : fld(std::make_shared<const rayoptics::specs::FieldSnapshot>(fld_)), geo_mtf(&geo_mtf_) {}

    GeoMTFPlot(const rayoptics::specs::Field &fld_,
               const rayoptics::analysis::MonochromaticGeometricMTF &geo_mtf_)
        : GeoMTFPlot(rayoptics::specs::FieldSnapshot(fld_), geo_mtf_) {}

    std::string plot() const;
};

class GeoMTFByFieldPlot {
public:
    std::vector<rayoptics::analysis::MTFResultByFreq> mtfs_by_freq;
    std::vector<double> fields;

    GeoMTFByFieldPlot(std::vector<rayoptics::analysis::MTFResultByFreq> mtfs_by_freq_,
                      const std::vector<double> &fields_)
        // Names of the above, in the same order, for describing the plot in text
        : mtfs_by_freq(std::move(mtfs_by_freq_)), fields(fields_) {}

    /**
     * Describes which color the plot gives each frequency, e.g. for the default
     * frequencies "10=red,30=blue,50=black". Kept next to FREQ_COLORS so a
     * report can never claim a color the plot did not actually use.
     */
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
