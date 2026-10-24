// C++ port of org.redukti.analysis.ContrastPupilMapTest.
//
// The Java checks identities (losses add up to the phasor, a pure shift loses
// nothing) and counts the plot's circles. Here the maps of the Nikkor model are
// also asserted against values from the Java classes (DumpPupilMap.java, a
// throwaway dumper), within the tolerance their trig allows.
#include "TestHarness.h"

#include "NikkorZ58Model.h"

#include "redukti/analysis/ContrastAnalysis.h"
#include "redukti/analysis/ContrastPupilMap.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/plotter/ContrastPupilMapPlot.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using redukti::analysis::ContrastAnalysis;
using redukti::analysis::ContrastAnalysisResult;
using redukti::analysis::ContrastOptions;
using redukti::analysis::ContrastPupilMap;
using redukti::mathlib::Vector2;
using redukti::plotter::ContrastPupilMapPlot;
namespace Orientation = redukti::rayoptics::util::Orientation;

ContrastAnalysisResult::WavelengthResult block(const std::vector<double> &differences,
                                               const std::vector<double> &weights,
                                               const std::vector<double> &references) {
    std::vector<ContrastAnalysisResult::Sample> samples;
    for (std::size_t i = 0; i < differences.size(); i++)
        samples.emplace_back(Vector2(0.1 * static_cast<double>(i), 0), differences[i],
                             differences[i], weights[i], true, std::nullopt, references[i]);
    return ContrastAnalysisResult::WavelengthResult(587.5618, 0.1, std::move(samples));
}

/** One map and its last point, as DumpPupilMap.java printed them. */
struct ExpectedMap {
    int n;
    double mean_phase, phasor, rms;
    double px, py, loss, average_phase, weight;
};

// Indexed [field][orientation], sagittal first.
const ExpectedMap EXPECTED_MAPS[3][2] = {
    {{32, 5.2726618333691995E-12, 0.9699321829011592, 0.039298783340513156,
      0.941526067178229, -0.02313353862199797, 0.02587720205137234, -3.7201880840813684,
      0.02174092782109086},
     {32, 2.7303092400924268E-11, 0.9699321828928835, 0.039298783345811564,
      0.941526067178229, -0.02313353862199797, 0.0, -3.88175399808225, 0.02174092782109086}},
    {{32, 0.007994411563234214, 0.582319184398872, 0.16770742379585094, 0.8595437621445654,
      -0.023133538621997865, 0.9519494493281953, -1.289002810736335, 0.018417322836358958},
     {32, 0.5187475957297895, 0.8272919986225218, 0.133672716064093, 0.8595437621445654,
      -0.023133538621997865, 0.36332843931671543, 0.9631359491066998, 0.018417322836358958}},
    {{32, 0.022469610568068062, 0.5060612608266882, 0.400290677052209, 0.7126955615659444,
      -0.023133538621997806, 0.37553963838573046, -34.544000417401115, 0.017565883144437286},
     {32, 0.20901400487049934, 0.7773860845067159, 0.12335940256736826, 0.7126955615659444,
      -0.023133538621997806, 0.11389367251927279, -30.99376027889083, 0.017565883144437286}}};

void checkMap(const ContrastPupilMap &map, const ExpectedMap &e, double atol, double rtol) {
    auto close = [&](double actual, double expected) {
        CHECK_CLOSE(actual, expected, atol + rtol * std::abs(expected));
    };
    CHECK_EQ(static_cast<int>(map.points.size()), e.n);
    close(map.mean_phase, e.mean_phase);
    close(map.phasor, e.phasor);
    close(map.rms, e.rms);
    const auto &p = map.points.back();
    close(p.pupil.x, e.px);
    close(p.pupil.y, e.py);
    close(p.loss, e.loss);
    close(p.average_phase, e.average_phase);
    close(p.weight, e.weight);
}

int count(const std::string &s, const std::string &what) {
    int n = 0;
    for (auto at = s.find(what); at != std::string::npos; at = s.find(what, at + 1)) n++;
    return n;
}

} // namespace

TEST(contrast_pupil_map_a_pure_image_shift_loses_no_contrast) {
    ContrastPupilMap map(block({0.3, 0.3, 0.3}, {0.2, 0.3, 0.5}, {0, 0, 0}), Orientation::SAGITTAL);
    CHECK_CLOSE(map.phasor, 1.0, 1e-15);
    for (const auto &p : map.points)
        CHECK_CLOSE(p.loss, 0.0, 1e-15);
}

TEST(contrast_pupil_map_half_a_wave_apart_is_full_loss) {
    // In phase with the mean, and half a wave from it
    ContrastPupilMap map(block({0.0, 0.0, 0.5}, {0.4, 0.4, 0.2}, {0, 0, 0}),
                         Orientation::TANGENTIAL);
    CHECK_CLOSE(map.points[0].loss, 0.0, 1e-15);
    CHECK_CLOSE(map.points[2].loss, 1.0, 1e-15);
}

TEST(contrast_pupil_map_weighted_losses_add_up_to_the_phasor) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto result =
        ContrastAnalysis::eval(model.get(), ContrastOptions(40.0).num_rings(4).num_spokes(8));
    for (std::size_t fi = 0; fi < result.fields.size(); fi++) {
        for (int orientation : {Orientation::SAGITTAL, Orientation::TANGENTIAL}) {
            const auto &blk = result.fields[fi].wavelengths[0];
            ContrastPupilMap map(blk, orientation);
            double loss = 0;
            for (const auto &p : map.points)
                loss += p.weight * p.loss;
            CHECK_CLOSE(1 - 2 * loss, map.phasor, 1e-12);
            CHECK(map.phasor > 0 && map.phasor <= 1);
            // contrast samples record the reference OPD
            for (const auto &p : map.points)
                CHECK(std::isfinite(p.average_phase));
            // The map itself is trig-fed (cos, sin, atan2), and so is the gaussian
            // quadrature pattern under it. At the 45 degree spoke the JVM's
            // sin(0.7853981633974483) is 0.7071067811865475, one ulp below its cos, while
            // MSVC's sin returns the cos value; the pupil y of that spoke then differs by an
            // ulp. On axis the OPD differences are small, so the ulp reaches ~1e-10 of rms
            // and ~2e-11 absolute in a mean phase near zero. Off axis it stays at 1e-12.
            if (fi == 0)
                checkMap(map, EXPECTED_MAPS[fi][orientation], 1e-10, 1e-9);
            else
                checkMap(map, EXPECTED_MAPS[fi][orientation], 1e-14, 1e-12);
        }
    }
}

TEST(contrast_pupil_map_average_phase_is_the_mean_of_the_two_wavefront_phases) {
    // W(p) = 0.25, W(p + s) = 0.25 + 0.5: mean 0.5 waves, so pi radians
    ContrastPupilMap map(block({0.5, 0.0}, {0.5, 0.5}, {0.25, 0.0}), Orientation::SAGITTAL);
    CHECK_CLOSE(map.points[0].average_phase, redukti::mathlib::M::PI, 1e-15);
    CHECK_CLOSE(map.points[1].average_phase, 0.0, 1e-15);
}

TEST(contrast_pupil_map_rejects_a_block_without_valid_samples) {
    std::vector<ContrastAnalysisResult::Sample> samples;
    samples.emplace_back(Vector2::vector2_0, 0, 0, 1, false, std::nullopt);
    ContrastAnalysisResult::WavelengthResult blk(587.5618, 0.1, samples);
    CHECK_THROWS(ContrastPupilMap(blk, Orientation::SAGITTAL), redukti::IllegalArgumentException);
}

TEST(contrast_pupil_map_plot_draws_a_circle_per_lossy_sample) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto result =
        ContrastAnalysis::eval(model.get(), ContrastOptions(40.0).num_rings(4).num_spokes(8));
    const auto &field = result.fields[1];
    const auto &blk = field.wavelengths[0];
    std::vector<ContrastPupilMapPlot::Row> maps(1);
    maps[0][0].emplace(blk, Orientation::SAGITTAL);
    maps[0][1].emplace(blk, Orientation::TANGENTIAL);
    long lossy = 0;
    for (const auto &row : maps)
        for (const auto &map : row)
            for (const auto &p : map->points)
                if (p.loss > 0) lossy++;
    std::string svg = ContrastPupilMapPlot(field.field, {40}, maps).plot();
    CHECK(svg.find("<svg") != std::string::npos);
    // RendererSvg writes circles as ellipses, and defines one more for its point marker
    CHECK_EQ(count(svg, "<ellipse"), static_cast<int>(lossy + 1));
    CHECK_EQ(count(svg, "<ellipse"), 65); // the Java's count
}
