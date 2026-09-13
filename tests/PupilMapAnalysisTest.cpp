// C++ port of org.redukti.rayoptics.analysis.PupilMapAnalysisTest.
//
// The pupil map on an ultra wide angle lens, where the vignetting factors expand the pupil
// rather than shrink it: the map is what shows that the bundle off axis reaches well
// outside the nominal pupil, which a sampling pattern confined to the unit circle misses.
#include "TestHarness.h"

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/plotter/Plotter.h"
#include "redukti/rayoptics/analysis/PupilMapAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/spec/Prescription.h"

#include <memory>
#include <string>
#include <vector>

namespace {

using redukti::importers::OpticalBenchDataImporter;
using redukti::plotter::PupilMapPlot;
using redukti::rayoptics::analysis::PupilMapAnalysis;
using redukti::rayoptics::optical::OpticalModel;
using redukti::spec::Prescription;
using redukti::spec::RayOpticsModelBuilder;
using redukti::spec::VigType;

/** 114 degrees; the factors reach 1.42 in x at full field, and 0.65 in -y. */
const char *const EF14 = REDUKTI_EXAMPLES_DIR "canon-ef14mm-f2.8L/JP1993-034592_Example02P.txt";

const std::vector<double> FIELDS{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
/** Coarse enough to be quick; a boundary is then placed to about a 20th of the pupil. */
const int SAMPLES = 41;

std::unique_ptr<OpticalModel> model() {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(EF14);
    Prescription prescription = Prescription::build_prescription(specs, true, false, false);
    return RayOpticsModelBuilder(prescription)
        .build_optical_model(true, FIELDS, false, VigType::SetPupil, true, 0);
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

TEST(pupilMap_measuresTheBundleTheLensActuallyPasses) {
    auto opm = model();
    auto result = PupilMapAnalysis::eval(opm.get(), SAMPLES, std::vector<int>{0, 10});
    CHECK_EQ(result.maps.size(), static_cast<std::size_t>(2));
    if (result.maps.size() != 2)
        return;
    const auto &axial = result.maps[0];
    const auto &full = result.maps[1];

    // On axis nothing is vignetted: the whole nominal pupil passes, and nothing outside
    // it does.
    CHECK_CLOSE(axial.nominal_passed, 1.0, 1e-9);
    double step = 2 * axial.reach / (SAMPLES - 1);
    CHECK_CLOSE(axial.passed_x, 1.0, step);
    CHECK_CLOSE(axial.passed_y, 1.0, step);

    // At full field the factors expand the pupil in x and shrink it below the axis, and
    // the traced bundle agrees with them: this is the measurement that says the nominal
    // pupil is the wrong bundle on a lens like this.
    CHECK_CLOSE(PupilMapAnalysis::scale(full.fld->vlx), 1.418, 1e-3);
    CHECK_CLOSE(PupilMapAnalysis::scale(full.fld->vly), 0.654, 1e-3);
    CHECK_CLOSE(full.passed_x, PupilMapAnalysis::scale(full.fld->vlx), step);
    // Part of the nominal pupil is blocked at full field.
    CHECK(full.nominal_passed < 0.95);
    // The bundle reaches outside the nominal pupil.
    CHECK(full.passed_x > 1.0 + step);
}

TEST(pupilMap_bothMappingsDescribeTheMeasuredRegion) {
    auto opm = model();
    auto result = PupilMapAnalysis::eval(opm.get(), SAMPLES, std::vector<int>{0, 6, 10});
    for (const auto &map : result.maps) {
        for (const auto *quality : {&map.piecewise_quality, &map.ellipse_quality}) {
            // Neither mapping is a fit to the traced boundary, so this is a sanity bound
            // rather than a measure of quality: the region is in the right place and
            // mostly usable.
            if (!(quality->sampled > 0.85 && quality->covered > 0.85))
                ::redukti::test::reportFailure(
                    __FILE__, __LINE__,
                    "field " + redukti::doubleToString(map.fld->yv()) + ": " +
                        quality->toString());
        }
    }
}

TEST(pupilMap_expandsToMeasureTheBundleWithClearedFactors) {
    auto opm = model();
    int n = PupilMapAnalysis::DEFAULT_NUM_SAMPLES;
    auto referenceResult = PupilMapAnalysis::eval(opm.get(), n, std::vector<int>{10});
    const auto &reference = referenceResult.maps.at(0);
    long passed = 0, nominalPassed = 0;
    for (const auto &s : reference.samples) {
        if (!s.passed)
            continue;
        passed++;
        if (s.x * s.x + s.y * s.y <= 1.0)
            nominalPassed++;
    }
    double expectedCoverage = static_cast<double>(nominalPassed) / passed;

    auto &fld = *opm->optical_spec->fov->fields[10];
    fld.clear_vignetting();
    double initialReach = PupilMapAnalysis::reach_for(fld);
    auto mapped = PupilMapAnalysis::eval(opm.get(), n, std::vector<int>{10});
    const auto &map = mapped.maps.at(0);

    CHECK(map.reach > initialReach);
    // Transmitted light outside the initial square is measured.
    CHECK(map.passed_x > 1.3);
    CHECK_CLOSE(map.passed_x, reference.passed_x, 2 * map.reach / (n - 1));
    CHECK_CLOSE(map.piecewise_quality.covered, expectedCoverage, 0.01);
    CHECK_CLOSE(map.ellipse_quality.covered, expectedCoverage, 0.01);
    CHECK_EQ(map.samples.size(), static_cast<std::size_t>(n * n));
    bool boundaryBlocked = true;
    for (int k = 0; k < n; k++)
        boundaryBlocked = boundaryBlocked && !map.sample(0, k).passed &&
                          !map.sample(n - 1, k).passed && !map.sample(k, 0).passed &&
                          !map.sample(k, n - 1).passed;
    CHECK(boundaryBlocked);
    CHECK_EQ(fld.vlx, 0.0);
    CHECK_EQ(fld.vux, 0.0);
    CHECK_EQ(fld.vly, 0.0);
    CHECK_EQ(fld.vuy, 0.0);
}

TEST(pupilMap_plotsTheMap) {
    auto opm = model();
    auto result = PupilMapAnalysis::eval(opm.get(), SAMPLES, std::vector<int>{10});
    std::string svg = PupilMapPlot(result.maps.at(0)).plot(320);
    CHECK(svg.rfind("<?xml", 0) == 0);
    CHECK(contains(svg, "<svg"));
    // The passed region, the two candidate outlines and the labels all made it in.
    CHECK(contains(svg, "polygon"));
    CHECK(contains(svg, "field 1.00"));
    CHECK(contains(svg, "piecewise"));
}

} // namespace
