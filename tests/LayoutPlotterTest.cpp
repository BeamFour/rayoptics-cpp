// End-to-end check of rayoptics.layout and the plotter package against the JVM.
//
// Builds the canon-rf70-200mm-f2.8LZ zoom from its Examples prescription and
// reproduces exactly what LensTool2.doLayoutDiagrams draws -- the ray-fan
// layout, the elements-only layout and the reference-ray layout -- then the
// four plot types: spot diagram, monochromatic geometric MTF, MTF by field
// (SVG and CSV) and the transverse/OPD ray aberration fans.
//
// Compared line by line with output dumped from JDK 25
// (scratchpad/DumpLayout.java).
#include "LayoutPlotterExpected.h"
#include "TestHarness.h"

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/plotter/Plotter.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/layout/Layout2D.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/spec/Prescription.h"

#include <memory>
#include <string>
#include <vector>

namespace {

using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;
using redukti::spec::RayOpticsModelBuilder;
using redukti::spec::VigType;
using namespace redukti::rayoptics;
using namespace redukti::plotter;

const char *const SPEC =
    REDUKTI_EXAMPLES_DIR "canon-rf70-200mm-f2.8LZ/US20250155694_Example01P.txt";

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        if (nl == std::string::npos) {
            out.push_back(text.substr(start));
            break;
        }
        out.push_back(text.substr(start, nl - start));
        start = nl + 1;
    }
    while (!out.empty() && out.back().empty())
        out.pop_back();
    return out;
}

Prescription buildPrescription() {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(SPEC);
    return Prescription::build_prescription(specs, true, false, false);
}

std::unique_ptr<optical::OpticalModel> layoutSystem(const Prescription &p, int config) {
    return RayOpticsModelBuilder(p).build_optical_model(
        true, std::vector<double>{0.0, 1.0}, false, VigType::SetPupil, true, config);
}

std::unique_ptr<optical::OpticalModel> analysisSystem(const Prescription &p) {
    return RayOpticsModelBuilder(p).build_optical_model(
        true, std::vector<double>{0.0, 0.3, 0.5, 0.707, 0.85, 1.0}, false,
        VigType::SetVig, false, 0);
}

} // namespace

#define CHECK_BLOCK(actual, expected)                                                    \
    do {                                                                                 \
        auto _al = lines(actual);                                                        \
        const std::size_t _en = sizeof(expected) / sizeof(expected[0]);                  \
        CHECK_EQ(_al.size(), _en);                                                       \
        for (std::size_t _i = 0; _i < _al.size() && _i < _en; _i++)                       \
            CHECK_STR_EQ(_al[_i], std::string(expected[_i]));                            \
    } while (0)

TEST(layout_svgs_match_jvm) {
    Prescription p = buildPrescription();
    auto opm = layoutSystem(p, 0);
    layout::Layout2D lay;

    layout::LayoutOptions fanOpts;
    fanOpts.drawReferenceRays_(false).fanRayCount_(9).clipRays_(true).useTraceFan_(true);
    CHECK_BLOCK(lay.renderSvg(opm.get(), 1000, 500, &fanOpts), EXPECTED_LAYOUT_FAN);

    layout::LayoutOptions elementsOpts;
    elementsOpts.drawReferenceRays_(false);
    CHECK_BLOCK(lay.renderSvg(opm.get(), 1000, 500, &elementsOpts), EXPECTED_LAYOUTONLY);

    layout::LayoutOptions defaultOpts;
    CHECK_BLOCK(lay.renderSvg(opm.get(), 1000, 500, &defaultOpts), EXPECTED_LAYOUT);
}

TEST(plotter_spot_and_mtf_match_jvm) {
    Prescription p = buildPrescription();
    auto model = analysisSystem(p);
    auto options = analysis::SpotOptions();
    options.use_hexapolar().num_rings(6);
    auto spot = analysis::SpotAnalysis::eval(model.get(), options);

    CHECK_BLOCK(SpotDiagram(spot.spot_results[0]).plot(std::nullopt), EXPECTED_SPOT_DIAGRAM);

    const auto &spotFld = spot.spot_results[1];
    auto cfg = spotFld.mtfHistogramConfig();
    analysis::MonochromaticGeometricMTF mono(spotFld.intercepts[0], cfg);
    CHECK_BLOCK(GeoMTFPlot(*spotFld.fld, mono).plot(), EXPECTED_GEO_MTF);

    std::vector<int> freqs{10, 20, 40};
    GeoMTFByFieldPlot byField(spot.computeMTFs(freqs), spot.fields());
    CHECK_BLOCK(byField.plot(), EXPECTED_MTF_BY_FIELD);
    CHECK_BLOCK(byField.toString(), EXPECTED_MTF_BY_FIELD_CSV);
    CHECK_BLOCK(GeoMTFByFieldPlot::freq_legend(freqs), EXPECTED_FREQ_LEGEND);
}

TEST(plotter_aberration_fans_match_jvm) {
    Prescription p = buildPrescription();
    auto model = analysisSystem(p);
    raytr::TraceOptions traceOptions;

    auto tra = analysis::TransverseRayAberrationAnalysis::eval(model.get(), 11, false,
                                                               traceOptions);
    CHECK_BLOCK(RayAberrationPlot(tra).plot(tra.results[4], 0), EXPECTED_RAY_ABERRATION);

    auto opd = analysis::WavefrontAberrationAnalysis::eval(model.get(), 11, false,
                                                           traceOptions);
    CHECK_BLOCK(RayAberrationPlot(opd).plot(opd.results[7], 0), EXPECTED_OPD_ABERRATION);
}
