// C++ port of org.redukti.rayoptics.analysis.MtfTest.
//
// Compares polychromatic geometric MTF from the two spot sampling patterns on
// the Nikkor Z 58mm f/0.95 S, against golden strings rounded to three decimals,
// and then the two patterns against each other to 0.01.
//
// This is the only test covering the Gaussian-quadrature MTF path.
// analysis_mtf_matches_jvm exercises computeMTFs too, but only from hexapolar
// sampling, and the two patterns reach the MTF by different routes -- the
// quadrature grid carries per-ray weights the hexapolar one does not.
#include "TestHarness.h"

#include "NikkorZ58Model.h"

#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace redukti::rayoptics;

/** Java: MtfTest.df, an M.decimal_format pinned to the root locale. */
redukti::DecimalFormat df() { return redukti::mathlib::M::decimal_format(); }


/** Java: MtfTest.generateMTFs. */
std::vector<analysis::MTFResultByFreq> generateMTFs(
    optical::OpticalModel *opm, const std::vector<int> &freqs,
    const std::map<double, double> &wv_wts, const analysis::SpotOptions &spotOptions) {
    auto spotAnalysis = analysis::SpotAnalysis::eval(opm, spotOptions);
    std::vector<analysis::PolyMTF> mtfs;
    for (std::size_t i = 0; i < spotAnalysis.spot_results.size(); i++) {
        const auto &spotFld = spotAnalysis.spot_results[i];
        auto cfg = spotFld.mtfHistogramConfig();
        // The Java starts from a null PolyMTF and builds it from the first
        // monochromatic MTF, so that it inherits that MTF's fft size and pixel
        // size. An optional stands in for the null.
        std::optional<analysis::PolyMTF> polyMtfForField;
        for (const auto &intercepts : spotFld.intercepts) {
            analysis::MonochromaticGeometricMTF mtf(intercepts, cfg);
            if (!polyMtfForField.has_value())
                polyMtfForField.emplace(mtf.mtf.fft_size, mtf.h2d.pixel_size);
            auto it = wv_wts.find(intercepts.wvl);
            double wt = it == wv_wts.end() ? 0.0 : it->second;
            if (wt != 0.0)
                polyMtfForField->add(mtf.mtf, wt);
        }
        if (polyMtfForField.has_value()) {
            polyMtfForField->compute();
            mtfs.push_back(std::move(*polyMtfForField));
        }
    }
    std::vector<analysis::MTFResultByFreq> mtfResults;
    for (auto freq : freqs)
        mtfResults.push_back(analysis::MTFResultByFreq(mtfs, freq));
    return mtfResults;
}

std::map<double, double> get_wvl_wts(const std::vector<double> &wvls,
                                     const std::vector<double> &wts) {
    std::map<double, double> map;
    for (std::size_t i = 0; i < wvls.size(); i++)
        map[wvls[i]] = wts[i];
    return map;
}

std::string toString(const std::vector<analysis::MTFResultByFreq> &mtfs_by_freq,
                     const std::vector<double> &fields) {
    auto fmt = df();
    std::string sb;
    sb += ",";
    for (std::size_t i = 0; i < mtfs_by_freq[0].tan_mtf_by_field.size(); i++) {
        if (i > 0)
            sb += ",";
        sb += fmt.format(fields[i]);
    }
    sb += "\n";
    for (std::size_t i = 0; i < mtfs_by_freq.size(); i++) {
        const auto &mtf = mtfs_by_freq[i];
        for (int xy = 0; xy < 2; xy++) {
            const auto &mtf_data =
                (xy == 0) ? mtf.sag_mtf_by_field : mtf.tan_mtf_by_field;
            sb += redukti::intToString(mtf.freq);
            sb += " ";
            sb += (xy == 0 ? "sag" : "tan");
            sb += ",";
            for (std::size_t j = 0; j < mtf_data.size(); j++) {
                if (j > 0)
                    sb += ",";
                sb += fmt.format(mtf_data[j]);
            }
            sb += "\n";
        }
    }
    return sb;
}

void assertMtfClose(const std::vector<analysis::MTFResultByFreq> &expected,
                    const std::vector<analysis::MTFResultByFreq> &actual,
                    double tolerance) {
    CHECK_EQ(static_cast<int>(expected.size()), static_cast<int>(actual.size()));
    for (std::size_t f = 0; f < expected.size() && f < actual.size(); f++) {
        CHECK_EQ(expected[f].freq, actual[f].freq);
        for (std::size_t i = 0; i < expected[f].sag_mtf_by_field.size(); i++)
            CHECK_CLOSE(actual[f].sag_mtf_by_field[i], expected[f].sag_mtf_by_field[i],
                        tolerance);
        for (std::size_t i = 0; i < expected[f].tan_mtf_by_field.size(); i++)
            CHECK_CLOSE(actual[f].tan_mtf_by_field[i], expected[f].tan_mtf_by_field[i],
                        tolerance);
    }
}

} // namespace

TEST(mtf_sampling_patterns_match_jvm) {
    auto opm = redukti::test::buildNikkorZ58TestModel();
    auto wavelengthWeights =
        get_wvl_wts({587.5618, 486.1327, 656.2725}, {1.0, 1.0, 1.0});

    analysis::SpotOptions hexOptions;
    hexOptions.use_hexapolar().num_rays(64);
    auto hexapolarMtf =
        generateMTFs(opm.get(), {10, 30, 50}, wavelengthWeights, hexOptions);

    analysis::SpotOptions gqOptions;
    gqOptions.use_gaussian_quadrature().num_rings(14).num_spokes(20);
    auto quadratureMtf =
        generateMTFs(opm.get(), {10, 30, 50}, wavelengthWeights, gqOptions);

    const std::string expectedHexapolar =
        ",0,0.7,1\n"
        "10 sag,0.967,0.955,0.763\n"
        "10 tan,0.967,0.949,0.939\n"
        "30 sag,0.754,0.709,0.506\n"
        "30 tan,0.754,0.648,0.589\n"
        "50 sag,0.519,0.47,0.414\n"
        "50 tan,0.519,0.385,0.319\n";
    const std::string expectedQuadrature =
        ",0,0.7,1\n"
        "10 sag,0.967,0.954,0.766\n"
        "10 tan,0.967,0.95,0.94\n"
        "30 sag,0.753,0.713,0.507\n"
        "30 tan,0.753,0.654,0.593\n"
        "50 sag,0.519,0.478,0.417\n"
        "50 tan,0.519,0.388,0.316\n";

    std::vector<double> fields{0.0, 0.7, 1.0};
    CHECK_STR_EQ(toString(hexapolarMtf, fields), expectedHexapolar);
    CHECK_STR_EQ(toString(quadratureMtf, fields), expectedQuadrature);
    assertMtfClose(hexapolarMtf, quadratureMtf, 0.01);
}
