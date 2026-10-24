// C++ port of org.redukti.analysis.ContrastReportTest and
// MonochromaticMTFReportTest.
//
// The Java tests check the table's shape and recompute a few cells. Here the
// whole table is also asserted against the Java output for the same model and
// options (DumpReports.java, a throwaway dumper run against the Java classes).
#include "TestHarness.h"

#include "NikkorZ58Model.h"

#include "redukti/analysis/ContrastAnalysis.h"
#include "redukti/analysis/ContrastReport.h"
#include "redukti/analysis/MonochromaticMTFReport.h"
#include "redukti/analysis/MTF.h"
#include "redukti/analysis/SpotAnalysis.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

using redukti::analysis::ContrastAnalysis;
using redukti::analysis::ContrastOptions;
using redukti::analysis::ContrastReport;
using redukti::analysis::MonochromaticGeometricMTF;
using redukti::analysis::MonochromaticMTFReport;
using redukti::analysis::MTFResultByFreq;
namespace Orientation = redukti::rayoptics::util::Orientation;
namespace M = redukti::mathlib::M;

/** Java's report.lines().filter(line -> line.startsWith(" ")). */
std::vector<std::string> tableLines(const std::string &report) {
    std::vector<std::string> lines;
    std::istringstream in(report);
    std::string line;
    while (std::getline(in, line))
        if (!line.empty() && line[0] == ' ') lines.push_back(line);
    return lines;
}

/** Java's line.trim().split("\\s+"). */
std::vector<std::string> split(const std::string &line) {
    std::vector<std::string> words;
    std::istringstream in(line);
    std::string word;
    while (in >> word) words.push_back(word);
    return words;
}

const char *const EXPECTED_CONTRAST_REPORT =
    "RMS OPD difference (waves), sampling 2 rings x 6 spokes, calibrate no, exit-pupil-aiming no, centering no\n"
    "\n"
    "     freq      wvl      dir        0      0.7        1\n"
    "       10 587.5618      tan   0.0048   0.0321   0.0273\n"
    "       10 587.5618      sag   0.0048   0.0544    0.119\n"
    "       10 486.1327      tan    0.031   0.0693   0.0426\n"
    "       10 486.1327      sag    0.031   0.0522   0.1153\n"
    "       10 656.2725      tan   0.0625    0.045   0.0617\n"
    "       10 656.2725      sag   0.0625   0.0767   0.1465\n"
    "       40 587.5618      tan   0.0197   0.1311   0.1092\n"
    "       40 587.5618      sag   0.0197   0.2174   0.4685\n"
    "       40 486.1327      tan   0.1244   0.2789   0.1744\n"
    "       40 486.1327      sag   0.1244   0.2066    0.451\n"
    "       40 656.2725      tan   0.2497   0.1802   0.2458\n"
    "       40 656.2725      sag   0.2497    0.313   0.5806\n";

const char *const EXPECTED_MONOCHROMATIC_MTF_REPORT =
    "Monochromatic geometric MTF\n"
    "\n"
    "     freq      wvl      dir        0      0.7        1\n"
    "       10 587.5618      tan    0.998    0.982    0.972\n"
    "       10 587.5618      sag    0.998    0.964    0.789\n"
    "       10 486.1327      tan    0.977    0.934    0.957\n"
    "       10 486.1327      sag    0.977    0.962    0.775\n"
    "       10 656.2725      tan    0.926    0.959    0.897\n"
    "       10 656.2725      sag    0.926    0.933    0.709\n"
    "       30 587.5618      tan    0.982    0.857    0.772\n"
    "       30 587.5618      sag    0.982    0.725    0.551\n"
    "       30 486.1327      tan    0.813    0.588    0.708\n"
    "       30 486.1327      sag    0.814    0.712    0.479\n"
    "       30 656.2725      tan    0.462    0.675    0.366\n"
    "       30 656.2725      sag    0.464    0.657    0.488\n";

} // namespace

TEST(contrast_report_tabulates_one_row_per_frequency_wavelength_and_orientation) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto options = ContrastOptions(0.0).num_rings(2).num_spokes(6);

    std::string report = ContrastReport::report(model.get(), {10, 40}, options);
    CHECK_STR_EQ(report, std::string(EXPECTED_CONTRAST_REPORT));

    auto lines = tableLines(report);
    auto header = split(lines[0]);
    int fields = static_cast<int>(model->optical_spec->fov->fields.size());
    int wavelengths = static_cast<int>(model->optical_spec->wvls->wavelengths.size());
    CHECK_EQ(static_cast<int>(header.size()), 3 + fields);
    CHECK_STR_EQ(header[0], std::string("freq"));
    CHECK_EQ(static_cast<int>(lines.size()), 1 + 2 * 2 * wavelengths);

    // First row: lowest frequency, first wavelength, tangential.
    auto row = split(lines[1]);
    CHECK_STR_EQ(row[0], std::string("10"));
    CHECK_STR_EQ(row[2], std::string("tan"));
    CHECK_STR_EQ(split(lines[2])[2], std::string("sag"));

    auto result = ContrastAnalysis::eval(model.get(), options.with_frequency(10));
    auto format = M::decimal_format(4);
    for (int fi = 0; fi < fields; fi++) {
        const auto &block = result.fields[static_cast<std::size_t>(fi)].wavelengths[0];
        CHECK_STR_EQ(format.format(ContrastReport::rms(block, Orientation::TANGENTIAL)),
                     row[static_cast<std::size_t>(3 + fi)]);
    }
}

TEST(contrast_report_rms_is_normalized_by_the_weight_of_valid_samples) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto result = ContrastAnalysis::eval(model.get(),
                                         ContrastOptions(40.0).num_rings(2).num_spokes(6));
    const auto &block = result.fields[1].wavelengths[0];

    double sum = 0.0;
    for (const auto &sample : block.samples)
        sum += sample.weight * sample.tangentialOpdDifference * sample.tangentialOpdDifference;

    double rms = ContrastReport::rms(block, Orientation::TANGENTIAL);
    CHECK(rms > 0.0);
    CHECK_CLOSE(rms, std::sqrt(sum), 1.0e-12);
}

TEST(monochromatic_mtf_report_tabulates_one_row_per_frequency_wavelength_and_orientation) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto spotAnalysis = redukti::analysis::SpotAnalysis::eval(
        model.get(), redukti::analysis::SpotOptions().use_hexapolar().num_rings(21));
    std::vector<std::vector<MonochromaticGeometricMTF>> mtfs;
    std::vector<double> fields(spotAnalysis.spot_results.size());
    for (std::size_t i = 0; i < fields.size(); i++) {
        const auto &spotFld = spotAnalysis.spot_results[i];
        fields[i] = spotFld.fld->y;
        auto cfg = spotFld.mtfHistogramConfig();
        auto &field = mtfs.emplace_back();
        for (const auto &intercepts : spotFld.intercepts)
            field.emplace_back(intercepts, cfg);
    }

    std::string report = MonochromaticMTFReport::report(fields, mtfs, {10, 30});
    // Three decimals hide the sin/cos ulps of the hexapolar pattern.
    CHECK_STR_EQ(report, std::string(EXPECTED_MONOCHROMATIC_MTF_REPORT));

    auto lines = tableLines(report);
    int wavelengths = static_cast<int>(mtfs[0].size());
    CHECK_EQ(static_cast<int>(split(lines[0]).size()), 3 + static_cast<int>(fields.size()));
    CHECK_EQ(static_cast<int>(lines.size()), 1 + 2 * 2 * wavelengths);

    // Second row: lowest frequency, first wavelength, sagittal.
    auto row = split(lines[2]);
    CHECK_STR_EQ(row[0], std::string("10"));
    CHECK_STR_EQ(row[2], std::string("sag"));
    auto format = M::decimal_format(3);
    for (std::size_t fi = 0; fi < fields.size(); fi++) {
        const auto &mtf = mtfs[fi][0].mtf;
        CHECK_STR_EQ(format.format(MTFResultByFreq::interpolate(mtf.freq, mtf.mag_x, 10)),
                     row[3 + fi]);
    }
}

TEST(monochromatic_mtf_report_rejects_a_field_label_count_mismatch) {
    std::vector<std::vector<MonochromaticGeometricMTF>> mtfs;
    CHECK_THROWS(MonochromaticMTFReport::report({0.0}, mtfs, {10}),
                 redukti::IllegalArgumentException);
}
