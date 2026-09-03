// End-to-end check of the analysis layer against the JVM.
//
// Builds the same Leica Summicron R 50mm f/2 prescription as RaytrTest, runs
// spot, MTF, ray-fan and contrast analyses, and compares the rendered results
// character-for-character with output dumped from the Java on JDK 25
// (scratchpad DumpAnalysis.java -> tests/AnalysisExpected.h).
//
// The second field vignettes: it traces 177 of 217 rays, so this also pins the
// path where the image filter answers null for a ray that failed.
#include "AnalysisExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"
#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <cmath>
#include <cctype>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace redukti::rayoptics;

/** Compare block by block so a failure points at the one row that moved. */
#define CHECK_BLOCK(actual, expected)                                         \
    do {                                                                      \
        auto _al = lines(actual);                                             \
        auto _el = lines(expected);                                           \
        CHECK_EQ(_al.size(), _el.size());                                     \
        for (std::size_t _i = 0; _i < _al.size() && _i < _el.size(); _i++)    \
            CHECK_STR_EQ(_al[_i], _el[_i]);                                   \
    } while (0)

/** Structure exact, trig-derived numbers within tolerance. See splitNums. */
#define CHECK_BLOCK_APPROX(actual, expected, atol, rtol)                                  \
    do {                                                                                  \
        auto _al = lines(actual);                                                         \
        auto _el = lines(expected);                                                       \
        CHECK_EQ(_al.size(), _el.size());                                                 \
        for (std::size_t _i = 0; _i < _al.size() && _i < _el.size(); _i++) {              \
            auto _sa = splitNums(_al[_i]);                                                 \
            auto _se = splitNums(_el[_i]);                                                 \
            CHECK_STR_EQ(_sa.skeleton, _se.skeleton);                                       \
            CHECK_EQ(_sa.nums.size(), _se.nums.size());                                     \
            for (std::size_t _k = 0; _k < _sa.nums.size() && _k < _se.nums.size(); _k++)    \
                CHECK_CLOSE(_sa.nums[_k], _se.nums[_k],                                     \
                            (atol) + (rtol) * std::abs(_se.nums[_k]));                     \
        }                                                                                 \
    } while (0)

using redukti::doubleToString;
using util::Pair;

namespace {

/** The Java writes `String.valueOf(double)`; this is the same text. */
std::string d(double v) { return doubleToString(v); }

std::string b(bool v) { return v ? "true" : "false"; }

constexpr char chr_nl = '\n';

std::unique_ptr<optical::OpticalModel> buildSummicron() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                               specs::ValueKey::Fnum),
        2.0);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp,
        Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                               specs::ValueKey::Angle),
        22.5, std::vector<double>{0., 1.}, true, true);
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0)}, 0);
    opm->system_spec->title = "Leica Summicron R 50mm f/2)";
    opm->system_spec->dimensions = "mm";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;

    auto add = [&](double curv, double thi, double index, double vd, double max_ap) {
        seq::SurfaceData sd(curv, thi);
        if (index > 0.0)
            sd.rindex(index, vd);
        sd.max_aperture(max_ap);
        sm->add_surface(sd);
    };
    add(42.71, 3.99, 1.73430, 28.19, 14.47);
    add(195.38, 0.2, 0.0, 0.0, 13.53);
    add(20.5, 7.18, 1.67133, 41.64, 12.01);
    add(0.0, 1.29, 1.79190, 25.55, 10.745);
    add(14.94, 5.35, 0.0, 0.0, 9.195);
    add(0.0, 7.61, 0.0, 0.0, 9.0295);
    sm->set_stop();
    add(-14.94, 1.0, 1.65222, 33.60, 8.75);
    add(0.0, 5.22, 1.79227, 47.15, 9.635);
    add(-20.5, 0.2, 0.0, 0.0, 10.19);
    add(0.0, 3.69, 1.79227, 47.15, 11.48);
    add(-42.71, 37.32, 0.0, 0.0, 11.985);
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_vig(opm.get());
    opm->update_model();
    return opm;
}

analysis::SpotAnalysisResult spotResult(optical::OpticalModel *opm) {
    auto options = analysis::SpotOptions();
    options.use_hexapolar().num_rings(8);
    return analysis::SpotAnalysis::eval(opm, options);
}

/**
 * Compare a rendered block where some numbers are trig-derived.
 *
 * Every non-numeric character must match exactly -- labels, counts, validity
 * flags, array indices -- so structure is still pinned. Numbers are compared
 * with a tolerance, because the hexapolar pupil generator calls sin/cos and
 * JDK 25 and MSVC disagree by 1-2 ulps on 18 of the 217 ring points (measured
 * with scratchpad/hexprobe.cpp). Those differences propagate through the trace
 * into the intercepts and on into the MTF. The fan and contrast analyses step
 * the pupil linearly, so they have no such divergence and stay exact.
 */
struct NumSplit {
    std::string skeleton;
    std::vector<double> nums;
};

NumSplit splitNums(const std::string &s) {
    NumSplit out;
    std::size_t i = 0;
    while (i < s.size()) {
        std::size_t j = i;
        if (s[j] == '-' || s[j] == '+')
            j++;
        std::size_t digits = j;
        while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j])))
            j++;
        bool any = j > digits;
        if (j < s.size() && s[j] == '.') {
            j++;
            while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) {
                j++;
                any = true;
            }
        }
        if (any && j < s.size() && (s[j] == 'E' || s[j] == 'e')) {
            std::size_t k = j + 1;
            if (k < s.size() && (s[k] == '-' || s[k] == '+'))
                k++;
            std::size_t ed = k;
            while (k < s.size() && std::isdigit(static_cast<unsigned char>(s[k])))
                k++;
            if (k > ed)
                j = k;
        }
        if (any) {
            out.nums.push_back(std::strtod(s.substr(i, j - i).c_str(), nullptr));
            out.skeleton += '#';
            i = j;
        } else {
            out.skeleton += s[i];
            i++;
        }
    }
    return out;
}

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find(chr_nl, start);
        if (nl == std::string::npos) {
            out.push_back(text.substr(start));
            break;
        }
        out.push_back(text.substr(start, nl - start));
        start = nl + 1;
    }
    return out;
}

} // namespace

TEST(analysis_spot_matches_jvm) {
    auto opm = buildSummicron();
    auto spot = spotResult(opm.get());

    std::string sb = spot.toString();
    auto flds = spot.fields();
    for (std::size_t i = 0; i < flds.size(); i++)
        sb += "field[" + std::to_string(i) + "]=" + d(flds[i]) + "\n";
    for (std::size_t i = 0; i < spot.spot_results.size(); i++) {
        const auto &r = spot.spot_results[i];
        sb += "spot[" + std::to_string(i) + "] max=" + d(r.max_radius) +
              " mean=" + d(r.mean_radius) +
              " ngrids=" + std::to_string(r.trace_results.size()) + "\n";
        for (std::size_t k = 0; k < r.intercepts.size(); k++) {
            const auto &ic = r.intercepts[k];
            sb += "  wvl=" + d(ic.wvl) + " n=" + std::to_string(ic.x.size()) + "\n";
            for (std::size_t j = 0; j < ic.x.size(); j += 17)
                sb += "   [" + std::to_string(j) + "] " + d(ic.x[j]) + " " + d(ic.y[j]) +
                      " w=" + d(ic.weights[j]) + " v=" + b(ic.valid[j] != 0) + "\n";
            auto c = ic.compute_centroid();
            sb += "  centroid=" + d(c.x) + "," + d(c.y) + "\n";
        }
        auto cfg = r.mtfHistogramConfig();
        sb += "  cfg bins=" + std::to_string(cfg.num_bins) + " px=" + d(cfg.pixel_size) +
              "\n";
    }
    CHECK_BLOCK_APPROX(sb, EXPECTED_SPOT, 1e-14, 1e-12);
}

TEST(analysis_mtf_matches_jvm) {
    auto opm = buildSummicron();
    auto spot = spotResult(opm.get());

    std::string sb;
    auto mtfs = spot.computeMTFs(std::vector<int>{5, 10, 20, 40});
    for (const auto &m : mtfs) {
        sb += "freq=" + std::to_string(m.freq) + "\n";
        for (std::size_t i = 0; i < m.sag_mtf_by_field.size(); i++)
            sb += "  sag[" + std::to_string(i) + "]=" + d(m.sag_mtf_by_field[i]) +
                  " tan[" + std::to_string(i) + "]=" + d(m.tan_mtf_by_field[i]) + "\n";
    }
    CHECK_BLOCK_APPROX(sb, EXPECTED_MTF, 1e-14, 1e-12);
}

TEST(analysis_transverse_fans_match_jvm) {
    auto opm = buildSummicron();
    raytr::TraceOptions trace_options;
    auto tra = analysis::TransverseRayAberrationAnalysis::eval(opm.get(), 11, false,
                                                               trace_options);
    CHECK_BLOCK(tra.list_ray_fans(), EXPECTED_TRA);
}

TEST(analysis_opd_fans_match_jvm) {
    auto opm = buildSummicron();
    raytr::TraceOptions trace_options;
    auto opd = analysis::WavefrontAberrationAnalysis::eval(opm.get(), 11, false,
                                                           trace_options);
    CHECK_BLOCK(opd.list_ray_fans(), EXPECTED_OPD);
}

TEST(analysis_contrast_matches_jvm) {
    auto opm = buildSummicron();
    auto options = analysis::ContrastOptions(30.0);
    options.num_rings(2).num_spokes(4).center_residuals(true);
    auto ca = analysis::ContrastAnalysis::eval(opm.get(), options);

    std::string sb = "spatialFrequency=" + d(ca.spatialFrequency) + "\n";
    for (std::size_t fi = 0; fi < ca.fields.size(); fi++) {
        const auto &f = ca.fields[fi];
        sb += "field " + std::to_string(fi) + " " + f.field->toString() + "\n";
        for (const auto &w : f.wavelengths) {
            sb += "  wvl=" + d(w.wavelength) + " shift=" + d(w.normalizedPupilShift) +
                  " sagOff=" + d(w.sagittalOffset) + " tanOff=" + d(w.tangentialOffset) +
                  " n=" + std::to_string(w.samples.size()) + "\n";
            for (std::size_t si = 0; si < w.samples.size(); si++) {
                const auto &s = w.samples[si];
                sb += "   [" + std::to_string(si) + "] pupil=" + d(s.pupil.x) + "," +
                      d(s.pupil.y) + " sag=" + d(s.sagittalDifference) +
                      " tan=" + d(s.tangentialDifference) + " w=" + d(s.weight) +
                      " valid=" + b(s.valid) + " fail=" +
                      (!s.failure.has_value()
                           ? std::string("null")
                           : s.failure->ray + "/" + s.failure->exceptionType + "/" +
                                 std::to_string(s.failure->surface)) +
                      " sres=" + d(w.sagittalResidual(static_cast<int>(si))) +
                      " tres=" + d(w.tangentialResidual(static_cast<int>(si))) + "\n";
            }
        }
    }
    CHECK_BLOCK(sb, EXPECTED_CONTRAST);
}
