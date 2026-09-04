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

#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/analysis/MTF.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/seq/SurfaceData.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/specs/WvlSpec.h"

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace redukti::rayoptics;
using redukti::rayoptics::util::Pair;
namespace profiles = redukti::rayoptics::elem::profiles;

/** Java: MtfTest.df, an M.decimal_format pinned to the root locale. */
redukti::DecimalFormat df() { return redukti::mathlib::M::decimal_format(); }

std::unique_ptr<optical::OpticalModel> buildTestModel() {
    auto opm = std::make_unique<optical::OpticalModel>();
    auto sm = opm->seq_model.get();
    auto osp = opm->optical_spec.get();
    osp->pupil = std::make_unique<specs::PupilSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Image,
                                                    specs::ValueKey::Fnum),
        0.98);
    osp->fov = std::make_unique<specs::FieldSpec>(
        osp, Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                    specs::ValueKey::Angle),
        std::vector<double>{0., 0.7, 1.0});
    osp->fov->value = 19.98;
    osp->fov->is_relative = true;
    osp->fov->is_wide_angle = true;
    osp->wvls = std::make_unique<specs::WvlSpec>(
        std::vector<specs::WvlWt>{specs::WvlWt(587.5618, 1.0),
                                  specs::WvlWt(486.1327, 1.0),
                                  specs::WvlWt(656.2725, 1.0)},
        0);
    opm->system_spec->title = "WO2019-229849 Example 1 (Nikkor Z 58mm f/0.95 S)";
    opm->system_spec->dimensions = "MM";
    opm->radius_mode = true;
    sm->gaps[0]->thi = 1e10;
    {
        seq::SurfaceData sd(108.488, 7.65);
        sd.rindex(1.90265, 35.77, "J-LASFH9A", "Hikari")->max_aperture(33.4);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(108.488)->setCc(0)->setCoefs(std::vector<double>{0.0, -3.82177e-07, -6.06486e-11, -3.80172e-15, -1.32266e-18, 0, 0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(-848.55, 2.8);
        sd.rindex(1.55298, 55.07, "J-KZFH4", "Hikari")->max_aperture(32.91);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(50.252, 18.12);
        sd.max_aperture(28.97);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-60.72, 2.8);
        sd.rindex(1.61266, 44.46, "J-KZFH1", "Hikari")->max_aperture(29.14);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(2497.5, 9.15);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(32.66);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-77.239, 0.4);
        sd.max_aperture(32.66);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(113.763, 10.95);
        sd.rindex(1.8485, 43.79, "J-LASFH22", "Hikari")->max_aperture(35.45);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-178.06, 0.4);
        sd.max_aperture(35.45);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(70.659, 9.74);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(32.5);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-1968.5, 0.2);
        sd.max_aperture(32.5);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(289.687, 8);
        sd.rindex(1.59319, 67.9, "J-PSKH1", "Hikari")->max_aperture(30.53);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-97.087, 2.8);
        sd.rindex(1.738, 32.33, "J-KZFH9", "Hikari")->max_aperture(29.71);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(47.074, 8.7);
        sd.max_aperture(25.12);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(0, 5.29);
        // 23.85 here, not the 23.959 the Nikkor58mm integration test uses for
        // the same lens: MtfTest builds its own model and its golden MTF values
        // depend on this stop semi-diameter.
        sd.max_aperture(23.85);
        sm->add_surface(sd);
    }
    sm->set_stop();
    {
        seq::SurfaceData sd(-95.23, 2.2);
        sd.rindex(1.61266, 44.46, "J-KZFH1", "Hikari")->max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(41.204, 11.55);
        sd.rindex(1.49782, 82.57, "J-FKH1", "Hikari")->max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-273.092, 0.2);
        sd.max_aperture(24.96);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(76.173, 9.5);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(25.56);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-101.575, 0.2);
        sd.max_aperture(25.56);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(176.128, 7.45);
        sd.rindex(1.95375, 32.33, "J-LASFH21", "Hikari")->max_aperture(23.4);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(176.128)->setCc(0)->setCoefs(std::vector<double>{0.0, -1.15028e-06, -4.51771e-10, 2.7267e-13, -7.66812e-17, 0, 0});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(-67.221, 1.8);
        sd.rindex(1.738, 32.33, "J-KZFH9", "Hikari")->max_aperture(22.68);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(55.51, 2.68);
        sd.max_aperture(19.92);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(71.413, 6.35);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-115.025, 1.81);
        sd.rindex(1.69895, 30.13, "J-SF15", "Hikari")->max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(46.943, 0.8);
        sd.max_aperture(19.73);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(55.281, 9.11);
        sd.rindex(1.883, 40.69, "J-LASF08A", "Hikari")->max_aperture(19.47);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(-144.041, 3);
        sd.rindex(1.76554, 46.76, "J-LASFH2", "Hikari")->max_aperture(19.14);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(52.858, 14.5);
        sd.max_aperture(19.14);
        sm->add_surface(sd);
    }
    {
        auto prof = std::make_shared<profiles::EvenPolynomial>();
        prof->r(52.858)->setCc(0)->setCoefs(std::vector<double>{0.0, 3.18645e-06, -1.14718e-08, 7.74567e-11, -2.24225e-13, 3.3479e-16, -1.7047e-19});
        sm->ifcs[static_cast<std::size_t>(*sm->cur_surface)]->profile = prof;
    }
    {
        seq::SurfaceData sd(0, 1.6);
        sd.rindex(1.5168, 64.14, "J-BK7A", "Hikari")->max_aperture(22.15);
        sm->add_surface(sd);
    }
    {
        seq::SurfaceData sd(0, 1);
        sd.max_aperture(22.15);
        sm->add_surface(sd);
    }
    sm->do_apertures = false;
    opm->update_model();
    raytr::VigCalc::set_pupil(opm.get());
    opm->update_model();
    return opm;
}

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
    auto opm = buildTestModel();
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
