// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.PolyDiffractionMTFTest.
//
// The Java tests are ported with their tolerances. In addition, the polychromatic MTF
// and MTFResultByFreq::from_diffraction of the Summicron prescription are asserted
// against the Java classes (DumpPolyDiffractionMTF.java, a throwaway dumper): the
// image points exactly, the MTF to 1e-12, since it goes through the PSF's cos and sin
// and the FFT wavetables (see DiffractionMTFTest).
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/JavaSemantics.h"
#include "redukti/analysis/DiffractionMTF.h"
#include "redukti/analysis/MTF.h"
#include "redukti/analysis/PolyDiffractionMTF.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/spec/Prescription.h"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace {

using redukti::analysis::DiffractionMTF;
using redukti::analysis::MTFResultByFreq;
using redukti::analysis::PolyDiffractionMTF;
using redukti::mathlib::M::PI;
using redukti::rayoptics::optical::OpticalModel;
using redukti::rayoptics::specs::Field;

const int N = 64;
const char *const SUMMICRON = REDUKTI_EXAMPLES_DIR "leica-summicron-50mm-f2/US004123144_Example08P.txt";

/** DiffractionMTFTest.polychromatic_model in the Java. */
std::unique_ptr<OpticalModel> polychromatic_model() {
    redukti::importers::OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(SUMMICRON);
    auto prescription = redukti::spec::Prescription::build_prescription(specs, true, false, false);
    return redukti::spec::RayOpticsModelBuilder(prescription)
        .build_optical_model(true, std::vector<double>{0.0, 1.0}, false,
                             redukti::spec::VigType::SetPupil, true, 0);
}

Field &field(OpticalModel *opm, std::size_t i) { return *opm->optical_spec->fov->fields[i]; }

std::string d(double v) { return redukti::doubleToString(v); }

/** The dumper's output for this port. */
std::string dump() {
    std::string sb;
    auto opm = polychromatic_model();
    const auto &wvls = *opm->optical_spec->wvls;
    for (double foc : {0.0, 0.01}) {
        std::vector<PolyDiffractionMTF> mtfs;
        mtfs.emplace_back(opm.get(), field(opm.get(), 0), wvls.wavelengths, wvls.spectral_wts, foc, 64, 128);
        mtfs.emplace_back(opm.get(), field(opm.get(), 1), wvls.wavelengths, wvls.spectral_wts, foc, 64, 128);
        for (std::size_t fi = 0; fi < 2; fi++) {
            const auto &m = mtfs[fi];
            sb += "foc=" + d(foc) + " field=" + std::to_string(fi) + " image_pt=" + d(m.image_pt_2d.x) +
                  "," + d(m.image_pt_2d.y) + " mono=" + std::to_string(m.mono.size()) + '\n';
            for (double nu : {0.0, 10.0, 20.0, 40.0, 80.0, 160.0})
                sb += "nu=" + d(nu) + " s=" + d(m.sagittal(nu)) + " t=" + d(m.tangential(nu)) + '\n';
        }
        for (int freq : {10, 40}) {
            auto r = MTFResultByFreq::from_diffraction(opm.get(), mtfs, freq);
            sb += "by_freq " + std::to_string(r.freq) + " s=" + d(r.sag_mtf_by_field[0]) + "," +
                  d(r.sag_mtf_by_field[1]) + " t=" + d(r.tan_mtf_by_field[0]) + "," +
                  d(r.tan_mtf_by_field[1]) + '\n';
        }
    }
    return sb;
}

std::vector<std::string> split(const std::string &s, char separator) {
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t at; (at = s.find(separator, start)) != std::string::npos; start = at + 1)
        parts.push_back(s.substr(start, at - start));
    parts.push_back(s.substr(start));
    return parts;
}

const char *const EXPECTED =
    "foc=0.0 field=0 image_pt=0.0,0.0 mono=3\n"
    "nu=0.0 s=1.0000000000000022 t=1.0000000000000024\n"
    "nu=10.0 s=0.8466515704595721 t=0.8466515704595721\n"
    "nu=20.0 s=0.5320210236110688 t=0.5320210236110687\n"
    "nu=40.0 s=0.2134577283337443 t=0.2134577283337443\n"
    "nu=80.0 s=0.10099404278911912 t=0.10099404278911912\n"
    "nu=160.0 s=0.09042958370848468 t=0.09042958370848468\n"
    "foc=0.0 field=1 image_pt=0.0,21.298442142955004 mono=3\n"
    "nu=0.0 s=1.0000000000000016 t=1.0000000000000016\n"
    "nu=10.0 s=0.6083177117853783 t=0.715397394122629\n"
    "nu=20.0 s=0.45798110345875603 t=0.2773072903533144\n"
    "nu=40.0 s=0.3414256630232342 t=0.11818011838002689\n"
    "nu=80.0 s=0.17937651908698796 t=0.034683705569081395\n"
    "nu=160.0 s=0.06852497588109269 t=0.012672979217163953\n"
    "by_freq 10 s=0.8466515704595721,0.6083177117853783 t=0.8466515704595721,0.715397394122629\n"
    "by_freq 40 s=0.2134577283337443,0.3414256630232342 t=0.2134577283337443,0.11818011838002689\n"
    "foc=0.01 field=0 image_pt=0.0,0.0 mono=3\n"
    "nu=0.0 s=1.0000000000000044 t=1.0000000000000044\n"
    "nu=10.0 s=0.8160915034183728 t=0.8160915034183728\n"
    "nu=20.0 s=0.4570561533470175 t=0.4570561533470175\n"
    "nu=40.0 s=0.1752626109282727 t=0.1752626109282727\n"
    "nu=80.0 s=0.1039974131240945 t=0.10399741312409448\n"
    "nu=160.0 s=0.040347124587261696 t=0.04034712458726173\n"
    "foc=0.01 field=1 image_pt=0.0,21.301698929447085 mono=3\n"
    "nu=0.0 s=1.000000000000001 t=1.0000000000000016\n"
    "nu=10.0 s=0.6206469830123512 t=0.7053261499259084\n"
    "nu=20.0 s=0.4707258495402146 t=0.2570544163897505\n"
    "nu=40.0 s=0.3416844793068645 t=0.11975015674383833\n"
    "nu=80.0 s=0.16132864252867643 t=0.025320469997528756\n"
    "nu=160.0 s=0.06901531234233087 t=0.012828414294924326\n"
    "by_freq 10 s=0.8160915034183728,0.6206469830123512 t=0.8160915034183728,0.7053261499259084\n"
    "by_freq 40 s=0.1752626109282727,0.3416844793068645 t=0.1752626109282727,0.11975015674383833\n";

} // namespace

TEST(poly_diffraction_mtf_matches_jvm) {
    auto got = split(dump(), '\n');
    auto want = split(EXPECTED, '\n');
    CHECK_EQ(got.size(), want.size());
    for (std::size_t i = 0; i < got.size() && i < want.size(); i++) {
        if (got[i].rfind("foc=", 0) == 0) {
            // the image point comes from the ray trace alone
            CHECK_STR_EQ(got[i], want[i]);
            continue;
        }
        auto a = split(got[i], ' ');
        auto e = split(want[i], ' ');
        CHECK_EQ(a.size(), e.size());
        for (std::size_t k = 0; k < a.size() && k < e.size(); k++) {
            auto at = a[k].find('=');
            CHECK_STR_EQ(a[k].substr(0, at + 1), e[k].substr(0, e[k].find('=') + 1));
            if (at == std::string::npos || k == 0) continue;
            auto av = split(a[k].substr(at + 1), ',');
            auto ev = split(e[k].substr(e[k].find('=') + 1), ',');
            CHECK_EQ(av.size(), ev.size());
            for (std::size_t j = 0; j < av.size() && j < ev.size(); j++)
                CHECK_CLOSE(std::stod(av[j]), std::stod(ev[j]), 1e-12);
        }
        // the frequency label of each row
        CHECK_STR_EQ(a[0], e[0]);
    }
}

TEST(poly_diffraction_mtf_one_wavelength_is_its_monochromatic_mtf) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 1);
    double ref = opm->optical_spec->wvls->central_wvl();
    DiffractionMTF mono(opm.get(), fld, ref, 0.0, N, 2 * N);
    // the weight only scales the sum, which is divided by it
    PolyDiffractionMTF poly(opm.get(), fld, {ref, 500.0}, {2.5, 0.0}, 0.0, N, 2 * N);
    CHECK_EQ(poly.mono.size(), static_cast<std::size_t>(1));
    for (double nu : {0.0, 10.0, 40.0, 80.0}) {
        auto s = mono.sagittal_otf(nu);
        auto t = mono.tangential_otf(nu);
        CHECK_CLOSE(poly.sagittal(nu), std::hypot(s[0], s[1]), 1e-12);
        CHECK_CLOSE(poly.tangential(nu), std::hypot(t[0], t[1]), 1e-12);
    }
}

TEST(poly_diffraction_mtf_lateral_colour_lowers_the_polychromatic_mtf) {
    auto opm = polychromatic_model();
    const auto &wvls = *opm->optical_spec->wvls;
    auto &fld = field(opm.get(), 1);
    PolyDiffractionMTF poly(opm.get(), fld, wvls.wavelengths, wvls.spectral_wts, 0.0, N, 2 * N);
    CHECK_EQ(poly.mono.size(), wvls.wavelengths.size());
    CHECK_CLOSE(poly.tangential(0.0), 1.0, 1e-12);
    // |sum w OTF| / sum w never exceeds the weighted mean of the moduli; with 5 microns of
    // lateral colour in F it falls clearly below it in the tangential direction
    double sum = 0;
    for (double w : poly.weights) sum += w;
    for (double nu : {10.0, 20.0, 40.0, 80.0}) {
        double mean_s = 0, mean_t = 0;
        for (std::size_t i = 0; i < poly.mono.size(); i++) {
            auto s = poly.mono[i].sagittal_otf(nu);
            auto t = poly.mono[i].tangential_otf(nu);
            mean_s += poly.weights[i] * std::hypot(s[0], s[1]) / sum;
            mean_t += poly.weights[i] * std::hypot(t[0], t[1]) / sum;
        }
        CHECK(poly.sagittal(nu) <= mean_s + 1e-12);
        CHECK(poly.tangential(nu) <= mean_t + 1e-12);
        if (nu >= 20)
            CHECK(poly.tangential(nu) < mean_t - 1e-3);
    }
}

TEST(poly_diffraction_mtf_combines_intrinsic_complex_otfs_with_analytic_lateral_colour_phase) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 1);
    const auto &wvls = *opm->optical_spec->wvls;
    double focus = 0.01;
    PolyDiffractionMTF poly(opm.get(), fld, wvls.wavelengths, wvls.spectral_wts, focus, N, 2 * N);
    std::vector<DiffractionMTF> own;
    for (double w : wvls.wavelengths)
        own.emplace_back(opm.get(), fld, w, focus, N, 2 * N);
    for (int axis : {0, 1})
        for (double nu : {-80.0, -40.0, 0.0, 40.0, 80.0}) {
            double re = 0, im = 0, sum = 0;
            for (std::size_t wi = 0; wi < own.size(); wi++) {
                auto point = DiffractionMTF::image_point(opm.get(), fld, wvls.wavelengths[wi], focus);
                double delta = axis == 0 ? point.x - poly.image_pt_2d.x : point.y - poly.image_pt_2d.y;
                double angle = -2 * PI * nu * delta;
                auto otf = axis == 0 ? own[wi].sagittal_otf(nu) : own[wi].tangential_otf(nu);
                double weight = wvls.spectral_wts[wi];
                re += weight * (otf[0] * std::cos(angle) - otf[1] * std::sin(angle));
                im += weight * (otf[0] * std::sin(angle) + otf[1] * std::cos(angle));
                sum += weight;
            }
            CHECK_CLOSE(axis == 0 ? poly.sagittal(nu) : poly.tangential(nu),
                        std::hypot(re, im) / sum, 1e-12);
        }
}

TEST(poly_diffraction_mtf_results_by_frequency_are_in_cycles_per_millimetre) {
    auto opm = polychromatic_model();
    const auto &wvls = *opm->optical_spec->wvls;
    std::vector<PolyDiffractionMTF> mtfs;
    mtfs.emplace_back(opm.get(), field(opm.get(), 0), wvls.wavelengths, wvls.spectral_wts, 0.0, N, 2 * N);
    mtfs.emplace_back(opm.get(), field(opm.get(), 1), wvls.wavelengths, wvls.spectral_wts, 0.0, N, 2 * N);
    auto result = MTFResultByFreq::from_diffraction(opm.get(), mtfs, 40);
    CHECK_EQ(result.freq, 40);
    for (std::size_t fi = 0; fi < 2; fi++) {
        CHECK_CLOSE(result.sag_mtf_by_field[fi], mtfs[fi].sagittal(40.0), 1e-15);
        CHECK_CLOSE(result.tan_mtf_by_field[fi], mtfs[fi].tangential(40.0), 1e-15);
    }
}

TEST(poly_diffraction_mtf_rejects_weights_that_do_not_match_or_are_all_zero) {
    auto opm = polychromatic_model();
    auto &fld = field(opm.get(), 0);
    double ref = opm->optical_spec->wvls->central_wvl();
    CHECK_THROWS(PolyDiffractionMTF(opm.get(), fld, {ref}, {1.0, 1.0}, 0.0, N, 2 * N),
                 redukti::IllegalArgumentException);
    CHECK_THROWS(PolyDiffractionMTF(opm.get(), fld, {ref}, {0.0}, 0.0, N, 2 * N),
                 redukti::IllegalArgumentException);
}
