// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.importers.obench.OpticalBenchAsphereTest.
//
// How AsphericalOddCount maps Optical Bench's packed columns onto the normalized
// coefficient array, and that the mapping survives a round trip through
// to_opt_bench_str. The Java writes each fixture to a temp file and parses it; this
// parses the same text from a buffer.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/elem/profiles/RadialPolynomial.h"
#include "redukti/spec/Prescription.h"
#include "redukti/spec/SurfaceType.h"
#include "redukti/tools/LensTool2.h"

#include <cmath>
#include <string>
#include <vector>

namespace {

using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;
using redukti::spec::SurfaceType;
using redukti::tools::LensTool2;
using AsphereType = OpticalBenchDataImporter::AsphereType;

OpticalBenchDataImporter::LensSpecifications parse(const std::string &text) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return specs;
}

OpticalBenchDataImporter::LensSpecifications lens(const std::string &constants,
                                                  const std::string &coefficients) {
    return parse("[constants]\n" + constants + "\n[variable distances]\n" +
                 "Focal Length\t36\nF-Number\t1.45\nAngle of View\t63.06\n" +
                 "[lens data]\n1\t-224.2964\t2.5\t1.7433\t29.46\t49.32\n" +
                 "[aspherical data]\n1\t-224.2964\t195\t" + coefficients + "\n");
}

void checkCoeffs(const std::vector<double> &actual, const std::vector<double> &expected,
                 double tol = 0.0) {
    CHECK_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); i++)
        CHECK_CLOSE(actual[i], expected[i], tol);
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

} // namespace

TEST(obench_asphere_nikkor35mm_coefficients_and_sag) {
    // Haruo Sato, US7663816B2, Table 1: https://patents.google.com/patent/US7663816B2/en
    auto specs = lens("AsphericalOddCount\t1",
                      "-2.0873E-07\t-1.2426E-05\t2.7998E-09\t-5.1736E-11\t1.7973E-13\t"
                      "-8.9748E-17");
    std::vector<double> expected{0,        0,         -2.0873e-7, -1.2426e-5,
                                 0,        2.7998e-9, 0,          -5.1736e-11,
                                 0,        1.7973e-13, 0,         -8.9748e-17};
    auto prescription = Prescription::build_prescription(specs, false);
    CHECK_EQ(prescription._aspherical_odd_count, 1);
    const auto &surface = prescription.get_surfaces()[0];
    checkCoeffs(surface.get_aspheric_coeffs(), expected);
    std::string readme = LensTool2::startREADME(specs);
    CHECK(contains(readme, "| ID  | Type | k   | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9"
                           " | P10 | P11 | P12 |\n"));
    CHECK(contains(readme,
                   "| 1| RADIAL | 195.0 | 0.0 | 0.0 | -2.0873E-7 | -1.2426E-5 | 0.0 | "
                   "2.7998E-9 | 0.0 | -5.1736E-11 | 0.0 | 1.7973E-13 | 0.0 | -8.9748E-17 |\n"));

    redukti::rayoptics::elem::profiles::RadialPolynomial profile;
    profile.r(surface._radius)->cc(surface._k)->setCoefs(*surface._coeffs);
    double r = 5;
    double c = 1 / surface._radius;
    double conic = c * r * r / (1 + std::sqrt(1 - 196 * c * c * r * r));
    double sag = conic - 2.0873e-7 * std::pow(r, 3) - 1.2426e-5 * std::pow(r, 4) +
                 2.7998e-9 * std::pow(r, 6) - 5.1736e-11 * std::pow(r, 8) +
                 1.7973e-13 * std::pow(r, 10) - 8.9748e-17 * std::pow(r, 12);
    CHECK_CLOSE(profile.sag(3, 4), sag, 1e-14);
    double h = 1e-5;
    CHECK_CLOSE(profile.df(redukti::mathlib::Vector3(3, 4, sag)).x,
                -(profile.sag(3 + h, 4) - profile.sag(3 - h, 4)) / (2 * h), 1e-10);

    std::string optBench;
    prescription.to_opt_bench_str(optBench);
    auto restored = Prescription::build_prescription(parse(optBench), false);
    CHECK_EQ(restored._aspherical_odd_count, 1);
    checkCoeffs(*restored.get_surfaces()[0]._coeffs, expected);
}

TEST(obench_asphere_two_odd_terms_are_interleaved_before_even_only_tail) {
    auto specs = lens("AsphericalOddCount\t2", "3\t4\t5\t6\t8\t10");
    checkCoeffs(specs.get_surfaces()[0].get_aspherical_data()->get_coeffs(),
                {0, 0, 3, 4, 5, 6, 0, 8, 0, 10});
    auto p = Prescription::build_prescription(specs, false);
    std::string optBench;
    p.to_opt_bench_str(optBench);
    auto restored = Prescription::build_prescription(parse(optBench), false);
    CHECK_EQ(restored._aspherical_odd_count, 2);
    checkCoeffs(*restored.get_surfaces()[0]._coeffs, *p.get_surfaces()[0]._coeffs);
}

TEST(obench_asphere_preserves_declared_count_even_when_odd_terms_are_zero) {
    auto p = Prescription::build_prescription(lens("AsphericalOddCount\t7", "0\t4\t0\t6"),
                                              false);
    std::string optBench;
    p.to_opt_bench_str(optBench);
    auto restored = Prescription::build_prescription(parse(optBench), false);
    CHECK_EQ(restored._aspherical_odd_count, 7);
    checkCoeffs(*restored.get_surfaces()[0]._coeffs, *p.get_surfaces()[0]._coeffs);
}

TEST(obench_asphere_zero_count_and_even_formats_keep_even_powers) {
    for (const char *constants : {"", "AsphericalOddCount\t0"}) {
        auto specs = lens(constants, "4\t6\t8");
        const auto *asphere = specs.get_surfaces()[0].get_aspherical_data();
        CHECK(asphere->get_asphere_type() == AsphereType::Even);
        checkCoeffs(asphere->get_coeffs(), {0, 4, 6, 8});
    }
    auto a2specs = lens("AsphericalA2", "2\t4\t6");
    checkCoeffs(a2specs.get_surfaces()[0].get_aspherical_data()->get_coeffs(), {2, 4, 6});
}

TEST(obench_asphere_export_retains_new_odd_terms_and_pads_even_surfaces) {
    auto p = Prescription::build_prescription(lens("AsphericalOddCount\t1", "3\t4\t6"),
                                              false);
    (*p.get_surfaces()[0]._coeffs)[4] = 5; // newly introduced A5
    p.surf(100, 1, 20).asph(SurfaceType::ASPH_EVEN, 0, {0, 4, 6}).build();
    std::string optBench;
    p.to_opt_bench_str(optBench);
    auto restored = Prescription::build_prescription(parse(optBench), false);
    CHECK_EQ(restored._aspherical_odd_count, 2);
    checkCoeffs(*restored.get_surfaces()[0]._coeffs, {0, 0, 3, 4, 5, 6});
    checkCoeffs(*restored.get_surfaces()[1]._coeffs, {0, 0, 0, 4, 0, 6});
}

TEST(obench_asphere_invalid_counts_are_rejected) {
    for (const char *count : {"-1", "1.5", "bad", ""})
        CHECK_THROWS(lens(std::string("AsphericalOddCount\t") + count, "3\t4"),
                     redukti::IllegalArgumentException);
}
