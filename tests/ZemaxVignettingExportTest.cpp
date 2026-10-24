// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.exporters.ZemaxVignettingExportTest.
//
// The Java tests are ported as they are. In addition, the whole .zmx text of each case
// is asserted against the Java ZemaxExporter (DumpZemaxVignetting.java, a throwaway
// dumper); nothing here involves a trace, so it matches byte for byte.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/exporters/ZemaxExporter.h"
#include "redukti/spec/Prescription.h"
#include "redukti/tools/LensTool2.h"

#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>

namespace {

using redukti::exporters::ZemaxExporter;
using redukti::spec::Prescription;
using redukti::tools::LensTool2;
using FieldVignetting = Prescription::FieldVignetting;

const std::string INPUT = "[descriptive data]\n"
                          "title\tTest lens\n"
                          "[variable distances]\n"
                          "Focal Length\t50\t60\n"
                          "F-Number\t4\t5\n"
                          "Angle of View\t40\t35\n"
                          "Bf\t45\t55\n"
                          "[lens data]\n"
                          "1\t50\t4\t1.5\t20\t60\n"
                          "2\t-50\tBf\t\t20\n"
                          "[report data]\n"
                          "lens name\tTest lens\n"
                          "scenarios\t0\t1\n"
                          "names\tWide\tLong\n";

std::string single() {
    std::string s = INPUT;
    const std::string scenarios = "scenarios\t0\t1\nnames\tWide\tLong\n";
    s.replace(s.find(scenarios), scenarios.size(), "");
    return s;
}

Prescription prescription(const std::string &input) {
    LensTool2::LensSpecifications specs;
    specs.parse_buffer(input);
    return LensTool2::createPrescription(specs, true, false);
}

std::vector<std::string> words(const std::string &line) {
    std::istringstream in(line);
    std::vector<std::string> out;
    for (std::string w; in >> w;) out.push_back(w);
    return out;
}

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    for (std::size_t at; (at = text.find('\n', start)) != std::string::npos; start = at + 1)
        out.push_back(text.substr(start, at - start));
    out.push_back(text.substr(start));
    return out;
}

/** The numbers following key on its line. */
std::vector<double> values(const std::string &zmx, const std::string &key) {
    for (const auto &line : lines(zmx))
        if (line.rfind(key + " ", 0) == 0) {
            auto w = words(line);
            std::vector<double> v;
            for (std::size_t i = 1; i < w.size(); i++) v.push_back(std::stod(w[i]));
            return v;
        }
    ::redukti::test::reportFailure(__FILE__, __LINE__, key + " not found");
    return {};
}

/** The value of the multi-configuration operand op for field and config, 1-based. */
double operand(const std::string &zmx, const std::string &op, int field, int config) {
    for (const auto &line : lines(zmx)) {
        auto w = words(line);
        if (w.size() > 3 && w[0] == op && std::stoi(w[1]) == field && std::stoi(w[2]) == config)
            return std::stod(w[3]);
    }
    ::redukti::test::reportFailure(__FILE__, __LINE__, op + " not found");
    return 0.0;
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

const char *const EXPECTED_PLAIN = R"zmx(VERS 161019 507 33785
MODE SEQ
NAME Test lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
FNUM 4.0
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT HOYA OHARA NIKON-HIKARI NIKON HIKARI SUMITA CDGM CORNING LACROIX
RAIM 0 2 1 1 0 1 0 0 0
SDMA 0 1 0
FTYP 3 0 9 3 0 0 0 9
ROPD 2
HYPR 0
PICB 1
XFLN 0 0 0 0 0 0 0 0 0 0 0 0
YFLN 0 2.16 4.33 6.5 8.65 10.8 12.98 15.14 21.600000 0 0 0 0 0 0 0 0
FWGN 10 9 9 8 8 7 7 7 3 1 1 1
VDXN 0 0 0 0 0 0 0 0 0 0 0 0
VDYN 0 0 0 0 0 0 0 0 0 0 0 0
VCXN 0 0 0 0 0 0 0 0 0 0 0 0
VCYN 0 0 0 0 0 0 0 0 0 0 0 0
VANN 0 0 0 0 0 0 0 0 0 0 0 0
WAVM 1 0.4861327 1
WAVM 2 0.5875618 1
WAVM 3 0.6562725 1
WAVM 4 0.550 0
WAVM 5 0.550 0
WAVM 6 0.550 0
WAVM 7 0.550 0
WAVM 8 0.550 0
WAVM 9 0.550 0
WAVM 10 0.550 0
WAVM 11 0.550 0
WAVM 12 0.550 0
WAVM 13 0.550 0
WAVM 14 0.550 0
WAVM 15 0.550 0
WAVM 16 0.550 0
WAVM 17 0.550 0
WAVM 18 0.550 0
WAVM 19 0.550 0
WAVM 20 0.550 0
WAVM 21 0.550 0
WAVM 22 0.550 0
WAVM 23 0.550 0
WAVM 24 0.550 0
PWAV 2
POLS 1 0 1 0 0 1 0
GSTD 0 100 100 100 100 100 100 0 1 1 0 0 1 1 1 1 1 1
NSCD 100 500 0 1.0E-3 5 1.0E-6 0 0 0 0 0 0 1000000 0 2
SURF 0
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ INFINITY
  DIAM 0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 1
  TYPE STANDARD
  CURV 0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 4.0
  GLAS ___BLANK 1 0 1.5 60.0 0 0 0 0 0 0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 2
  TYPE STANDARD
  CURV -0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 45.0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 3
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 0
  DIAM 21.600000 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
TOL TOFF   0   0              0              0   0 0 0 0
)zmx";

const char *const EXPECTED_SINGLE = R"zmx(VERS 161019 507 33785
MODE SEQ
NAME Test lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
FNUM 4.0
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT HOYA OHARA NIKON-HIKARI NIKON HIKARI SUMITA CDGM CORNING LACROIX
RAIM 0 2 1 1 0 1 0 0 0
SDMA 0 1 0
FTYP 0 0 3 3 0 0 0 3
ROPD 2
HYPR 0
PICB 1
XFLN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
YFLN 0.0 10.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
FWGN 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VDXN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VDYN 0.0 -0.010000000000000002 -0.04999999999999999 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VCXN 0.0 0.02 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VCYN 0.0 0.09 0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VANN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
WAVM 1 0.4861327 1
WAVM 2 0.5875618 1
WAVM 3 0.6562725 1
WAVM 4 0.550 0
WAVM 5 0.550 0
WAVM 6 0.550 0
WAVM 7 0.550 0
WAVM 8 0.550 0
WAVM 9 0.550 0
WAVM 10 0.550 0
WAVM 11 0.550 0
WAVM 12 0.550 0
WAVM 13 0.550 0
WAVM 14 0.550 0
WAVM 15 0.550 0
WAVM 16 0.550 0
WAVM 17 0.550 0
WAVM 18 0.550 0
WAVM 19 0.550 0
WAVM 20 0.550 0
WAVM 21 0.550 0
WAVM 22 0.550 0
WAVM 23 0.550 0
WAVM 24 0.550 0
PWAV 2
POLS 1 0 1 0 0 1 0
GSTD 0 100 100 100 100 100 100 0 1 1 0 0 1 1 1 1 1 1
NSCD 100 500 0 1.0E-3 5 1.0E-6 0 0 0 0 0 0 1000000 0 2
SURF 0
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ INFINITY
  DIAM 0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 1
  TYPE STANDARD
  CURV 0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 4.0
  GLAS ___BLANK 1 0 1.5 60.0 0 0 0 0 0 0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 2
  TYPE STANDARD
  CURV -0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 45.0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 3
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 0
  DIAM 21.600000 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
TOL TOFF   0   0              0              0   0 0 0 0
)zmx";

const char *const EXPECTED_MULTI = R"zmx(VERS 161019 507 33785
MODE SEQ
NAME Test lens
PFIL 0 0 0
LANG 0
UNIT MM X W X CM MR CPMM
FLOA
ENVD 20 1 0
GFAC 0 0
GCAT SCHOTT HOYA OHARA NIKON-HIKARI NIKON HIKARI SUMITA CDGM CORNING LACROIX
RAIM 0 2 1 1 0 1 0 0 0
SDMA 0 1 0
FTYP 0 0 2 3 0 0 0 2
ROPD 2
HYPR 0
PICB 1
XFLN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
YFLN 0.0 20.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
FWGN 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VDXN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VDYN 0.0 -0.04999999999999999 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VCXN 0.0 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VCYN 0.0 0.25 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
VANN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
WAVM 1 0.5875618 1
WAVM 2 0.550 0
WAVM 3 0.550 0
WAVM 4 0.550 0
WAVM 5 0.550 0
WAVM 6 0.550 0
WAVM 7 0.550 0
WAVM 8 0.550 0
WAVM 9 0.550 0
WAVM 10 0.550 0
WAVM 11 0.550 0
WAVM 12 0.550 0
WAVM 13 0.550 0
WAVM 14 0.550 0
WAVM 15 0.550 0
WAVM 16 0.550 0
WAVM 17 0.550 0
WAVM 18 0.550 0
WAVM 19 0.550 0
WAVM 20 0.550 0
WAVM 21 0.550 0
WAVM 22 0.550 0
WAVM 23 0.550 0
WAVM 24 0.550 0
PWAV 2
POLS 1 0 1 0 0 1 0
GSTD 0 100 100 100 100 100 100 0 1 1 0 0 1 1 1 1 1 1
NSCD 100 500 0 1.0E-3 5 1.0E-6 0 0 0 0 0 0 1000000 0 2
SURF 0
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ INFINITY
  DIAM 0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 1
  TYPE STANDARD
  CURV 0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 4.0
  GLAS ___BLANK 1 0 1.5 60.0 0 0 0 0 0 0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 2
  TYPE STANDARD
  CURV -0.02 0 0 0 0
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 45.0
  DIAM 10.0 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
SURF 3
  TYPE STANDARD
  CURV 0.0 0 0 0 0 ""
  HIDE 0 0 0 0 0 0 0 0 0 0
  MIRR 2 0
  DISZ 0
  DIAM 21.600000 1 0 0 1 ""
  POPS 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0
TOL TOFF   0   0              0              0   0 0 0 0
MNUM 2 1
THIC    2   1 45.0  0 0 0 1 1 1 0 0
THIC    2   2 55.0  0 0 0 1 1 1 0 0
YFIE    1   1 0.0  0 0 0 1 1 1 0 0
YFIE    1   2 0.0  0 0 0 1 1 1 0 0
FVDX    1   1 0.0  0 0 0 1 1 1 0 0
FVDX    1   2 0.0  0 0 0 1 1 1 0 0
FVDY    1   1 0.0  0 0 0 1 1 1 0 0
FVDY    1   2 0.0  0 0 0 1 1 1 0 0
FVCX    1   1 0.0  0 0 0 1 1 1 0 0
FVCX    1   2 0.0  0 0 0 1 1 1 0 0
FVCY    1   1 0.0  0 0 0 1 1 1 0 0
FVCY    1   2 0.0  0 0 0 1 1 1 0 0
YFIE    2   1 20.0  0 0 0 1 1 1 0 0
YFIE    2   2 17.5  0 0 0 1 1 1 0 0
FVDX    2   1 0.0  0 0 0 1 1 1 0 0
FVDX    2   2 0.0  0 0 0 1 1 1 0 0
FVDY    2   1 -0.04999999999999999  0 0 0 1 1 1 0 0
FVDY    2   2 -0.15000000000000002  0 0 0 1 1 1 0 0
FVCX    2   1 0.1  0 0 0 1 1 1 0 0
FVCX    2   2 0.05  0 0 0 1 1 1 0 0
FVCY    2   1 0.25  0 0 0 1 1 1 0 0
FVCY    2   2 0.25  0 0 0 1 1 1 0 0
)zmx";

} // namespace

TEST(zemax_vignetting_without_vignetting_the_fields_are_unchanged) {
    std::string zmx = ZemaxExporter().generate(prescription(single()), false);
    CHECK(contains(zmx, "SDMA 0 1 0\nFTYP 3 0 9 3 0 0 0 9\nROPD 2\nHYPR 0\nPICB 1\n"
                        "XFLN 0 0 0 0 0 0 0 0 0 0 0 0\nYFLN 0 2.16 4.33 6.5 8.65 10.8 12.98 15.14 "));
    CHECK(contains(zmx, "VCYN 0 0 0 0 0 0 0 0 0 0 0 0\n"));
    CHECK(!contains(zmx, "FVDY"));
    CHECK_STR_EQ(zmx, EXPECTED_PLAIN);
}

TEST(zemax_vignetting_angles_and_factors_of_one_configuration) {
    auto p = prescription(single());
    p._vignetting = std::vector<FieldVignetting>{
        {0, 0.0, 0, 0, 0, 0}, {0, 0.5, 0.02, 0.02, 0.1, 0.08}, {0, 1.0, 0.1, 0.1, 0.3, 0.2}};
    std::string zmx = ZemaxExporter().generate(p, false);
    CHECK(contains(zmx, "\nFTYP 0 0 3 3 0 0 0 3\n"));
    auto y = values(zmx, "YFLN");
    CHECK_EQ(y.size(), static_cast<std::size_t>(12));
    if (y.size() == 12) {
        CHECK_CLOSE(y[0], 0.0, 1e-12);
        CHECK_CLOSE(y[1], 10.0, 1e-12);
        CHECK_CLOSE(y[2], 20.0, 1e-12);
        for (std::size_t i = 3; i < 12; i++) CHECK_EQ(y[i], 0.0);
    }
    auto vdx = values(zmx, "VDXN");
    auto vdy = values(zmx, "VDYN");
    auto vcx = values(zmx, "VCXN");
    auto vcy = values(zmx, "VCYN");
    // Zemax maps P to VD + P(1 - VC): the pupil edges land where our factors put them
    for (std::size_t i = 0; i < 3; i++) {
        const auto &v = (*p._vignetting)[i];
        CHECK_CLOSE(vdy[i] + (1 - vcy[i]), 1 - v.vuy, 1e-12);
        CHECK_CLOSE(vdy[i] - (1 - vcy[i]), -1 + v.vly, 1e-12);
        CHECK_CLOSE(vdx[i] + (1 - vcx[i]), 1 - v.vux, 1e-12);
        CHECK_CLOSE(vdx[i] - (1 - vcx[i]), -1 + v.vlx, 1e-12);
    }
    CHECK(!contains(zmx, "MNUM"));
    CHECK_STR_EQ(zmx, EXPECTED_SINGLE);
}

TEST(zemax_vignetting_each_configuration_has_its_own_angles_and_factors) {
    auto p = prescription(INPUT);
    p._vignetting = std::vector<FieldVignetting>{{0, 0.0, 0, 0, 0, 0},
                                                 {0, 1.0, 0.1, 0.1, 0.3, 0.2},
                                                 {1, 0.0, 0, 0, 0, 0},
                                                 {1, 1.0, 0.05, 0.05, 0.4, 0.1}};
    std::string zmx = ZemaxExporter().generate(p, false);
    CHECK(contains(zmx, "\nMNUM 2 1\n"));
    // half of 40 and of 35 degrees
    CHECK_CLOSE(operand(zmx, "YFIE", 2, 1), 20.0, 1e-12);
    CHECK_CLOSE(operand(zmx, "YFIE", 2, 2), 17.5, 1e-12);
    CHECK_CLOSE(operand(zmx, "YFIE", 1, 2), 0.0, 1e-12);
    CHECK_CLOSE(operand(zmx, "FVDY", 2, 1), (0.2 - 0.3) / 2, 1e-12);
    CHECK_CLOSE(operand(zmx, "FVDY", 2, 2), (0.1 - 0.4) / 2, 1e-12);
    CHECK_CLOSE(operand(zmx, "FVCY", 2, 2), (0.4 + 0.1) / 2, 1e-12);
    CHECK_CLOSE(operand(zmx, "FVCX", 2, 2), 0.05, 1e-12);
    CHECK_CLOSE(operand(zmx, "FVDX", 2, 2), 0.0, 1e-12);
    // one operand row holds all its configurations, in order
    CHECK(contains(zmx, "YFIE    2   1 20.0  0 0 0 1 1 1 0 0\nYFIE    2   2 17.5  0 0 0 1 1 1 0 0\n"));
    // the dumper writes this case with d_line_only, to cover the other wavelength block
    CHECK_STR_EQ(ZemaxExporter().generate(p, true), EXPECTED_MULTI);
}

TEST(zemax_vignetting_rejects_configurations_with_different_fields) {
    auto p = prescription(INPUT);
    std::vector<FieldVignetting> rows{{0, 0.0, 0, 0, 0, 0}, {0, 1.0, 0, 0, 0, 0}, {1, 0.0, 0, 0, 0, 0}};
    p._vignetting = rows;
    CHECK_THROWS(ZemaxExporter().generate(p, false), redukti::IllegalArgumentException);
    rows.push_back({1, 0.7, 0, 0, 0, 0});
    p._vignetting = rows;
    CHECK_THROWS(ZemaxExporter().generate(p, false), redukti::IllegalArgumentException);
}
