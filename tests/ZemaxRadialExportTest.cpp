// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.exporters.ZemaxRadialExportTest.
//
// A radial asphere exports as XOSPHERE with its terms in the Extra Data section, so the
// omitted odd powers keep their slots; an even asphere still uses PARM.
#include "TestHarness.h"

#include "redukti/Exceptions.h"
#include "redukti/exporters/ZemaxExporter.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/spec/Prescription.h"
#include "redukti/spec/SurfaceType.h"

#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

using redukti::exporters::ZemaxExporter;
using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;
using redukti::spec::SurfaceType;

const char *const SIGMA =
    REDUKTI_EXAMPLES_DIR "sigma-14-24mm-f2.8-ml/JP2020-042221_Example01P.txt";

std::string surface(const std::string &zmx, int number) {
    std::string marker = "SURF " + std::to_string(number) + "\n";
    std::size_t start = zmx.find(marker);
    CHECK(start != std::string::npos);
    if (start == std::string::npos)
        return "";
    std::size_t end = zmx.find("\nSURF ", start);
    return zmx.substr(start, end == std::string::npos ? std::string::npos : end - start);
}

/** The block's XDAT lines, keyed by their index; a repeated index is a failure. */
std::map<int, double> extraData(const std::string &block) {
    std::map<int, double> result;
    std::istringstream lines(block);
    std::string line;
    while (std::getline(lines, line)) {
        std::istringstream fields(line);
        std::string name;
        if (!(fields >> name) || name != "XDAT")
            continue;
        int index = 0;
        double value = 0;
        fields >> index >> value;
        CHECK(result.find(index) == result.end());
        result[index] = value;
    }
    return result;
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

} // namespace

TEST(zemax_radial_sigma_patent_has_nine_odd_powers_and_twenty_radial_slots) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(SIGMA);
    auto prescription = Prescription::build_prescription(specs, false);
    CHECK_EQ(prescription._aspherical_odd_count, 9);
    // JP2020-042221, Example 1, supplied patent table: A3 through A20.
    const int ids[] = {1, 5, 6, 31, 32};
    const std::vector<std::vector<double>> patent = {
        {0, 8.58209e-6, 0, -1.40764e-8, 0, 3.05748e-11, 0, -5.97803e-14, 0, 9.08590e-17, 0,
         -9.58737e-20, 0, 6.40051e-23, 0, -2.39147e-26, 0, 3.78519e-30},
        {-5.23111e-5, -1.26716e-5, -1.13040e-5, 1.95245e-6, -9.38134e-8, -7.82976e-10,
         1.22496e-10, 1.97968e-12, -1.33295e-14, -1.10265e-14, 5.32582e-17, 7.28532e-18,
         1.89244e-19, -1.81192e-20, 1.13899e-21, -2.99255e-23, 1.48595e-25, -3.31214e-27},
        {-3.48559e-5, -1.50563e-5, -1.10043e-5, 1.72222e-6, -4.71099e-8, -3.15483e-9,
         -5.86163e-11, 1.76453e-11, -1.10783e-13, -5.28181e-15, -8.35047e-16, 5.89441e-17,
         -9.54814e-18, 3.21284e-19, 6.43253e-21, -4.91029e-23, 1.37261e-24, -5.09803e-25},
        {0, -2.45667e-5, 0, -8.17092e-8, 0, 2.81370e-9, 0, -7.26008e-11, 0, 1.11778e-12, 0,
         -9.83681e-15, 0, 4.86452e-17, 0, -1.24975e-19, 0, 1.29336e-22},
        {0, 5.66654e-6, 0, 5.42201e-8, 0, -2.06458e-9, 0, 3.25090e-11, 0, -2.56410e-13, 0,
         1.03356e-15, 0, -1.78037e-18, 0, -8.07610e-22, 0, 5.22718e-24}};
    std::string exported = ZemaxExporter().generate(prescription, true);
    CHECK(!contains(exported, "ODDASPHE"));
    for (std::size_t i = 0; i < 5; i++) {
        const auto &s = prescription.get_surfaces()[static_cast<std::size_t>(ids[i] - 1)];
        std::vector<double> expected(20, 0.0);
        for (std::size_t k = 0; k < 18; k++)
            expected[k + 2] = patent[i][k];
        CHECK_EQ(s._asph_type, SurfaceType::ASPH_RADIAL);
        CHECK_EQ(s._coeffs->size(), expected.size());
        for (std::size_t k = 0; k < expected.size() && k < s._coeffs->size(); k++)
            CHECK_EQ((*s._coeffs)[k], expected[k]);
        std::string block = surface(exported, ids[i]);
        CHECK(contains(block, "TYPE XOSPHERE\n"));
        CHECK(!contains(block, "PARM "));
        std::string conic = ids[i] == 6 ? "-0.0364842" : "0.0";
        CHECK(contains(block, "CONI " + conic + "\n"));
        auto xd = extraData(block);
        CHECK_EQ(xd.size(), static_cast<std::size_t>(22));
        CHECK_EQ(xd[1], 20.0);
        CHECK_EQ(xd[2], 1.0);
        for (int power = 1; power <= 20; power++)
            CHECK_CLOSE(xd[power + 2], expected[static_cast<std::size_t>(power - 1)],
                        1e-30);
    }
    OpticalBenchDataImporter::LensSpecifications restored;
    std::string optBench;
    prescription.to_opt_bench_str(optBench);
    restored.parse_buffer(optBench);
    CHECK_EQ(restored.get_aspherical_odd_count(), 9);
    std::string readme;
    prescription.to_markdown_str(readme);
    CHECK(contains(readme, "| RADIAL |"));
    CHECK(!contains(readme, "| ODD |"));
}

TEST(zemax_radial_sparse_export_retains_zeros_and_conic_while_even_keeps_parameters) {
    auto p = Prescription(50, 4, 20, 30, true)
                 .surf(50, 5, 20)
                 .asph(SurfaceType::ASPH_RADIAL, -0.25,
                       {0, 0, 1e-6, 2e-7, 0, 3e-9, 0, 4e-11, 0, 5e-13, 0, 6e-15})
                 .surf(-50, 40, 20)
                 .asph(SurfaceType::ASPH_EVEN, -0.5, {0, 1e-6})
                 .build();
    std::string zmx = ZemaxExporter().generate(p, true);
    std::string radial = surface(zmx, 1);
    CHECK(contains(radial, "CONI -0.25\n"));
    auto xd = extraData(radial);
    CHECK_EQ(xd[1], 12.0);
    for (int power : {1, 2, 5, 7, 9, 11})
        CHECK_EQ(xd[power + 2], 0.0);
    CHECK_CLOSE(xd[14], 6e-15, 1e-30);
    std::string even = surface(zmx, 2);
    CHECK(contains(even, "TYPE EVENASPH\n"));
    CHECK(contains(even, "CONI -0.5\n"));
    CHECK(contains(even, "PARM 1 0.0\n"));
    CHECK(contains(even, "PARM 2 1.0E-6\n"));
    CHECK(!contains(even, "XDAT"));
}

TEST(zemax_radial_rejects_orders_beyond_zemax_limit) {
    auto p = Prescription(50, 4, 20, 30, true)
                 .surf(50, 5, 20)
                 .asph(SurfaceType::ASPH_RADIAL, 0, std::vector<double>(241, 0.0))
                 .build();
    CHECK_THROWS(ZemaxExporter().generate(p, true), redukti::IllegalArgumentException);
}
