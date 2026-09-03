// Checks the Zemax exporter against the .zmx file committed under Examples/.
//
// The comparison target is the repository file itself rather than a dump from
// the JVM: the Java was confirmed to reproduce that file byte for byte, so this
// pins the C++ against an artifact already in the tree. It is the strongest
// check available -- the whole pipeline from prescription file to exported lens
// runs, and the answer is a file a person has already looked at.
//
// canon-rf70-200mm-f2.8LZ is a recent two-configuration zoom, so it exercises
// the MNUM/SDIA/THIC configuration block as well as aspheres and named glasses.
#include "TestHarness.h"

#include "redukti/exporters/ZemaxExporter.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/spec/Prescription.h"

#include <fstream>
#include <string>
#include <vector>

namespace {

using redukti::exporters::ZemaxExporter;
using redukti::importers::OpticalBenchDataImporter;
using redukti::spec::Prescription;

const char *const DIR = REDUKTI_EXAMPLES_DIR "canon-rf70-200mm-f2.8LZ/";

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        std::string line;
        if (nl == std::string::npos) {
            line = text.substr(start);
            start = text.size();
        } else {
            line = text.substr(start, nl - start);
            start = nl + 1;
        }
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (!line.empty())
            out.push_back(line);
    }
    return out;
}

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    CHECK(in.good());
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

} // namespace

TEST(zemax_exporter_matches_committed_zmx) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(std::string(DIR) + "US20250155694_Example01P.txt");
    Prescription p = Prescription::build_prescription(specs, true, false, false);

    std::string generated = ZemaxExporter().generate(p, false);
    std::string committed = readFile(std::string(DIR) + "US20250155694_Example01P.zmx");

    auto actual = lines(generated);
    auto expected = lines(committed);
    CHECK_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size() && i < expected.size(); i++)
        CHECK_STR_EQ(actual[i], expected[i]);
}

TEST(zemax_exporter_d_line_only_switches_wavelengths) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_file(std::string(DIR) + "US20250155694_Example01P.txt");
    Prescription p = Prescription::build_prescription(specs, true, false, false);

    std::string generated = ZemaxExporter().generate(p, true);
    auto rows = lines(generated);
    // The d-line branch replaces the first three WAVM rows and leaves the rest.
    bool found1 = false, found2 = false, found3 = false;
    for (const auto &r : rows) {
        if (r == "WAVM 1 0.5875618 1")
            found1 = true;
        if (r == "WAVM 2 0.550 0")
            found2 = true;
        if (r == "WAVM 3 0.550 0")
            found3 = true;
        CHECK(r != "WAVM 1 0.4861327 1");
    }
    CHECK(found1);
    CHECK(found2);
    CHECK(found3);
}
