// C++ port of org.redukti.examples.OptimizationTrialExamplesTest.
//
// The Java checks the four trials Documentation/OPTIMIZER.md prints against the setups its
// examples package builds in code. Only the Otus is ported here: this repository carries no
// examples package, and of the four lenses only the Otus prescription is under Examples/.
// The setup below is ZeissOtusML50mm.createSetup, written out.
#include "TestHarness.h"

#include "SetupAssertions.h"

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/spec/Prescription.h"

#include <fstream>
#include <iterator>
#include <string>

namespace {

using redukti::importers::OpticalBenchDataImporter;
using redukti::optim::OptimizationBuilder;
using redukti::optim::OptimizationTrial;
using redukti::spec::Prescription;
using redukti::test::assertSameSetup;

const char *const OTUS =
    REDUKTI_EXAMPLES_DIR "cosina-otus-ml-50mm-f1.4/JP2026-105585_Example01.txt";

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    CHECK(in.good());
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::string out;
    for (std::size_t i = 0; i < s.size(); i++) {
        if (s[i] == '\r' && i + 1 < s.size() && s[i + 1] == '\n')
            continue;
        out.push_back(s[i]);
    }
    return out;
}

Prescription prescription(const std::string &text, bool weighted, bool dLineOnly) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return Prescription::build_prescription(specs, true, weighted, dLineOnly);
}

/** ZeissOtusML50mm.createSetup. */
OptimizationBuilder otusBuilder(Prescription *prescription_, bool weighted, bool dLineOnly) {
    return OptimizationBuilder::builder(prescription_)
        .fields({0.0, 0.3, 0.7, 1.0})
        .mtfFrequencies({10, 20, 40})
        .varyCurvatures({0, 1, 2, 3, 4, 5, 6, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21,
                         22, 23})
        .varyThicknesses({25})
        .varyExistingAspherics()
        .weighted(weighted)
        .dLineOnly(dLineOnly)
        .rayAberrationGoals()
        .mtfGoals({OptimizationBuilder::mtf(10, {93, 93, 94, 93}, {93, 93, 90, 82}),
                   OptimizationBuilder::mtf(20, {85, 85, 85, 80}, {85, 85, 78, 62}),
                   OptimizationBuilder::mtf(40, {65, 65, 64, 58}, {65, 62, 45, 38})});
}

/** The trial as Documentation/OPTIMIZER.md prints it. */
const char *const TRIAL = R"([trial 1]
description       MTF targets, selected curvatures and the back focus
fields            0 0.3 0.7 1.0
frequencies       10 20 40

vary curvatures   all except 7 10 24
vary thicknesses  25
vary aspherics    existing

goal mtf   10 sag   93 93 94 93
goal mtf   10 tan   93 93 90 82
goal mtf   20 sag   85 85 85 80
goal mtf   20 tan   85 85 78 62
goal mtf   40 sag   65 65 64 58
goal mtf   40 tan   65 62 45 38
goal ray-aberrations  yes
)";

TEST(trialExamples_zeissOtusMtf) {
    std::string lens = readFile(OTUS);
    auto owned = prescription(lens, true, false);
    auto expectedBuilder = otusBuilder(&owned, true, false);
    auto expected = expectedBuilder.build();

    // The trial builds the example's setup, before and after a round trip through toTrial.
    auto trial = OptimizationTrial::read(lens + "\n" + TRIAL, 1, true);
    std::string written = trial.builder.toTrial(1);
    auto actual = trial.builder.build();
    assertSameSetup("otus trial", expected, actual);

    auto reread = OptimizationTrial::read(lens + "\n" + written, 1, true);
    CHECK_STR_EQ(reread.builder.toTrial(1), written);
    auto rereadSetup = reread.builder.build();
    assertSameSetup("otus trial, written and read back", expected, rereadSetup);
}

} // namespace
