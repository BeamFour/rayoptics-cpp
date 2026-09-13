// C++ port of org.redukti.tools.OptimizationTrialRunTest.
//
// Running a [trial n] from LensTool2, and the --optimize argument that asks for it. JUnit's
// @TempDir becomes a directory under the system temp directory, removed at the end of each
// test.
#include "TestHarness.h"

#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/spec/Prescription.h"
#include "redukti/tools/LensTool2.h"
#include "redukti/util/Args.h"

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

using redukti::importers::OpticalBenchDataImporter;
using redukti::optim::OptimizationTrial;
using redukti::optim::VarThickness;
using redukti::spec::Prescription;
using redukti::tools::LensTool2;
using redukti::util::Args;

const char *const OTUS =
    REDUKTI_EXAMPLES_DIR "cosina-otus-ml-50mm-f1.4/JP2026-105585_Example01.txt";
/** A zoom of two configurations, 0 and 1, whose surface 8 varies between them. */
const char *const ZOOM =
    REDUKTI_EXAMPLES_DIR "canon-rf70-200mm-f2.8LZ/US20250155694_Example01P.txt";

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

/** Java's String.replace: every occurrence, left to right. */
std::string replaceAll(std::string text, const std::string &from, const std::string &to) {
    std::size_t at = 0;
    while ((at = text.find(from, at)) != std::string::npos) {
        text.replace(at, from.size(), to);
        at += to.size();
    }
    return text;
}

void writeFile(const fs::path &path, const std::string &text) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << text;
}

/** JUnit's @TempDir: a directory of its own per test, removed when the test ends. */
class TempDir {
public:
    explicit TempDir(const char *name)
        : _path(fs::temp_directory_path() / ("rayoptics-" + std::string(name))) {
        fs::remove_all(_path);
        fs::create_directories(_path);
    }

    ~TempDir() {
        std::error_code ignored;
        fs::remove_all(_path, ignored);
    }

    const fs::path &path() const { return _path; }
    fs::path resolve(const std::string &child) const { return _path / child; }

private:
    fs::path _path;
};

const char *const TRIALS = R"(
[trial 3]
description   Back focus for contrast
fields        0
frequencies   20
vary thicknesses  25
goal contrast     20
goal contrast     sampling 3 6

[trial 4]
description   The same, reported separately
outdir        trials/back focus
fields        0
frequencies   20
vary thicknesses  25
goal contrast     20
goal contrast     sampling 3 6
)";

fs::path otus(const TempDir &dir) {
    fs::path spec = dir.resolve("otus.txt");
    writeFile(spec, readFile(OTUS) + TRIALS);
    return spec;
}

Prescription prescription(const std::string &text) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return Prescription::build_prescription(specs, true, false, false);
}

Args arguments(const std::vector<std::string> &args) {
    return Args::parseArguments(args);
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

bool sameFile(const fs::path &a, const fs::path &b) {
    return fs::absolute(a) == fs::absolute(b);
}

TEST(trialRun_optimizeTakesAnOptionalTrialNumber) {
    auto trial = arguments({"--specfile", "lens.txt", "--optimize", "2"});
    CHECK(trial.optimize_trial.has_value());
    CHECK_EQ(*trial.optimize_trial, 2);
    CHECK(!trial.optimize);

    auto routine = arguments({"--specfile", "lens.txt", "--optimize", "--mtf", "10,20"});
    CHECK(routine.optimize);
    CHECK(!routine.optimize_trial.has_value());
    CHECK(routine.mtf_freqs == (std::vector<int>{10, 20}));

    auto last = arguments({"--specfile", "lens.txt", "--optimize"});
    CHECK(last.optimize);

    try {
        arguments({"--optimize", "two"});
        CHECK(false);
    } catch (const redukti::IllegalArgumentException &e) {
        CHECK(contains(e.getMessage(), "[trial n]"));
    }
}

TEST(trialRun_writesToTheTrialsOwnDirectory) {
    TempDir dir("trialRun-own-directory");
    fs::path spec = otus(dir);
    Args args = arguments({"--specfile", spec.string(), "--optimize", "4"});
    LensTool2::runOptimizationTrial(readFile(spec.string()), args);
    fs::path output = dir.resolve("trials") / "back focus" / "otus-trial4.txt";
    CHECK(fs::exists(output));
    CHECK(sameFile(*args.specfile, output));
    CHECK(!args.outdir.has_value());
}

TEST(trialRun_commandLineOutdirTakesPrecedence) {
    TempDir dir("trialRun-outdir");
    fs::path spec = otus(dir);
    fs::path elsewhere = dir.resolve("elsewhere");
    Args args = arguments({"--specfile", spec.string(), "--outdir", elsewhere.string(),
                           "--optimize", "4"});
    LensTool2::runOptimizationTrial(readFile(spec.string()), args);
    CHECK(fs::exists(elsewhere / "otus-trial4.txt"));
    CHECK(!fs::exists(dir.resolve("trials")));
}

TEST(trialRun_writesThePrescriptionFollowedByTheTrial) {
    TempDir dir("trialRun-prescription");
    fs::path spec = otus(dir);
    std::string specText = readFile(spec.string());
    Args args = arguments({"--specfile", spec.string(), "--optimize", "3"});

    std::string optimized = LensTool2::runOptimizationTrial(specText, args);

    // Without an outdir the optimized prescription is written next to the input, and the
    // rest of the run is pointed at it.
    fs::path output = dir.resolve("otus-trial3.txt");
    CHECK(sameFile(*args.specfile, output));
    CHECK_STR_EQ(readFile(output.string()), optimized);

    // Beam42's own prescription, then the trial that ran as the builder writes it, and no
    // other trial.
    CHECK(optimized.rfind("[descriptive data]\n", 0) == 0);
    CHECK(contains(optimized, "Generated by Beam42"));
    std::string trailer = "\n" + OptimizationTrial::read(specText, 3, true).builder.toTrial(3);
    CHECK(optimized.size() >= trailer.size() &&
          optimized.compare(optimized.size() - trailer.size(), trailer.size(), trailer) == 0);
    CHECK(!contains(optimized, "[trial 4]"));

    // The back focus moved, and the trial carried over still names it: it can be run again
    // on the output as it stands.
    auto moved = prescription(optimized);
    CHECK(moved._surface_list[25]._thickness != 18.665);
    auto again = OptimizationTrial::read(optimized, 3, true);
    auto setup = again.builder.build();
    const auto *variable = dynamic_cast<VarThickness *>(setup.variables()[0].get());
    CHECK(variable != nullptr);
    CHECK_EQ(variable->_surface_id, 25);
    const_cast<VarThickness *>(variable)->read_from_prescription();
    CHECK_EQ(variable->get_unscaled_value(), moved._surface_list[25]._thickness);
}

/** One configuration at a time: the tele stage starts from the wide stage's result. */
const char *const PIPELINE = R"(
[trial 5]
description       Wide end
configuration     0
fields            0
frequencies       20
vary thicknesses  8
goal paraxial     bfl 40

[trial 6]
description       Tele end
configuration     1
fields            0
frequencies       20
vary thicknesses  8
goal paraxial     bfl 40

[pipeline 7]
description       Wide, then tele
outdir            trials/zoom
trials            5 6
)";

TEST(trialRun_runsAPipelineStageByStage) {
    TempDir dir("trialRun-pipeline");
    fs::path spec = dir.resolve("zoom.txt");
    writeFile(spec, readFile(ZOOM) + replaceAll(PIPELINE, "trials            5 6",
                                                "trials            5 6 5"));
    Args args = arguments({"--specfile", spec.string(), "--optimize", "7"});

    std::string optimized = LensTool2::runOptimizationTrial(readFile(spec.string()), args);

    // One result, in the pipeline's own directory, and the rest of the run uses it.
    fs::path output = dir.resolve("trials") / "zoom" / "zoom-pipeline7.txt";
    CHECK(fs::exists(output));
    CHECK(sameFile(*args.specfile, output));
    CHECK_STR_EQ(readFile(output.string()), optimized);

    // Each stage moved its own configuration's column, so both changed.
    auto moved = prescription(optimized);
    CHECK((*moved._surface_list[8]._thickness_by_scenario)[0] != 8.46);
    CHECK((*moved._surface_list[8]._thickness_by_scenario)[1] != 13.14);

    // The result carries the pipeline and both its trials, so it can run again as it stands.
    auto again = OptimizationTrial::readPipeline(optimized, 7);
    CHECK(again.has_value());
    CHECK(again->trials() == (std::vector<int>{5, 6, 5}));
    CHECK_STR_EQ(OptimizationTrial::readPipeline(optimized + "\n", 7)->toPipeline(),
                 again->toPipeline());
    CHECK_STR_EQ(*again->outdir(), "trials/zoom");
    for (int stage : again->trials()) {
        auto trial = OptimizationTrial::read(optimized, stage, true);
        auto setup = trial.builder.build();
        const auto *variable = dynamic_cast<VarThickness *>(setup.variables()[0].get());
        CHECK(variable != nullptr);
        CHECK_EQ(variable->_surface_id, 8);
    }
    for (int stage : again->distinctTrials()) {
        auto original = OptimizationTrial::parse(readFile(spec.string()), stage);
        auto restored = OptimizationTrial::parse(optimized, stage);
        CHECK_STR_EQ(restored.toTrial(&moved), original.toTrial(&moved));
    }
}

} // namespace
