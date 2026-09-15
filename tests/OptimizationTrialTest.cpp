// C++ port of org.redukti.optim.OptimizationTrialTest.
//
// The Java holds the builder alone and lets the collector keep the prescription alive;
// OptimizationTrial::read returns the two together, so each test here keeps the Trial in a
// named local that outlives the setup built from it.
#include "TestHarness.h"

#include "SetupAssertions.h"

#include "redukti/Text.h"
#include "redukti/importers/OpticalBenchDataImporter.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/optim/Analysis.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/spec/Prescription.h"
#include "redukti/spec/SurfaceType.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <iterator>
#include <string>
#include <vector>

namespace {

using redukti::doubleToString;
using redukti::IllegalStateException;
using redukti::importers::OpticalBenchDataImporter;
using redukti::optim::GoalParax;
using redukti::optim::OptimizationBuilder;
using redukti::optim::OptimizationTrial;
using redukti::optim::ParaxHelper;
using redukti::optim::Var;
using redukti::optim::VarAsphCoeff;
using redukti::optim::VarAsphK;
using redukti::optim::VarThickness;
using redukti::spec::Prescription;
using redukti::spec::SurfaceType;
using redukti::test::assertSameSetup;
using redukti::IllegalArgumentException;
using redukti::rayoptics::analysis::ContrastOptions;
using redukti::rayoptics::analysis::SpotOptions;
using TrialException = OptimizationTrial::TrialException;

/** Eleven surfaces, 0 to 10: the stop is 5, surfaces 3, 7 and 9 are flat, and 10 the last. */
const char *const SUMMICRON =
    REDUKTI_EXAMPLES_DIR "leica-summicron-50mm-f2/US004123144_Example08P.txt";
/** Ten surfaces whose ids run 1 to 6, then the stop 6AS, then 7 to 9. */
const char *const FD300 =
    REDUKTI_EXAMPLES_DIR "canon-fd300mm-f2.8-fluorite/US003868174_Example01P.txt";
/** A zoom of three scenarios, configured as 0 and 2. */
const char *const ZOOM =
    REDUKTI_EXAMPLES_DIR "canon-rf70-200mm-f2.8LZ/US20250155694_Example01P.txt";

std::string readFile(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    CHECK(in.good());
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    // The committed files use bare newlines; normalise in case of checkout translation, so
    // the line numbers a trial reports are the ones the test counts.
    std::string out;
    for (std::size_t i = 0; i < s.size(); i++) {
        if (s[i] == '\r' && i + 1 < s.size() && s[i + 1] == '\n')
            continue;
        out.push_back(s[i]);
    }
    return out;
}

std::string withTrials(const char *example, const std::string &trials) {
    return readFile(example) + "\n" + trials;
}

OptimizationTrial::Trial read(const std::string &text) {
    return OptimizationTrial::read(text, 1, true);
}

Prescription prescription(const std::string &text, bool weighted, bool dLineOnly) {
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    return Prescription::build_prescription(specs, true, weighted, dLineOnly);
}

std::vector<std::string> names(const std::vector<std::shared_ptr<Var>> &variables) {
    std::vector<std::string> result;
    for (const auto &variable : variables)
        result.push_back(OptimizationTrial::describe(*variable));
    return result;
}

void checkNames(const std::vector<std::shared_ptr<Var>> &variables,
                const std::vector<std::string> &expected) {
    auto actual = names(variables);
    if (actual == expected)
        return;
    std::string got, want;
    for (const auto &name : actual)
        got += "[" + name + "]";
    for (const auto &name : expected)
        want += "[" + name + "]";
    ::redukti::test::reportFailure(__FILE__, __LINE__, "got " + got + " want " + want);
}

/** The 1-based number of the line holding exactly `wanted`, as the Java's indexOf + 1. */
int lineOf(const std::string &text, const std::string &wanted) {
    int line = 1;
    std::size_t start = 0;
    while (start <= text.size()) {
        std::size_t end = text.find('\n', start);
        std::string current =
            text.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (current == wanted)
            return line;
        if (end == std::string::npos)
            break;
        start = end + 1;
        line++;
    }
    return 0;
}

bool contains(const std::string &text, const std::string &part) {
    return text.find(part) != std::string::npos;
}

/** CHECK that `message` mentions `part`, reporting the message itself when it does not. */
void checkMentions(const std::string &message, const std::string &part) {
    if (!contains(message, part))
        ::redukti::test::reportFailure(__FILE__, __LINE__,
                                       "message \"" + message + "\" does not mention \"" +
                                           part + "\"");
}

/** Java's replaceAll("\\s+", " "): runs of whitespace collapsed to one space. */
std::string normalizeSpaces(const std::string &text) {
    std::string out;
    bool space = false;
    for (char ch : text) {
        bool isSpace = ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r' || ch == '\f' ||
                       ch == '\v';
        if (isSpace) {
            if (!space)
                out.push_back(' ');
        } else
            out.push_back(ch);
        space = isSpace;
    }
    return out;
}

/** Java's text.split("\n", -1).length: the lines, counting a trailing empty one. */
int lineCount(const std::string &text) {
    return static_cast<int>(std::count(text.begin(), text.end(), '\n')) + 1;
}

/** The message of the TrialException `call` throws; a failure when it throws nothing. */
template <typename Call> std::string trialErrorOf(Call call) {
    try {
        call();
    } catch (const TrialException &e) {
        return e.getMessage();
    }
    ::redukti::test::reportFailure(__FILE__, __LINE__, "expected a TrialException");
    return "";
}

const char *const MTF_FIVE_FIELDS = R"(fields        0 to 1 step 0.25
frequencies   20
goal mtf      20 sag   50 50 50 50 50
goal mtf      20 tan   50 50 50 50 50
)";

std::string trialOf(const std::string &settings) {
    return std::string("[trial 1]\n") + MTF_FIVE_FIELDS + settings;
}

TEST(trial_reusesParsedDefinitionAndRoundTripsItsSettings) {
    std::string text =
        withTrials(SUMMICRON, trialOf("vary thicknesses 10\nweighted no\nd-line-only yes\n"));
    auto definition = OptimizationTrial::parse(text, 1);
    auto first = definition.createBuilder(text, true);
    std::string canonical = definition.toTrial(first.prescription.get());
    CHECK_STR_EQ(definition.toTrial(), canonical);
    CHECK_STR_EQ(first.builder.toTrial(1), canonical);
    CHECK_EQ(definition.number(), 1);
    CHECK(contains(canonical, "weighted"));
    {
        auto expectedTrial = read(text);
        auto expected = expectedTrial.builder.build();
        auto actual = first.builder.build();
        assertSameSetup("first stage", expected, actual);
    }

    // A stage hands on its changed design, without needing to carry trial text for the
    // already-parsed definition to construct a fresh stage.
    first.prescription->_surface_list[10]._thickness += 0.25;
    std::string design;
    first.prescription->to_opt_bench_str(design);
    auto next = definition.createBuilder(design, true);
    CHECK(first.prescription.get() != next.prescription.get());
    CHECK_CLOSE(next.prescription->_surface_list[10]._thickness,
                first.prescription->_surface_list[10]._thickness, 1e-8);
    CHECK_STR_EQ(definition.toTrial(next.prescription.get()), canonical);
    CHECK(first.prescription->_wvls == next.prescription->_wvls);

    auto reread = OptimizationTrial::parse(design + "\n" + canonical, 1);
    CHECK_STR_EQ(reread.toTrial(next.prescription.get()), canonical);
    {
        auto expected = next.builder.build();
        auto rebuilt = reread.createBuilder(design, true);
        auto actual = rebuilt.builder.build();
        assertSameSetup("next stage", expected, actual);
    }
    // Changing one returned builder must not change the reusable definition.
    next.builder.fields({0.0});
    next.builder.contrastGoals({OptimizationBuilder::contrast(30, std::vector<double>{1.0})});
    next.builder.varyCurvatures({0});
    next.builder.contrastBalanceGoals(std::vector<bool>{true});
    CHECK_STR_EQ(definition.toTrial(first.prescription.get()), canonical);
    CHECK_STR_EQ(definition.toTrial(), canonical);
    CHECK_STR_EQ(definition.createBuilder(design, true).builder.toTrial(1), canonical);
}

TEST(trial_effectiveAnalysesAndSamplingSurviveRoundTrip) {
    struct Case {
        const char *goals;
        bool spots, rays, mtf, hexapolar;
    };
    const Case cases[] = {
        {"", false, false, false, false},
        {"goal spot-rms 10\n", true, false, false, false},
        {"goal spot-max-radius 10\n", true, false, false, true},
        {"goal spot-deviation 1\n", true, false, false, false},
        {"goal mtf 20 sag 50\ngoal mtf 20 tan 50\n", true, false, true, false},
        {"goal ray-aberrations yes\n", false, true, false, false},
        {"goal contrast 20\n", false, false, false, false},
        {"goal spot sampling hexapolar 32\n", false, false, false, true}};
    for (const Case &c : cases) {
        std::string text =
            withTrials(SUMMICRON, std::string("[trial 1]\nfields 0\nfrequencies 20\n") + c.goals);
        auto trial = read(text);
        auto setup = trial.builder.build();
        auto *analysis = setup.analysis();
        CHECK_EQ(analysis->_compute_spots, c.spots);
        CHECK_EQ(analysis->_compute_ray_aberrations, c.rays);
        CHECK_EQ(analysis->_compute_mtf, c.mtf);
        CHECK_EQ(analysis->_spot_pattern, c.hexapolar ? SpotOptions::PATTERN_HEXAPOLAR
                                                      : SpotOptions::PATTERN_GAUSS_QUADRATURE);
        std::string written = trial.builder.toTrial(1);
        CHECK_STR_EQ(OptimizationTrial::parse(text, 1).toTrial(), written);
        std::string normalized = normalizeSpaces(written);
        CHECK_EQ(contains(normalized, "sampling hexapolar"), c.hexapolar);
        CHECK_EQ(contains(normalized, "sampling gaussian"), c.spots && !c.hexapolar);
        auto restored = read(withTrials(SUMMICRON, written));
        CHECK_STR_EQ(restored.builder.toTrial(1), written);
        auto restoredSetup = restored.builder.build();
        assertSameSetup(c.goals, setup, restoredSetup);
    }
}

TEST(trial_solverTolerancesRoundTripAndAreAbsentUntilAskedFor) {
    std::string settings = R"(solver ftol            1.0E-6
solver xtol            1.0E-5
solver gtol            0
solver max-evaluations 250
)";
    auto trial = read(withTrials(SUMMICRON, trialOf(settings)));
    std::string written = trial.builder.toTrial(1);
    std::string normalized;
    for (char c : written) {
        if (c == ' ' || c == '\t' || c == '\n') {
            if (!normalized.empty() && normalized.back() != ' ')
                normalized.push_back(' ');
        } else
            normalized.push_back(c);
    }
    CHECK(contains(normalized, "solver ftol 1.0E-6"));
    CHECK(contains(normalized, "solver xtol 1.0E-5"));
    CHECK(contains(normalized, "solver gtol 0"));
    CHECK(contains(normalized, "solver max-evaluations 250"));
    auto restored = read(withTrials(SUMMICRON, written));
    CHECK_STR_EQ(restored.builder.toTrial(1), written);

    // A trial that says nothing about the solver writes nothing back, so the solver
    // keeps the defaults it has always used.
    auto silentTrial = read(withTrials(SUMMICRON, trialOf("")));
    CHECK(!contains(silentTrial.builder.toTrial(1), "solver "));

    try {
        read(withTrials(SUMMICRON, trialOf("solver wibble 1\n")));
        CHECK(false);
    } catch (const TrialException &e) {
        checkMentions(e.getMessage(), "unknown solver setting 'wibble'");
    }
}

TEST(trial_numbersSurfacesByPosition) {
    std::string text = withTrials(SUMMICRON, trialOf(R"(vary curvatures   all except 2 6
vary thicknesses  10 4
)"));
    auto trial = read(text);
    auto setup = trial.builder.build();
    // Surfaces 3, 7 and 9 are flat and 5 is the stop, so "all" leaves them out without
    // being told.
    checkNames(setup.variables(),
               {"surface 0 radius", "surface 1 radius", "surface 4 radius",
                "surface 8 radius", "surface 10 radius", "surface 10 thickness",
                "surface 4 thickness"});
}

TEST(trial_rejectsFractionalCurvatureConstraintsOnFlatSurfaces) {
    for (const char *settings : {"vary curvatures 3\nconstrain curvatures\n",
                                 "constrain curvatures\nvary curvatures 3\n"}) {
        std::string text = withTrials(SUMMICRON, trialOf(settings));
        int constraintLine = lineOf(text, "constrain curvatures");
        try {
            read(text);
            CHECK(false);
        } catch (const TrialException &e) {
            checkMentions(e.getMessage(), "line " + redukti::intToString(constraintLine) + ":");
            checkMentions(e.getMessage(), "flat surface 3");
        }
    }
    // Listing a flat surface without a fractional constraint remains supported.
    auto trial = read(withTrials(SUMMICRON, trialOf("vary curvatures 3\n")));
    auto setup = trial.builder.build();
    checkNames(setup.variables(), {"surface 3 radius"});
}

/** The Java's assertRejected: the message of the TrialException the trial raises. */
std::string rejection(const char *example, const std::string &trialText) {
    std::string text = withTrials(example, trialText);
    try {
        auto trial = read(text);
        trial.builder.build();
    } catch (const TrialException &e) {
        return e.getMessage();
    }
    ::redukti::test::reportFailure(__FILE__, __LINE__, "expected a TrialException");
    return "";
}

std::string rejection(const std::string &trialText) {
    return rejection(SUMMICRON, trialText);
}

TEST(trial_ignoresTheIdsInTheFile) {
    // Position 6 is the stop, although the row with id 6 is an ordinary surface.
    checkMentions(rejection(FD300, trialOf("vary curvatures 6\n")), "surface 6 is a stop");
    checkMentions(rejection(FD300, trialOf("vary thicknesses 6AS\n")),
                  "expected a surface number, found '6AS'");

    auto trial = read(withTrials(FD300, trialOf(R"(vary curvatures   5 7
vary thicknesses  6 9
)")));
    auto setup = trial.builder.build();
    checkNames(setup.variables(), {"surface 5 radius", "surface 7 radius",
                                   "surface 6 thickness", "surface 9 thickness"});
    // Surface 5 is the row with id 6, and surface 6 the stop, whose gap is 30.10.
    CHECK_EQ(trial.prescription->_surface_list[5]._radius, -524.3616);
    CHECK_EQ(trial.prescription->_surface_list[6]._thickness, 30.10);
}

TEST(trial_fieldShorthandGivesTheDecimalValues) {
    auto trial = read(withTrials(SUMMICRON, R"([trial 1]
fields        0 to 1 step 0.1
frequencies   20
vary thicknesses 10
)"));
    auto setup = trial.builder.build();
    std::vector<double> expected{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    CHECK_EQ(setup.analysis()->_fields.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); i++)
        CHECK_EQ(setup.analysis()->_fields[i], expected[i]);
}

TEST(trial_focalLengthGoalReplacesTheAutomaticOne) {
    auto trial = read(withTrials(SUMMICRON, R"([trial 1]
fields        0
frequencies   20
vary thicknesses 10
goal paraxial efl 42 weight 2
goal paraxial bfl 30
)"));
    auto setup = trial.builder.build();
    auto goals = setup.goals();
    CHECK_EQ(goals.size(), static_cast<std::size_t>(3));
    std::vector<const GoalParax *> paraxial;
    for (const auto &goal : goals)
        paraxial.push_back(dynamic_cast<const GoalParax *>(goal.get()));
    for (const auto *goal : paraxial)
        CHECK(goal != nullptr);
    CHECK_EQ(paraxial[0]->_parax_id, ParaxHelper::Effective_focal_length);
    CHECK_EQ(paraxial[0]->_target, 42.0);
    CHECK_EQ(paraxial[0]->_weight, 2.0);
    CHECK_EQ(paraxial[1]->_parax_id, ParaxHelper::Fno);
    CHECK_EQ(paraxial[1]->_target, 2.0);
    CHECK_EQ(paraxial[1]->_weight, 1.0);
    CHECK_EQ(paraxial[2]->_parax_id, ParaxHelper::Back_focal_length);
    CHECK_EQ(paraxial[2]->_target, 30.0);
    CHECK_EQ(paraxial[2]->_weight, 1.0);
}

TEST(trial_makesASphereAsphericWithScaledCoefficients) {
    auto trial = read(withTrials(SUMMICRON, trialOf("vary aspherics 0 K 1 2:1e9\n")));
    auto setup = trial.builder.build();
    auto variables = setup.variables();
    const SurfaceType &surface = trial.prescription->_surface_list[0];
    CHECK_EQ(surface._asph_type, SurfaceType::ASPH_EVEN);
    CHECK_EQ(surface._coeffs->size(), static_cast<std::size_t>(3));
    checkNames(variables,
               {"surface 0 K", "surface 0 coefficient 1", "surface 0 coefficient 2"});
    const auto *a4 = dynamic_cast<const VarAsphCoeff *>(variables[1].get());
    const auto *a6 = dynamic_cast<const VarAsphCoeff *>(variables[2].get());
    CHECK(a4 != nullptr && a6 != nullptr);
    // Coefficient 1 of an even asphere is A4. Half the 28.94 diameter, to the fourth
    // power, is 43800: nearest decade 1e5.
    CHECK_EQ(a4->_index, 1);
    CHECK_EQ(a4->_scaling_factor, 1e5);
    CHECK_EQ(a6->_index, 2);
    CHECK_EQ(a6->_scaling_factor, 1e9);
}

/** A trial using every setting, in the form the builder writes it back. */
const char *const EVERYTHING = R"([trial 1]
description           Everything a trial can say
outdir                trials/one
configuration         0
fields                0 0.25 0.5 0.75 1
frequencies           10 20
weighted              no
d-line-only           yes
vignetting            set-vig frozen
check-spot-apertures  no
vary curvatures       all except 2 6
vary thicknesses      0 10
vary aspherics        existing
vary aspherics        0 K 1:100000 2
constrain curvatures  2
constrain thicknesses 1
constrain edges       0.5
goal contrast         10 20
goal contrast         sag 3 3 2 2 1
goal contrast         20 tan 1 1 1 0.5 0.5
goal contrast         balance all except 0.25 weight 0.4
goal contrast         sampling 4 8
goal contrast         calibrate yes
goal contrast         exit-pupil-aiming no
goal contrast         centering no
goal mtf              10 sag 90 85 80 70 60
goal mtf              10 tan 90 85 80 70 60
goal mtf              10 tan weights 1 1 2 2 4
goal spot-rms         10 20 30 40 50
goal spot sampling    gaussian 6 12 0.2
goal spot sampling    hexapolar 32
goal ray-aberrations  yes
goal paraxial         efl 51 weight 2
goal paraxial         bfl 37
)";

TEST(trial_readsAndWritesTheSameTrial) {
    // Written differently - shorthand, other spacing, other order, defaults spelt out -
    // and read, the trial is written back in the builder's form. Goals of one kind keep
    // the order they were given in, since that is the order of the residuals.
    const char *const loose = R"([trial 1]
goal paraxial efl 51 weight 2
goal paraxial bfl 37
description   Everything a trial can say   # a comment is not kept
outdir  trials/one
fields 0 to 1 step 0.25
frequencies 10 20
weighted no
d-line-only yes
vignetting set-vig frozen
check-spot-apertures no
vary curvatures all except 2 6
vary thicknesses 0 10
vary aspherics existing
vary aspherics 0 K 1:1e5 2
constrain curvatures 2
constrain thicknesses
constrain edges 0.5
goal contrast 10 20
goal contrast sag 3 3 2 2 1
goal contrast 20 tan 1 1 1 0.5 0.5
goal contrast balance all except 0.25 weight 0.4
goal contrast sampling 4 8
goal contrast calibrate yes
goal contrast centering no
goal mtf 10 sag 90 85 80 70 60
goal mtf 10 tan 90 85 80 70 60
goal mtf 10 tan weights 1 1 2 2 4
goal spot-rms 10 20 30 40 50
goal spot sampling gaussian 6 12 0.2
goal spot sampling hexapolar 32
goal ray-aberrations yes
)";
    auto trial = read(withTrials(SUMMICRON, loose));
    CHECK_STR_EQ(trial.builder.toTrial(1), EVERYTHING);
    // Parsing and canonical writing need no prescription or solver construction.
    auto definition = OptimizationTrial::parse(withTrials(SUMMICRON, loose), 1);
    CHECK_STR_EQ(definition.toTrial(), EVERYTHING);
    CHECK_STR_EQ(
        OptimizationTrial::parse(withTrials(SUMMICRON, definition.toTrial()), 1).toTrial(),
        EVERYTHING);

    // Reading what was written gives the same setup, and writes the same text again.
    auto again = read(withTrials(SUMMICRON, EVERYTHING));
    CHECK_STR_EQ(again.builder.toTrial(1), EVERYTHING);
    auto expected = trial.builder.build();
    auto actual = again.builder.build();
    assertSameSetup("everything", expected, actual);
}

TEST(trial_writesASetupBuiltInCode) {
    std::string lens = withTrials(SUMMICRON, "");
    auto owned = prescription(lens, true, false);
    auto builder = OptimizationBuilder::builder(&owned)
                       .fields({0.0, 0.5, 1.0})
                       .mtfFrequencies({10})
                       .varyThicknesses({10, 4})
                       .varyAsphericCoefficient(0, 1, 1e5)
                       .spotDeviationGoals({1, 2, 3}, {1, 1, 1})
                       .gaussianQuadratureSampling(3, 6)
                       .paraxialGoal(ParaxHelper::Back_focal_length, 37.3);
    std::string written = builder.toTrial(7);
    CHECK_STR_EQ(written, R"([trial 7]
configuration         0
fields                0 0.5 1
frequencies           10
weighted              yes
d-line-only           no
vignetting            set-pupil
check-spot-apertures  yes
vary thicknesses      10 4
vary aspherics        0 1:100000
goal spot-deviation   x 1 2 3
goal spot-deviation   y 1 1 1
goal spot sampling    gaussian 3 6
goal ray-aberrations  no
goal paraxial         bfl 37.3
)");
    auto reread = OptimizationTrial::read(lens + written, 7, true);
    CHECK_STR_EQ(reread.builder.toTrial(7), written);
    auto expected = builder.build();
    auto actual = reread.builder.build();
    assertSameSetup("code-built", expected, actual);
}

TEST(trial_cannotWriteVariablesOrGoalsGivenAsCode) {
    auto owned = prescription(withTrials(SUMMICRON, ""), true, false);
    auto builder = OptimizationBuilder::builder(&owned)
                       .fields({0.0})
                       .mtfFrequencies({10})
                       .additionalGoals({[](redukti::optim::Analysis *analysis) {
                           return std::shared_ptr<redukti::optim::Goal>(
                               std::make_shared<GoalParax>(
                                   analysis, ParaxHelper::Back_focal_length, 37.3, 1.0));
                       }});
    try {
        builder.toTrial(1);
        CHECK(false);
    } catch (const IllegalStateException &e) {
        checkMentions(e.getMessage(), "added as code");
    }
}

/** Moves every variable, as a solve would, and writes the result to the prescription. */
void move(const std::vector<std::shared_ptr<Var>> &variables) {
    for (const auto &variable : variables) {
        variable->read_from_prescription();
        double value = dynamic_cast<VarAsphCoeff *>(variable.get()) != nullptr ? 1.25e-5
                       : dynamic_cast<VarAsphK *>(variable.get()) != nullptr
                           ? -0.75
                           : variable->get_unscaled_value() * 1.01;
        variable->set_unscaled_value(value);
        variable->write_to_prescription();
    }
}

/** What LensTool2 saves: the prescription as Beam42 writes it, then the trial. */
std::string optimized(OptimizationBuilder &builder) {
    std::string sb;
    builder.prescription()->to_opt_bench_str(sb);
    return sb + "\n" + builder.toTrial(1);
}

TEST(trial_thePrescriptionAndTrialReadBackAsOptimized) {
    auto trial = read(withTrials(SUMMICRON, trialOf(R"(vary curvatures   0 2
vary thicknesses  0 10
vary aspherics    0 K 1
)")));
    auto setup = trial.builder.build();
    move(setup.variables());
    std::string written = optimized(trial.builder);

    // Reads back as the optimized prescription, exactly, and the trial carried over
    // refers to the same surfaces.
    auto again = read(written);
    const auto &surfaces = trial.prescription->_surface_list;
    const auto &reread = again.prescription->_surface_list;
    CHECK_EQ(reread.size(), surfaces.size());
    for (std::size_t i = 0; i < surfaces.size() && i < reread.size(); i++) {
        CHECK_EQ(reread[i]._radius, surfaces[i]._radius);
        CHECK_EQ(reread[i]._thickness, surfaces[i]._thickness);
        CHECK_EQ(reread[i]._k, surfaces[i]._k);
        CHECK_EQ(reread[i]._asph_type, surfaces[i]._asph_type);
        CHECK_EQ(reread[i]._coeffs.has_value(), surfaces[i]._coeffs.has_value());
        if (reread[i]._coeffs.has_value() && surfaces[i]._coeffs.has_value())
            CHECK(*reread[i]._coeffs == *surfaces[i]._coeffs);
    }
    auto rereadSetup = again.builder.build();
    checkNames(rereadSetup.variables(),
               {"surface 0 radius", "surface 2 radius", "surface 0 thickness",
                "surface 10 thickness", "surface 0 K", "surface 0 coefficient 1"});
}

TEST(trial_writesAZoomWithItsConfiguredScenariosOnly) {
    auto trial = read(withTrials(ZOOM, R"([trial 1]
configuration     1
fields            0
frequencies       10
vary thicknesses  8 10
goal mtf          10 sag 50
goal mtf          10 tan 50
)"));
    auto setup = trial.builder.build();
    for (const auto &variable : setup.variables()) {
        variable->read_from_prescription();
        variable->set_unscaled_value(variable->get_unscaled_value() + 0.5);
        variable->write_to_prescription();
    }
    std::string written = optimized(trial.builder);

    // The 120mm scenario is not configured, so it is not written; the configurations are
    // renumbered from 0 in the same order.
    CHECK(!contains(written, "120.07"));
    CHECK(contains(written, "scenarios\t0\t1\n"));
    auto again = read(written);
    const auto &reread = again.prescription->_surface_list;
    CHECK_CLOSE((*reread[8]._thickness_by_scenario)[1], 13.64, 1e-12);
    CHECK_EQ((*reread[8]._thickness_by_scenario)[0], 8.46);
    CHECK_EQ((*reread[10]._thickness_by_scenario)[0], 20.31);
    // The trial carried over still means configuration 1 and the same surfaces.
    auto rereadSetup = again.builder.build();
    const auto *moved = dynamic_cast<const VarThickness *>(rereadSetup.variables()[0].get());
    CHECK(moved != nullptr);
    CHECK_EQ(moved->_scenario, 1);
    CHECK_EQ(moved->_surface_id, 8);
}

TEST(trial_reportsProblemsWithTheirLine) {
    std::string text = withTrials(SUMMICRON, "[trial 1]\nfields 0\nbogus 1\n");
    int line = lineOf(text, "bogus 1");
    try {
        read(text);
        CHECK(false);
    } catch (const TrialException &e) {
        CHECK_STR_EQ(e.getMessage(),
                     "trial 1, line " + redukti::intToString(line) +
                         ": unknown keyword 'bogus'");
    }

    // A problem only the builder can see, against the prescription, still names its line.
    std::string aspheric = withTrials(SUMMICRON, trialOf("vary aspherics 0 0\n"));
    int asphericLine = lineOf(aspheric, "vary aspherics 0 0");
    try {
        read(aspheric);
        CHECK(false);
    } catch (const TrialException &e) {
        CHECK_STR_EQ(e.getMessage(),
                     "trial 1, line " + redukti::intToString(asphericLine) +
                         ": coefficient 0 is not a term of an even asphere, whose terms "
                         "start at index 1, the A4 term");
    }
}

TEST(trial_sharedDomainRulesRejectBothEntryPointsWithSourceContext) {
    struct Invalid {
        std::string rows;
        std::string offendingRow;
        std::function<void(OptimizationBuilder &)> configure;
    };
    const std::vector<Invalid> invalid{
        {"goal spot-rms -1\n", "goal spot-rms -1",
         [](OptimizationBuilder &b) { b.spotRmsGoals(std::vector<double>{-1.0}); }},
        {"goal spot-deviation -1\n", "goal spot-deviation -1",
         [](OptimizationBuilder &b) { b.spotDeviationGoals(std::vector<double>{-1.0}); }},
        {"goal spot-rms 1 2\n", "goal spot-rms 1 2",
         [](OptimizationBuilder &b) { b.spotRmsGoals(std::vector<double>{1.0, 2.0}); }},
        {"goal mtf 20 sag 101\ngoal mtf 20 tan 50\n", "goal mtf 20 sag 101",
         [](OptimizationBuilder &b) {
             b.mtfGoals({OptimizationBuilder::mtf(20, std::vector<double>{101.0},
                                                  std::vector<double>{50.0})});
         }},
        {"goal mtf 30 sag 50\ngoal mtf 30 tan 50\n", "goal mtf 30 sag 50",
         [](OptimizationBuilder &b) {
             b.mtfGoals({OptimizationBuilder::mtf(30, std::vector<double>{50.0},
                                                  std::vector<double>{50.0})});
         }},
        {"goal contrast 20 20\n", "goal contrast 20 20", [](OptimizationBuilder &b) {
             b.contrastGoals({OptimizationBuilder::contrast(20, std::vector<double>{1.0}),
                              OptimizationBuilder::contrast(20, std::vector<double>{1.0})});
         }}};
    std::string lens = withTrials(SUMMICRON, "");
    for (const auto &c : invalid) {
        auto owned = prescription(lens, true, false);
        auto builder = OptimizationBuilder::builder(&owned).fields({0.0}).mtfFrequencies({20});
        c.configure(builder);
        CHECK_THROWS(builder.build(), IllegalArgumentException);
        std::string text = lens + "\n[trial 1]\nfields 0\nfrequencies 20\n" + c.rows;
        std::string message = trialErrorOf([&] { OptimizationTrial::parse(text, 1); });
        int line = lineOf(text, c.offendingRow);
        if (message.rfind("trial 1, line " + redukti::intToString(line) + ":", 0) != 0)
            ::redukti::test::reportFailure(__FILE__, __LINE__,
                                           c.offendingRow + " reported as \"" + message + "\"");
    }
}

TEST(trial_writingAnInvalidBuilderPreservesItsHexapolarSetting) {
    std::string lens = withTrials(SUMMICRON, "");
    auto owned = prescription(lens, true, false);
    auto builder = OptimizationBuilder::builder(&owned)
                       .fields({0.0})
                       .mtfFrequencies({20})
                       .spotDeviationGoals(std::vector<double>{1.0})
                       .hexapolarSampling(32);
    CHECK_THROWS(builder.build(), IllegalArgumentException);
    std::string written = builder.toTrial(1);
    CHECK(contains(normalizeSpaces(written), "goal spot sampling hexapolar 32"));
    CHECK_THROWS(OptimizationTrial::parse(lens + "\n" + written, 1), TrialException);
}

TEST(trial_incompatibleSpotSamplingReportsTheConflictingLineInEitherOrder) {
    for (const char *deviation :
         {"goal spot-deviation 1", "goal spot-deviation x 1\ngoal spot-deviation y 1"}) {
        for (const char *conflict : {"goal spot sampling hexapolar 32", "goal spot-max-radius 10"}) {
            for (bool reverse : {false, true}) {
                std::string rows = reverse ? std::string(conflict) + "\n" + deviation
                                           : std::string(deviation) + "\n" + conflict;
                std::string text =
                    withTrials(SUMMICRON, "[trial 1]\nfields 0\nfrequencies 20\n" + rows);
                std::string message = trialErrorOf([&] { OptimizationTrial::parse(text, 1); });
                CHECK_STR_EQ(message, "trial 1, line " + redukti::intToString(lineCount(text)) +
                                          ": spot deviation goals require Gaussian-quadrature "
                                          "spot sampling");
            }
        }
    }
}

TEST(trial_trialSamplingMatchesTheAvailableTracerPatterns) {
    std::string prefix = "[trial 1]\nfields 0\nfrequencies 20\n";
    std::string lens = withTrials(SUMMICRON, "");
    for (int spokes : {1, 2}) {
        std::string text = lens + prefix + "goal contrast 20\ngoal contrast sampling 1 " +
                           redukti::intToString(spokes);
        std::string message = trialErrorOf([&] { OptimizationTrial::parse(text, 1); });
        checkMentions(message, "line " + redukti::intToString(lineCount(text)) + ":");
        auto owned = prescription(lens, true, false);
        CHECK_THROWS(OptimizationBuilder::builder(&owned).contrastSampling(1, spokes),
                     IllegalArgumentException);
        CHECK_THROWS(ContrastOptions(20).num_spokes(spokes), IllegalArgumentException);
    }
    CHECK_THROWS(OptimizationTrial::parse(lens + prefix + "goal spot sampling grid 9\n", 1),
                 TrialException);
    // Independent spot and contrast settings are valid together, at the lower bound.
    std::string text = lens + prefix +
                       "goal spot-rms 10\ngoal spot sampling hexapolar 2\n"
                       "goal contrast 20\ngoal contrast sampling 1 3\n";
    auto trial = read(text);
    auto setup = trial.builder.build();
    setup.analysis()->compute();
    CHECK(setup.analysis()->_spots.has_value());
    CHECK_EQ(setup.analysis()->_contrasts.at(0).fields.at(0).wavelengths.at(0).samples.size(),
             static_cast<std::size_t>(3));
    CHECK_EQ(setup.analysis()->_spot_pattern, SpotOptions::PATTERN_HEXAPOLAR);
    auto restored = read(lens + trial.builder.toTrial(1));
    auto expected = trial.builder.build();
    auto actual = restored.builder.build();
    assertSameSetup("sampling", expected, actual);
    CHECK_STR_EQ(trial.builder.toTrial(1), restored.builder.toTrial(1));
}

TEST(trial_rejectsMistakes) {
    checkMentions(rejection("[trial 1]\nfields 0\nfields 0\n"),
                  "'fields' is given more than once");
    checkMentions(rejection("[trial 1]\nfields 0\n"), "'frequencies' is required");
    checkMentions(rejection("[trial 1]\nfields 0 to 1 step 0.25\nfrequencies 20\n"
                            "goal mtf 20 sag 50 50 50\ngoal mtf 20 tan 50 50 50 50 50\n"),
                  "expected 5 targets, one per field, but found 3");
    checkMentions(rejection(trialOf("vary curvatures 11\n")),
                  "there is no surface 11; [lens data] has surfaces 0 to 10");
    checkMentions(rejection(trialOf("vary curvatures 5\n")), "surface 5 is a stop");
    checkMentions(rejection(trialOf("vary thicknesses Bf\n")),
                  "expected a surface number, found 'Bf'");
    checkMentions(rejection(trialOf("vary thicknesses 4 04\n")), "surface 4 is listed twice");
    checkMentions(rejection(trialOf("vary aspherics 0 K:10\n")), "K takes no scale");
    checkMentions(rejection(trialOf("vary aspherics 0 A4\n")),
                  "expected K or a coefficient index, found 'A4'");
    checkMentions(rejection(trialOf("vary aspherics 5 K\n")),
                  "surface 5 is a stop; it cannot be aspheric");
    checkMentions(rejection("[trial 1]\nfields 0 0.5\nfrequencies 20\ngoal contrast sag 1 1\n"),
                  "contrast settings need a 'goal contrast <frequencies>' line");
    checkMentions(rejection("[trial 1]\nfields 0 0.5\nfrequencies 20\nvary thicknesses 10\n"
                            "goal contrast 20\ngoal contrast balance all except 0.3\n"),
                  "there is no field 0.3 in 'fields'");
    checkMentions(rejection("[trial 1]\nfields 0 0.5\nfrequencies 20\ngoal spot-deviation 1 1\n"
                            "goal spot-deviation x 1 1\n"),
                  "either as one row or as x and y rows");
    checkMentions(rejection("[trial 1]\nfields 0\nfrequencies 20\ngoal paraxial focus 3\n"),
                  "unknown paraxial quantity 'focus'");
}

TEST(trial_namesTheTrialsAFileDefines) {
    std::string text = withTrials(
        SUMMICRON, "[trial 1]\nfields 0\nfrequencies 20\n[trial 4]\nfields 0\nfrequencies 20\n");
    try {
        OptimizationTrial::read(text, 2, true);
        CHECK(false);
    } catch (const TrialException &e) {
        CHECK_STR_EQ(e.getMessage(),
                     "there is no [trial 2] in this prescription; it defines trials 1, 4");
    }
    // The prescription reader is unaffected by the trials.
    OpticalBenchDataImporter::LensSpecifications specs;
    specs.parse_buffer(text);
    CHECK_EQ(specs.get_surfaces().size(), static_cast<std::size_t>(11));
}

} // namespace
