// C++ port of org.redukti.optim.TestOptimLMder.
//
// A full Levenberg-Marquardt solve of an 18-variable Noct-Nikkor-like design
// against a fixed final RMS. That number is the strongest single check on the
// optimizer port: it depends on every Var, every Goal, the finite-difference
// Jacobian and MinPack's lmder all agreeing with the Java at once.
//
// Disabled by default, as the Java is: it gates this behind
// -Doptimization.runSlowTests=true. Set RAYOPTICS_RUN_SLOW_TESTS to a non-empty
// value other than 0 to run it.
//
// NOT YET PASSING. The C++ raises "Bad optional access" partway through the
// solve, from one of the sites that deliberately mirrors a Java
// NullPointerException with a throwing optional. Whether the Java reaches the
// same point is unverified -- being gated off, this test may not have been run
// against the current Java either. Left here, disabled, as the next thing to
// look at rather than deleted.
#include "TestHarness.h"

#include "redukti/optim/LMDer.h"
#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/spec/Prescription.h"
#include "redukti/spec/SurfaceType.h"

#include <cstdlib>
#include <memory>
#include <vector>

namespace {

using redukti::optim::Analysis;
using redukti::optim::Goal;
using redukti::optim::GoalParax;
using redukti::optim::GoalRayAberration;
using redukti::optim::GoalSpotMaxRadius;
using redukti::optim::GoalSpotRMS;
using redukti::optim::LMDerMeritFunction;
using redukti::optim::ParaxHelper;
using redukti::optim::Var;
using redukti::optim::VarAsphCoeff;
using redukti::optim::VarAsphK;
using redukti::optim::VarRadius;
using redukti::spec::Prescription;
using redukti::spec::SurfaceType;

bool runSlowTests() {
#ifdef _WIN32
    char *value = nullptr;
    std::size_t size = 0;
    _dupenv_s(&value, &size, "RAYOPTICS_RUN_SLOW_TESTS");
    const bool run = value != nullptr && value[0] != '\0' && value[0] != '0';
    std::free(value);
#else
    const char *value = std::getenv("RAYOPTICS_RUN_SLOW_TESTS");
    const bool run = value != nullptr && value[0] != '\0' && value[0] != '0';
#endif
    return run;
}

/** Prescription based off measurements from Nikkor tale. */
Prescription getPrescription() {
    return Prescription(58.0, 1.2, 40.9, 43.28, false)
        .surf(79.9975, 6.885, 50.4875, 1.795, 45.31, "J-LASF017", "Hikari")
        .asph(SurfaceType::ASPH_EVEN, 0, {0.0, 0.0, 0.0, 0.0, 0.0})
        .surf(0, 0.1, 50.4875)
        .surf(33.737, 9.75, 44.832, 1.8485, 43.79, "J-LASFH22", "Hikari")
        .surf(70.18675, 1.56, 44.832)
        .surf(134.505, 2.87, 42.169, 1.74, 28.3, "S-TIH3", "Ohara")
        .surf(22.3687, 8.44, 32.12841)
        .stop(7.95, 31.227)
        .surf(-23.02418, 1.64, 31.445, 1.74077, 27.79, "S-TIH13", "Ohara")
        .surf(306.553, 8.196, 40.2, 1.788, 47.37, "TAF4", "Hoya")
        .surf(-37.555, 0.15, 40.2)
        .surf(-396.94, 6.147, 39.5, 1.7725, 46.62, "J-LASF016", "Hikari")
        .surf(-52.56789, 0.0, 39.5)
        .surf(223.8426, 4.016, 38.275, 1.795, 45.31, "J-LASF017", "Hikari")
        .surf(-94.08052, 37.78, 38.275)
        .build();
}

} // namespace

TEST(optim_nikkor58mm_noct_lmder_solve) {
    if (!runSlowTests())
        return;

    auto prescription = getPrescription();
    Analysis analysis(&prescription, {0.0, 0.1, 0.3, 0.5, 0.7, 1.0}, {30});

    std::vector<std::shared_ptr<Var>> vars{
        std::make_shared<VarRadius>(&prescription, 0),
        std::make_shared<VarAsphK>(&prescription, 0),
        std::make_shared<VarAsphCoeff>(&prescription, 0, 1, 1E6),
        std::make_shared<VarAsphCoeff>(&prescription, 0, 2, 1E9),
        std::make_shared<VarAsphCoeff>(&prescription, 0, 3, 1E11),
        std::make_shared<VarAsphCoeff>(&prescription, 0, 4, 1E14),
        std::make_shared<VarRadius>(&prescription, 2),
        std::make_shared<VarRadius>(&prescription, 3),
        std::make_shared<VarRadius>(&prescription, 4),
        std::make_shared<VarRadius>(&prescription, 5),
        std::make_shared<VarRadius>(&prescription, 7),
        std::make_shared<VarRadius>(&prescription, 8),
        std::make_shared<VarRadius>(&prescription, 9),
        std::make_shared<VarRadius>(&prescription, 10),
        std::make_shared<VarRadius>(&prescription, 11),
        std::make_shared<VarRadius>(&prescription, 12),
        std::make_shared<VarRadius>(&prescription, 13),
    };

    std::vector<std::shared_ptr<Goal>> goals{
        std::make_shared<GoalSpotRMS>(&analysis, 1, 10.0, 7.0),
        std::make_shared<GoalSpotRMS>(&analysis, 2, 12.0, 5.0),
        std::make_shared<GoalSpotRMS>(&analysis, 3, 20.0, 2.0),
        std::make_shared<GoalSpotRMS>(&analysis, 4, 30.0, 2.0),
        std::make_shared<GoalSpotRMS>(&analysis, 5, 40.0, 2.0),
        std::make_shared<GoalSpotRMS>(&analysis, 6, 50.0, 3.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 1, 30.0, 2.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 2, 35.0, 2.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 3, 40.0, 2.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 4, 80.0, 2.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 5, 120.0, 2.0),
        std::make_shared<GoalSpotMaxRadius>(&analysis, 6, 200.0, 2.0),
        std::make_shared<GoalRayAberration>(&analysis, 6, 0, 0, 587.5618, 0, 1),
        std::make_shared<GoalRayAberration>(&analysis, 6, 0, -1, 587.5618, 0, 1),
        std::make_shared<GoalRayAberration>(&analysis, 6, 0, 0, 486.1327, 0, 1),
        std::make_shared<GoalRayAberration>(&analysis, 6, 0, -1, 486.1327, 0, 1),
        std::make_shared<GoalParax>(&analysis, ParaxHelper::Effective_focal_length, 58.0, 1.0),
        std::make_shared<GoalParax>(&analysis, ParaxHelper::Fno, 1.2, 1.0),
        std::make_shared<GoalParax>(&analysis, ParaxHelper::Back_focal_length, 37.78, 1.0),
    };

    LMDerMeritFunction f(&analysis, vars, goals, false);
    analysis.compute();
    auto lm = f.getSolver();
    double initialRMS = f.getRMS();
    lm->solve();
    double finalRMS = f.getRMS();

    CHECK(finalRMS < initialRMS);
    CHECK_CLOSE(finalRMS, 22.392, 1e-3);
}
