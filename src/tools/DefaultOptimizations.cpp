// C++ port of org.redukti.tools.DefaultOptimizations
#include "redukti/tools/DefaultOptimizations.h"

#include "redukti/Exceptions.h"
#include "redukti/optim/LMDer.h"
#include "redukti/optim/OptimizationBuilder.h"

namespace redukti::tools {

namespace {

using optim::OptimizationBuilder;

/** Central field, the single field these optimizations are judged on. */
const std::vector<double> CENTRAL_FIELD{0.0};
/** Goals aim at perfect contrast, so the solver simply maximizes it. */
constexpr double PERFECT_MTF = 1.0;
/** Pupil sampling for the contrast objective. */
constexpr int CONTRAST_RINGS = 6;
constexpr int CONTRAST_SPOKES = 12;

} // namespace

int DefaultOptimizations::findBackFocusSurface(const spec::Prescription &prescription) {
    const auto &surfaces = prescription.get_surfaces();
    if (surfaces.empty())
        return -1;
    // Walk back over a trailing cover glass block, if there is one.
    int firstCoverGlass = -1;
    for (int i = static_cast<int>(surfaces.size()) - 1; i >= 0; i--) {
        if (surfaces[static_cast<std::size_t>(i)].is_cover_glass())
            firstCoverGlass = i;
        else if (firstCoverGlass >= 0)
            break;
    }
    int candidate =
        firstCoverGlass >= 0 ? firstCoverGlass - 1 : static_cast<int>(surfaces.size()) - 1;
    if (candidate < 0)
        return -1;
    // The back focus is an airspace; a glass thickness here means the
    // prescription is not shaped the way this rule assumes.
    if (surfaces[static_cast<std::size_t>(candidate)].get_refractive_index() != 0.0)
        return -1;
    return candidate;
}

std::vector<int> DefaultOptimizations::findVariableThicknesses(
    const spec::Prescription &prescription, int backFocusSurface) {
    const auto &surfaces = prescription.get_surfaces();
    std::vector<int> result;
    for (int i = 0; i < static_cast<int>(surfaces.size()); i++)
        if (isVariableAirspace(surfaces[static_cast<std::size_t>(i)], i, backFocusSurface))
            result.push_back(i);
    return result;
}

bool DefaultOptimizations::isVariableAirspace(const spec::SurfaceType &surface, int index,
                                              int backFocusSurface) {
    return index != backFocusSurface && surface._thickness_by_scenario.has_value() &&
           surface.get_refractive_index() == 0.0;
}

DefaultOptimizations::Result DefaultOptimizations::optimizeThicknesses(
    spec::Prescription *prescription, const std::vector<int> &surfaces,
    const std::vector<int> &mtfFrequencies, int configuration, spec::VigType vigType,
    bool dLineOnly, Objective objective) {
    if (surfaces.empty())
        throw IllegalArgumentException("no surfaces to optimize");
    auto builder = OptimizationBuilder::builder(prescription);
    builder.fields(CENTRAL_FIELD)
        .mtfFrequencies(mtfFrequencies)
        .scenario(configuration)
        .vignetting(vigType)
        .dLineOnly(dLineOnly)
        .varyThicknesses(surfaces)
        .applyThicknessConstraints();
    // Goals at the central field, both meridians, driving contrast up. Without
    // them the only residuals are the automatic focal length and f-number
    // anchors, which do not depend on an airspace at all, and the solve sits
    // still.
    if (objective == Objective::CONTRAST) {
        std::vector<OptimizationBuilder::ContrastGoals> goals;
        for (int frequency : mtfFrequencies)
            goals.push_back(OptimizationBuilder::contrast(frequency, {1.0}));
        builder.contrastSampling(CONTRAST_RINGS, CONTRAST_SPOKES)
            .calibrateContrastFrequency(true)
            .contrastGoals(goals);
    } else {
        std::vector<OptimizationBuilder::MtfGoals> goals;
        for (int frequency : mtfFrequencies)
            goals.push_back(OptimizationBuilder::mtf(frequency, {PERFECT_MTF}, {PERFECT_MTF}));
        builder.mtfGoals(goals);
    }
    auto setup = builder.build();
    auto *analysis = setup.analysis();
    auto merit = setup.meritFunction(false);
    analysis->compute();
    double before = merit.getRMS();
    int status = merit.getSolver()->solve();
    double after = merit.getRMS();
    return Result{surfaces, status, before, after};
}

} // namespace redukti::tools
