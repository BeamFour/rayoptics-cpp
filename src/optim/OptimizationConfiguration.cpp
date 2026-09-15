// C++ port of org.redukti.optim.OptimizationConfiguration
#include "redukti/optim/OptimizationConfiguration.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/optim/Analysis.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/util/Args.h"

#include <algorithm>
#include <cmath>

namespace redukti::optim {

namespace {

std::string yesNo(bool value) {
    return value ? "yes" : "no";
}

/** OptimizationTrial::line, reached through the same short name the Java uses. */
void line(std::string &sb, const std::string &key, const std::string &values) {
    OptimizationTrial::line(sb, key, values);
}

} // namespace

void OptimizationConfiguration::gaussianSampling(int rings, int spokes,
                                                 double innerPupilRadius) {
    if (rings < 1 || spokes < 3)
        throw IllegalArgumentException(
            "Gaussian quadrature requires at least 1 ring and 3 spokes");
    if (!std::isfinite(innerPupilRadius) || innerPupilRadius < 0.0 ||
        innerPupilRadius >= 1.0)
        throw IllegalArgumentException("Inner pupil radius must be finite and in [0, 1)");
    this->gaussianQuadratureRings = rings;
    this->gaussianQuadratureSpokes = spokes;
    this->gaussianQuadratureInnerRadius = innerPupilRadius;
}

std::string OptimizationConfiguration::toTrial(int number) const {
    auto effective = effectiveAnalysis(false);
    // Preserve configured intent, even when the combination is invalid. Using
    // execution's Gaussian override here would silently repair a rejected setup
    // into a different, valid one when its trial is read back.
    bool configuredHexapolar = useHexapolarSpotPattern || spotMaxRadiusGoals.has_value();
    std::string sb;
    sb += "[trial " + intToString(number) + "]\n";
    if (description.has_value())
        line(sb, "description", *description);
    if (outdir.has_value())
        line(sb, "outdir", *outdir);
    line(sb, "configuration", intToString(scenario));
    if (fields.has_value())
        line(sb, "fields", OptimizationTrial::format(*fields));
    if (mtfFrequencies.has_value())
        line(sb, "frequencies", OptimizationTrial::format(*mtfFrequencies));
    line(sb, "weighted", yesNo(weighted));
    line(sb, "d-line-only", yesNo(dLineOnly));
    line(sb, "vignetting", OptimizationTrial::kebab(util::Args::vig_type_name(vigType)) +
                               (freezeVignetting ? " frozen" : ""));
    if (!checkSpotApertures || (effective.spots && !configuredHexapolar))
        line(sb, "check-spot-apertures", yesNo(checkSpotApertures));
    // Only what this trial changed: the values came from the same constants, so
    // an untouched trial compares equal and writes nothing.
    if (solverTolerances.ftol() != SolverTolerances::defaultFtol())
        line(sb, "solver ftol", OptimizationTrial::format(solverTolerances.ftol()));
    if (solverTolerances.xtol() != SolverTolerances::DEFAULT_XTOL)
        line(sb, "solver xtol", OptimizationTrial::format(solverTolerances.xtol()));
    if (solverTolerances.gtol() != SolverTolerances::defaultGtol())
        line(sb, "solver gtol", OptimizationTrial::format(solverTolerances.gtol()));
    if (solverTolerances.maxEvaluations() != SolverTolerances::FROM_VARIABLE_COUNT)
        line(sb, "solver max-evaluations", intToString(solverTolerances.maxEvaluations()));

    if (allCurvatureSurfaces)
        line(sb, "vary curvatures", allExcept(curvatureExclusions));
    else if (!curvatureSurfaces.empty())
        line(sb, "vary curvatures", OptimizationTrial::format(curvatureSurfaces));
    if (allThicknessSurfaces)
        line(sb, "vary thicknesses", allExcept(thicknessExclusions));
    else if (!thicknessSurfaces.empty())
        line(sb, "vary thicknesses", OptimizationTrial::format(thicknessSurfaces));
    if (includeExistingAspherics)
        line(sb, "vary aspherics", "existing");
    // A LinkedHashMap in the Java: one row per surface, surfaces in the order their first
    // term was given.
    std::vector<int> termSurfaces;
    std::vector<std::string> termRows;
    for (const auto &term : asphericTerms) {
        std::string written =
            term.index < 0 ? std::string("K")
                           : intToString(term.index) +
                                 (term.scale.has_value()
                                      ? ":" + OptimizationTrial::format(*term.scale)
                                      : "");
        auto found = std::find(termSurfaces.begin(), termSurfaces.end(), term.surface);
        if (found == termSurfaces.end()) {
            termSurfaces.push_back(term.surface);
            termRows.push_back(written);
        } else
            termRows[static_cast<std::size_t>(found - termSurfaces.begin())] += " " + written;
    }
    for (std::size_t i = 0; i < termSurfaces.size(); i++)
        line(sb, "vary aspherics", intToString(termSurfaces[i]) + " " + termRows[i]);

    if (curvatureConstraintWeight.has_value())
        line(sb, "constrain curvatures", OptimizationTrial::format(*curvatureConstraintWeight));
    if (thicknessConstraintWeight.has_value())
        line(sb, "constrain thicknesses", OptimizationTrial::format(*thicknessConstraintWeight));
    if (edgeThicknessConstraintWeight.has_value())
        line(sb, "constrain edges", OptimizationTrial::format(*edgeThicknessConstraintWeight));

    if (!contrastGoals.empty()) {
        std::vector<int> frequencies;
        for (const auto &goal : contrastGoals)
            frequencies.push_back(goal.frequency);
        line(sb, "goal contrast", OptimizationTrial::format(frequencies));
        contrastWeights(sb, true);
        contrastWeights(sb, false);
        if (contrastBalanceFields.has_value())
            line(sb, "goal contrast",
                 "balance " + balance() + " weight " +
                     OptimizationTrial::format(contrastBalanceWeight));
        line(sb, "goal contrast",
             "sampling " + intToString(contrastRings) + " " + intToString(contrastSpokes));
        line(sb, "goal contrast", "calibrate " + yesNo(calibrateContrastFrequency));
        line(sb, "goal contrast", "exit-pupil-aiming " + yesNo(aimContrastAtExitPupil));
        line(sb, "goal contrast", "centering " + yesNo(centerContrastResiduals));
    }
    for (const auto &goal : mtfGoals) {
        std::string frequency = intToString(goal.frequency);
        line(sb, "goal mtf", frequency + " sag " + OptimizationTrial::format(goal.sagittal));
        line(sb, "goal mtf", frequency + " tan " + OptimizationTrial::format(goal.tangential));
        if (goal.sagittalWeights == goal.tangentialWeights) {
            if (!allOnes(goal.sagittalWeights))
                line(sb, "goal mtf",
                     frequency + " weights " + OptimizationTrial::format(goal.sagittalWeights));
        } else {
            if (!allOnes(goal.sagittalWeights))
                line(sb, "goal mtf", frequency + " sag weights " +
                                         OptimizationTrial::format(goal.sagittalWeights));
            if (!allOnes(goal.tangentialWeights))
                line(sb, "goal mtf", frequency + " tan weights " +
                                         OptimizationTrial::format(goal.tangentialWeights));
        }
    }
    spotGoals(sb, "goal spot-rms", spotRmsGoals);
    spotGoals(sb, "goal spot-max-radius", spotMaxRadiusGoals);
    if (addSpotDeviationGoals) {
        if (spotDeviationXWeights == spotDeviationYWeights)
            line(sb, "goal spot-deviation", OptimizationTrial::format(*spotDeviationXWeights));
        else {
            line(sb, "goal spot-deviation",
                 "x " + OptimizationTrial::format(*spotDeviationXWeights));
            line(sb, "goal spot-deviation",
                 "y " + OptimizationTrial::format(*spotDeviationYWeights));
        }
    }
    if (gaussianQuadratureRings != DEFAULT_GAUSSIAN_QUADRATURE_RINGS ||
        gaussianQuadratureSpokes != DEFAULT_GAUSSIAN_QUADRATURE_SPOKES ||
        gaussianQuadratureInnerRadius != 0.0 || (effective.spots && !configuredHexapolar))
        line(sb, "goal spot sampling",
             "gaussian " + intToString(gaussianQuadratureRings) + " " +
                 intToString(gaussianQuadratureSpokes) +
                 (gaussianQuadratureInnerRadius != 0.0
                      ? " " + OptimizationTrial::format(gaussianQuadratureInnerRadius)
                      : ""));
    if (configuredHexapolar)
        line(sb, "goal spot sampling", "hexapolar " + intToString(hexapolarSpotRays));
    line(sb, "goal ray-aberrations", yesNo(addRayAberrationGoals));
    for (const auto &goal : paraxialGoals)
        line(sb, "goal paraxial",
             OptimizationTrial::paraxialName(goal.paraxId) + " " +
                 OptimizationTrial::format(goal.target) +
                 (goal.weight != 1.0 ? " weight " + OptimizationTrial::format(goal.weight)
                                     : ""));
    return sb;
}

std::string OptimizationConfiguration::allExcept(const std::vector<int> &exclusions) {
    return exclusions.empty() ? "all" : "all except " + OptimizationTrial::format(exclusions);
}

void OptimizationConfiguration::contrastWeights(std::string &sb, bool sagittal) const {
    std::string direction = sagittal ? "sag" : "tan";
    const std::vector<double> &first =
        sagittal ? contrastGoals[0].sagittalWeights : contrastGoals[0].tangentialWeights;
    bool shared = true;
    for (const auto &goal : contrastGoals)
        if ((sagittal ? goal.sagittalWeights : goal.tangentialWeights) != first)
            shared = false;
    if (shared) {
        if (!allOnes(first))
            line(sb, "goal contrast", direction + " " + OptimizationTrial::format(first));
        return;
    }
    for (const auto &goal : contrastGoals) {
        const std::vector<double> &weights =
            sagittal ? goal.sagittalWeights : goal.tangentialWeights;
        if (!allOnes(weights))
            line(sb, "goal contrast", intToString(goal.frequency) + " " + direction + " " +
                                          OptimizationTrial::format(weights));
    }
}

std::string OptimizationConfiguration::balance() const {
    bool all = true, none = true;
    for (bool flag : *contrastBalanceFields) {
        all = all && flag;
        none = none && !flag;
    }
    if (all)
        return "all";
    if (none || !fields.has_value() || fields->size() != contrastBalanceFields->size()) {
        std::string flags;
        for (bool flag : *contrastBalanceFields) {
            if (!flags.empty())
                flags += " ";
            flags += yesNo(flag);
        }
        return flags;
    }
    std::string except;
    for (std::size_t i = 0; i < fields->size(); i++)
        if (!(*contrastBalanceFields)[i]) {
            if (!except.empty())
                except += " ";
            except += OptimizationTrial::format((*fields)[i]);
        }
    return "all except " + except;
}

void OptimizationConfiguration::spotGoals(
    std::string &sb, const char *key, const std::optional<OptimizationBuilder::SpotGoals> &goals) {
    if (!goals.has_value())
        return;
    line(sb, key, OptimizationTrial::format(goals->targets));
    if (!allOnes(goals->weights))
        line(sb, key, "weights " + OptimizationTrial::format(goals->weights));
}

OptimizationConfiguration::EffectiveAnalysis OptimizationConfiguration::effectiveAnalysis(
    bool customMaximumRadius) const {
    bool mtf = !mtfGoals.empty();
    bool spots = spotRmsGoals.has_value() || spotMaxRadiusGoals.has_value() ||
                 addSpotDeviationGoals || mtf;
    bool hexapolar = !addSpotDeviationGoals &&
                     (useHexapolarSpotPattern || spotMaxRadiusGoals.has_value() ||
                      customMaximumRadius);
    return EffectiveAnalysis{spots, addRayAberrationGoals, mtf, hexapolar,
                             addSpotDeviationGoals};
}

void OptimizationConfiguration::configureAnalysis(Analysis &analysis,
                                                  const EffectiveAnalysis &effective,
                                                  bool customGoals) const {
    analysis.vignetting(vigType)
        .freezing_vignetting(freezeVignetting)
        .checking_spot_apertures(checkSpotApertures);
    if (effective.hexapolar)
        analysis.using_hexapolar_pattern(hexapolarSpotRays);
    else
        analysis.using_gauss_quadrature_pattern(gaussianQuadratureRings,
                                                gaussianQuadratureSpokes,
                                                gaussianQuadratureInnerRadius);
    if (effective.retainFailedRays)
        analysis.retaining_failed_spot_rays(true);
    if (!contrastGoals.empty()) {
        if (calibrateContrastFrequency && aimContrastAtExitPupil)
            throw IllegalArgumentException("Contrast frequency calibration and exit-pupil "
                                           "aiming are mutually exclusive");
        std::vector<int> frequencies;
        frequencies.reserve(contrastGoals.size());
        for (const auto &goal : contrastGoals)
            frequencies.push_back(goal.frequency);
        analysis.using_contrast_analysis(frequencies, contrastRings, contrastSpokes);
        analysis.calibrating_contrast_frequency(calibrateContrastFrequency);
        analysis.aiming_contrast_at_exit_pupil(aimContrastAtExitPupil);
        analysis.centering_contrast_residuals(centerContrastResiduals);
    }
    // Unknown factories may need any analysis. Keep the Analysis defaults (or
    // the factory's explicit choices) instead of disabling work based on built-in goals.
    if (!customGoals)
        analysis.required_analyses(effective.spots, effective.rayAberrations, effective.mtf);
}

bool OptimizationConfiguration::allOnes(const std::vector<double> &values) {
    for (double value : values)
        if (value != 1.0)
            return false;
    return true;
}

} // namespace redukti::optim
