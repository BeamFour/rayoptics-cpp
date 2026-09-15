// C++ port of org.redukti.optim.OptimizationBuilder
#include "redukti/optim/OptimizationBuilder.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/optim/OptimizationConfiguration.h"
#include "redukti/optim/OptimizationTrial.h"
#include "redukti/optim/OptimizationValidation.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/util/Orientation.h"
#include "redukti/util/Args.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace redukti::optim {

namespace Orientation = rayoptics::util::Orientation;
using rayoptics::seq::Glass;

namespace {

/** Java's `weights == null ? unitWeights(targets) : copy(weights)`. */
std::vector<double> unitWeightsFor(const std::vector<double> &targets) {
    return std::vector<double>(targets.size(), 1.0);
}

/** Java's OptimizationBuilder.contains(int[], int). */
bool contains(const std::vector<int> &values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

/** True when any goal in the list is a T; the Java uses `instanceof`. */
template <typename T> bool anyGoalIs(const std::vector<std::shared_ptr<Goal>> &goals) {
    for (const auto &goal : goals)
        if (dynamic_cast<const T *>(goal.get()) != nullptr)
            return true;
    return false;
}

} // namespace

OptimizationBuilder::OptimizationBuilder(spec::Prescription *prescription)
    : OptimizationBuilder(prescription, OptimizationConfiguration()) {}

OptimizationBuilder::OptimizationBuilder(spec::Prescription *prescription,
                                         const OptimizationConfiguration &configuration_)
    : prescription_(prescription),
      configuration(std::make_unique<OptimizationConfiguration>(configuration_.copy())) {
    if (prescription == nullptr)
        throw IllegalArgumentException("prescription must not be null");
    // The Java tests `prescription._surfaces == null`, the array build() fills;
    // here the flag records the same thing without a second copy of the list.
    if (!prescription->_built)
        throw IllegalArgumentException("prescription must be built before optimization");
}

OptimizationBuilder::OptimizationBuilder(const OptimizationBuilder &other)
    : prescription_(other.prescription_),
      configuration(std::make_unique<OptimizationConfiguration>(*other.configuration)),
      additionalVariables_(other.additionalVariables_),
      additionalGoalFactories(other.additionalGoalFactories) {}

OptimizationBuilder &OptimizationBuilder::operator=(const OptimizationBuilder &other) {
    if (this != &other) {
        prescription_ = other.prescription_;
        configuration = std::make_unique<OptimizationConfiguration>(*other.configuration);
        additionalVariables_ = other.additionalVariables_;
        additionalGoalFactories = other.additionalGoalFactories;
    }
    return *this;
}

OptimizationBuilder::OptimizationBuilder(OptimizationBuilder &&other) noexcept = default;
OptimizationBuilder &OptimizationBuilder::operator=(OptimizationBuilder &&other) noexcept = default;
OptimizationBuilder::~OptimizationBuilder() = default;

const std::optional<std::string> &OptimizationBuilder::description() const {
    return configuration->description;
}

const std::optional<std::string> &OptimizationBuilder::outdir() const {
    return configuration->outdir;
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

OptimizationBuilder &OptimizationBuilder::description(
    const std::optional<std::string> &description) {
    configuration->description = description;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::outdir(const std::optional<std::string> &outdir) {
    configuration->outdir = outdir;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::fields(const std::vector<double> &fields_) {
    configuration->fields = fields_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::mtfFrequencies(
    const std::vector<int> &frequencies) {
    configuration->mtfFrequencies = frequencies;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::weighted(bool weighted_) {
    configuration->weighted = weighted_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::dLineOnly(bool dLineOnly_) {
    configuration->dLineOnly = dLineOnly_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::scenario(int scenario_) {
    if (scenario_ < 0)
        throw IllegalArgumentException("scenario must be non-negative, got " +
                                       intToString(scenario_));
    configuration->scenario = scenario_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::vignetting(spec::VigType vigType_) {
    configuration->vigType = vigType_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::freezeVignetting(bool freeze) {
    configuration->freezeVignetting = freeze;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::gaussianQuadratureSampling(
    int rings, int spokes, double innerPupilRadius) {
    configuration->gaussianSampling(rings, spokes, innerPupilRadius);
    return *this;
}

OptimizationBuilder &OptimizationBuilder::checkSpotApertures(bool check) {
    configuration->checkSpotApertures = check;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::hexapolarSampling(int numRays) {
    if (numRays < 1)
        throw IllegalArgumentException("hexapolar spot rings must be at least 1");
    configuration->useHexapolarSpotPattern = true;
    configuration->hexapolarSpotRays = numRays;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastSampling(int rings, int spokes) {
    if (rings < 1 || spokes < 3)
        throw IllegalArgumentException(
            "contrast sampling requires at least 1 ring and 3 spokes");
    configuration->contrastRings = rings;
    configuration->contrastSpokes = spokes;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::calibrateContrastFrequency(bool value) {
    configuration->calibrateContrastFrequency = value;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::aimContrastAtExitPupil(bool value) {
    configuration->aimContrastAtExitPupil = value;
    return *this;
}

/**
 * Subtract the constant part of each contrast block, so the residuals carry the
 * variance the OTF modulus depends on rather than the un-centred second moment.
 *
 * A constant wavefront difference across the pupil is tilt, which displaces the
 * image and costs no MTF - but it is reducible, so leaving it in offers the solver
 * merit reduction that buys nothing. It is identically zero in the sagittal direction
 * by symmetry and reaches 57% of an outer-field tangential block on the Leica 75/2,
 * which biases the astigmatic focus split toward tangential.
 *
 * Off by default because it changes every contrast residual. See
 * ContrastAnalysis#center_residuals(ContrastAnalysisResult, int).
 */
OptimizationBuilder &OptimizationBuilder::centerContrastResiduals(bool value) {
    configuration->centerContrastResiduals = value;
    return *this;
}

// ---------------------------------------------------------------------------
// Variables
// ---------------------------------------------------------------------------

OptimizationBuilder &OptimizationBuilder::varyCurvatures(
    const std::vector<int> &surfaces) {
    configuration->curvatureSurfaces = surfaces;
    configuration->allCurvatureSurfaces = false;
    configuration->curvatureExclusions.clear();
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyAllCurvatures() {
    return varyAllCurvaturesExcept({});
}

OptimizationBuilder &OptimizationBuilder::varyAllCurvaturesExcept(
    const std::vector<int> &surfaces) {
    configuration->curvatureSurfaces.clear();
    configuration->allCurvatureSurfaces = true;
    configuration->curvatureExclusions = surfaces;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyThicknesses(
    const std::vector<int> &surfaces) {
    configuration->thicknessSurfaces = surfaces;
    configuration->allThicknessSurfaces = false;
    configuration->thicknessExclusions.clear();
    return *this;
}

/**
 * Vary every thickness, air spaces and element thicknesses alike. Surfaces with zero
 * thickness are excluded, being coincident rather than a space to open up.
 *
 * The counterpart to #varyAllCurvatures(), and best paired with
 * #applyThicknessConstraints(double) - with every space free and nothing holding the
 * layout, the solver will collapse gaps and drive elements through one another.
 */
OptimizationBuilder &OptimizationBuilder::varyAllThicknesses() {
    return varyAllThicknessesExcept({});
}

OptimizationBuilder &OptimizationBuilder::varyAllThicknessesExcept(
    const std::vector<int> &surfaces) {
    configuration->thicknessSurfaces.clear();
    configuration->allThicknessSurfaces = true;
    configuration->thicknessExclusions = surfaces;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyExistingAspherics(bool include) {
    configuration->includeExistingAspherics = include;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyConic(int surface) {
    return addAsphericTerm(surface, -1, std::nullopt);
}

OptimizationBuilder &OptimizationBuilder::varyAsphericCoefficient(int surface, int index) {
    return addAsphericTerm(surface, index, std::nullopt);
}

OptimizationBuilder &OptimizationBuilder::varyAsphericCoefficient(int surface, int index,
                                                                  double scale) {
    if (!std::isfinite(scale) || scale <= 0.0)
        throw IllegalArgumentException(
            "the scale of an aspheric coefficient must be finite and positive");
    return addAsphericTerm(surface, index, scale);
}

OptimizationBuilder &OptimizationBuilder::addAsphericTerm(
    int surface, int index, const std::optional<double> &scale) {
    if (surface < 0 || surface >= static_cast<int>(prescription_->_surface_list.size()))
        throw IllegalArgumentException("aspheric surface is out of range: " +
                                       intToString(surface));
    const auto &definition = prescription_->_surface_list[static_cast<std::size_t>(surface)];
    if (definition.is_aperture_stop() || definition.is_field_stop())
        throw IllegalArgumentException("surface " + intToString(surface) +
                                       " is a stop; it cannot be aspheric");
    for (const AsphericTerm &term : configuration->asphericTerms)
        if (term.surface == surface && term.index == index)
            throw IllegalArgumentException(
                (index < 0 ? std::string("the conic constant")
                           : "coefficient " + intToString(index)) +
                " of surface " + intToString(surface) + " is varied twice");
    if (index >= 0) {
        powerOf(asphereTypeOf(surface), index);
        if (!scale.has_value() && coefficientOf(surface, index) == 0.0 &&
            !(definition._diameter > 0.0))
            throw IllegalArgumentException(
                "surface " + intToString(surface) +
                " has no diameter to derive a scale for coefficient " + intToString(index) +
                " from; give the coefficient a scale");
    }
    configuration->asphericTerms.push_back(AsphericTerm{surface, index, scale});
    return *this;
}

int OptimizationBuilder::asphereTypeOf(int surface) const {
    const auto &definition = prescription_->_surface_list[static_cast<std::size_t>(surface)];
    if (definition.is_aspheric())
        return definition._asph_type;
    if (prescription_->has_odd_aspheric())
        return spec::SurfaceType::ASPH_ODD;
    if (prescription_->has_even_a2_aspheric())
        return spec::SurfaceType::ASPH_EVEN_A2;
    return spec::SurfaceType::ASPH_EVEN;
}

int OptimizationBuilder::powerOf(int asphereType, int index) {
    switch (asphereType) {
    case spec::SurfaceType::ASPH_ODD:
        if (index >= 2)
            return index + 1;
        throw IllegalArgumentException(
            "coefficient " + intToString(index) +
            " is not a term of an odd asphere, whose terms start at index 2, the A3 term");
    case spec::SurfaceType::ASPH_EVEN_A2:
        return 2 * (index + 1);
    default:
        if (index >= 1)
            return 2 * (index + 1);
        throw IllegalArgumentException(
            "coefficient " + intToString(index) +
            " is not a term of an even asphere, whose terms start at index 1, the A4 term");
    }
}

double OptimizationBuilder::coefficientOf(int surface, int index) const {
    const auto &coefficients =
        prescription_->_surface_list[static_cast<std::size_t>(surface)]._coeffs;
    return coefficients.has_value() && index < static_cast<int>(coefficients->size())
               ? (*coefficients)[static_cast<std::size_t>(index)]
               : 0.0;
}

bool OptimizationBuilder::hasExplicitAsphericTerms(int surface) const {
    for (const AsphericTerm &term : configuration->asphericTerms)
        if (term.surface == surface)
            return true;
    return false;
}

OptimizationBuilder &OptimizationBuilder::additionalVariables(
    const std::vector<std::shared_ptr<Var>> &variables) {
    for (const auto &variable : variables) {
        if (variable == nullptr)
            throw IllegalArgumentException("additional variables must not contain null");
        if (variable->_prescription != prescription_)
            throw IllegalArgumentException(
                "additional variables must use this builder's prescription");
        additionalVariables_.push_back(variable);
    }
    return *this;
}

double OptimizationBuilder::thicknessOf(int surface) const {
    const auto &definition =
        prescription_->_surface_list[static_cast<std::size_t>(surface)];
    return definition._thickness_by_scenario.has_value()
               ? (*definition._thickness_by_scenario)[static_cast<std::size_t>(configuration->scenario)]
               : definition._thickness;
}

double OptimizationBuilder::focalLengthOf() const {
    // The Java field is a nullable array; here an empty vector is the same state.
    return !prescription_->_focal_length_by_scenario.empty()
               ? prescription_
                     ->_focal_length_by_scenario[static_cast<std::size_t>(configuration->scenario)]
               : prescription_->_focal_length;
}

double OptimizationBuilder::fNumberOf() const {
    return !prescription_->_f_number_by_scenario.empty()
               ? prescription_->_f_number_by_scenario[static_cast<std::size_t>(configuration->scenario)]
               : prescription_->_fno;
}

void OptimizationBuilder::validateScenario() const {
    if (configuration->scenario == 0)
        return;
    int available = scenarioCount();
    if (configuration->scenario >= available)
        throw IllegalArgumentException(
            "scenario " + intToString(configuration->scenario) +
            " requested but the prescription defines " + intToString(available) +
            (available == 1 ? " (it is not multi-configuration)" : ""));
}

int OptimizationBuilder::scenarioCount() const {
    int count = 1;
    if (!prescription_->_focal_length_by_scenario.empty())
        count = std::max(
            count, static_cast<int>(prescription_->_focal_length_by_scenario.size()));
    for (const auto &surface : prescription_->_surface_list)
        if (surface._thickness_by_scenario.has_value())
            count =
                std::max(count, static_cast<int>(surface._thickness_by_scenario->size()));
    return count;
}

// ---------------------------------------------------------------------------
// Goals
// ---------------------------------------------------------------------------

OptimizationBuilder::MtfGoals OptimizationBuilder::mtf(
    int frequency, const std::vector<double> &sagittal,
    const std::vector<double> &tangential) {
    return MtfGoals(frequency, sagittal, tangential, std::nullopt, std::nullopt);
}

OptimizationBuilder::MtfGoals OptimizationBuilder::mtf(
    int frequency, const std::vector<double> &sagittal,
    const std::vector<double> &tangential, const std::vector<double> &weights) {
    return MtfGoals(frequency, sagittal, tangential, weights, weights);
}

OptimizationBuilder::MtfGoals OptimizationBuilder::mtf(
    int frequency, const std::vector<double> &sagittal,
    const std::vector<double> &tangential, const std::vector<double> &sagittalWeights,
    const std::vector<double> &tangentialWeights) {
    return MtfGoals(frequency, sagittal, tangential, sagittalWeights, tangentialWeights);
}

OptimizationBuilder::ContrastGoals OptimizationBuilder::contrast(
    int frequency, const std::vector<double> &sagittalWeights,
    const std::vector<double> &tangentialWeights) {
    return ContrastGoals(frequency, sagittalWeights, tangentialWeights);
}

OptimizationBuilder::ContrastGoals OptimizationBuilder::contrast(
    int frequency, const std::vector<double> &weights) {
    return ContrastGoals(frequency, weights, weights);
}

OptimizationBuilder &OptimizationBuilder::mtfGoals(const std::vector<MtfGoals> &goals) {
    configuration->mtfGoals.insert(configuration->mtfGoals.end(), goals.begin(), goals.end());
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastGoals(
    const std::vector<ContrastGoals> &goals) {
    configuration->contrastGoals.insert(configuration->contrastGoals.end(), goals.begin(), goals.end());
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastBalanceGoals(
    const std::vector<bool> &fields_, double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "contrast balance weight must be finite and non-negative");
    configuration->contrastBalanceFields = fields_;
    configuration->contrastBalanceWeight = weight;
    return *this;
}

/**
 * One aggregate GoalSpotRMS per field, each aiming at a target RMS spot
 * radius in microns. To minimize spot size rather than hit a number, prefer
 * #spotDeviationGoals(double...), which takes weights instead.
 */
OptimizationBuilder &OptimizationBuilder::spotRmsGoals(
    const std::vector<double> &targets) {
    configuration->spotRmsGoals = SpotGoals{targets, unitWeightsFor(targets)};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotRmsGoals(
    const std::vector<double> &targets, const std::vector<double> &weights) {
    configuration->spotRmsGoals = SpotGoals{targets, weights};
    return *this;
}

/**
 * Minimize RMS spot radius through the individual signed X/Y ray deviations that
 * make it up, one GoalSpotDeviation per orientation per sampled ray.
 * Differentiating those exposes far more to the solver than one square-rooted
 * aggregate does.
 *
 * These take <em>weights</em>, not targets - every residual aims at zero. One
 * weight per field is applied to every wavelength, sample and orientation.
 */
OptimizationBuilder &OptimizationBuilder::spotDeviationGoals(
    const std::vector<double> &fieldWeights) {
    configuration->addSpotDeviationGoals = true;
    configuration->spotDeviationXWeights = fieldWeights;
    configuration->spotDeviationYWeights = fieldWeights;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotDeviationGoals(
    const std::vector<double> &xWeights, const std::vector<double> &yWeights) {
    configuration->addSpotDeviationGoals = true;
    configuration->spotDeviationXWeights = xWeights;
    configuration->spotDeviationYWeights = yWeights;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotMaxRadiusGoals(
    const std::vector<double> &targets) {
    configuration->spotMaxRadiusGoals = SpotGoals{targets, unitWeightsFor(targets)};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotMaxRadiusGoals(
    const std::vector<double> &targets, const std::vector<double> &weights) {
    configuration->spotMaxRadiusGoals = SpotGoals{targets, weights};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::rayAberrationGoals(bool enabled) {
    configuration->addRayAberrationGoals = enabled;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::paraxialGoal(int paraxId, double target,
                                                       double weight) {
    if (paraxId < 0 || paraxId >= static_cast<int>(std::size(ParaxHelper::Names)))
        throw IllegalArgumentException("unknown paraxial quantity: " + intToString(paraxId));
    if (!std::isfinite(target))
        throw IllegalArgumentException("paraxial target must be finite");
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException("paraxial weight must be finite and non-negative");
    for (const ParaxialGoal &goal : configuration->paraxialGoals)
        if (goal.paraxId == paraxId)
            throw IllegalArgumentException(std::string("there is already a goal for ") +
                                           ParaxHelper::Names[paraxId]);
    configuration->paraxialGoals.push_back(ParaxialGoal{paraxId, target, weight});
    return *this;
}

OptimizationBuilder &OptimizationBuilder::additionalGoals(
    const std::vector<GoalFactory> &factories) {
    for (const auto &factory : factories) {
        if (!factory)
            throw IllegalArgumentException(
                "additional goal factories must not contain null");
        additionalGoalFactories.push_back(factory);
    }
    return *this;
}

// ---------------------------------------------------------------------------
// Constraints
// ---------------------------------------------------------------------------

/**
 * Hold the varied thicknesses near their starting values, at a chosen weight.
 *
 * @param weight relative strength; see #NOMINAL_CONSTRAINT_WEIGHT for why the
 *               nominal value is usually the right one
 */
OptimizationBuilder &OptimizationBuilder::applyThicknessConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "thickness constraint weight must be finite and non-negative");
    configuration->thicknessConstraintWeight = weight;
    return *this;
}

/**
 * Hold the varied gaps near their starting edge separations, at a chosen weight.
 *
 * Gaps whose starting edge separation is not positive and finite are skipped: a
 * fractional constraint cannot be formed around zero, and a design that already starts
 * with coincident or crossed surfaces has nothing useful to anchor to. See
 * ConstraintEdgeThickness#is_constrainable(Analysis, int).
 *
 * @param weight relative strength; see #NOMINAL_CONSTRAINT_WEIGHT
 */
OptimizationBuilder &OptimizationBuilder::applyEdgeThicknessConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "edge thickness constraint weight must be finite and non-negative");
    configuration->edgeThicknessConstraintWeight = weight;
    return *this;
}

/**
 * Hold the varied surfaces near their starting curvatures, at a chosen weight.
 *
 * @param weight relative strength; see #NOMINAL_CONSTRAINT_WEIGHT
 */
OptimizationBuilder &OptimizationBuilder::applyCurvatureConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "curvature constraint weight must be finite and non-negative");
    configuration->curvatureConstraintWeight = weight;
    return *this;
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------

OptimizationBuilder::OptimizationSetup OptimizationBuilder::build() {
    validate();
    auto analysis = std::make_shared<Analysis>(prescription_, *configuration->fields,
                                               *configuration->mtfFrequencies,
                                               configuration->scenario);
    auto variables = buildVariables();
    auto goals = buildGoals(analysis.get(), variables);
    if (goals.size() < variables.size())
        throw IllegalArgumentException(
            "optimization requires at least as many goals as variables: " +
            intToString(static_cast<int>(goals.size())) + " goals for " +
            intToString(static_cast<int>(variables.size())) +
            " variables; add optical goals or enable rayAberrationGoals()");
    bool customMaximumRadius =
        !additionalGoalFactories.empty() && anyGoalIs<GoalSpotMaxRadius>(goals);
    auto effective = configuration->effectiveAnalysis(customMaximumRadius);
    configuration->configureAnalysis(*analysis, effective, !additionalGoalFactories.empty());
    return OptimizationSetup(std::move(analysis), std::move(variables), std::move(goals),
                             configuration->solverTolerances);
}

/**
 * The gaps whose edge separation some varied parameter can move, sorted and
 * deduplicated: the gap a varied thickness <em>is</em>, and <em>both</em> gaps beside
 * a surface whose shape is varied.
 *
 * The second half is easy to miss and was missed here originally, in both the
 * penalty and the bound form. The separation is
 * gap(h) = t + sag_next(h) - sag_this(h), so moving a radius, conic constant
 * or aspheric coefficient closes the gap on either side of that surface with no
 * thickness variable involved anywhere. A setup that varies curvatures and aspherics
 * but no thicknesses therefore got <em>no</em> edge protection at all, silently -
 * which is precisely the configuration in which curvature is the only freedom, and
 * curvature-driven crossing is the failure the edge constraint exists to catch.
 *
 * Out-of-range gap indices produced at either end of the surface list are left in
 * and rejected by the caller's is_constrainable / is_boundable check.
 */
std::set<int> OptimizationBuilder::edgeAffectedGaps(
    const std::vector<std::shared_ptr<Var>> &variables) {
    std::set<int> gaps;
    for (const auto &variable : variables) {
        if (const auto *thickness = dynamic_cast<const VarThickness *>(variable.get())) {
            gaps.insert(thickness->_surface_id);
            continue;
        }
        int surface;
        if (const auto *radius = dynamic_cast<const VarRadius *>(variable.get()))
            surface = radius->_surface_id;
        else if (const auto *conic = dynamic_cast<const VarAsphK *>(variable.get()))
            surface = conic->_surface_id;
        else if (const auto *coefficient =
                     dynamic_cast<const VarAsphCoeff *>(variable.get()))
            surface = coefficient->_surface_id;
        else
            continue;
        gaps.insert(surface - 1);
        gaps.insert(surface);
    }
    return gaps;
}

std::vector<std::shared_ptr<Var>> OptimizationBuilder::buildVariables() const {
    std::vector<std::shared_ptr<Var>> result;
    const auto &surfaces = prescription_->_surface_list;
    if (configuration->allCurvatureSurfaces) {
        for (int surface = 0; surface < static_cast<int>(surfaces.size()); surface++) {
            const auto &definition = surfaces[static_cast<std::size_t>(surface)];
            if (!definition.is_aperture_stop() && !definition.is_field_stop() &&
                definition._radius != 0.0 && !contains(configuration->curvatureExclusions, surface))
                result.push_back(std::make_shared<VarRadius>(prescription_, surface));
        }
    } else {
        for (int surface : configuration->curvatureSurfaces)
            result.push_back(std::make_shared<VarRadius>(prescription_, surface));
    }
    if (configuration->allThicknessSurfaces) {
        for (int surface = 0; surface < static_cast<int>(surfaces.size()); surface++) {
            // A zero thickness is a coincident surface, not a space to open up,
            // and it gives the fractional ConstraintThickness no base to work from.
            if (thicknessOf(surface) != 0.0 && !contains(configuration->thicknessExclusions, surface))
                result.push_back(
                    std::make_shared<VarThickness>(prescription_, surface, configuration->scenario));
        }
    } else {
        for (int surface : configuration->thicknessSurfaces)
            result.push_back(
                std::make_shared<VarThickness>(prescription_, surface, configuration->scenario));
    }
    if (configuration->includeExistingAspherics) {
        for (int surfaceId = 0; surfaceId < static_cast<int>(surfaces.size());
             surfaceId++) {
            if (hasExplicitAsphericTerms(surfaceId))
                continue;
            const auto &surface = surfaces[static_cast<std::size_t>(surfaceId)];
            if (surface._k != 0.0)
                result.push_back(std::make_shared<VarAsphK>(prescription_, surfaceId));
            if (!surface._coeffs.has_value())
                continue;
            for (int coefficient = 0;
                 coefficient < static_cast<int>(surface._coeffs->size()); coefficient++) {
                double value = (*surface._coeffs)[static_cast<std::size_t>(coefficient)];
                if (value != 0.0)
                    result.push_back(std::make_shared<VarAsphCoeff>(
                        prescription_, surfaceId, coefficient, scalingFor(value)));
            }
        }
    }
    auto explicitTerms = explicitAsphericVariables();
    result.insert(result.end(), explicitTerms.begin(), explicitTerms.end());
    result.insert(result.end(), additionalVariables_.begin(), additionalVariables_.end());
    return result;
}

std::vector<std::shared_ptr<Var>> OptimizationBuilder::explicitAsphericVariables() const {
    std::vector<std::shared_ptr<Var>> result;
    for (const AsphericTerm &term : configuration->asphericTerms) {
        auto &surface = prescription_->_surface_list[static_cast<std::size_t>(term.surface)];
        if (!surface.is_aspheric())
            surface._asph_type = asphereTypeOf(term.surface);
        if (!surface._coeffs.has_value())
            surface._coeffs = std::vector<double>();
        if (term.index < 0) {
            result.push_back(std::make_shared<VarAsphK>(prescription_, term.surface));
            continue;
        }
        auto index = static_cast<std::size_t>(term.index);
        if (surface._coeffs->size() <= index)
            surface._coeffs->resize(index + 1, 0.0);
        double value = (*surface._coeffs)[index];
        double scale;
        if (term.scale.has_value())
            scale = *term.scale;
        else if (value != 0.0)
            scale = scalingFor(value);
        else
            // Java's Math.round(double), which is floor(x + 0.5) rather than
            // std::round's half-away-from-zero.
            scale = std::pow(10.0, std::floor(powerOf(surface._asph_type, term.index) *
                                                  std::log10(surface._diameter / 2.0) +
                                              0.5));
        result.push_back(std::make_shared<VarAsphCoeff>(prescription_, term.surface,
                                                        term.index, scale));
    }
    return result;
}

std::vector<std::shared_ptr<Goal>> OptimizationBuilder::buildGoals(
    Analysis *analysis, const std::vector<std::shared_ptr<Var>> &variables) const {
    std::vector<std::shared_ptr<Goal>> result;
    const auto &fields_ = *configuration->fields;
    // Anchor the varied parameters to where they started. Built from the variable
    // list so the goals attach to exactly what is free to move, and built here
    // while the prescription still holds its original values.
    if (configuration->thicknessConstraintWeight.has_value()) {
        for (const auto &variable : variables)
            if (const auto *thickness =
                    dynamic_cast<const VarThickness *>(variable.get()))
                result.push_back(std::make_shared<ConstraintThickness>(
                    analysis, thickness->_surface_id, *configuration->thicknessConstraintWeight));
    }
    if (configuration->edgeThicknessConstraintWeight.has_value()) {
        for (int gap : edgeAffectedGaps(variables))
            if (ConstraintEdgeThickness::is_constrainable(analysis, gap))
                result.push_back(std::make_shared<ConstraintEdgeThickness>(
                    analysis, gap, *configuration->edgeThicknessConstraintWeight));
    }
    if (configuration->curvatureConstraintWeight.has_value()) {
        for (const auto &variable : variables)
            if (const auto *radius = dynamic_cast<const VarRadius *>(variable.get()))
                result.push_back(std::make_shared<ConstraintCurvature>(
                    analysis, radius->_surface_id, *configuration->curvatureConstraintWeight));
    }
    for (const auto &curve : configuration->mtfGoals) {
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            result.push_back(std::make_shared<GoalGeoMTF>(
                analysis, field + 1, Orientation::SAGITTAL, curve.frequency,
                curve.sagittal[f] / 100.0, curve.sagittalWeights[f]));
            result.push_back(std::make_shared<GoalGeoMTF>(
                analysis, field + 1, Orientation::TANGENTIAL, curve.frequency,
                curve.tangential[f] / 100.0, curve.tangentialWeights[f]));
        }
    }

    int contrastSamples = configuration->contrastRings * configuration->contrastSpokes;
    const auto &wvls = prescription_->_wvls;
    const auto &wts = prescription_->_wts;
    for (int contrast_index = 0; contrast_index < static_cast<int>(configuration->contrastGoals.size());
         contrast_index++) {
        const auto &curve = configuration->contrastGoals[static_cast<std::size_t>(contrast_index)];
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                 wavelength++) {
                double wavelengthWeight =
                    configuration->weighted ? wts[static_cast<std::size_t>(wavelength)] : 1.0;
                for (int sample = 0; sample < contrastSamples; sample++) {
                    result.push_back(std::make_shared<GoalContrast>(
                        analysis, contrast_index, curve.frequency, field + 1, wavelength,
                        sample, Orientation::SAGITTAL,
                        wavelengthWeight * curve.sagittalWeights[f]));
                    result.push_back(std::make_shared<GoalContrast>(
                        analysis, contrast_index, curve.frequency, field + 1, wavelength,
                        sample, Orientation::TANGENTIAL,
                        wavelengthWeight * curve.tangentialWeights[f]));
                }
            }
        }
    }

    if (configuration->contrastBalanceFields.has_value()) {
        std::vector<double> wavelengthWeights(wvls.size(), 0.0);
        for (std::size_t w = 0; w < wavelengthWeights.size(); w++)
            wavelengthWeights[w] = configuration->weighted ? wts[w] : 1.0;
        for (int contrast_index = 0;
             contrast_index < static_cast<int>(configuration->contrastGoals.size()); contrast_index++) {
            const auto &curve = configuration->contrastGoals[static_cast<std::size_t>(contrast_index)];
            for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
                auto f = static_cast<std::size_t>(field);
                if (!(*configuration->contrastBalanceFields)[f])
                    continue;
                result.push_back(std::make_shared<GoalContrastBalance>(
                    analysis, contrast_index, curve.frequency, field + 1,
                    wavelengthWeights, curve.sagittalWeights[f],
                    curve.tangentialWeights[f], configuration->contrastBalanceWeight));
            }
        }
    }

    if (configuration->spotRmsGoals.has_value()) {
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            result.push_back(std::make_shared<GoalSpotRMS>(
                analysis, field + 1, configuration->spotRmsGoals->targets[f], configuration->spotRmsGoals->weights[f]));
        }
    }
    if (configuration->addSpotDeviationGoals) {
        int samples = configuration->gaussianQuadratureRings * configuration->gaussianQuadratureSpokes;
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                 wavelength++) {
                double wavelengthWeight =
                    configuration->weighted ? wts[static_cast<std::size_t>(wavelength)] : 1.0;
                for (int sample = 0; sample < samples; sample++) {
                    result.push_back(std::make_shared<GoalSpotDeviation>(
                        analysis, field + 1, wavelength, sample, Orientation::X,
                        wavelengthWeight * (*configuration->spotDeviationXWeights)[f]));
                    result.push_back(std::make_shared<GoalSpotDeviation>(
                        analysis, field + 1, wavelength, sample, Orientation::Y,
                        wavelengthWeight * (*configuration->spotDeviationYWeights)[f]));
                }
            }
        }
    }
    if (configuration->spotMaxRadiusGoals.has_value()) {
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            result.push_back(std::make_shared<GoalSpotMaxRadius>(
                analysis, field + 1, configuration->spotMaxRadiusGoals->targets[f],
                configuration->spotMaxRadiusGoals->weights[f]));
        }
    }

    // Anchor first-order properties to the requested prescription values, unless the
    // caller has set targets of its own.
    result.push_back(anchor(analysis, ParaxHelper::Effective_focal_length, focalLengthOf()));
    result.push_back(anchor(analysis, ParaxHelper::Fno, fNumberOf()));

    if (configuration->addRayAberrationGoals) {
        for (int field = 1; field <= static_cast<int>(fields_.size()); field++) {
            for (int orientation = Orientation::SAGITTAL;
                 orientation <= Orientation::TANGENTIAL; orientation++) {
                for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                     wavelength++) {
                    auto w = static_cast<std::size_t>(wavelength);
                    if (configuration->dLineOnly && !sameWavelength(wvls[w], Glass::d))
                        continue;
                    double weight = configuration->weighted ? wts[w] : 1.0;
                    for (int sample = 0; sample < RAY_FAN_SAMPLES; sample++)
                        result.push_back(std::make_shared<GoalRayAberration>(
                            analysis, field, orientation, sample, wvls[w], 0.0, weight));
                }
            }
        }
    }
    for (const ParaxialGoal &goal : configuration->paraxialGoals)
        if (goal.paraxId != ParaxHelper::Effective_focal_length &&
            goal.paraxId != ParaxHelper::Fno)
            result.push_back(std::make_shared<GoalParax>(analysis, goal.paraxId, goal.target,
                                                         goal.weight));
    for (const auto &factory : additionalGoalFactories) {
        auto goal = factory(analysis);
        if (goal == nullptr)
            throw IllegalArgumentException("an additional goal factory returned null");
        if (goal->_analysis != analysis)
            throw IllegalArgumentException(
                "additional goals must use the Analysis supplied to their factory");
        result.push_back(std::move(goal));
    }
    return result;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

void OptimizationBuilder::validate() const {
    if (!configuration->fields.has_value() || configuration->fields->empty())
        throw IllegalArgumentException("at least one field is required");
    validateScenario();
    const auto &fields_ = *configuration->fields;
    int fieldCount = static_cast<int>(fields_.size());
    if (configuration->contrastBalanceFields.has_value()) {
        if (configuration->contrastBalanceFields->size() != fields_.size())
            throw IllegalArgumentException(
                "contrast balance needs one flag per field: " + intToString(fieldCount) +
                " fields but " +
                intToString(static_cast<int>(configuration->contrastBalanceFields->size())) +
                " flags");
        if (configuration->contrastGoals.empty())
            throw IllegalArgumentException(
                "contrast balance goals require contrast goals to balance");
    }
    OptimizationValidation::range(fields_, 1.0, [] {
        return IllegalArgumentException("fields must be finite values between 0 and 1");
    });
    if (fields_[0] != 0.0)
        throw IllegalArgumentException("the first field must be 0.0");

    if (!configuration->mtfFrequencies.has_value() || configuration->mtfFrequencies->empty())
        throw IllegalArgumentException("at least one MTF frequency is required");
    std::set<int> frequencies;
    for (int frequency : *configuration->mtfFrequencies)
        OptimizationValidation::positiveUnique(frequency, frequencies, [] {
            return IllegalArgumentException("MTF frequencies must be positive and unique");
        });
    std::set<int> goalFrequencies;
    for (const auto &curve : configuration->mtfGoals) {
        OptimizationValidation::frequency(curve.frequency, *configuration->mtfFrequencies, [&] {
            return IllegalArgumentException(
                "MTF goal frequency was not requested for measurement: " +
                intToString(curve.frequency));
        });
        if (!goalFrequencies.insert(curve.frequency).second)
            throw IllegalArgumentException("duplicate MTF goal frequency: " +
                                           intToString(curve.frequency));
        curve.validate(fieldCount);
    }
    std::set<int> contrastFrequencies;
    for (const auto &curve : configuration->contrastGoals) {
        OptimizationValidation::positiveUnique(curve.frequency, contrastFrequencies, [] {
            return IllegalArgumentException("contrast frequencies must be positive and unique");
        });
        curve.validate(fieldCount);
    }
    if (configuration->spotRmsGoals.has_value())
        configuration->spotRmsGoals->validate(fieldCount, "spot RMS");
    if (configuration->addSpotDeviationGoals) {
        MtfGoals::validateWeights(
            configuration->spotDeviationXWeights.value_or(std::vector<double>()), fieldCount,
            "spot deviation X weights");
        MtfGoals::validateWeights(
            configuration->spotDeviationYWeights.value_or(std::vector<double>()), fieldCount,
            "spot deviation Y weights");
        if (configuration->spotRmsGoals.has_value())
            throw IllegalArgumentException("aggregate spot RMS goals and per-ray spot "
                                           "deviation goals cannot both be enabled");
        if (configuration->spotMaxRadiusGoals.has_value() ||
            configuration->useHexapolarSpotPattern)
            throw IllegalArgumentException(
                "spot deviation goals require Gaussian-quadrature spot sampling");
    }
    if (configuration->spotMaxRadiusGoals.has_value())
        configuration->spotMaxRadiusGoals->validate(fieldCount, "spot maximum radius");
    validateSurfaces(configuration->curvatureSurfaces, "curvature");
    validateSurfaces(configuration->thicknessSurfaces, "thickness");
    validateSurfaces(configuration->curvatureExclusions, "excluded curvature");
    validateSurfaces(configuration->thicknessExclusions, "excluded thickness");
    if (configuration->addRayAberrationGoals && configuration->dLineOnly) {
        bool any = false;
        for (double w : prescription_->_wvls)
            if (sameWavelength(w, Glass::d))
                any = true;
        if (!any)
            throw IllegalArgumentException("d-line optimization requires the "
                                           "prescription to contain the d-line wavelength");
    }
}

void OptimizationBuilder::validateSurfaces(const std::vector<int> &surfaces,
                                           const char *kind) const {
    std::set<int> seen;
    for (int surface : surfaces) {
        if (surface < 0 || surface >= static_cast<int>(prescription_->_surface_list.size()))
            throw IllegalArgumentException(std::string(kind) +
                                           " surface is out of range: " +
                                           intToString(surface));
        if (!seen.insert(surface).second)
            throw IllegalArgumentException("duplicate " + std::string(kind) +
                                           " surface: " + intToString(surface));
    }
}

std::shared_ptr<Goal> OptimizationBuilder::anchor(Analysis *analysis, int paraxId,
                                                 double prescribed) const {
    for (const ParaxialGoal &goal : configuration->paraxialGoals)
        if (goal.paraxId == paraxId)
            return std::make_shared<GoalParax>(analysis, paraxId, goal.target, goal.weight);
    return std::make_shared<GoalParax>(analysis, paraxId, prescribed, 1.0);
}

double OptimizationBuilder::scalingFor(double value) {
    double magnitude = std::abs(value);
    return std::pow(10.0, -std::floor(std::log10(magnitude)));
}

bool OptimizationBuilder::sameWavelength(double a, double b) {
    return std::abs(a - b) < 1.0e-3;
}

// ---------------------------------------------------------------------------
// Nested value types
// ---------------------------------------------------------------------------

OptimizationBuilder::MtfGoals::MtfGoals(
    int frequency_, const std::vector<double> &sagittal_,
    const std::vector<double> &tangential_,
    const std::optional<std::vector<double>> &sagittalWeights_,
    const std::optional<std::vector<double>> &tangentialWeights_)
    : frequency(frequency_), sagittal(sagittal_), tangential(tangential_),
      sagittalWeights(sagittalWeights_.has_value() ? *sagittalWeights_
                                                   : unitWeightsFor(sagittal_)),
      tangentialWeights(tangentialWeights_.has_value() ? *tangentialWeights_
                                                       : unitWeightsFor(tangential_)) {}

void OptimizationBuilder::MtfGoals::validate(int fieldCount) const {
    validateTargets(sagittal, fieldCount, "sagittal targets");
    validateTargets(tangential, fieldCount, "tangential targets");
    validateWeights(sagittalWeights, fieldCount, "sagittal weights");
    validateWeights(tangentialWeights, fieldCount, "tangential weights");
}

void OptimizationBuilder::MtfGoals::validateTargets(const std::vector<double> &values,
                                                    int count, const char *name) {
    OptimizationValidation::fieldCount(values, count, [&] {
        return IllegalArgumentException(std::string(name) +
                                        " must contain one value per field");
    });
    OptimizationValidation::range(values, 100.0, [&] {
        return IllegalArgumentException(std::string(name) +
                                        " must be percentages between 0 and 100");
    });
}

void OptimizationBuilder::MtfGoals::validateWeights(const std::vector<double> &values,
                                                    int count, const char *name) {
    OptimizationValidation::fieldCount(values, count, [&] {
        return IllegalArgumentException(std::string(name) +
                                        " must contain one value per field");
    });
    OptimizationValidation::range(values, std::numeric_limits<double>::infinity(), [&] {
        return IllegalArgumentException(std::string(name) +
                                        " must be finite and non-negative");
    });
}

void OptimizationBuilder::ContrastGoals::validate(int fieldCount) const {
    MtfGoals::validateWeights(sagittalWeights, fieldCount, "sagittal contrast weights");
    MtfGoals::validateWeights(tangentialWeights, fieldCount,
                              "tangential contrast weights");
}

void OptimizationBuilder::SpotGoals::validate(int fieldCount, const char *name) const {
    OptimizationValidation::fieldCount(targets, fieldCount, [&] {
        return IllegalArgumentException(std::string(name) +
                                        " targets must contain one value per field");
    });
    OptimizationValidation::fieldCount(weights, fieldCount, [&] {
        return IllegalArgumentException(std::string(name) +
                                        " weights must contain one value per field");
    });
    OptimizationValidation::range(targets, std::numeric_limits<double>::infinity(), [&] {
        return IllegalArgumentException(std::string(name) +
                                        " targets must be finite and non-negative");
    });
    OptimizationValidation::range(weights, std::numeric_limits<double>::infinity(), [&] {
        return IllegalArgumentException(std::string(name) +
                                        " weights must be finite and non-negative");
    });
}

// ---------------------------------------------------------------------------
// Writing - the setup as a [trial n] section
// ---------------------------------------------------------------------------

std::string OptimizationBuilder::toTrial(int number) const {
    if (!additionalVariables_.empty() || !additionalGoalFactories.empty())
        throw IllegalStateException("variables and goals added as code have no written form, "
                                    "so this setup cannot be written as a trial");
    return configuration->toTrial(number);
}

} // namespace redukti::optim
