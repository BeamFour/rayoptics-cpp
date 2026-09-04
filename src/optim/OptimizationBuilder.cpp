// C++ port of org.redukti.optim.OptimizationBuilder
#include "redukti/optim/OptimizationBuilder.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace redukti::optim {

namespace Orientation = rayoptics::util::Orientation;
using rayoptics::seq::Glass;

namespace {

/** Java's `weights == null ? unitWeights(targets) : copy(weights)`. */
std::vector<double> unitWeightsFor(const std::vector<double> &targets) {
    return std::vector<double>(targets.size(), 1.0);
}

/** True when any goal in the list is a T; the Java uses `instanceof`. */
template <typename T> bool anyGoalIs(const std::vector<std::shared_ptr<Goal>> &goals) {
    for (const auto &goal : goals)
        if (dynamic_cast<const T *>(goal.get()) != nullptr)
            return true;
    return false;
}

} // namespace

OptimizationBuilder::OptimizationBuilder(spec::Prescription *prescription_)
    : prescription(prescription_) {
    if (prescription_ == nullptr)
        throw IllegalArgumentException("prescription must not be null");
    // The Java tests `prescription._surfaces == null`, the array build() fills;
    // here the flag records the same thing without a second copy of the list.
    if (!prescription_->_built)
        throw IllegalArgumentException("prescription must be built before optimization");
}

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

OptimizationBuilder &OptimizationBuilder::fields(const std::vector<double> &fields_) {
    this->_fields = fields_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::mtfFrequencies(
    const std::vector<int> &frequencies) {
    this->_mtfFrequencies = frequencies;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::weighted(bool weighted_) {
    this->_weighted = weighted_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::dLineOnly(bool dLineOnly_) {
    this->_dLineOnly = dLineOnly_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::scenario(int scenario_) {
    if (scenario_ < 0)
        throw IllegalArgumentException("scenario must be non-negative, got " +
                                       intToString(scenario_));
    this->_scenario = scenario_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::vignetting(spec::VigType vigType_) {
    this->vigType = vigType_;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::freezeVignetting(bool freeze) {
    this->freezeVignetting_ = freeze;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::gaussianQuadratureSampling(
    int rings, int spokes, double innerPupilRadius) {
    if (rings < 1 || spokes < 3)
        throw IllegalArgumentException(
            "Gaussian quadrature requires at least 1 ring and 3 spokes");
    if (!std::isfinite(innerPupilRadius) || innerPupilRadius < 0.0 ||
        innerPupilRadius >= 1.0)
        throw IllegalArgumentException("Inner pupil radius must be finite and in [0, 1)");
    this->gaussianQuadratureRings = rings;
    this->gaussianQuadratureSpokes = spokes;
    this->gaussianQuadratureInnerRadius = innerPupilRadius;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::checkSpotApertures(bool check) {
    this->_checkSpotApertures = check;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::hexapolarSampling(int numRays) {
    if (numRays < 1)
        throw IllegalArgumentException("hexapolar spot rays must be at least 1");
    this->useHexapolarSpotPattern = true;
    this->hexapolarSpotRays = numRays;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastSampling(int rings, int spokes) {
    if (rings < 1 || spokes < 1)
        throw IllegalArgumentException("contrast rings and spokes must be at least 1");
    contrastRings = rings;
    contrastSpokes = spokes;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::calibrateContrastFrequency(bool value) {
    calibrateContrastFrequency_ = value;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::aimContrastAtExitPupil(bool value) {
    aimContrastAtExitPupil_ = value;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::centerContrastResiduals(bool value) {
    centerContrastResiduals_ = value;
    return *this;
}

// ---------------------------------------------------------------------------
// Variables
// ---------------------------------------------------------------------------

OptimizationBuilder &OptimizationBuilder::varyCurvatures(
    const std::vector<int> &surfaces) {
    this->curvatureSurfaces = surfaces;
    this->allCurvatureSurfaces = false;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyAllCurvatures() {
    this->allCurvatureSurfaces = true;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyThicknesses(
    const std::vector<int> &surfaces) {
    this->thicknessSurfaces = surfaces;
    this->allThicknessSurfaces = false;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyAllThicknesses() {
    this->allThicknessSurfaces = true;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::varyExistingAspherics(bool include) {
    this->includeExistingAspherics = include;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::additionalVariables(
    const std::vector<std::shared_ptr<Var>> &variables) {
    for (const auto &variable : variables) {
        if (variable == nullptr)
            throw IllegalArgumentException("additional variables must not contain null");
        if (variable->_prescription != prescription)
            throw IllegalArgumentException(
                "additional variables must use this builder's prescription");
        additionalVariables_.push_back(variable);
    }
    return *this;
}

double OptimizationBuilder::thicknessOf(int surface) const {
    const auto &definition =
        prescription->_surface_list[static_cast<std::size_t>(surface)];
    return definition._thickness_by_scenario.has_value()
               ? (*definition._thickness_by_scenario)[static_cast<std::size_t>(_scenario)]
               : definition._thickness;
}

double OptimizationBuilder::focalLengthOf() const {
    // The Java field is a nullable array; here an empty vector is the same state.
    return !prescription->_focal_length_by_scenario.empty()
               ? prescription
                     ->_focal_length_by_scenario[static_cast<std::size_t>(_scenario)]
               : prescription->_focal_length;
}

double OptimizationBuilder::fNumberOf() const {
    return !prescription->_f_number_by_scenario.empty()
               ? prescription->_f_number_by_scenario[static_cast<std::size_t>(_scenario)]
               : prescription->_fno;
}

void OptimizationBuilder::validateScenario() const {
    if (_scenario == 0)
        return;
    int available = scenarioCount();
    if (_scenario >= available)
        throw IllegalArgumentException(
            "scenario " + intToString(_scenario) +
            " requested but the prescription defines " + intToString(available) +
            (available == 1 ? " (it is not multi-configuration)" : ""));
}

int OptimizationBuilder::scenarioCount() const {
    int count = 1;
    if (!prescription->_focal_length_by_scenario.empty())
        count = std::max(
            count, static_cast<int>(prescription->_focal_length_by_scenario.size()));
    for (const auto &surface : prescription->_surface_list)
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
    _mtfGoals.insert(_mtfGoals.end(), goals.begin(), goals.end());
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastGoals(
    const std::vector<ContrastGoals> &goals) {
    _contrastGoals.insert(_contrastGoals.end(), goals.begin(), goals.end());
    return *this;
}

OptimizationBuilder &OptimizationBuilder::contrastBalanceGoals(
    const std::vector<bool> &fields_, double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "contrast balance weight must be finite and non-negative");
    this->contrastBalanceFields = fields_;
    this->contrastBalanceWeight = weight;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotRmsGoals(
    const std::vector<double> &targets) {
    spotRmsGoals_ = SpotGoals{targets, unitWeightsFor(targets)};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotRmsGoals(
    const std::vector<double> &targets, const std::vector<double> &weights) {
    spotRmsGoals_ = SpotGoals{targets, weights};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotDeviationGoals(
    const std::vector<double> &fieldWeights) {
    this->addSpotDeviationGoals = true;
    this->spotDeviationXWeights = fieldWeights;
    this->spotDeviationYWeights = fieldWeights;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotDeviationGoals(
    const std::vector<double> &xWeights, const std::vector<double> &yWeights) {
    this->addSpotDeviationGoals = true;
    this->spotDeviationXWeights = xWeights;
    this->spotDeviationYWeights = yWeights;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotMaxRadiusGoals(
    const std::vector<double> &targets) {
    spotMaxRadiusGoals_ = SpotGoals{targets, unitWeightsFor(targets)};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::spotMaxRadiusGoals(
    const std::vector<double> &targets, const std::vector<double> &weights) {
    spotMaxRadiusGoals_ = SpotGoals{targets, weights};
    return *this;
}

OptimizationBuilder &OptimizationBuilder::rayAberrationGoals(bool enabled) {
    this->addRayAberrationGoals = enabled;
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

OptimizationBuilder &OptimizationBuilder::applyThicknessConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "thickness constraint weight must be finite and non-negative");
    this->thicknessConstraintWeight = weight;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::applyEdgeThicknessConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "edge thickness constraint weight must be finite and non-negative");
    this->edgeThicknessConstraintWeight = weight;
    return *this;
}

OptimizationBuilder &OptimizationBuilder::applyCurvatureConstraints(double weight) {
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException(
            "curvature constraint weight must be finite and non-negative");
    this->curvatureConstraintWeight = weight;
    return *this;
}

// ---------------------------------------------------------------------------
// Build
// ---------------------------------------------------------------------------

OptimizationBuilder::OptimizationSetup OptimizationBuilder::build() {
    validate();
    auto analysis =
        std::make_shared<Analysis>(prescription, *_fields, *_mtfFrequencies, _scenario);
    auto variables = buildVariables();
    auto goals = buildGoals(analysis.get(), variables);
    if (goals.size() < variables.size())
        throw IllegalArgumentException(
            "optimization requires at least as many goals as variables: " +
            intToString(static_cast<int>(goals.size())) + " goals for " +
            intToString(static_cast<int>(variables.size())) +
            " variables; add optical goals or enable rayAberrationGoals()");
    analysis->vignetting(vigType)
        .freezing_vignetting(freezeVignetting_)
        .checking_spot_apertures(_checkSpotApertures);
    configureSpotPattern(*analysis, goals);
    configureContrastAnalysis(*analysis);
    configureRequiredAnalyses(*analysis, goals);
    return OptimizationSetup(std::move(analysis), std::move(variables), std::move(goals));
}

void OptimizationBuilder::configureContrastAnalysis(Analysis &analysis) const {
    if (_contrastGoals.empty())
        return;
    if (calibrateContrastFrequency_ && aimContrastAtExitPupil_)
        throw IllegalArgumentException("Contrast frequency calibration and exit-pupil "
                                       "aiming are mutually exclusive");
    std::vector<int> frequencies;
    frequencies.reserve(_contrastGoals.size());
    for (const auto &goal : _contrastGoals)
        frequencies.push_back(goal.frequency);
    analysis.using_contrast_analysis(frequencies, contrastRings, contrastSpokes);
    analysis.calibrating_contrast_frequency(calibrateContrastFrequency_);
    analysis.aiming_contrast_at_exit_pupil(aimContrastAtExitPupil_);
    analysis.centering_contrast_residuals(centerContrastResiduals_);
}

void OptimizationBuilder::configureRequiredAnalyses(
    Analysis &analysis, const std::vector<std::shared_ptr<Goal>> &goals) const {
    // Additional goal factories are conservatively assumed to require all analyses.
    if (additionalGoalFactories.empty()) {
        bool spots = anyGoalIs<GoalSpotRMS>(goals) ||
                     anyGoalIs<GoalSpotDeviation>(goals) ||
                     anyGoalIs<GoalSpotMaxRadius>(goals) || anyGoalIs<GoalGeoMTF>(goals);
        bool mtf = anyGoalIs<GoalGeoMTF>(goals);
        bool rayAberrations =
            anyGoalIs<GoalRayAberration>(goals) || anyGoalIs<GoalMTFProxy>(goals);
        analysis.required_analyses(spots, rayAberrations, mtf);
    }
}

void OptimizationBuilder::configureSpotPattern(
    Analysis &analysis, const std::vector<std::shared_ptr<Goal>> &goals) const {
    bool hasSpotMaxRadiusGoal = anyGoalIs<GoalSpotMaxRadius>(goals);
    if (addSpotDeviationGoals) {
        analysis
            .using_gauss_quadrature_pattern(gaussianQuadratureRings,
                                            gaussianQuadratureSpokes,
                                            gaussianQuadratureInnerRadius)
            .retaining_failed_spot_rays(true);
    } else if (useHexapolarSpotPattern || hasSpotMaxRadiusGoal) {
        analysis.using_hexapolar_pattern(hexapolarSpotRays);
    } else {
        analysis.using_gauss_quadrature_pattern(gaussianQuadratureRings,
                                                gaussianQuadratureSpokes,
                                                gaussianQuadratureInnerRadius);
    }
}

std::vector<std::shared_ptr<Var>> OptimizationBuilder::buildVariables() const {
    std::vector<std::shared_ptr<Var>> result;
    const auto &surfaces = prescription->_surface_list;
    if (allCurvatureSurfaces) {
        for (int surface = 0; surface < static_cast<int>(surfaces.size()); surface++) {
            const auto &definition = surfaces[static_cast<std::size_t>(surface)];
            if (!definition.is_aperture_stop() && !definition.is_field_stop() &&
                definition._radius != 0.0)
                result.push_back(std::make_shared<VarRadius>(prescription, surface));
        }
    } else {
        for (int surface : curvatureSurfaces)
            result.push_back(std::make_shared<VarRadius>(prescription, surface));
    }
    if (allThicknessSurfaces) {
        for (int surface = 0; surface < static_cast<int>(surfaces.size()); surface++) {
            // A zero thickness is a coincident surface, not a space to open up,
            // and it gives the fractional ConstraintThickness no base to work from.
            if (thicknessOf(surface) != 0.0)
                result.push_back(
                    std::make_shared<VarThickness>(prescription, surface, _scenario));
        }
    } else {
        for (int surface : thicknessSurfaces)
            result.push_back(
                std::make_shared<VarThickness>(prescription, surface, _scenario));
    }
    if (includeExistingAspherics) {
        for (int surfaceId = 0; surfaceId < static_cast<int>(surfaces.size());
             surfaceId++) {
            const auto &surface = surfaces[static_cast<std::size_t>(surfaceId)];
            if (surface._k != 0.0)
                result.push_back(std::make_shared<VarAsphK>(prescription, surfaceId));
            if (!surface._coeffs.has_value())
                continue;
            for (int coefficient = 0;
                 coefficient < static_cast<int>(surface._coeffs->size()); coefficient++) {
                double value = (*surface._coeffs)[static_cast<std::size_t>(coefficient)];
                if (value != 0.0)
                    result.push_back(std::make_shared<VarAsphCoeff>(
                        prescription, surfaceId, coefficient, scalingFor(value)));
            }
        }
    }
    result.insert(result.end(), additionalVariables_.begin(), additionalVariables_.end());
    return result;
}

std::vector<std::shared_ptr<Goal>> OptimizationBuilder::buildGoals(
    Analysis *analysis, const std::vector<std::shared_ptr<Var>> &variables) const {
    std::vector<std::shared_ptr<Goal>> result;
    const auto &fields_ = *_fields;
    // Anchor the varied parameters to where they started. Built from the variable
    // list so the goals attach to exactly what is free to move, and built here
    // while the prescription still holds its original values.
    if (thicknessConstraintWeight.has_value()) {
        for (const auto &variable : variables)
            if (const auto *thickness =
                    dynamic_cast<const VarThickness *>(variable.get()))
                result.push_back(std::make_shared<ConstraintThickness>(
                    analysis, thickness->_surface_id, *thicknessConstraintWeight));
    }
    if (edgeThicknessConstraintWeight.has_value()) {
        for (const auto &variable : variables) {
            const auto *thickness = dynamic_cast<const VarThickness *>(variable.get());
            if (thickness != nullptr &&
                ConstraintEdgeThickness::is_constrainable(analysis,
                                                          thickness->_surface_id))
                result.push_back(std::make_shared<ConstraintEdgeThickness>(
                    analysis, thickness->_surface_id, *edgeThicknessConstraintWeight));
        }
    }
    if (curvatureConstraintWeight.has_value()) {
        for (const auto &variable : variables)
            if (const auto *radius = dynamic_cast<const VarRadius *>(variable.get()))
                result.push_back(std::make_shared<ConstraintCurvature>(
                    analysis, radius->_surface_id, *curvatureConstraintWeight));
    }
    for (const auto &curve : _mtfGoals) {
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

    int contrastSamples = contrastRings * contrastSpokes;
    const auto &wvls = prescription->_wvls;
    const auto &wts = prescription->_wts;
    for (int contrast_index = 0; contrast_index < static_cast<int>(_contrastGoals.size());
         contrast_index++) {
        const auto &curve = _contrastGoals[static_cast<std::size_t>(contrast_index)];
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                 wavelength++) {
                double wavelengthWeight =
                    _weighted ? wts[static_cast<std::size_t>(wavelength)] : 1.0;
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

    if (contrastBalanceFields.has_value()) {
        std::vector<double> wavelengthWeights(wvls.size(), 0.0);
        for (std::size_t w = 0; w < wavelengthWeights.size(); w++)
            wavelengthWeights[w] = _weighted ? wts[w] : 1.0;
        for (int contrast_index = 0;
             contrast_index < static_cast<int>(_contrastGoals.size()); contrast_index++) {
            const auto &curve = _contrastGoals[static_cast<std::size_t>(contrast_index)];
            for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
                auto f = static_cast<std::size_t>(field);
                if (!(*contrastBalanceFields)[f])
                    continue;
                result.push_back(std::make_shared<GoalContrastBalance>(
                    analysis, contrast_index, curve.frequency, field + 1,
                    wavelengthWeights, curve.sagittalWeights[f],
                    curve.tangentialWeights[f], contrastBalanceWeight));
            }
        }
    }

    if (spotRmsGoals_.has_value()) {
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            result.push_back(std::make_shared<GoalSpotRMS>(
                analysis, field + 1, spotRmsGoals_->targets[f], spotRmsGoals_->weights[f]));
        }
    }
    if (addSpotDeviationGoals) {
        int samples = gaussianQuadratureRings * gaussianQuadratureSpokes;
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                 wavelength++) {
                double wavelengthWeight =
                    _weighted ? wts[static_cast<std::size_t>(wavelength)] : 1.0;
                for (int sample = 0; sample < samples; sample++) {
                    result.push_back(std::make_shared<GoalSpotDeviation>(
                        analysis, field + 1, wavelength, sample, Orientation::X,
                        wavelengthWeight * (*spotDeviationXWeights)[f]));
                    result.push_back(std::make_shared<GoalSpotDeviation>(
                        analysis, field + 1, wavelength, sample, Orientation::Y,
                        wavelengthWeight * (*spotDeviationYWeights)[f]));
                }
            }
        }
    }
    if (spotMaxRadiusGoals_.has_value()) {
        for (int field = 0; field < static_cast<int>(fields_.size()); field++) {
            auto f = static_cast<std::size_t>(field);
            result.push_back(std::make_shared<GoalSpotMaxRadius>(
                analysis, field + 1, spotMaxRadiusGoals_->targets[f],
                spotMaxRadiusGoals_->weights[f]));
        }
    }

    // Anchor first-order properties to the requested prescription values.
    result.push_back(std::make_shared<GoalParax>(
        analysis, ParaxHelper::Effective_focal_length, focalLengthOf(), 1.0));
    result.push_back(
        std::make_shared<GoalParax>(analysis, ParaxHelper::Fno, fNumberOf(), 1.0));

    if (addRayAberrationGoals) {
        for (int field = 1; field <= static_cast<int>(fields_.size()); field++) {
            for (int orientation = Orientation::SAGITTAL;
                 orientation <= Orientation::TANGENTIAL; orientation++) {
                for (int wavelength = 0; wavelength < static_cast<int>(wvls.size());
                     wavelength++) {
                    auto w = static_cast<std::size_t>(wavelength);
                    if (_dLineOnly && !sameWavelength(wvls[w], Glass::d))
                        continue;
                    double weight = _weighted ? wts[w] : 1.0;
                    for (int sample = 0; sample < RAY_FAN_SAMPLES; sample++)
                        result.push_back(std::make_shared<GoalRayAberration>(
                            analysis, field, orientation, sample, wvls[w], 0.0, weight));
                }
            }
        }
    }
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
    if (!_fields.has_value() || _fields->empty())
        throw IllegalArgumentException("at least one field is required");
    validateScenario();
    const auto &fields_ = *_fields;
    if (contrastBalanceFields.has_value()) {
        if (contrastBalanceFields->size() != fields_.size())
            throw IllegalArgumentException(
                "contrast balance needs one flag per field: " +
                intToString(static_cast<int>(fields_.size())) + " fields but " +
                intToString(static_cast<int>(contrastBalanceFields->size())) + " flags");
        if (_contrastGoals.empty())
            throw IllegalArgumentException(
                "contrast balance goals require contrast goals to balance");
    }
    for (double field : fields_)
        if (!std::isfinite(field) || field < 0.0 || field > 1.0)
            throw IllegalArgumentException("fields must be finite values between 0 and 1");
    if (fields_[0] != 0.0)
        throw IllegalArgumentException("the first field must be 0.0");

    if (!_mtfFrequencies.has_value() || _mtfFrequencies->empty())
        throw IllegalArgumentException("at least one MTF frequency is required");
    std::set<int> frequencies;
    for (int frequency : *_mtfFrequencies) {
        if (frequency <= 0 || !frequencies.insert(frequency).second)
            throw IllegalArgumentException("MTF frequencies must be positive and unique");
    }
    std::set<int> goalFrequencies;
    for (const auto &curve : _mtfGoals) {
        if (frequencies.count(curve.frequency) == 0)
            throw IllegalArgumentException(
                "MTF goal frequency was not requested for measurement: " +
                intToString(curve.frequency));
        if (!goalFrequencies.insert(curve.frequency).second)
            throw IllegalArgumentException("duplicate MTF goal frequency: " +
                                           intToString(curve.frequency));
        curve.validate(static_cast<int>(fields_.size()));
    }
    std::set<int> contrastFrequencies;
    for (const auto &curve : _contrastGoals) {
        if (curve.frequency <= 0 || !contrastFrequencies.insert(curve.frequency).second)
            throw IllegalArgumentException(
                "contrast frequencies must be positive and unique");
        curve.validate(static_cast<int>(fields_.size()));
    }
    if (spotRmsGoals_.has_value())
        spotRmsGoals_->validate(static_cast<int>(fields_.size()), "spot RMS");
    if (addSpotDeviationGoals) {
        MtfGoals::validateWeights(spotDeviationXWeights.value_or(std::vector<double>()),
                                  static_cast<int>(fields_.size()),
                                  "spot deviation X weights");
        MtfGoals::validateWeights(spotDeviationYWeights.value_or(std::vector<double>()),
                                  static_cast<int>(fields_.size()),
                                  "spot deviation Y weights");
        if (spotRmsGoals_.has_value())
            throw IllegalArgumentException("aggregate spot RMS goals and per-ray spot "
                                           "deviation goals cannot both be enabled");
        if (spotMaxRadiusGoals_.has_value() || useHexapolarSpotPattern)
            throw IllegalArgumentException(
                "spot deviation goals require Gaussian-quadrature spot sampling");
    }
    if (spotMaxRadiusGoals_.has_value())
        spotMaxRadiusGoals_->validate(static_cast<int>(fields_.size()),
                                      "spot maximum radius");
    validateSurfaces(curvatureSurfaces, "curvature");
    validateSurfaces(thicknessSurfaces, "thickness");
    if (addRayAberrationGoals && _dLineOnly) {
        bool any = false;
        for (double w : prescription->_wvls)
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
        if (surface < 0 || surface >= static_cast<int>(prescription->_surface_list.size()))
            throw IllegalArgumentException(std::string(kind) +
                                           " surface is out of range: " +
                                           intToString(surface));
        if (!seen.insert(surface).second)
            throw IllegalArgumentException("duplicate " + std::string(kind) +
                                           " surface: " + intToString(surface));
    }
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
    if (static_cast<int>(values.size()) != count)
        throw IllegalArgumentException(std::string(name) +
                                       " must contain one value per field");
    for (double value : values)
        if (!std::isfinite(value) || value < 0.0 || value > 100.0)
            throw IllegalArgumentException(std::string(name) +
                                           " must be percentages between 0 and 100");
}

void OptimizationBuilder::MtfGoals::validateWeights(const std::vector<double> &values,
                                                    int count, const char *name) {
    if (static_cast<int>(values.size()) != count)
        throw IllegalArgumentException(std::string(name) +
                                       " must contain one value per field");
    for (double value : values)
        if (!std::isfinite(value) || value < 0.0)
            throw IllegalArgumentException(std::string(name) +
                                           " must be finite and non-negative");
}

void OptimizationBuilder::ContrastGoals::validate(int fieldCount) const {
    MtfGoals::validateWeights(sagittalWeights, fieldCount, "sagittal contrast weights");
    MtfGoals::validateWeights(tangentialWeights, fieldCount,
                              "tangential contrast weights");
}

void OptimizationBuilder::SpotGoals::validate(int fieldCount, const char *name) const {
    if (static_cast<int>(targets.size()) != fieldCount)
        throw IllegalArgumentException(std::string(name) +
                                       " targets must contain one value per field");
    if (static_cast<int>(weights.size()) != fieldCount)
        throw IllegalArgumentException(std::string(name) +
                                       " weights must contain one value per field");
    for (double target : targets)
        if (!std::isfinite(target) || target < 0.0)
            throw IllegalArgumentException(std::string(name) +
                                           " targets must be finite and non-negative");
    for (double weight : weights)
        if (!std::isfinite(weight) || weight < 0.0)
            throw IllegalArgumentException(std::string(name) +
                                           " weights must be finite and non-negative");
}

} // namespace redukti::optim
