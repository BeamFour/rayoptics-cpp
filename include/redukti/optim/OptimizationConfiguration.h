// C++ port of org.redukti.optim.OptimizationConfiguration
#ifndef REDUKTI_OPTIM_OPTIMIZATIONCONFIGURATION_H
#define REDUKTI_OPTIM_OPTIMIZATIONCONFIGURATION_H

#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/spec/Prescription.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::optim {

class Analysis;

/**
 * Serializable optimization settings, independent of a prescription or solver.
 * Each builder owns its settings; trial definitions keep a private copy and hand
 * out fresh copies so changing one stage cannot change another stage's settings.
 *
 * Package-private in the Java. Here it is a plain value type, so copying it is the
 * deep copy the Java's copy() makes field by field.
 */
class OptimizationConfiguration {
public:
    static constexpr int DEFAULT_HEXAPOLAR_RAYS = 64;
    static constexpr int DEFAULT_GAUSSIAN_QUADRATURE_RINGS = 14;
    static constexpr int DEFAULT_GAUSSIAN_QUADRATURE_SPOKES = 20;
    static constexpr int DEFAULT_CONTRAST_RINGS = 6;
    static constexpr int DEFAULT_CONTRAST_SPOKES = 12;

    /** Null until fields() is called; validate() rejects that. */
    std::optional<std::vector<double>> fields;
    /** Null until mtfFrequencies() is called; validate() rejects that. */
    std::optional<std::vector<int>> mtfFrequencies;
    std::vector<int> curvatureSurfaces;
    bool allCurvatureSurfaces = false;
    std::vector<int> thicknessSurfaces;
    bool allThicknessSurfaces = false;
    bool includeExistingAspherics = false;
    bool weighted = true;
    bool dLineOnly = false;
    bool addRayAberrationGoals = false;
    bool useHexapolarSpotPattern = false;
    int hexapolarSpotRays = DEFAULT_HEXAPOLAR_RAYS;
    int gaussianQuadratureRings = DEFAULT_GAUSSIAN_QUADRATURE_RINGS;
    int gaussianQuadratureSpokes = DEFAULT_GAUSSIAN_QUADRATURE_SPOKES;
    double gaussianQuadratureInnerRadius = 0.0;
    bool checkSpotApertures = true;
    std::optional<std::vector<double>> spotDeviationXWeights;
    std::optional<std::vector<double>> spotDeviationYWeights;
    bool addSpotDeviationGoals = false;
    // 3x6 is enough to measure a fixed design but not to optimize against: the
    // solver drives the 18 sampled points further than the wavefront between
    // them, so the merit reads better than the lens is. 6x12 is converged - 8x16
    // reproduces it - and 12 spokes samples the x and y axes alike, so sagittal
    // and tangential residuals stay comparable.
    int contrastRings = DEFAULT_CONTRAST_RINGS;
    int contrastSpokes = DEFAULT_CONTRAST_SPOKES;
    bool calibrateContrastFrequency = false;
    bool aimContrastAtExitPupil = false;
    bool centerContrastResiduals = false;
    /** Null until contrastBalanceGoals() is called. */
    std::optional<std::vector<bool>> contrastBalanceFields;
    double contrastBalanceWeight = OptimizationBuilder::NOMINAL_BALANCE_WEIGHT;
    int scenario = 0;
    spec::VigType vigType = spec::VigType::SetPupil;
    bool freezeVignetting = false;
    /** Null unless the matching applyXConstraints() was called. */
    std::optional<double> thicknessConstraintWeight;
    std::optional<double> edgeThicknessConstraintWeight;
    std::optional<double> curvatureConstraintWeight;
    std::vector<OptimizationBuilder::MtfGoals> mtfGoals;
    std::vector<OptimizationBuilder::ContrastGoals> contrastGoals;
    std::optional<OptimizationBuilder::SpotGoals> spotRmsGoals;
    std::optional<OptimizationBuilder::SpotGoals> spotMaxRadiusGoals;
    std::vector<int> curvatureExclusions;
    std::vector<int> thicknessExclusions;
    /** Aspheric terms varied explicitly, in the order given. */
    std::vector<OptimizationBuilder::AsphericTerm> asphericTerms;
    /** First-order goals, in the order given; efl and fno replace the automatic ones. */
    std::vector<OptimizationBuilder::ParaxialGoal> paraxialGoals;
    std::optional<std::string> description;
    std::optional<std::string> outdir;

    void gaussianSampling(int rings, int spokes, double innerPupilRadius);

    OptimizationConfiguration copy() const { return *this; }

    std::string toTrial(int number) const;

    /** The decisions shared by setup construction and canonical writing. */
    struct EffectiveAnalysis {
        bool spots;
        bool rayAberrations;
        bool mtf;
        bool hexapolar;
        bool retainFailedRays;
    };

    /**
     * Resolve the decisions shared by setup construction and canonical writing.
     * A custom maximum-radius goal also requests hexapolar sampling, but per-ray
     * deviation goals retain Gaussian sampling so their residual count stays fixed.
     * Custom goals have no written form, so the writer passes false.
     */
    EffectiveAnalysis effectiveAnalysis(bool customMaximumRadius) const;

    void configureAnalysis(Analysis &analysis, const EffectiveAnalysis &effective,
                           bool customGoals) const;

private:
    static std::string allExcept(const std::vector<int> &exclusions);

    /** Contrast weights: one row when every frequency shares them, else a row per frequency. */
    void contrastWeights(std::string &sb, bool sagittal) const;

    /** The balanced fields: all, all except the listed field values, or yes/no for each. */
    std::string balance() const;

    static void spotGoals(std::string &sb, const char *key,
                          const std::optional<OptimizationBuilder::SpotGoals> &goals);

    static bool allOnes(const std::vector<double> &values);
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_OPTIMIZATIONCONFIGURATION_H
