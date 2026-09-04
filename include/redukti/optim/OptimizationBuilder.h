// C++ port of org.redukti.optim.OptimizationBuilder
#ifndef REDUKTI_OPTIM_OPTIMIZATIONBUILDER_H
#define REDUKTI_OPTIM_OPTIMIZATIONBUILDER_H

#include "redukti/optim/Constraint.h"
#include "redukti/optim/Goals.h"
#include "redukti/optim/LMDer.h"
#include "redukti/optim/Var.h"
#include "redukti/spec/Prescription.h"

#include <functional>
#include <memory>
#include <optional>
#include <vector>

namespace redukti::optim {

/**
 * Builds the repetitive variables and goals used by an optical optimization.
 * Surface numbers are zero based, consistently with VarRadius and the other
 * optimization variable classes.
 */
class OptimizationBuilder {
public:
    /**
     * Default strength for applyThicknessConstraints() and
     * applyCurvatureConstraints().
     *
     * A single number suffices because Constraint residuals are fractions of
     * each parameter's starting value, so the per-parameter scaling is already
     * handled: a 0.1mm air gap and a 39mm back focus resist the same
     * PROPORTIONAL change equally. What this weight sets is only the global
     * trade between optical performance and preserving the layout.
     *
     * Raising it does tighten the design - on a 15-element f/2 with every space
     * free, the worst thickness excursion fell from 39% to 11% to 3% at weights
     * of 1, 10 and 100. But the useful range is narrow: past the nominal value
     * the optical cost outruns the benefit and the constraints start to dominate
     * the Jacobian, stalling the solver. Both a badly aberrated starting design
     * and a well corrected one behaved best here.
     *
     * For a one-off override on a particular surface, construct the constraint
     * directly through additionalGoals() rather than shifting the global weight.
     */
    static constexpr double NOMINAL_CONSTRAINT_WEIGHT = 1.0;

    /**
     * Default strength for contrastBalanceGoals().
     *
     * Much smaller than NOMINAL_CONSTRAINT_WEIGHT, and for a concrete reason: a
     * balance residual is a difference of sums of squares, so it is large where
     * a per-sample contrast residual is small. Measured on the Leica 75/2
     * starting design at 10/30/50 cyc/mm over 11 fields, the balance block at
     * weight 1.0 came to 43.6 against the contrast block's 52.6 - 83% of the
     * optical merit from 33 residuals against 14256. It would have run the solve.
     *
     * 0.1 puts it near 8% there, which is visible without dominating. It is a
     * starting point, not a normalization: unlike the constraints, nothing here
     * adapts to the design. Check the actual share on your own case before
     * trusting it.
     */
    static constexpr double NOMINAL_BALANCE_WEIGHT = 0.1;

    class MtfGoals;
    class ContrastGoals;
    class OptimizationSetup;

    /** Java's @FunctionalInterface GoalFactory. */
    using GoalFactory = std::function<std::shared_ptr<Goal>(Analysis *)>;

    static OptimizationBuilder builder(spec::Prescription *prescription) {
        return OptimizationBuilder(prescription);
    }

    // ------------------------------------------------------------------
    // Configuration - what gets evaluated, and how finely
    // ------------------------------------------------------------------

    OptimizationBuilder &fields(const std::vector<double> &fields_);
    OptimizationBuilder &mtfFrequencies(const std::vector<int> &frequencies);

    /** Use the prescription's wavelength weights; false assigns every wavelength weight 1.0. */
    OptimizationBuilder &weighted(bool weighted_);

    /** Restrict enabled ray-aberration goals to the Fraunhofer d-line. */
    OptimizationBuilder &dLineOnly(bool dLineOnly_);

    /**
     * Which configuration of a multi-configuration prescription to optimize;
     * zero based, and zero for a single-configuration lens.
     *
     * A zoom prescription carries a thickness per configuration on the varying
     * spaces, and its own focal length, f-number and angle of view. This selects
     * all of them together: the analysis builds the model for that
     * configuration, thickness variables read and write that configuration's
     * value, and the paraxial anchors target that configuration's focal length
     * and f-number.
     *
     * Configurations are optimized one at a time. Nothing here couples them, so
     * a space moved for one configuration is not reconciled against the others -
     * that has to be checked separately.
     */
    OptimizationBuilder &scenario(int scenario_);

    /**
     * How each rebuilt model establishes vignetting. See
     * Analysis::vignetting(VigType) for the trade-offs; VigType::Paraxial in
     * particular leaves the sagittal pupil unvignetted and breaks on-axis
     * rotational symmetry.
     */
    OptimizationBuilder &vignetting(spec::VigType vigType_);

    /**
     * Measure vignetting once at the start and hold it fixed for the run, so
     * every iteration is compared on the same pupil. See
     * Analysis::freezing_vignetting(bool) for the trade-off.
     */
    OptimizationBuilder &freezeVignetting() { return freezeVignetting(true); }
    OptimizationBuilder &freezeVignetting(bool freeze);

    /**
     * Configure the ordinary Gaussian-quadrature spot pattern shared by spot and
     * geometric-MTF analyses. Contrast uses its separate sheared-pupil pattern.
     */
    OptimizationBuilder &gaussianQuadratureSampling(int rings, int spokes) {
        return gaussianQuadratureSampling(rings, spokes, 0.0);
    }

    /** Configure Gaussian quadrature for a concentric annular pupil. */
    OptimizationBuilder &gaussianQuadratureSampling(int rings, int spokes,
                                                    double innerPupilRadius);

    /**
     * Whether Gaussian-quadrature spot rays are rejected when they cross a
     * physical surface aperture. Disable this with frozen vignetting to optimize
     * a fixed factor-defined pupil, matching the usual Zemax GQ merit-function
     * behaviour. Grid and hexapolar sampling always retain physical aperture
     * checking.
     */
    OptimizationBuilder &checkSpotApertures(bool check);

    /**
     * Use hexapolar sampling for spot analysis even when no maximum-radius goal
     * requires it. The default is Gaussian quadrature.
     */
    OptimizationBuilder &hexapolarSampling() { return hexapolarSampling(64); }

    /**
     * Use hexapolar sampling with the requested number of pupil rings for spot
     * analysis even when no maximum-radius goal requires it.
     */
    OptimizationBuilder &hexapolarSampling(int numRays);

    OptimizationBuilder &contrastSampling(int rings, int spokes);

    /**
     * Correct the contrast pupil shift so each sample realises the requested
     * spatial frequency in image space.
     *
     * The shift is applied in entrance-pupil coordinates but derives from an
     * exit-pupil relation, so pupil aberration makes the realised frequency fall
     * short, increasingly with field - measured 8.5% low at full field
     * tangential on an f/2 lens. Enabling this measures the shortfall per field,
     * wavelength and direction and scales the shift to compensate, which brought
     * that case to within 0.1%.
     *
     * Off by default because it changes every contrast residual.
     */
    OptimizationBuilder &calibrateContrastFrequency(bool value);

    /**
     * Iteratively aim each sheared contrast ray so its separation from the
     * reference ray is realised on the exit-pupil reference sphere. This
     * directly accounts for pupil aberration, including cross-axis displacement.
     *
     * This is mutually exclusive with calibrateContrastFrequency(bool). Off by
     * default because it changes the sampled rays and costs additional traces.
     */
    OptimizationBuilder &aimContrastAtExitPupil() { return aimContrastAtExitPupil(true); }
    OptimizationBuilder &aimContrastAtExitPupil(bool value);

    /**
     * Subtract the constant part of each contrast block, so the residuals carry
     * the variance the OTF modulus depends on rather than the un-centred second
     * moment.
     *
     * A constant wavefront difference across the pupil is tilt, which displaces
     * the image and costs no MTF - but it is reducible, so leaving it in offers
     * the solver merit reduction that buys nothing. It is identically zero in
     * the sagittal direction by symmetry and reaches 57% of an outer-field
     * tangential block on the Leica 75/2, which biases the astigmatic focus
     * split toward tangential.
     *
     * Off by default because it changes every contrast residual.
     */
    OptimizationBuilder &centerContrastResiduals(bool value);

    // ------------------------------------------------------------------
    // Variables - what the solver is allowed to change
    // ------------------------------------------------------------------

    OptimizationBuilder &varyCurvatures(const std::vector<int> &surfaces);

    /**
     * Vary every curved optical surface. Aperture/field stops and surfaces whose
     * radius is zero (and therefore intentionally flat) are excluded.
     */
    OptimizationBuilder &varyAllCurvatures();

    OptimizationBuilder &varyThicknesses(const std::vector<int> &surfaces);

    /**
     * Vary every thickness, air spaces and element thicknesses alike. Surfaces
     * with zero thickness are excluded, being coincident rather than a space to
     * open up.
     *
     * The counterpart to varyAllCurvatures(), and best paired with
     * applyThicknessConstraints() - with every space free and nothing holding
     * the layout, the solver will collapse gaps and drive elements through one
     * another.
     */
    OptimizationBuilder &varyAllThicknesses();

    /**
     * Vary the conic constants and polynomial coefficients already present in
     * the prescription. Only nonzero terms become variables, so a spherical
     * surface stays spherical and an asphere does not gain orders it did not
     * have.
     */
    OptimizationBuilder &varyExistingAspherics() { return varyExistingAspherics(true); }
    OptimizationBuilder &varyExistingAspherics(bool include);

    /** Adds caller-defined variables after the automatically generated variables. */
    OptimizationBuilder &additionalVariables(
        const std::vector<std::shared_ptr<Var>> &variables);

    // ------------------------------------------------------------------
    // Goals - what the solver optimizes towards
    // ------------------------------------------------------------------

    static MtfGoals mtf(int frequency, const std::vector<double> &sagittal,
                        const std::vector<double> &tangential);

    /** Applies one weight per field to both sagittal and tangential goals. */
    static MtfGoals mtf(int frequency, const std::vector<double> &sagittal,
                        const std::vector<double> &tangential,
                        const std::vector<double> &weights);

    static MtfGoals mtf(int frequency, const std::vector<double> &sagittal,
                        const std::vector<double> &tangential,
                        const std::vector<double> &sagittalWeights,
                        const std::vector<double> &tangentialWeights);

    /** Contrast optimization at one frequency, with one directional weight per field. */
    static ContrastGoals contrast(int frequency,
                                  const std::vector<double> &sagittalWeights,
                                  const std::vector<double> &tangentialWeights);

    /** Applies one weight per field to both contrast directions. */
    static ContrastGoals contrast(int frequency, const std::vector<double> &weights);

    OptimizationBuilder &mtfGoals(const std::vector<MtfGoals> &goals);
    OptimizationBuilder &contrastGoals(const std::vector<ContrastGoals> &goals);

    /**
     * Hold sagittal and tangential contrast in balance at the selected fields, at
     * the nominal balance weight. One flag per configured field, in field order;
     * a false adds no explicit balance constraint there, which is usually what
     * the outermost field wants. The ordinary contrast residuals still constrain
     * both meridians.
     *
     * Applies to every configured contrast frequency, so this adds one residual
     * per enabled field per frequency. See GoalContrastBalance for what it
     * measures and why the contrast merit does not already care.
     *
     * Two settings that surprise: the sagittal and tangential contrast weights
     * set the RATIO this goal targets rather than merely prioritising a
     * meridian, and on axis it cannot balance anything - with unequal weights it
     * quietly becomes an axial contrast goal instead. Both are covered on
     * GoalContrastBalance.
     */
    OptimizationBuilder &contrastBalanceGoals(const std::vector<bool> &fields_) {
        return contrastBalanceGoals(fields_, NOMINAL_BALANCE_WEIGHT);
    }

    /**
     * Hold sagittal and tangential contrast in balance at the selected fields, at
     * a chosen weight.
     *
     * The weight needs setting deliberately, because the scale is nothing like
     * the per-sample contrast residuals: a balance residual is a difference of
     * sums of squares, of order 0.1 to 3 waves squared on the test lenses, so a
     * handful of them can outweigh thousands of contrast residuals. See
     * NOMINAL_BALANCE_WEIGHT for the measured example. Compare the two blocks'
     * sum-of-squares contributions on your own case rather than assuming the
     * default is proportionate.
     *
     * @param fields_ one flag per configured field, in field order
     * @param weight  relative strength
     */
    OptimizationBuilder &contrastBalanceGoals(const std::vector<bool> &fields_,
                                              double weight);

    /**
     * One aggregate GoalSpotRMS per field, each aiming at a target RMS spot
     * radius in microns. To minimize spot size rather than hit a number, prefer
     * spotDeviationGoals(), which takes weights instead.
     */
    OptimizationBuilder &spotRmsGoals(const std::vector<double> &targets);
    OptimizationBuilder &spotRmsGoals(const std::vector<double> &targets,
                                      const std::vector<double> &weights);

    /**
     * Minimize RMS spot radius through the individual signed X/Y ray deviations
     * that make it up, one GoalSpotDeviation per orientation per sampled ray.
     * Differentiating those exposes far more to the solver than one
     * square-rooted aggregate does.
     *
     * These take WEIGHTS, not targets - every residual aims at zero. One weight
     * per field is applied to every wavelength, sample and orientation.
     */
    OptimizationBuilder &spotDeviationGoals(const std::vector<double> &fieldWeights);

    /** Assign separate per-field weights to the signed X and Y spot deviations. */
    OptimizationBuilder &spotDeviationGoals(const std::vector<double> &xWeights,
                                            const std::vector<double> &yWeights);

    OptimizationBuilder &spotMaxRadiusGoals(const std::vector<double> &targets);
    OptimizationBuilder &spotMaxRadiusGoals(const std::vector<double> &targets,
                                            const std::vector<double> &weights);

    /**
     * Add the legacy ten-sample sagittal and tangential ray-aberration fans to
     * the merit function. They are disabled by default: dense contrast or spot
     * goals already provide enough residuals, and a failed fan edge ray should
     * not invalidate an otherwise usable merit function.
     */
    OptimizationBuilder &rayAberrationGoals() { return rayAberrationGoals(true); }
    OptimizationBuilder &rayAberrationGoals(bool enabled);

    /**
     * Adds caller-defined goals after the automatically generated goals.
     * Factories receive the exact Analysis owned by the resulting setup.
     */
    OptimizationBuilder &additionalGoals(const std::vector<GoalFactory> &factories);

    // ------------------------------------------------------------------
    // Constraints - what holds the starting design together
    // ------------------------------------------------------------------

    /**
     * Hold the varied thicknesses near their starting values, at the nominal
     * weight.
     *
     * An optical merit function has no opinion about mechanical layout, so left
     * alone the solver will collapse air spaces and push elements through each
     * other and through the stop. This anchors each varied thickness to where it
     * began: it is still free to move, it just costs merit to do so.
     */
    OptimizationBuilder &applyThicknessConstraints() {
        return applyThicknessConstraints(NOMINAL_CONSTRAINT_WEIGHT);
    }
    OptimizationBuilder &applyThicknessConstraints(double weight);

    /**
     * Hold the EDGE separation of each varied gap near its starting value, at
     * the nominal weight.
     *
     * Complements applyThicknessConstraints() rather than replacing it. That one
     * holds axial centre thickness, which two surfaces can honour while still
     * crossing away from the axis once curvature moves - the failure mode that
     * produced overlapping first and second surfaces on the Leica 75/2 with
     * thickness constraints already in place. Use both when curvatures and
     * thicknesses are varied together.
     *
     * Gaps whose starting edge separation is not positive and finite are
     * skipped: a fractional constraint cannot be formed around zero, and a
     * design that already starts with coincident or crossed surfaces has nothing
     * useful to anchor to. See ConstraintEdgeThickness::is_constrainable.
     */
    OptimizationBuilder &applyEdgeThicknessConstraints() {
        return applyEdgeThicknessConstraints(NOMINAL_CONSTRAINT_WEIGHT);
    }
    OptimizationBuilder &applyEdgeThicknessConstraints(double weight);

    /**
     * Hold the varied surfaces near their starting CURVATURES, at the nominal
     * weight.
     *
     * Curvature, not radius: radius runs away towards infinity on a near-flat
     * surface for a negligible optical change, so a fractional radius constraint
     * would barely restrain it there while over-restraining a strongly curved
     * one. See ConstraintCurvature.
     */
    OptimizationBuilder &applyCurvatureConstraints() {
        return applyCurvatureConstraints(NOMINAL_CONSTRAINT_WEIGHT);
    }
    OptimizationBuilder &applyCurvatureConstraints(double weight);

    // ------------------------------------------------------------------
    // Build
    // ------------------------------------------------------------------

    OptimizationSetup build();

    /** Normalizes an existing coefficient to a scaled value in [1, 10). */
    static double scalingFor(double value);

    class MtfGoals {
    public:
        int frequency;
        std::vector<double> sagittal;
        std::vector<double> tangential;
        std::vector<double> sagittalWeights;
        std::vector<double> tangentialWeights;

        MtfGoals(int frequency_, const std::vector<double> &sagittal_,
                 const std::vector<double> &tangential_,
                 const std::optional<std::vector<double>> &sagittalWeights_,
                 const std::optional<std::vector<double>> &tangentialWeights_);

        void validate(int fieldCount) const;

        static void validateTargets(const std::vector<double> &values, int count,
                                    const char *name);
        static void validateWeights(const std::vector<double> &values, int count,
                                    const char *name);
    };

    class ContrastGoals {
    public:
        int frequency;
        std::vector<double> sagittalWeights;
        std::vector<double> tangentialWeights;

        ContrastGoals(int frequency_, const std::vector<double> &sagittalWeights_,
                      const std::vector<double> &tangentialWeights_)
            : frequency(frequency_), sagittalWeights(sagittalWeights_),
              tangentialWeights(tangentialWeights_) {}

        void validate(int fieldCount) const;
    };

    class OptimizationSetup {
    public:
        /** The Analysis is shared with merit functions and solvers built here. */
        Analysis *analysis() const { return _analysis.get(); }
        std::vector<std::shared_ptr<Var>> variables() const { return _variables; }
        std::vector<std::shared_ptr<Goal>> goals() const { return _goals; }

        LMDerMeritFunction meritFunction(bool useNative) const {
            return LMDerMeritFunction(_analysis, _variables, _goals, useNative);
        }

    private:
        friend class OptimizationBuilder;

        OptimizationSetup(std::shared_ptr<Analysis> analysis_,
                          std::vector<std::shared_ptr<Var>> variables_,
                          std::vector<std::shared_ptr<Goal>> goals_)
            : _analysis(std::move(analysis_)), _variables(std::move(variables_)),
              _goals(std::move(goals_)) {}

        std::shared_ptr<Analysis> _analysis;
        std::vector<std::shared_ptr<Var>> _variables;
        std::vector<std::shared_ptr<Goal>> _goals;
    };

private:
    explicit OptimizationBuilder(spec::Prescription *prescription);

    /** Thickness of a surface for the scenario this builder targets. */
    double thicknessOf(int surface) const;
    /** Effective focal length this scenario is anchored to. */
    double focalLengthOf() const;
    /** F-number this scenario is anchored to. */
    double fNumberOf() const;

    /**
     * Reject a scenario the prescription does not define, rather than letting it
     * surface as an array index failure from somewhere inside the solve.
     */
    void validateScenario() const;

    /** Number of configurations the prescription defines; 1 when it is not a zoom. */
    int scenarioCount() const;

    void configureContrastAnalysis(Analysis &analysis) const;
    void configureRequiredAnalyses(Analysis &analysis,
                                   const std::vector<std::shared_ptr<Goal>> &goals) const;
    void configureSpotPattern(Analysis &analysis,
                              const std::vector<std::shared_ptr<Goal>> &goals) const;

    std::vector<std::shared_ptr<Var>> buildVariables() const;
    std::vector<std::shared_ptr<Goal>> buildGoals(
        Analysis *analysis, const std::vector<std::shared_ptr<Var>> &variables) const;

    void validate() const;
    void validateSurfaces(const std::vector<int> &surfaces, const char *kind) const;

    static bool sameWavelength(double a, double b);

    static constexpr int RAY_FAN_SAMPLES = 10;

    spec::Prescription *prescription;
    /** Null until fields() is called; validate() rejects that. */
    std::optional<std::vector<double>> _fields;
    /** Null until mtfFrequencies() is called; validate() rejects that. */
    std::optional<std::vector<int>> _mtfFrequencies;
    std::vector<int> curvatureSurfaces;
    bool allCurvatureSurfaces = false;
    std::vector<int> thicknessSurfaces;
    bool allThicknessSurfaces = false;
    bool includeExistingAspherics = false;
    bool _weighted = true;
    bool _dLineOnly = false;
    bool addRayAberrationGoals = false;
    bool useHexapolarSpotPattern = false;
    int hexapolarSpotRays = 64;
    int gaussianQuadratureRings = 14;
    int gaussianQuadratureSpokes = 20;
    double gaussianQuadratureInnerRadius = 0.0;
    bool _checkSpotApertures = true;
    std::optional<std::vector<double>> spotDeviationXWeights;
    std::optional<std::vector<double>> spotDeviationYWeights;
    bool addSpotDeviationGoals = false;
    // 3x6 is enough to measure a fixed design but not to optimize against: the
    // solver drives the 18 sampled points further than the wavefront between
    // them, so the merit reads better than the lens is. 6x12 is converged - 8x16
    // reproduces it - and 12 spokes samples the x and y axes alike, so sagittal
    // and tangential residuals stay comparable.
    int contrastRings = 6;
    int contrastSpokes = 12;
    bool calibrateContrastFrequency_ = false;
    bool aimContrastAtExitPupil_ = false;
    bool centerContrastResiduals_ = false;
    /** Null until contrastBalanceGoals() is called. */
    std::optional<std::vector<bool>> contrastBalanceFields;
    double contrastBalanceWeight = NOMINAL_BALANCE_WEIGHT;
    int _scenario = 0;
    spec::VigType vigType = spec::VigType::SetPupil;
    bool freezeVignetting_ = false;
    /** Null unless the matching applyXConstraints() was called. */
    std::optional<double> thicknessConstraintWeight;
    std::optional<double> edgeThicknessConstraintWeight;
    std::optional<double> curvatureConstraintWeight;
    std::vector<MtfGoals> _mtfGoals;
    std::vector<ContrastGoals> _contrastGoals;
    std::vector<std::shared_ptr<Var>> additionalVariables_;
    std::vector<GoalFactory> additionalGoalFactories;
    /** Java's nullable SpotGoals; targets and weights, one per field. */
    struct SpotGoals {
        std::vector<double> targets;
        std::vector<double> weights;

        void validate(int fieldCount, const char *name) const;
    };
    std::optional<SpotGoals> spotRmsGoals_;
    std::optional<SpotGoals> spotMaxRadiusGoals_;
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_OPTIMIZATIONBUILDER_H
