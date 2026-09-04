// C++ port of org.redukti.optim.DirectionalGoalTest, LMDerMeritFunctionTest and
// OptimizationBuilderTest.
//
// One porting note runs through the whole file. The Java writes
// `OptimizationBuilder.builder(prescription())` against a freshly allocated
// prescription the collector keeps alive; here the builder, the Analysis and
// every Var borrow that object, so each test has to hold it in a named local
// that outlives the setup.
#include "TestHarness.h"

#include "redukti/mathlib/LMLSolver.h"
#include "redukti/optim/ConfigurationReport.h"
#include "redukti/optim/OptimizationBuilder.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/util/Orientation.h"
#include "redukti/rayoptics/analysis/SpotAnalysis.h"
#include "redukti/rayoptics/seq/Glass.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/spec/SurfaceType.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

namespace {

using redukti::optim::Analysis;
using redukti::optim::ConstraintEdgeThickness;
using redukti::optim::ConstraintThickness;
using redukti::optim::Goal;
using redukti::optim::GoalContrast;
using redukti::optim::GoalContrastBalance;
using redukti::optim::GoalGeoMTF;
using redukti::optim::GoalMTFProxy;
using redukti::optim::GoalParax;
using redukti::optim::GoalRayAberration;
using redukti::optim::GoalSpotDeviation;
using redukti::optim::GoalSpotMaxRadius;
using redukti::optim::GoalSpotRMS;
using redukti::optim::LMDerMeritFunction;
using redukti::optim::OptimizationBuilder;
using redukti::optim::ParaxHelper;
using redukti::optim::Var;
using redukti::optim::VarAsphCoeff;
using redukti::optim::VarAsphK;
using redukti::optim::VarRadius;
using redukti::optim::VarThickness;
using redukti::rayoptics::analysis::SpotOptions;
using redukti::rayoptics::seq::Glass;
using redukti::spec::Prescription;
using redukti::spec::SurfaceType;
using redukti::spec::VigType;
namespace Orientation = redukti::rayoptics::util::Orientation;

Prescription prescription() {
    return Prescription(50.0, 1.4, 40.0, 43.28, std::vector<double>{Glass::d, Glass::F, Glass::C},
                        std::vector<double>{1.0, 0.5, 0.25})
        .surf(50.0, 5.0, 30.0, 1.5, 50.0)
        .asph(SurfaceType::ASPH_EVEN, -0.5, {0.0, 2.5e-6, 0.0, -4.0e-10})
        .surf(-50.0, 20.0, 30.0)
        .stop(1.0, 20.0)
        .surf(0.0, 5.0, 30.0)
        .build();
}

/** Every goal in `goals` that is a T, in order. */
template <typename T> std::vector<T *> goalsOfType(const std::vector<std::shared_ptr<Goal>> &goals) {
    std::vector<T *> out;
    for (const auto &goal : goals)
        if (auto *typed = dynamic_cast<T *>(goal.get()))
            out.push_back(typed);
    return out;
}

void assertMtf(Goal *goal, int frequency, int field, int orientation, double target,
               double weight) {
    auto *mtf = dynamic_cast<GoalGeoMTF *>(goal);
    CHECK(mtf != nullptr);
    if (mtf == nullptr)
        return;
    CHECK_EQ(mtf->_freq, frequency);
    CHECK_EQ(mtf->_field, field);
    CHECK_EQ(mtf->_orientation, orientation);
    CHECK_EQ(mtf->_target, target);
    CHECK_EQ(mtf->_weight, weight);
}

std::vector<double> vignettingOf(Analysis &analysis, int fieldIndex) {
    const auto &field =
        analysis._opt_model->optical_spec->fov->fields[static_cast<std::size_t>(fieldIndex)];
    return {field->vux, field->vlx, field->vuy, field->vly};
}

/** Sum of squares of the merit's own contrast residuals for one field and orientation. */
double blockSumOfSquares(const OptimizationBuilder::OptimizationSetup &setup, int orientation,
                         int field) {
    double sum = 0.0;
    for (auto *goal : goalsOfType<GoalContrast>(setup.goals()))
        if (goal->_field == field && goal->_orientation == orientation)
            sum += goal->value() * goal->value();
    return sum;
}

double wavelengthSumOfSquares(const OptimizationBuilder::OptimizationSetup &setup, int orientation,
                              int field, int wavelength) {
    double sum = 0.0;
    for (auto *goal : goalsOfType<GoalContrast>(setup.goals()))
        if (goal->_field == field && goal->_orientation == orientation &&
            goal->_wavelength_index == wavelength)
            sum += goal->value() * goal->value();
    return sum;
}

} // namespace

// ===========================================================================
// DirectionalGoalTest
// ===========================================================================

/**
 * The goals that come in sagittal/tangential pairs all validate their
 * orientation against Orientation. Three of these used to accept any int and
 * read the sagittal fan for it, so a typo produced a plausible but wrong
 * residual rather than an error.
 */
TEST(optim_directional_goals_reject_a_non_meridian_orientation) {
    using redukti::IllegalArgumentException;
    CHECK_THROWS(GoalGeoMTF(nullptr, 1, 2, 10, 0.5, 1.0), IllegalArgumentException);
    CHECK_THROWS(GoalRayAberration(nullptr, 1, 2, 0, 0.5876, 0.0, 1.0), IllegalArgumentException);
    CHECK_THROWS(GoalMTFProxy(nullptr, 1, 2, 0, 0.5876, 10, 0.0, 1.0), IllegalArgumentException);
    CHECK_THROWS(GoalContrast(nullptr, 0, 10, 1, 0, 0, 2, 1.0), IllegalArgumentException);
    CHECK_THROWS(GoalSpotDeviation(nullptr, 1, 0, 0, 2, 1.0), IllegalArgumentException);
}

// ===========================================================================
// LMDerMeritFunctionTest
// ===========================================================================

namespace {

using Validity = bool (*)(double);

/** Analysis stub whose compute() succeeds only for values the fixture allows. */
class TestAnalysis : public Analysis {
public:
    Validity validity;
    double value = 1.0;
    bool valid = true;

    explicit TestAnalysis(Validity validity_)
        : Analysis(nullptr, {0.0}, {1}), validity(validity_) {}

    void compute() override {
        valid = validity(value);
        if (!valid)
            throw redukti::IllegalStateException("synthetic killed ray");
    }
};

class TestVar : public Var {
public:
    TestAnalysis *analysis;

    explicit TestVar(TestAnalysis *analysis_) : Var(nullptr), analysis(analysis_) {
        _d_delta = 1.0e-4;
    }

    double read_from_prescription() override {
        set_unscaled_value(analysis->value);
        return get_scaled_value();
    }

    void write_to_prescription() override { analysis->value = get_unscaled_value(); }
};

class TestGoal : public Goal {
public:
    TestAnalysis *analysis;

    explicit TestGoal(TestAnalysis *analysis_)
        : Goal(analysis_, 0.0, 1.0), analysis(analysis_) {}

    double value() override {
        return analysis->valid ? analysis->value * analysis->value
                               : redukti::mathlib::LMLSolver::BIGVAL;
    }
};

struct Fixture {
    TestAnalysis analysis;
    std::shared_ptr<TestVar> variable;
    std::shared_ptr<TestGoal> goal;
    LMDerMeritFunction merit;

    explicit Fixture(Validity validity)
        : analysis(validity), variable(std::make_shared<TestVar>(&analysis)),
          goal(std::make_shared<TestGoal>(&analysis)),
          merit(&analysis, {variable}, {goal}, false) {}
};

} // namespace

TEST(optim_uses_central_difference_when_both_perturbations_are_valid) {
    Fixture fixture([](double) { return true; });
    std::vector<double> jacobian(1, 0.0);
    std::vector<double> x{1.0};

    CHECK(fixture.merit.buildJacobian(x, jacobian, 1));

    CHECK_CLOSE(jacobian[0], 2.0, 1.0e-8);
    CHECK_EQ(fixture.analysis.value, 1.0);
}

TEST(optim_reduces_step_and_uses_one_sided_difference_at_validity_boundary) {
    Fixture fixture([](double value) { return value <= 1.0; });
    std::vector<double> jacobian(1, 0.0);
    std::vector<double> x{1.0};

    CHECK(fixture.merit.buildJacobian(x, jacobian, 1));

    CHECK_CLOSE(jacobian[0], 2.0, 1.0e-5);
    CHECK_EQ(fixture.analysis.value, 1.0);
}

TEST(optim_rejects_unknown_derivative_and_restores_base_point) {
    Fixture fixture([](double value) { return value == 1.0; });
    std::vector<double> jacobian(1, 0.0);
    std::vector<double> x{1.0};

    CHECK(!fixture.merit.buildJacobian(x, jacobian, 1));

    CHECK_EQ(fixture.analysis.value, 1.0);
    CHECK(fixture.analysis.valid);
    CHECK_EQ(fixture.goal->value(), 1.0);
}

// ===========================================================================
// OptimizationBuilderTest
// ===========================================================================

TEST(optim_builds_variables_and_curve_oriented_mtf_goals) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10, 20})
                     .varyCurvatures({0, 1})
                     .varyThicknesses({1})
                     .varyExistingAspherics()
                     .rayAberrationGoals()
                     .mtfGoals({OptimizationBuilder::mtf(10, {90, 80}, {85, 65}, {2.0, 0.5})})
                     .build();

    auto variables = setup.variables();
    CHECK_EQ(static_cast<int>(variables.size()), 6);
    CHECK(dynamic_cast<VarRadius *>(variables[0].get()) != nullptr);
    CHECK(dynamic_cast<VarThickness *>(variables[2].get()) != nullptr);
    CHECK(dynamic_cast<VarAsphK *>(variables[3].get()) != nullptr);
    auto *coeff = dynamic_cast<VarAsphCoeff *>(variables[4].get());
    CHECK(coeff != nullptr);
    if (coeff != nullptr)
        CHECK_EQ(coeff->_scaling_factor, 1.0e6);

    auto goals = setup.goals();
    CHECK_EQ(static_cast<int>(goals.size()), 126);
    assertMtf(goals[0].get(), 10, 1, Orientation::SAGITTAL, 0.90, 2.0);
    assertMtf(goals[1].get(), 10, 1, Orientation::TANGENTIAL, 0.85, 2.0);
    assertMtf(goals[2].get(), 10, 2, Orientation::SAGITTAL, 0.80, 0.5);
    assertMtf(goals[3].get(), 10, 2, Orientation::TANGENTIAL, 0.65, 0.5);

    auto *focalLength = dynamic_cast<GoalParax *>(goals[4].get());
    auto *fNumber = dynamic_cast<GoalParax *>(goals[5].get());
    CHECK(focalLength != nullptr && fNumber != nullptr);
    if (focalLength != nullptr && fNumber != nullptr) {
        CHECK_EQ(focalLength->_parax_id, ParaxHelper::Effective_focal_length);
        CHECK_EQ(focalLength->_target, 50.0);
        CHECK_EQ(fNumber->_parax_id, ParaxHelper::Fno);
        CHECK_EQ(fNumber->_target, 1.4);
    }

    std::set<double> rayWeights;
    for (auto *goal : goalsOfType<GoalRayAberration>(goals))
        rayWeights.insert(goal->_weight);
    std::vector<double> sortedWeights(rayWeights.begin(), rayWeights.end());
    CHECK_EQ(static_cast<int>(sortedWeights.size()), 3);
    if (sortedWeights.size() == 3) {
        CHECK_EQ(sortedWeights[0], 0.25);
        CHECK_EQ(sortedWeights[1], 0.5);
        CHECK_EQ(sortedWeights[2], 1.0);
    }
}

TEST(optim_supports_d_line_only_and_unweighted_goals) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .weighted(false)
                     .dLineOnly(true)
                     .rayAberrationGoals()
                     .build();

    // two paraxial + 2 fields * 2 fans * 10 samples
    CHECK_EQ(static_cast<int>(setup.goals().size()), 42);
    for (auto *goal : goalsOfType<GoalRayAberration>(setup.goals())) {
        CHECK_EQ(goal->_wvl, Glass::d);
        CHECK_EQ(goal->_weight, 1.0);
    }
}

TEST(optim_ray_aberration_goals_are_optional_and_disabled_by_default) {
    auto p = prescription();
    auto setup =
        OptimizationBuilder::builder(&p).fields({0.0, 1.0}).mtfFrequencies({10}).build();

    // only the paraxial anchors are automatic
    CHECK_EQ(static_cast<int>(setup.goals().size()), 2);
    CHECK(goalsOfType<GoalRayAberration>(setup.goals()).empty());
    CHECK(!setup.analysis()->_compute_ray_aberrations);
}

TEST(optim_rejects_underdetermined_merit_instead_of_adding_hidden_ray_goals) {
    auto p = prescription();
    bool threw = false;
    std::string message;
    try {
        OptimizationBuilder::builder(&p)
            .fields({0.0})
            .mtfFrequencies({10})
            .varyCurvatures({0, 1})
            .varyThicknesses({0, 1})
            .build();
    } catch (const redukti::IllegalArgumentException &e) {
        threw = true;
        message = e.what();
    }
    CHECK(threw);
    CHECK(message.find("rayAberrationGoals()") != std::string::npos);
}

TEST(optim_builds_spot_goals_with_optional_field_weights) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .spotRmsGoals({8.0, 20.0}, {3.0, 1.0})
                     .spotMaxRadiusGoals({15.0, 50.0})
                     .build();

    auto goals = setup.goals();
    auto *rms0 = dynamic_cast<GoalSpotRMS *>(goals[0].get());
    auto *rms1 = dynamic_cast<GoalSpotRMS *>(goals[1].get());
    auto *max0 = dynamic_cast<GoalSpotMaxRadius *>(goals[2].get());
    auto *max1 = dynamic_cast<GoalSpotMaxRadius *>(goals[3].get());
    CHECK(rms0 != nullptr && rms1 != nullptr && max0 != nullptr && max1 != nullptr);
    if (rms0 == nullptr || rms1 == nullptr || max0 == nullptr || max1 == nullptr)
        return;
    CHECK_EQ(rms0->_target, 8.0);
    CHECK_EQ(rms0->_weight, 3.0);
    CHECK_EQ(rms1->_target, 20.0);
    CHECK_EQ(rms1->_weight, 1.0);
    CHECK_EQ(max0->_target, 15.0);
    CHECK_EQ(max0->_weight, 1.0);
    CHECK_EQ(max1->_target, 50.0);
    CHECK_EQ(max1->_weight, 1.0);
}

TEST(optim_selects_all_curved_non_stop_surfaces) {
    auto p = prescription();
    auto variables = OptimizationBuilder::builder(&p)
                         .fields({0.0})
                         .mtfFrequencies({10})
                         .varyAllCurvatures()
                         .build()
                         .variables();

    CHECK_EQ(static_cast<int>(variables.size()), 2);
    CHECK_EQ(dynamic_cast<VarRadius *>(variables[0].get())->_surface_id, 0);
    CHECK_EQ(dynamic_cast<VarRadius *>(variables[1].get())->_surface_id, 1);
}

TEST(optim_appends_user_supplied_variables_and_analysis_bound_goals) {
    auto p = prescription();
    auto customVariable = std::make_shared<VarThickness>(&p, 0);
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .dLineOnly(true)
                     .additionalVariables({customVariable})
                     .additionalGoals({[](Analysis *analysis) -> std::shared_ptr<Goal> {
                         return std::make_shared<GoalSpotRMS>(analysis, 1, 7.5, 4.0);
                     }})
                     .build();

    CHECK_EQ(static_cast<int>(setup.variables().size()), 1);
    CHECK(setup.variables()[0].get() == customVariable.get());
    auto goals = setup.goals();
    auto *spotGoal = dynamic_cast<GoalSpotRMS *>(goals.back().get());
    CHECK(spotGoal != nullptr);
    if (spotGoal == nullptr)
        return;
    CHECK(spotGoal->_analysis == setup.analysis());
    CHECK_EQ(spotGoal->_target, 7.5);
    CHECK_EQ(spotGoal->_weight, 4.0);
}

/**
 * SetPupil live is what every committed regression value was generated under.
 * Both halves of this matter: a changed default silently moves every golden
 * number, and a default that froze would quietly stop tracking the design.
 */
TEST(optim_vignetting_defaults_to_live_set_pupil_and_is_configurable) {
    auto p = prescription();
    CHECK(OptimizationBuilder::builder(&p)
              .fields({0.0})
              .mtfFrequencies({10})
              .build()
              .analysis()
              ->_vig_type == VigType::SetPupil);
    CHECK(!OptimizationBuilder::builder(&p)
               .fields({0.0})
               .mtfFrequencies({10})
               .build()
               .analysis()
               ->_freeze_vignetting);

    auto configured = OptimizationBuilder::builder(&p)
                          .fields({0.0})
                          .mtfFrequencies({10})
                          .vignetting(VigType::Paraxial)
                          .freezeVignetting()
                          .build();
    CHECK(configured.analysis()->_vig_type == VigType::Paraxial);
    CHECK(configured.analysis()->_freeze_vignetting);
}

TEST(optim_spot_aperture_checking_defaults_on_and_is_configurable) {
    auto p = prescription();
    auto defaults =
        OptimizationBuilder::builder(&p).fields({0.0}).mtfFrequencies({10}).build();
    CHECK(defaults.analysis()->_check_spot_apertures);

    auto fixedPupil = OptimizationBuilder::builder(&p)
                          .fields({0.0})
                          .mtfFrequencies({10})
                          .freezeVignetting()
                          .checkSpotApertures(false)
                          .build();
    CHECK(!fixedPupil.analysis()->_check_spot_apertures);
}

/**
 * The point of freezing: the factors must not follow the design. Nothing else in
 * the merit would notice if they did - the residuals would simply be evaluated
 * on a pupil that quietly moved, which is the drift this option exists to remove.
 */
TEST(optim_frozen_vignetting_survives_a_prescription_change) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .freezeVignetting()
                     .build();
    Analysis *analysis = setup.analysis();

    analysis->compute();
    auto captured = analysis->frozen_vignetting();
    CHECK(captured.has_value());
    if (!captured.has_value())
        return;
    CHECK_EQ(static_cast<int>(captured->size()), 2);
    auto modelFactors = vignettingOf(*analysis, 1);
    for (std::size_t i = 0; i < 4; i++)
        CHECK_EQ((*captured)[1][i], modelFactors[i]);

    // Move a surface far enough that live factors would certainly follow.
    p._surface_list[0]._radius *= 0.5;
    analysis->compute();
    auto after = vignettingOf(*analysis, 1);
    for (std::size_t i = 0; i < 4; i++)
        CHECK_EQ((*captured)[1][i], after[i]); // frozen vignetting must not follow
    auto stillFrozen = analysis->frozen_vignetting();
    CHECK(stillFrozen.has_value());
    if (stillFrozen.has_value())
        for (std::size_t i = 0; i < 4; i++)
            CHECK_EQ((*captured)[1][i], (*stillFrozen)[1][i]);

    // Unfreezing drops the capture, so the next compute measures the design again.
    analysis->freezing_vignetting(false);
    CHECK(!analysis->frozen_vignetting().has_value());
}

/**
 * The case ConstraintThickness cannot see: axial thickness untouched, but the
 * surfaces bent until they cross away from the axis. If this ever stops failing
 * for the thickness constraint, the edge constraint has lost its reason to exist.
 */
TEST(optim_edge_constraint_sees_curvature_driven_crossing) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .varyThicknesses({0})
                     .applyThicknessConstraints()
                     .applyEdgeThicknessConstraints()
                     .rayAberrationGoals()
                     .build();

    auto axials = goalsOfType<ConstraintThickness>(setup.goals());
    auto edges = goalsOfType<ConstraintEdgeThickness>(setup.goals());
    CHECK(!axials.empty() && !edges.empty());
    if (axials.empty() || edges.empty())
        return;
    ConstraintThickness *axial = axials[0];
    ConstraintEdgeThickness *edge = edges[0];

    // Surface 0 is r=50 with a 30mm diameter, so the gap is measured at h=15.
    CHECK_CLOSE(edge->_height, 15.0, 1.0e-12);
    CHECK_CLOSE(axial->fractional_deviation(), 0.0, 1.0e-12);
    CHECK_CLOSE(edge->fractional_deviation(), 0.0, 1.0e-12);

    // Bend surface 0 hard toward surface 1 without touching any thickness.
    p._surface_list[0]._radius = 16.0;

    // axial thickness is unchanged, which is exactly why it cannot see this
    CHECK_CLOSE(axial->fractional_deviation(), 0.0, 1.0e-12);
    CHECK(edge->value() < 0.0);
    CHECK(edge->fractional_deviation() < -1.0);
}

TEST(optim_edge_constraints_skip_gaps_that_cannot_be_anchored) {
    auto p = prescription();
    // The last surface has no following surface to form a gap with.
    int last = static_cast<int>(p._surface_list.size()) - 1;
    Analysis probe(&p, {0.0}, {10});
    CHECK(!ConstraintEdgeThickness::is_constrainable(&probe, last));

    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .varyThicknesses({last})
                     .applyEdgeThicknessConstraints()
                     .rayAberrationGoals()
                     .build();

    CHECK(goalsOfType<ConstraintEdgeThickness>(setup.goals()).empty());
}

TEST(optim_contrast_balance_goals_are_built_per_enabled_field_and_frequency) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 0.5, 1.0})
                     .mtfFrequencies({10, 20})
                     .contrastSampling(2, 4)
                     .contrastGoals({OptimizationBuilder::contrast(10, {1, 1, 1}),
                                     OptimizationBuilder::contrast(20, {1, 1, 1})})
                     // Outer field left unconstrained, as leniency there is usual.
                     .contrastBalanceGoals({true, true, false}, 4.0)
                     .build();

    auto balance = goalsOfType<GoalContrastBalance>(setup.goals());

    // Two enabled fields times two frequencies, and nothing for the disabled field.
    CHECK_EQ(static_cast<int>(balance.size()), 4);
    if (balance.size() != 4)
        return;
    const int expectedFields[] = {0, 1, 0, 1};
    const int expectedFreqs[] = {10, 10, 20, 20};
    for (std::size_t i = 0; i < 4; i++) {
        CHECK_EQ(balance[i]->_field, expectedFields[i]);
        CHECK_EQ(balance[i]->_frequency, expectedFreqs[i]);
        CHECK_EQ(balance[i]->_target, 0.0);
        CHECK_EQ(balance[i]->_weight, 4.0);
    }
}

TEST(optim_contrast_balance_rejects_a_flag_per_field_mismatch) {
    using redukti::IllegalArgumentException;
    auto p = prescription();
    bool threw = false;
    std::string message;
    try {
        OptimizationBuilder::builder(&p)
            .fields({0.0, 0.5, 1.0})
            .mtfFrequencies({10})
            .contrastSampling(2, 4)
            .contrastGoals({OptimizationBuilder::contrast(10, {1, 1, 1})})
            .contrastBalanceGoals({true, true})
            .build();
    } catch (const IllegalArgumentException &e) {
        threw = true;
        message = e.what();
    }
    CHECK(threw);
    CHECK(message.find("one flag per field") != std::string::npos);

    // balance without contrast goals has nothing to balance
    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .contrastBalanceGoals({true})
                     .build(),
                 IllegalArgumentException);
}

/**
 * The point of the goal: zero when the two meridians contribute equally, and
 * signed by which one is worse. It must not care about the overall aberration
 * level, only the split - that is what the contrast merit already handles.
 */
TEST(optim_contrast_balance_reads_zero_when_meridians_are_equal) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     // Unweighted so the goal's wavelength pooling and the sum
                     // below coincide; weighted pooling is covered separately.
                     .weighted(false)
                     .contrastSampling(2, 4)
                     .contrastGoals({OptimizationBuilder::contrast(10, {1, 1})})
                     .contrastBalanceGoals({true, true})
                     .build();
    setup.analysis()->compute();

    auto balance = goalsOfType<GoalContrastBalance>(setup.goals());
    GoalContrastBalance *onAxis = nullptr;
    GoalContrastBalance *offAxis = nullptr;
    for (auto *goal : balance) {
        if (goal->_field == 0 && onAxis == nullptr)
            onAxis = goal;
        if (goal->_field == 1 && offAxis == nullptr)
            offAxis = goal;
    }
    CHECK(onAxis != nullptr && offAxis != nullptr);
    if (onAxis == nullptr || offAxis == nullptr)
        return;

    // On axis the two meridians are identical by rotational symmetry, so a
    // balance goal there is satisfied however aberrated the lens is.
    CHECK_CLOSE(onAxis->value(), 0.0, 1.0e-12);

    double sagittal = blockSumOfSquares(setup, Orientation::SAGITTAL, 1);
    double tangential = blockSumOfSquares(setup, Orientation::TANGENTIAL, 1);
    // balance must equal the difference the two meridians contribute to the merit
    CHECK_CLOSE(offAxis->value(), sagittal - tangential, 1.0e-9);
}

/**
 * The goal pools per-wavelength blocks using the same wavelength weights the
 * merit applies to the contrast residuals, so an unweighted setup is the case
 * where this plain sum matches.
 */
TEST(optim_contrast_balance_pools_wavelengths_with_the_configured_weights) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .weighted(true) // prescription weights are 1.0, 0.5, 0.25
                     .contrastSampling(2, 4)
                     .contrastGoals({OptimizationBuilder::contrast(10, {1, 1})})
                     .contrastBalanceGoals({false, true})
                     .build();
    setup.analysis()->compute();

    auto balance = goalsOfType<GoalContrastBalance>(setup.goals());
    CHECK(!balance.empty());
    if (balance.empty())
        return;

    const double weights[] = {1.0, 0.5, 0.25};
    double expected = 0.0;
    for (int w = 0; w < 3; w++) {
        double sagittal = wavelengthSumOfSquares(setup, Orientation::SAGITTAL, 1, w);
        double tangential = wavelengthSumOfSquares(setup, Orientation::TANGENTIAL, 1, w);
        expected += weights[w] * (sagittal - tangential);
    }
    CHECK_CLOSE(balance[0]->value(), expected, 1.0e-9);
}

TEST(optim_contrast_balance_includes_sagittal_and_tangential_field_weights) {
    auto p = prescription();
    auto setup =
        OptimizationBuilder::builder(&p)
            .fields({0.0, 1.0})
            .mtfFrequencies({10})
            .weighted(false)
            .contrastSampling(2, 4)
            .contrastGoals({OptimizationBuilder::contrast(10, {1.0, 2.5}, {1.0, 0.4})})
            .contrastBalanceGoals({false, true})
            .build();
    setup.analysis()->compute();

    auto balance = goalsOfType<GoalContrastBalance>(setup.goals());
    CHECK(!balance.empty());
    if (balance.empty())
        return;

    double sagittal = blockSumOfSquares(setup, Orientation::SAGITTAL, 1);
    double tangential = blockSumOfSquares(setup, Orientation::TANGENTIAL, 1);
    // balance must include the contrast merit's field/orientation weights
    CHECK_CLOSE(balance[0]->value(), 2.5 * sagittal - 0.4 * tangential, 1.0e-9);
}

TEST(optim_uses_gaussian_quadrature_by_default) {
    auto p = prescription();
    auto setup =
        OptimizationBuilder::builder(&p).fields({0.0}).mtfFrequencies({10}).build();
    Analysis *analysis = setup.analysis();

    CHECK_EQ(analysis->_spot_pattern, SpotOptions::PATTERN_GAUSS_QUADRATURE);
    CHECK_EQ(analysis->_num_rings, 14);
    CHECK_EQ(analysis->_num_spokes, 20);
    CHECK_EQ(analysis->_inner_pupil_radius, 0.0);

    auto configuredSetup = OptimizationBuilder::builder(&p)
                               .fields({0.0})
                               .mtfFrequencies({10})
                               .gaussianQuadratureSampling(3, 8, 0.5)
                               .build();
    CHECK_EQ(configuredSetup.analysis()->_num_rings, 3);
    CHECK_EQ(configuredSetup.analysis()->_num_spokes, 8);
    CHECK_EQ(configuredSetup.analysis()->_inner_pupil_radius, 0.5);
}

TEST(optim_permits_explicit_hexapolar_spot_sampling) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .hexapolarSampling(32)
                     .build();

    CHECK_EQ(setup.analysis()->_spot_pattern, SpotOptions::PATTERN_HEXAPOLAR);
    CHECK_EQ(setup.analysis()->_num_rays, 32);
}

TEST(optim_maximum_spot_radius_goal_forces_hexapolar_sampling) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .spotMaxRadiusGoals({20.0})
                     .build();

    CHECK_EQ(setup.analysis()->_spot_pattern, SpotOptions::PATTERN_HEXAPOLAR);
    CHECK_EQ(setup.analysis()->_num_rays, 64);
}

TEST(optim_additional_maximum_spot_radius_goal_forces_hexapolar_sampling) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .additionalGoals({[](Analysis *a) -> std::shared_ptr<Goal> {
                         return std::make_shared<GoalSpotMaxRadius>(a, 1, 20.0, 1.0);
                     }})
                     .build();

    CHECK_EQ(setup.analysis()->_spot_pattern, SpotOptions::PATTERN_HEXAPOLAR);
}

TEST(optim_rejects_invalid_hexapolar_sample_count) {
    auto p = prescription();
    CHECK_THROWS(OptimizationBuilder::builder(&p).hexapolarSampling(0),
                 redukti::IllegalArgumentException);
}

TEST(optim_builds_stable_per_sample_contrast_goals_and_skips_unused_mtf) {
    auto p = prescription();
    auto setup =
        OptimizationBuilder::builder(&p)
            .fields({0.0, 1.0})
            .mtfFrequencies({10, 20})
            .contrastSampling(2, 4)
            .contrastGoals({OptimizationBuilder::contrast(20, {2.0, 0.5}, {3.0, 0.25})})
            .build();

    auto contrastGoals = goalsOfType<GoalContrast>(setup.goals());
    CHECK_EQ(static_cast<int>(contrastGoals.size()), 2 * 3 * 2 * 4 * 2);
    if (contrastGoals.size() < 2)
        return;
    CHECK_EQ(contrastGoals[0]->_contrast_index, 0);
    CHECK_EQ(contrastGoals[0]->_frequency, 20);
    CHECK_EQ(contrastGoals[0]->_field, 0);
    CHECK_EQ(contrastGoals[0]->_wavelength_index, 0);
    CHECK_EQ(contrastGoals[0]->_sample_index, 0);
    CHECK_EQ(contrastGoals[0]->_orientation, Orientation::SAGITTAL);
    CHECK_EQ(contrastGoals[0]->_weight, 2.0);
    CHECK_EQ(contrastGoals[1]->_orientation, Orientation::TANGENTIAL);
    CHECK_EQ(contrastGoals[1]->_weight, 3.0);

    Analysis *analysis = setup.analysis();
    CHECK_EQ(static_cast<int>(analysis->_contrast_freqs.size()), 1);
    CHECK_EQ(analysis->_contrast_freqs[0], 20);
    CHECK_EQ(analysis->_contrast_num_rings, 2);
    CHECK_EQ(analysis->_contrast_num_spokes, 4);
    CHECK(!analysis->_compute_spots);
    CHECK(!analysis->_compute_mtf);
    CHECK(!analysis->_compute_ray_aberrations);
}

TEST(optim_configures_exit_pupil_contrast_aiming_and_rejects_calibration_combination) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .aimContrastAtExitPupil()
                     .contrastGoals({OptimizationBuilder::contrast(20, {1.0})})
                     .build();

    CHECK(setup.analysis()->_contrast_aim_exit_pupil);
    CHECK(!setup.analysis()->_contrast_calibrate_frequency);

    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .aimContrastAtExitPupil()
                     .calibrateContrastFrequency(true)
                     .contrastGoals({OptimizationBuilder::contrast(20, {1.0})})
                     .build(),
                 redukti::IllegalArgumentException);
}

TEST(optim_validates_contrast_configuration) {
    using redukti::IllegalArgumentException;
    auto p = prescription();
    CHECK_THROWS(OptimizationBuilder::builder(&p).contrastSampling(0, 6),
                 IllegalArgumentException);
    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({20})
                     .contrastGoals({OptimizationBuilder::contrast(20, {1.0})})
                     .build(),
                 IllegalArgumentException);
    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .contrastGoals({OptimizationBuilder::contrast(20, {1.0}),
                                     OptimizationBuilder::contrast(20, {1.0})})
                     .build(),
                 IllegalArgumentException);
}

TEST(optim_computes_contrast_goals_without_computing_spot_mtf) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .dLineOnly(true)
                     .contrastSampling(2, 4)
                     .contrastGoals({OptimizationBuilder::contrast(10, {1.0}),
                                     OptimizationBuilder::contrast(20, {1.0})})
                     .build();

    setup.analysis()->compute();

    CHECK(!setup.analysis()->_spots.has_value());
    CHECK(!setup.analysis()->_mtfs.has_value());
    CHECK_EQ(static_cast<int>(setup.analysis()->_contrasts.size()), 2);
    auto goals = goalsOfType<GoalContrast>(setup.goals());
    bool sawFirst = false, sawSecond = false;
    for (auto *goal : goals) {
        if (goal->_contrast_index == 0 && goal->_frequency == 10)
            sawFirst = true;
        if (goal->_contrast_index == 1 && goal->_frequency == 20)
            sawSecond = true;
        CHECK(std::isfinite(goal->value()));
    }
    CHECK(sawFirst);
    CHECK(sawSecond);
}

TEST(optim_builds_per_ray_spot_goals_whose_sum_matches_rms_spot_radius) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .weighted(false)
                     .gaussianQuadratureSampling(2, 4)
                     .spotDeviationGoals({1.0}, {1.0})
                     .build();

    CHECK_EQ(setup.analysis()->_spot_pattern, SpotOptions::PATTERN_GAUSS_QUADRATURE);
    CHECK_EQ(setup.analysis()->_num_rings, 2);
    CHECK_EQ(setup.analysis()->_num_spokes, 4);
    CHECK(setup.analysis()->_append_failed_spot_rays);
    setup.analysis()->compute();

    auto goals = goalsOfType<GoalSpotDeviation>(setup.goals());
    CHECK_EQ(static_cast<int>(goals.size()), 3 * 2 * 4 * 2);
    if (goals.size() < 2)
        return;
    CHECK_EQ(goals[0]->_orientation, Orientation::X);
    CHECK_EQ(goals[1]->_orientation, Orientation::Y);
    CHECK_EQ(goals[0]->_weight, 1.0);
    CHECK_EQ(goals[1]->_weight, 1.0);
    // Fields are one based on the way in, as for every other field-addressed
    // goal, and stored zero based. An off-by-one here would silently attach the
    // goal to the wrong field rather than fail, so pin it.
    CHECK_EQ(goals[0]->_field, 0);
    CHECK_THROWS(GoalSpotDeviation(setup.analysis(), 0, 0, 0, Orientation::X, 1.0),
                 redukti::IllegalArgumentException);

    double sumOfSquares = 0.0;
    for (auto *goal : goals)
        sumOfSquares += goal->value() * goal->value();
    double rmsFromGoals = std::sqrt(sumOfSquares / 3.0);
    CHECK_CLOSE((*setup.analysis()->_spots)[0].get_mean_radius(), rmsFromGoals, 1.0e-9);
}

TEST(optim_per_ray_spot_goals_improve_spot_merit_during_optimization) {
    auto p = prescription();
    auto setup = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .varyCurvatures({0})
                     .weighted(false)
                     .gaussianQuadratureSampling(2, 4)
                     .spotDeviationGoals({1.0})
                     .build();
    setup.analysis()->compute();
    auto merit = setup.meritFunction(false);
    double initial = merit.getRMS();

    int status = merit.getSolver()->solve();

    CHECK(status > 0);
    CHECK(merit.getRMS() < initial);
}

TEST(optim_validates_mtf_array_lengths_and_measured_frequencies) {
    using redukti::IllegalArgumentException;
    auto p = prescription();
    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .mtfGoals({OptimizationBuilder::mtf(20, {90, 80}, {90, 80})})
                     .build(),
                 IllegalArgumentException);

    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .mtfGoals({OptimizationBuilder::mtf(10, {90}, {90, 80})})
                     .build(),
                 IllegalArgumentException);

    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .spotRmsGoals({10.0})
                     .build(),
                 IllegalArgumentException);

    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0, 1.0})
                     .mtfFrequencies({10})
                     .spotDeviationGoals({1.0})
                     .build(),
                 IllegalArgumentException);

    CHECK_THROWS(OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({10})
                     .spotDeviationGoals({1.0})
                     .hexapolarSampling()
                     .build(),
                 IllegalArgumentException);

    auto other = prescription();
    CHECK_THROWS(OptimizationBuilder::builder(&p).additionalVariables(
                     {std::make_shared<VarThickness>(&other, 0)}),
                 IllegalArgumentException);
}

TEST(optim_builder_products_keep_analysis_alive_after_temporary_setup_dies) {
    auto p = prescription();
    auto merit = OptimizationBuilder::builder(&p)
                     .fields({0.0})
                     .mtfFrequencies({20})
                     .varyCurvatures({0})
                     .weighted(false)
                     .gaussianQuadratureSampling(2, 4)
                     .spotDeviationGoals({1.0})
                     .build()
                     .meritFunction(false);

    auto solver = merit.getSolver();
    CHECK(std::isfinite(merit.getRMS()));
    CHECK(solver != nullptr);

    std::unique_ptr<redukti::optim::Solver> detachedSolver;
    {
        auto shortLivedMerit = OptimizationBuilder::builder(&p)
                                   .fields({0.0})
                                   .mtfFrequencies({20})
                                   .varyCurvatures({0})
                                   .weighted(false)
                                   .gaussianQuadratureSampling(2, 4)
                                   .spotDeviationGoals({1.0})
                                   .build()
                                   .meritFunction(false);
        detachedSolver = shortLivedMerit.getSolver();
    }
    CHECK(detachedSolver->solve() > 0);
}

TEST(optim_variables_and_constraints_reject_invalid_indices) {
    using redukti::IllegalArgumentException;
    auto p = prescription();

    CHECK_THROWS(VarRadius(nullptr, 0), IllegalArgumentException);
    CHECK_THROWS(VarRadius(&p, -1), IllegalArgumentException);
    CHECK_THROWS(VarRadius(&p, static_cast<int>(p._surface_list.size())),
                 IllegalArgumentException);
    CHECK_THROWS(VarThickness(&p, -1), IllegalArgumentException);
    CHECK_THROWS(VarAsphK(&p, static_cast<int>(p._surface_list.size())),
                 IllegalArgumentException);
    CHECK_THROWS(VarAsphCoeff(&p, 0, -1, 1.0), IllegalArgumentException);
    CHECK_THROWS(VarAsphCoeff(&p, 0, 100, 1.0), IllegalArgumentException);
    CHECK_THROWS(VarAsphCoeff(&p, 1, 0, 1.0), IllegalArgumentException);

    p._surface_list[0]._thickness_by_scenario = std::vector<double>{5.0};
    CHECK_THROWS(VarThickness(&p, 0, 1), IllegalArgumentException);

    Analysis analysis(&p, {0.0}, {20});
    CHECK_THROWS(redukti::optim::ConstraintCurvature(&analysis, -1, 1.0),
                 IllegalArgumentException);
    CHECK_THROWS(ConstraintThickness(&analysis,
                                     static_cast<int>(p._surface_list.size()), 1.0),
                 IllegalArgumentException);
    Analysis invalidScenarioAnalysis(&p, {0.0}, {20}, 1);
    CHECK_THROWS(ConstraintThickness(&invalidScenarioAnalysis, 0, 1.0),
                 IllegalArgumentException);
}
