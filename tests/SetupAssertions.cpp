// C++ port of org.redukti.optim.SetupAssertions
#include "SetupAssertions.h"

#include "TestHarness.h"

#include "redukti/Text.h"
#include "redukti/optim/Constraint.h"
#include "redukti/optim/Goals.h"
#include "redukti/util/Args.h"

namespace redukti::test {

using namespace redukti::optim;

namespace {

std::string field(const char *name, const std::string &value) {
    return std::string(name) + "=" + value + " ";
}

std::string field(const char *name, int value) {
    return field(name, intToString(value));
}

std::string field(const char *name, double value) {
    return field(name, doubleToString(value));
}

std::string field(const char *name, bool value) {
    return field(name, std::string(value ? "true" : "false"));
}

std::string field(const char *name, const std::vector<double> &values) {
    std::string text = "[";
    for (std::size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            text += ", ";
        text += doubleToString(values[i]);
    }
    return field(name, text + "]");
}

std::string field(const char *name, const std::vector<int> &values) {
    std::string text = "[";
    for (std::size_t i = 0; i < values.size(); i++) {
        if (i > 0)
            text += ", ";
        text += intToString(values[i]);
    }
    return field(name, text + "]");
}

/** The values every Var holds, whatever its type. */
std::string base(const Var &variable) {
    return field("_unscaled_value", variable.get_unscaled_value()) +
           field("_scaled_value", variable.get_scaled_value());
}

/** The values every Goal holds, whatever its type. */
std::string base(const Goal &goal) {
    return field("_weight", goal._weight) + field("_target", goal._target);
}

} // namespace

std::string describe(const Var &variable) {
    if (const auto *v = dynamic_cast<const VarRadius *>(&variable))
        return "VarRadius{" + base(*v) + field("_surface_id", v->_surface_id) + "}";
    if (const auto *v = dynamic_cast<const VarThickness *>(&variable))
        return "VarThickness{" + base(*v) + field("_surface_id", v->_surface_id) +
               field("_scenario", v->_scenario) + "}";
    if (const auto *v = dynamic_cast<const VarAsphK *>(&variable))
        return "VarAsphK{" + base(*v) + field("_surface_id", v->_surface_id) + "}";
    if (const auto *v = dynamic_cast<const VarAsphCoeff *>(&variable))
        return "VarAsphCoeff{" + base(*v) + field("_surface_id", v->_surface_id) +
               field("_index", v->_index) + field("_scaling_factor", v->_scaling_factor) + "}";
    return "Var{" + base(variable) + "}";
}

std::string describe(Goal &goal) {
    if (const auto *g = dynamic_cast<const ConstraintCurvature *>(&goal))
        return "ConstraintCurvature{" + base(*g) + field("_surface_id", g->_surface_id) + "}";
    if (const auto *g = dynamic_cast<const ConstraintThickness *>(&goal))
        return "ConstraintThickness{" + base(*g) + field("_surface_id", g->_surface_id) + "}";
    if (const auto *g = dynamic_cast<const ConstraintEdgeThickness *>(&goal))
        return "ConstraintEdgeThickness{" + base(*g) + field("_surface_id", g->_surface_id) +
               field("_height", g->_height) + "}";
    if (const auto *g = dynamic_cast<const GoalParax *>(&goal))
        return "GoalParax{" + base(*g) + field("_parax_id", g->_parax_id) + "}";
    if (const auto *g = dynamic_cast<const GoalSpotRMS *>(&goal))
        return "GoalSpotRMS{" + base(*g) + field("_field", g->_field) + "}";
    if (const auto *g = dynamic_cast<const GoalSpotMaxRadius *>(&goal))
        return "GoalSpotMaxRadius{" + base(*g) + field("_field", g->_field) + "}";
    if (const auto *g = dynamic_cast<const GoalSpotDeviation *>(&goal))
        return "GoalSpotDeviation{" + base(*g) + field("_field", g->_field) +
               field("_wavelength_index", g->_wavelength_index) +
               field("_sample_index", g->_sample_index) +
               field("_orientation", g->_orientation) + "}";
    if (const auto *g = dynamic_cast<const GoalRayAberration *>(&goal))
        return "GoalRayAberration{" + base(*g) + field("_field", g->_field) +
               field("_orientation", g->_orientation) + field("_pos", g->_pos) +
               field("_wvl", g->_wvl) + "}";
    if (const auto *g = dynamic_cast<const GoalGeoMTF *>(&goal))
        return "GoalGeoMTF{" + base(*g) + field("_freq", g->_freq) +
               field("_orientation", g->_orientation) + field("_field", g->_field) + "}";
    if (const auto *g = dynamic_cast<const GoalContrast *>(&goal))
        return "GoalContrast{" + base(*g) + field("_contrast_index", g->_contrast_index) +
               field("_frequency", g->_frequency) + field("_field", g->_field) +
               field("_wavelength_index", g->_wavelength_index) +
               field("_sample_index", g->_sample_index) +
               field("_orientation", g->_orientation) + "}";
    if (const auto *g = dynamic_cast<const GoalContrastBalance *>(&goal))
        return "GoalContrastBalance{" + base(*g) +
               field("_contrast_index", g->_contrast_index) +
               field("_frequency", g->_frequency) + field("_field", g->_field) + "}";
    return "Goal{" + base(goal) + "}";
}

std::string describe(const Analysis &analysis) {
    return "Analysis{" + field("_fields", analysis._fields) + field("_freqs", analysis._freqs) +
           field("_contrast_freqs", analysis._contrast_freqs) +
           field("_scenario", analysis._scenario) +
           field("_spot_pattern", analysis._spot_pattern) +
           field("_num_rays", analysis._num_rays) + field("_num_rings", analysis._num_rings) +
           field("_num_spokes", analysis._num_spokes) +
           field("_inner_pupil_radius", analysis._inner_pupil_radius) +
           field("_append_failed_spot_rays", analysis._append_failed_spot_rays) +
           field("_check_spot_apertures", analysis._check_spot_apertures) +
           field("_contrast_num_rings", analysis._contrast_num_rings) +
           field("_contrast_num_spokes", analysis._contrast_num_spokes) +
           field("_contrast_calibrate_frequency", analysis._contrast_calibrate_frequency) +
           field("_contrast_aim_exit_pupil", analysis._contrast_aim_exit_pupil) +
           field("_contrast_center_residuals", analysis._contrast_center_residuals) +
           field("_compute_spots", analysis._compute_spots) +
           field("_compute_ray_aberrations", analysis._compute_ray_aberrations) +
           field("_compute_mtf", analysis._compute_mtf) +
           field("_freeze_vignetting", analysis._freeze_vignetting) +
           field("_vig_type", util::Args::vig_type_name(analysis._vig_type)) + "}";
}

void assertSameSetup(const char *what, OptimizationBuilder::OptimizationSetup &expected,
                     OptimizationBuilder::OptimizationSetup &actual) {
    auto expectedVariables = expected.variables();
    auto actualVariables = actual.variables();
    if (expectedVariables.size() != actualVariables.size())
        ::redukti::test::reportFailure(__FILE__, __LINE__,
                                       std::string(what) + ": " +
                                           std::to_string(actualVariables.size()) +
                                           " variables, want " +
                                           std::to_string(expectedVariables.size()));
    for (std::size_t i = 0; i < expectedVariables.size() && i < actualVariables.size(); i++)
        if (describe(*expectedVariables[i]) != describe(*actualVariables[i]))
            ::redukti::test::reportFailure(__FILE__, __LINE__,
                                           std::string(what) + ": variable " +
                                               std::to_string(i) + " got " +
                                               describe(*actualVariables[i]) + " want " +
                                               describe(*expectedVariables[i]));
    auto expectedGoals = expected.goals();
    auto actualGoals = actual.goals();
    if (expectedGoals.size() != actualGoals.size())
        ::redukti::test::reportFailure(__FILE__, __LINE__,
                                       std::string(what) + ": " +
                                           std::to_string(actualGoals.size()) +
                                           " goals, want " +
                                           std::to_string(expectedGoals.size()));
    for (std::size_t i = 0; i < expectedGoals.size() && i < actualGoals.size(); i++)
        if (describe(*expectedGoals[i]) != describe(*actualGoals[i]))
            ::redukti::test::reportFailure(__FILE__, __LINE__,
                                           std::string(what) + ": goal " + std::to_string(i) +
                                               " got " + describe(*actualGoals[i]) +
                                               " want " + describe(*expectedGoals[i]));
    if (describe(*expected.analysis()) != describe(*actual.analysis()))
        ::redukti::test::reportFailure(__FILE__, __LINE__,
                                       std::string(what) + ": analysis got " +
                                           describe(*actual.analysis()) + " want " +
                                           describe(*expected.analysis()));
}

} // namespace redukti::test
