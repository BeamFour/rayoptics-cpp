// C++ port of org.redukti.optim.Constraint and its subclasses
#include "redukti/optim/Constraint.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/rayoptics/elem/profiles/EvenPolynomial.h"
#include "redukti/rayoptics/elem/profiles/RadialPolynomial.h"
#include "redukti/rayoptics/elem/profiles/Spherical.h"
#include "redukti/rayoptics/elem/profiles/SurfaceProfile.h"
#include "redukti/spec/SurfaceType.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace redukti::optim {

namespace profiles = rayoptics::elem::profiles;

namespace {

const spec::SurfaceType &checked_surface(Analysis *analysis, int surfaceId) {
    if (analysis == nullptr || analysis->_prescription == nullptr)
        throw IllegalArgumentException("analysis and prescription must not be null");
    const auto &surfaces = analysis->_prescription->_surface_list;
    if (surfaceId < 0 || surfaceId >= static_cast<int>(surfaces.size()))
        throw IllegalArgumentException("surface index out of range: " + intToString(surfaceId));
    return surfaces[static_cast<std::size_t>(surfaceId)];
}

void check_thickness_scenario(const spec::SurfaceType &surface, int scenario) {
    if (scenario < 0)
        throw IllegalArgumentException("scenario must be non-negative");
    if (surface._thickness_by_scenario.has_value() &&
        scenario >= static_cast<int>(surface._thickness_by_scenario->size()))
        throw IllegalArgumentException("scenario index out of range: " + intToString(scenario));
}

} // namespace

Constraint::Constraint(Analysis *analysis, int surfaceId, double base, double weight)
    : Goal(analysis, base, normalized_weight(base, weight)), _surface_id(surfaceId) {
    if (surfaceId < 0)
        throw IllegalArgumentException("surface must be non-negative");
}

double Constraint::normalized_weight(double base, double weight) {
    if (!std::isfinite(base) || base == 0.0)
        throw IllegalArgumentException(
            "a fractional constraint needs a finite non-zero starting value, got " +
            doubleToString(base));
    if (!std::isfinite(weight) || weight < 0.0)
        throw IllegalArgumentException("weight must be finite and non-negative");
    return weight / (base * base);
}

double Constraint::value() {
    double current = current_value();
    if (std::isfinite(current))
        return current;
    // An undefined parameter (a radius driven to zero, say) must still read as a
    // large finite miss rather than BIGVAL, which would abort the solve outright.
    return _target + std::copysign(UNDEFINED_DEVIATION * std::abs(_target), current);
}

std::string Constraint::toString() {
    return std::string(kind()) + " surface=" + intToString(_surface_id) +
           ", value=" + doubleToString(value()) + ", start=" + doubleToString(_target) +
           ", change=" + doubleToString(fractional_deviation()) +
           ", weight=" + doubleToString(_weight);
}

// ---------------------------------------------------------------------------
// ConstraintCurvature
// ---------------------------------------------------------------------------

double ConstraintCurvature::curvature(Analysis *analysis, int surfaceId) {
    return 1.0 / checked_surface(analysis, surfaceId)._radius;
}

// ---------------------------------------------------------------------------
// ConstraintThickness
// ---------------------------------------------------------------------------

double ConstraintThickness::thickness(Analysis *analysis, int surfaceId) {
    const auto &surface = checked_surface(analysis, surfaceId);
    check_thickness_scenario(surface, analysis->_scenario);
    return surface._thickness_by_scenario.has_value()
               ? (*surface._thickness_by_scenario)[static_cast<std::size_t>(
                     analysis->_scenario)]
               : surface._thickness;
}

// ---------------------------------------------------------------------------
// ConstraintEdgeThickness
// ---------------------------------------------------------------------------

namespace {

/**
 * The same profile the model builder would construct for this surface, so the
 * sag here and the sag the ray tracer sees cannot drift apart.
 */
std::unique_ptr<profiles::SurfaceProfile> make_profile(const spec::SurfaceType &surface) {
    double radius = surface.get_radius_of_curvature();
    if (!surface.is_aspheric())
        return std::make_unique<profiles::Spherical>(radius == 0.0 ? 0.0 : 1.0 / radius);
    if (surface.is_odd_asphere()) {
        auto prof = std::make_unique<profiles::RadialPolynomial>();
        prof->r(radius)->cc(surface.get_cc())->setCoefs(surface.get_aspheric_coeffs());
        return prof;
    }
    auto prof = std::make_unique<profiles::EvenPolynomial>();
    prof->r(radius)->setCc(surface.get_cc())->setCoefs(surface.get_aspheric_coeffs());
    return prof;
}

} // namespace

bool ConstraintEdgeThickness::is_constrainable(Analysis *analysis, int surfaceId) {
    const auto &surfaces = analysis->_prescription->_surface_list;
    if (surfaceId < 0 || surfaceId >= static_cast<int>(surfaces.size()) - 1)
        return false;
    double height = default_height(analysis, surfaceId);
    if (!(height > 0.0))
        return false;
    double gap = edge_gap(analysis, surfaceId, height);
    return std::isfinite(gap) && gap > 0.0;
}

double ConstraintEdgeThickness::default_height(Analysis *analysis, int surfaceId) {
    const auto &surfaces = analysis->_prescription->_surface_list;
    if (surfaceId < 0 || surfaceId >= static_cast<int>(surfaces.size()) - 1)
        throw IllegalArgumentException(
            "an edge gap needs a following surface, but surface " +
            intToString(surfaceId) + " is the last of " +
            intToString(static_cast<int>(surfaces.size())));
    return 0.5 * std::min(surfaces[static_cast<std::size_t>(surfaceId)]
                              .get_diameter_by_scenario(analysis->_scenario),
                          surfaces[static_cast<std::size_t>(surfaceId) + 1]
                              .get_diameter_by_scenario(analysis->_scenario));
}

double ConstraintEdgeThickness::checked_height(double height) {
    if (!std::isfinite(height) || height <= 0.0)
        throw IllegalArgumentException("edge height must be finite and positive, got " +
                                       doubleToString(height));
    return height;
}

double ConstraintEdgeThickness::edge_gap(Analysis *analysis, int surfaceId,
                                         double height) {
    const auto &surfaces = analysis->_prescription->_surface_list;
    double thickness = surfaces[static_cast<std::size_t>(surfaceId)]
                           .get_thickness_by_scenario(analysis->_scenario);
    try {
        return thickness +
               make_profile(surfaces[static_cast<std::size_t>(surfaceId) + 1])
                   ->sag(0.0, height) -
               make_profile(surfaces[static_cast<std::size_t>(surfaceId)])
                   ->sag(0.0, height);
    } catch (const RuntimeException &) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

std::string ConstraintEdgeThickness::toString() {
    return Constraint::toString() + ", height=" + doubleToString(_height);
}

} // namespace redukti::optim
