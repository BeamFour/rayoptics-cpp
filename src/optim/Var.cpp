// C++ port of org.redukti.optim.Var and its subclasses
#include "redukti/optim/Var.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"

#include <algorithm>
#include <typeinfo>
#include <cmath>

namespace redukti::optim {

namespace {

spec::SurfaceType &checked_surface(spec::Prescription *prescription, int surfaceId) {
    if (prescription == nullptr)
        throw IllegalArgumentException("prescription must not be null");
    if (surfaceId < 0 || surfaceId >= static_cast<int>(prescription->_surface_list.size()))
        throw IllegalArgumentException("surface index out of range: " + intToString(surfaceId));
    return prescription->_surface_list[static_cast<std::size_t>(surfaceId)];
}

void check_scenario(const spec::SurfaceType &surface, int scenario) {
    if (scenario < 0)
        throw IllegalArgumentException("scenario must be non-negative");
    if (surface._thickness_by_scenario.has_value() &&
        scenario >= static_cast<int>(surface._thickness_by_scenario->size()))
        throw IllegalArgumentException("scenario index out of range: " + intToString(scenario));
}

} // namespace

std::string Var::toString() const { return typeid(*this).name(); }

void Var::set_unscaled_value(double d) {
    if (std::isnan(d))
        throw IllegalArgumentException("NaN value supplied");
    _unscaled_value = d;
    _scaled_value = d * get_scaling_factor();
}

void Var::set_scaled_value(double d) {
    if (std::isnan(d))
        throw IllegalArgumentException("NaN value supplied");
    _scaled_value = d;
    _unscaled_value = mathlib::M::isZero(d) ? 0.0 : (d / get_scaling_factor());
}

VarRadius::VarRadius(spec::Prescription *prescription, int surfaceId)
    : Var(prescription), _surface_id(surfaceId) {
    const auto &surface = checked_surface(prescription, surfaceId);
    // same step as optimr's VarRadius, with a floor for flat/near-flat surfaces
    _d_delta = std::max(std::abs(surface._radius) * 0.001, 1.0e-3);
}

double VarRadius::read_from_prescription() {
    set_unscaled_value(
        _prescription->_surface_list[static_cast<std::size_t>(_surface_id)]._radius);
    return get_scaled_value();
}

void VarRadius::write_to_prescription() {
    _prescription->_surface_list[static_cast<std::size_t>(_surface_id)]._radius =
        get_unscaled_value();
}

std::string VarRadius::toString() const {
    return "Surface ID: " + intToString(_surface_id) +
           " Radius: " + doubleToString(get_unscaled_value());
}

VarThickness::VarThickness(spec::Prescription *prescription, int surfaceId, int scenario)
    : Var(prescription), _surface_id(surfaceId), _scenario(scenario) {
    check_scenario(checked_surface(prescription, surfaceId), scenario);
    _d_delta = 1.0e-3; // matches optimr's VarThickness
}

VarAsphK::VarAsphK(spec::Prescription *prescription, int surfaceId)
    : Var(prescription), _surface_id(surfaceId) {
    (void)checked_surface(prescription, surfaceId);
}

VarAsphCoeff::VarAsphCoeff(spec::Prescription *prescription, int surfaceId, int index,
                           double scalingFactor)
    : Var(prescription), _surface_id(surfaceId), _index(index),
      _scaling_factor(scalingFactor != 0.0 ? scalingFactor : 1.0) {
    const auto &surface = checked_surface(prescription, surfaceId);
    if (index < 0 || !surface._coeffs.has_value() ||
        index >= static_cast<int>(surface._coeffs->size()))
        throw IllegalArgumentException("aspheric coefficient index out of range: " +
                                       intToString(index));
}

double VarThickness::read_from_prescription() {
    auto &surface = _prescription->_surface_list[static_cast<std::size_t>(_surface_id)];
    if (surface._thickness_by_scenario.has_value())
        set_unscaled_value(
            (*surface._thickness_by_scenario)[static_cast<std::size_t>(_scenario)]);
    else
        set_unscaled_value(surface._thickness);
    return get_scaled_value();
}

void VarThickness::write_to_prescription() {
    auto &surface = _prescription->_surface_list[static_cast<std::size_t>(_surface_id)];
    if (surface._thickness_by_scenario.has_value())
        (*surface._thickness_by_scenario)[static_cast<std::size_t>(_scenario)] =
            get_unscaled_value();
    else
        surface._thickness = get_unscaled_value();
}

std::string VarThickness::toString() const {
    return "Surface ID: " + intToString(_surface_id) +
           " Thickness: " + doubleToString(get_unscaled_value());
}

double VarAsphK::read_from_prescription() {
    set_unscaled_value(
        _prescription->_surface_list[static_cast<std::size_t>(_surface_id)]._k);
    return get_scaled_value();
}

void VarAsphK::write_to_prescription() {
    _prescription->_surface_list[static_cast<std::size_t>(_surface_id)]._k =
        get_unscaled_value();
}

std::string VarAsphK::toString() const {
    return "Surface ID: " + intToString(_surface_id) +
           " Asph k: " + doubleToString(get_unscaled_value());
}

double VarAsphCoeff::read_from_prescription() {
    auto &surface = _prescription->_surface_list[static_cast<std::size_t>(_surface_id)];
    set_unscaled_value((*surface._coeffs)[static_cast<std::size_t>(_index)]);
    return get_scaled_value();
}

void VarAsphCoeff::write_to_prescription() {
    auto &surface = _prescription->_surface_list[static_cast<std::size_t>(_surface_id)];
    (*surface._coeffs)[static_cast<std::size_t>(_index)] = get_unscaled_value();
}

std::string VarAsphCoeff::toString() const {
    return "Surface ID: " + intToString(_surface_id) + " Asph Coeff [" +
           intToString(_index) + "]: " + doubleToString(get_unscaled_value()) +
           " scaling factor " + doubleToString(_scaling_factor);
}

} // namespace redukti::optim
