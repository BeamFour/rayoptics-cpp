// C++ port of the org.redukti.optim.Goal subclasses
#include "redukti/optim/Goals.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/LMLSolver.h"
#include "redukti/mathlib/M.h"
#include "redukti/optim/ParaxHelper.h"
#include "redukti/rayoptics/util/Lists.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>
#include <typeinfo>

namespace redukti::optim {

namespace {

using mathlib::LMLSolver;
namespace Orientation = rayoptics::util::Orientation;

/**
 * The Java dereferences a null array here and raises NullPointerException,
 * which LMDerMeritFunction catches along with every other RuntimeException.
 * There is no NullPointerException in the ported hierarchy, so this raises the
 * nearest equivalent: still a RuntimeException, so every catch site behaves the
 * same. It can only fire if a goal was built without its analysis enabled,
 * which OptimizationBuilder does not do.
 */
[[noreturn]] void analysis_not_computed(const char *what) {
    throw IllegalStateException(std::string(what) +
                                " were not computed for this analysis");
}

/**
 * Java indexes arrays here and raises ArrayIndexOutOfBoundsException when the
 * goal addresses a field the analysis did not produce - a RuntimeException that
 * LMDerMeritFunction catches and turns into a BIGVAL residual. A bare
 * vector::operator[] would instead be undefined behaviour, so the bound is
 * checked and the nearest ported exception raised.
 */
std::size_t checked(std::size_t size, int index, const char *what) {
    if (index < 0 || static_cast<std::size_t>(index) >= size)
        throw IndexOutOfBoundsException(std::string(what) + " index " +
                                        intToString(index) + " out of range for " +
                                        intToString(static_cast<int>(size)));
    return static_cast<std::size_t>(index);
}

} // namespace

std::string Goal::toString() { return typeid(*this).name(); }

// ---------------------------------------------------------------------------
// GoalParax
// ---------------------------------------------------------------------------

double GoalParax::value() {
    // _pfo is null in the Java until compute() has run, and empty here.
    return _analysis->_pfo[checked(_analysis->_pfo.size(), _parax_id, "paraxial")];
}

std::string GoalParax::toString() {
    return std::string(ParaxHelper::Names[_parax_id]) + " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalSpotRMS / GoalSpotMaxRadius
// ---------------------------------------------------------------------------

double GoalSpotRMS::value() {
    if (!_analysis->_spots.has_value())
        analysis_not_computed("spots");
    return (*_analysis->_spots)[checked(_analysis->_spots->size(), _field - 1, "spot field")]
        .get_mean_radius();
}

std::string GoalSpotRMS::toString() {
    return "SpotRMS field=" + intToString(_field) + ", target=" + doubleToString(_target) +
           ", weight=" + doubleToString(_weight) + " = " + doubleToString(value());
}

double GoalSpotMaxRadius::value() {
    if (!_analysis->_spots.has_value())
        analysis_not_computed("spots");
    return (*_analysis->_spots)[checked(_analysis->_spots->size(), _field - 1, "spot field")]
        .get_max_radius();
}

std::string GoalSpotMaxRadius::toString() {
    return "SpotMaxRadius field=" + intToString(_field) +
           ", target=" + doubleToString(_target) + ", weight=" + doubleToString(_weight) +
           " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalSpotDeviation
// ---------------------------------------------------------------------------

GoalSpotDeviation::GoalSpotDeviation(Analysis *analysis, int field, int wavelength_index,
                                     int sample_index, int orientation, double weight)
    : Goal(analysis, 0.0, weight), _field(field - 1), _wavelength_index(wavelength_index),
      _sample_index(sample_index), _orientation(orientation) {
    if (field < 1)
        throw IllegalArgumentException("field is one based and must be positive");
    if (wavelength_index < 0 || sample_index < 0)
        throw IllegalArgumentException("indices must be non-negative");
    Orientation::checked(orientation);
}

double GoalSpotDeviation::value() {
    if (!_analysis->_spots.has_value() ||
        _field >= static_cast<int>(_analysis->_spots->size()))
        return LMLSolver::BIGVAL;
    // The Java also tests the element for null; a Java array of objects can hold
    // one, a std::vector element cannot, and SpotAnalysisResult never stores one.
    const auto &field = (*_analysis->_spots)[static_cast<std::size_t>(_field)];
    if (_wavelength_index >= static_cast<int>(field.intercepts.size()))
        return LMLSolver::BIGVAL;
    const auto &intercepts =
        field.intercepts[static_cast<std::size_t>(_wavelength_index)];
    if (_sample_index >= static_cast<int>(intercepts.x.size()) ||
        !intercepts.valid[static_cast<std::size_t>(_sample_index)])
        return LMLSolver::BIGVAL;
    double deviation = _orientation == Orientation::X
                           ? intercepts.x[static_cast<std::size_t>(_sample_index)]
                           : intercepts.y[static_cast<std::size_t>(_sample_index)];
    // SpotAnalysis stores system units (normally mm); public spot radii and
    // optimization targets use microns.
    return 1000.0 *
           std::sqrt(intercepts.weights[static_cast<std::size_t>(_sample_index)]) *
           deviation;
}

std::string GoalSpotDeviation::toString() {
    return "SpotDeviation field=" + intToString(_field) +
           ", wavelength=" + intToString(_wavelength_index) +
           ", sample=" + intToString(_sample_index) +
           ", orientation=" + (_orientation == Orientation::X ? "x" : "y") +
           ", target=" + doubleToString(_target) + ", weight=" + doubleToString(_weight) +
           " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalRayAberration
// ---------------------------------------------------------------------------

GoalRayAberration::GoalRayAberration(Analysis *analysis, int field, int orientation,
                                     int pos, double wvl, double target, double weight)
    : Goal(analysis, target, weight), _field(field - 1),
      _orientation(Orientation::checked(orientation)), _pos(pos), _wvl(wvl) {
    // The Java stores pos before adjusting it, so a negative _pos survives into
    // value() where Lists::get resolves it Python-style from the end. The local
    // adjustment below exists only to range-check.
    int checked = pos;
    if (checked < 0)
        checked += Analysis::NUM_TRANSVERSE_RAYS;
    if (checked < 0 || checked >= Analysis::NUM_TRANSVERSE_RAYS)
        throw IllegalArgumentException(
            "position out of range, max number of rays is " +
            intToString(Analysis::NUM_TRANSVERSE_RAYS));
}

double GoalRayAberration::value() {
    if (!_analysis->_ray_aberrations.has_value())
        analysis_not_computed("ray aberrations");
    const auto *fans = _analysis->_ray_aberrations->get_fans(_field, _orientation, _wvl);
    if (fans != nullptr && _pos < static_cast<int>(fans->fan_x.size())) {
        const auto &result = rayoptics::util::Lists::get(fans->fan_y, _pos);
        return result.has_value() && std::isfinite(*result) ? *result : LMLSolver::BIGVAL;
    }
    return LMLSolver::BIGVAL;
}

std::string GoalRayAberration::toString() {
    return "RayAberration field=" + intToString(_field) +
           ", orientation=" + Orientation::name(_orientation) +
           ", pos=" + intToString(_pos) + ", target=" + doubleToString(_target) +
           ", weight=" + doubleToString(_weight) + " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalGeoMTF
// ---------------------------------------------------------------------------

GoalGeoMTF::GoalGeoMTF(Analysis *analysis, int field, int orientation, int freq,
                       double target, double weight)
    : Goal(analysis, target, weight), _freq(freq),
      _orientation(Orientation::checked(orientation)), _field(field) {}

double GoalGeoMTF::value() {
    if (!_analysis->_mtfs.has_value())
        analysis_not_computed("MTFs");
    const auto &mtfs = *_analysis->_mtfs;
    for (std::size_t i = 0; i < mtfs.size(); i++) {
        if (mtfs[i].freq == _freq) {
            const auto &byField = _orientation == Orientation::TANGENTIAL
                                      ? mtfs[i].tan_mtf_by_field
                                      : mtfs[i].sag_mtf_by_field;
            return byField[checked(byField.size(), _field - 1, "MTF field")];
        }
    }
    throw IllegalArgumentException("");
}

std::string GoalGeoMTF::toString() {
    return "GeoMTF field=" + intToString(_field) +
           ", orientation=" + Orientation::name(_orientation) +
           ", freq=" + intToString(_freq) + ", target=" + doubleToString(_target) +
           ", weight=" + doubleToString(_weight) + ", value=" + doubleToString(value()) +
           "}";
}

// ---------------------------------------------------------------------------
// GoalMTFProxy
// ---------------------------------------------------------------------------

GoalMTFProxy::GoalMTFProxy(Analysis *analysis, int field, int orientation, int pos,
                           double wvl, double freq, double target, double weight)
    : Goal(analysis, target, weight), _field(field - 1),
      _orientation(Orientation::checked(orientation)), _pos(pos), _wvl(wvl), _freq(freq) {}

double GoalMTFProxy::value() {
    if (!_analysis->_ray_aberrations.has_value())
        analysis_not_computed("ray aberrations");
    const auto *fans = _analysis->_ray_aberrations->get_fans(_field, _orientation, _wvl);
    if (fans != nullptr) {
        // The Java wraps Lists.get in a catch-all and returns 1.0. Two distinct
        // Java failures land there: an out-of-range index throws
        // IndexOutOfBoundsException, and a null element throws
        // NullPointerException when it is unboxed to a double. Lists::get here
        // neither range-checks nor unboxes, so both have to be tested for.
        int index = _pos;
        if (index < 0)
            index += static_cast<int>(fans->fan_y.size());
        if (index < 0 || index >= static_cast<int>(fans->fan_y.size()))
            return 1.0;
        const auto &entry = fans->fan_y[static_cast<std::size_t>(index)];
        if (!entry.has_value())
            return 1.0;
        return std::sin(mathlib::M::PI * _freq * *entry);
    }
    throw IllegalArgumentException("Invalid field, orientation or position");
}

std::string GoalMTFProxy::toString() {
    return "MTFProxy field=" + intToString(_field) +
           ", orientation=" + Orientation::name(_orientation) +
           ", pos=" + intToString(_pos) + ", target=" + doubleToString(_target) +
           ", weight=" + doubleToString(_weight) + " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalContrast
// ---------------------------------------------------------------------------

GoalContrast::GoalContrast(Analysis *analysis, int contrast_index, int frequency,
                           int field, int wavelength_index, int sample_index,
                           int orientation, double weight)
    : Goal(analysis, 0.0, weight), _contrast_index(contrast_index), _frequency(frequency),
      _field(field - 1), _wavelength_index(wavelength_index),
      _sample_index(sample_index), _orientation(orientation) {
    if (contrast_index < 0)
        throw IllegalArgumentException("contrast index must be non-negative");
    if (field < 1)
        throw IllegalArgumentException("field is one based and must be positive");
    if (wavelength_index < 0 || sample_index < 0)
        throw IllegalArgumentException("indices must be non-negative");
    Orientation::checked(orientation);
}

double GoalContrast::value() {
    if (_contrast_index >= static_cast<int>(_analysis->_contrasts.size()))
        return LMLSolver::BIGVAL;
    const auto &contrast =
        _analysis->_contrasts[static_cast<std::size_t>(_contrast_index)];
    if (_field >= static_cast<int>(contrast.fields.size()))
        return LMLSolver::BIGVAL;
    const auto &wavelengths =
        contrast.fields[static_cast<std::size_t>(_field)].wavelengths;
    if (_wavelength_index >= static_cast<int>(wavelengths.size()))
        return LMLSolver::BIGVAL;
    const auto &wavelength = wavelengths[static_cast<std::size_t>(_wavelength_index)];
    if (_sample_index >= static_cast<int>(wavelength.samples.size()))
        return LMLSolver::BIGVAL;
    if (!wavelength.samples[static_cast<std::size_t>(_sample_index)].valid)
        return LMLSolver::BIGVAL;
    // Read through the block rather than the sample: the residual carries the
    // block's constant offset, which is zero unless residual centering is enabled.
    return _orientation == Orientation::SAGITTAL
               ? wavelength.sagittalResidual(_sample_index)
               : wavelength.tangentialResidual(_sample_index);
}

std::string GoalContrast::toString() {
    return "Contrast index=" + intToString(_contrast_index) +
           ", frequency=" + intToString(_frequency) + ", field=" + intToString(_field) +
           ", wavelength=" + intToString(_wavelength_index) +
           ", sample=" + intToString(_sample_index) +
           ", orientation=" + Orientation::name(_orientation) +
           ", weight=" + doubleToString(_weight) + " = " + doubleToString(value());
}

// ---------------------------------------------------------------------------
// GoalContrastBalance
// ---------------------------------------------------------------------------

GoalContrastBalance::GoalContrastBalance(Analysis *analysis, int contrast_index,
                                         int frequency, int field,
                                         const std::vector<double> &wavelengthWeights,
                                         double sagittalWeight, double tangentialWeight,
                                         double weight)
    : Goal(analysis, 0.0, weight), _contrast_index(contrast_index), _frequency(frequency),
      _field(field - 1), _wavelength_weights(wavelengthWeights),
      _sagittal_weight(sagittalWeight), _tangential_weight(tangentialWeight) {
    if (contrast_index < 0)
        throw IllegalArgumentException("contrast index must be non-negative");
    if (field < 1)
        throw IllegalArgumentException("field is one based and must be positive");
    if (wavelengthWeights.empty())
        throw IllegalArgumentException("balance needs at least one wavelength weight");
    if (!std::isfinite(sagittalWeight) || sagittalWeight < 0.0 ||
        !std::isfinite(tangentialWeight) || tangentialWeight < 0.0)
        throw IllegalArgumentException(
            "contrast orientation weights must be finite and non-negative");
}

double GoalContrastBalance::value() {
    if (_contrast_index >= static_cast<int>(_analysis->_contrasts.size()))
        return LMLSolver::BIGVAL;
    const auto &contrast =
        _analysis->_contrasts[static_cast<std::size_t>(_contrast_index)];
    if (_field >= static_cast<int>(contrast.fields.size()))
        return LMLSolver::BIGVAL;
    const auto &wavelengths =
        contrast.fields[static_cast<std::size_t>(_field)].wavelengths;
    if (wavelengths.size() > _wavelength_weights.size())
        return LMLSolver::BIGVAL;

    double difference = 0.0;
    bool sampled = false;
    for (std::size_t wi = 0; wi < wavelengths.size(); wi++) {
        const auto &block = wavelengths[wi];
        double sagittal = 0.0;
        double tangential = 0.0;
        for (std::size_t i = 0; i < block.samples.size(); i++) {
            if (!block.samples[i].valid)
                continue;
            double rs = block.sagittalResidual(static_cast<int>(i));
            double rt = block.tangentialResidual(static_cast<int>(i));
            if (!std::isfinite(rs) || !std::isfinite(rt))
                continue;
            sagittal += rs * rs;
            tangential += rt * rt;
            sampled = true;
        }
        difference += _wavelength_weights[wi] *
                      (_sagittal_weight * sagittal - _tangential_weight * tangential);
    }
    if (!sampled || !std::isfinite(difference))
        return LMLSolver::BIGVAL;
    return difference;
}

std::string GoalContrastBalance::toString() {
    return "ContrastBalance frequency=" + intToString(_frequency) +
           ", field=" + intToString(_field) + ", weight=" + doubleToString(_weight) +
           " = " + doubleToString(value());
}

} // namespace redukti::optim
