// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of ContrastOptions, ContrastAnalysisResult and ContrastAnalysis,
// including the contrast trace that was in rayoptics::raytr::Trace before the Java
// refactor moved it here.
#include "redukti/analysis/ContrastAnalysis.h"

#include "redukti/Exceptions.h"
#include "redukti/analysis/ExitPupilAiming.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>

namespace redukti::analysis {

namespace Orientation = rayoptics::util::Orientation;

using mathlib::Vector2;
using rayoptics::raytr::PupilType;
using rayoptics::raytr::RayPkg;
using rayoptics::raytr::Trace;

// ---------------------------------------------------------------------------
// ContrastOptions
// ---------------------------------------------------------------------------

ContrastOptions::ContrastOptions(double spatialFrequency_) {
    if (!std::isfinite(spatialFrequency_) || spatialFrequency_ < 0.0) {
        throw IllegalArgumentException(
            "Spatial frequency must be finite and non-negative");
    }
    this->spatialFrequency = spatialFrequency_;
    // Contrast samples already occupy the common vignetted-pupil overlap.
    // Do not additionally reject them against surface apertures by default.
    this->traceOptions.check_apertures = false;
}

ContrastOptions ContrastOptions::with_frequency(double spatialFrequency_) const {
    ContrastOptions result(spatialFrequency_);
    result.numRings = numRings;
    result.numSpokes = numSpokes;
    result.traceOptions = traceOptions.copy();
    result.calibrateFrequency = calibrateFrequency;
    result.aimExitPupil = aimExitPupil;
    result.centerResiduals = centerResiduals;
    return result;
}

ContrastOptions &ContrastOptions::num_rings(int value) {
    if (value < 1)
        throw IllegalArgumentException("Number of rings must be at least 1");
    numRings = value;
    return *this;
}

ContrastOptions &ContrastOptions::num_spokes(std::optional<int> value) {
    if (value.has_value() && *value < 3)
        throw IllegalArgumentException("Number of spokes must be at least 3");
    numSpokes = value;
    return *this;
}

ContrastOptions &ContrastOptions::trace_options(const rayoptics::raytr::TraceOptions &value) {
    // The Java rejects a null here; a reference cannot be null, so the check
    // has nothing left to test.
    traceOptions = value;
    return *this;
}

/**
 * Subtract the constant part of each wavefront-difference block, so the residuals
 * carry the variance the OTF modulus depends on rather than the un-centred second
 * moment. Off by default because it changes every contrast residual. See
 * ContrastAnalysis#center_residuals(ContrastAnalysisResult, int).
 */
ContrastOptions &ContrastOptions::center_residuals(bool value) {
    centerResiduals = value;
    return *this;
}

/** Whether traced contrast rays are also rejected by surface apertures. */
ContrastOptions &ContrastOptions::check_apertures(bool value) {
    traceOptions.check_apertures = value;
    return *this;
}

/**
 * Correct the pupil shift so the sampled pair realises the requested spatial
 * frequency in image space.
 *
 * The shift is applied in entrance-pupil coordinates but is derived from an
 * exit-pupil relation, so pupil aberration makes the realised frequency fall short -
 * measured around 8% low at full field on an f/2 lens, and worsening with field.
 * Enabling this measures the shortfall with one probe pair per field, wavelength and
 * direction, and scales the shift to compensate. It costs four extra rays per field
 * and wavelength.
 *
 * Off by default: it changes the sampled frequency and therefore every contrast
 * residual, so it is opt-in until the correct exit-pupil treatment is settled.
 */
ContrastOptions &ContrastOptions::calibrate_frequency(bool value) {
    calibrateFrequency = value;
    return *this;
}

/**
 * Inverse-aim every displaced contrast ray so its separation from the reference ray
 * is the requested vector on the exit-pupil reference sphere.
 *
 * This is the physically direct alternative to block calibration. It retains the
 * entrance-pupil quadrature for the reference rays, but does not approximate the
 * partner ray with a rigid entrance-pupil displacement. It is unavailable for afocal
 * systems and cannot be combined with #calibrate_frequency(boolean).
 */
ContrastOptions &ContrastOptions::aim_exit_pupil(bool value) {
    aimExitPupil = value;
    return *this;
}

// ---------------------------------------------------------------------------
// ContrastAnalysisResult
// ---------------------------------------------------------------------------

double ContrastAnalysisResult::Sample::sagittalResidual() const {
    return std::sqrt(weight) * sagittalOpdDifference;
}

double ContrastAnalysisResult::Sample::tangentialResidual() const {
    return std::sqrt(weight) * tangentialOpdDifference;
}

ContrastAnalysisResult::WavelengthResult
ContrastAnalysisResult::WavelengthResult::withOffsets(double sagittal,
                                                      double tangential) const {
    return WavelengthResult(wavelength, normalizedPupilShift, samples, sagittal,
                            tangential);
}

/**
 * The optimizer residual for a sample: the sample's own residual less the
 * constant part, which is a pure image displacement and costs no MTF.
 *
 * Note this is sqrt(w) * (dW - offset), not
 * sqrt(w) * dW - offset.
 */
double ContrastAnalysisResult::WavelengthResult::sagittalResidual(int index) const {
    const auto &sample = samples[static_cast<std::size_t>(index)];
    return std::sqrt(sample.weight) * (sample.sagittalOpdDifference - sagittalOffset);
}

double ContrastAnalysisResult::WavelengthResult::tangentialResidual(int index) const {
    const auto &sample = samples[static_cast<std::size_t>(index)];
    return std::sqrt(sample.weight) * (sample.tangentialOpdDifference - tangentialOffset);
}

// ---------------------------------------------------------------------------
// ContrastAnalysis
// ---------------------------------------------------------------------------

ContrastAnalysisResult ContrastAnalysis::eval(rayoptics::optical::OpticalModel *opticalModel,
                                              const ContrastOptions &options) {
    if (options.aimExitPupil && options.calibrateFrequency) {
        throw IllegalArgumentException(
            "Exit-pupil aiming cannot be combined with contrast frequency calibration");
    }
    ContrastAnalysisResult result(options.spatialFrequency);
    auto &fields = opticalModel->optical_spec->fov->fields;
    const auto &wavelengths = opticalModel->optical_spec->wvls->wavelengths;
    for (std::size_t fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++) {
        std::vector<ContrastAnalysisResult::WavelengthResult> wavelengthResults;
        for (std::size_t wavelengthIndex = 0; wavelengthIndex < wavelengths.size();
             wavelengthIndex++) {
            double wavelength = wavelengths[wavelengthIndex];
            // Convert frequency to a nominal shift in normalized entrance-pupil radius units.
            // The paper requires the shift to be calculated at exit pupil, so this
            // strategy is only good when there is little pupil aberration.
            // In practice, it seems to work well except for outer-most fields / wide lenses
            // For those scenarios we have additional calibration below but these are optional.
            double shift = normalized_entry_pupil_shift(opticalModel, wavelength,
                                                        options.spatialFrequency);
            double sagittalShift = shift;
            double tangentialShift = shift;
            if (options.calibrateFrequency) {
                // This is the first calibration:
                // Block calibration applies one scale per orientation for this field and wavelength.
                // This traces one probe pair per field, wavelength and direction at ±shift/2,
                // measures the realised direction-cosine difference, and scales the entrance-pupil
                // shift by λν / realised.
                // Four extra rays per field and wavelength.
                // This is a relatively cheap entrance to exit pupil correction that gets us
                // 90% of the required correction.
                sagittalShift *= exit_pupil_frequency_calibration(
                    opticalModel, *fields[fieldIndex], wavelength, shift, Orientation::X,
                    options);
                tangentialShift *= exit_pupil_frequency_calibration(
                    opticalModel, *fields[fieldIndex], wavelength, shift, Orientation::Y,
                    options);
            }
            // For greater accuracy in calibration,
            // Direct exit-pupil aiming is the alternative: when enabled, trace_triplets() adjusts each
            // partner separately. This is an expensive iterative adjustment of the shift amount.
            // It is enabled by options.aimExitPupil. By default, aimExitPupil is off.
            ContrastTraceCallback<ContrastAnalysisResult::Sample> callback =
                [opticalModel, &options](const ContrastRayTriplet &rays,
                                         rayoptics::specs::Field &field, double tracedWavelength,
                                         double focus) {
                    return sample(opticalModel, rays, field, tracedWavelength, focus,
                                  options);
                };
            auto traced =
                trace_by_wavelength<ContrastAnalysisResult::Sample>(
                    opticalModel,
                    callback, static_cast<int>(fieldIndex),
                    static_cast<int>(wavelengthIndex), options.numRings,
                    options.numSpokes, Vector2(sagittalShift, 0.0),
                    Vector2(0.0, tangentialShift), options.spatialFrequency,
                    options.traceOptions, options.aimExitPupil);
            wavelengthResults.push_back(ContrastAnalysisResult::WavelengthResult(
                wavelength, shift, traced[0].samples));
        }
        result.fields.push_back(ContrastAnalysisResult::FieldResult(
            fields[fieldIndex].get(), std::move(wavelengthResults)));
    }
    if (options.centerResiduals)
        center_residuals(result, opticalModel->optical_spec->wvls->reference_wvl);
    return result;
}

void ContrastAnalysis::center_residuals(ContrastAnalysisResult &result,
                                        int referenceWavelengthIndex) {
    for (auto &field : result.fields) {
        auto &wavelengths = field.wavelengths;
        int reference = referenceWavelengthIndex;
        if (reference < 0 || reference >= static_cast<int>(wavelengths.size()))
            continue;
        double sagittalOffset =
            weighted_mean(wavelengths[static_cast<std::size_t>(reference)],
                          Orientation::SAGITTAL);
        double tangentialOffset =
            weighted_mean(wavelengths[static_cast<std::size_t>(reference)],
                          Orientation::TANGENTIAL);
        if (!std::isfinite(sagittalOffset) || !std::isfinite(tangentialOffset))
            continue;
        for (std::size_t i = 0; i < wavelengths.size(); i++)
            wavelengths[i] = wavelengths[i].withOffsets(sagittalOffset, tangentialOffset);
    }
}

double ContrastAnalysis::weighted_mean(
    const ContrastAnalysisResult::WavelengthResult &wavelength, int orientation) {
    double weightedSum = 0.0;
    double weightSum = 0.0;
    for (const auto &sample : wavelength.samples) {
        if (!sample.valid)
            continue;
        double difference = orientation == Orientation::SAGITTAL
                                ? sample.sagittalOpdDifference
                                : sample.tangentialOpdDifference;
        if (!std::isfinite(difference))
            continue;
        weightedSum += sample.weight * difference;
        weightSum += sample.weight;
    }
    return weightSum > 0.0 ? weightedSum / weightSum : 0.0;
}

/**
 * Convert frequency to a displacement in normalized <em>entrance</em> pupil-radius
 * units. The displacement is 2 at the incoherent diffraction cutoff 1/(lambda F#).
 *
 * Note the mismatch the name records: the magnitude comes from an exit-pupil
 * relation, since that is where the OTF autocorrelation is defined, but it is applied
 * in entrance-pupil coordinates because that is what REL_PUPIL means. Those
 * coincide only where the pupil imaging is aberration free. See
 * #exit_pupil_frequency_calibration, which measures and corrects the
 * difference when ContrastOptions#calibrate_frequency(boolean) is enabled.
 */
double ContrastAnalysis::normalized_entry_pupil_shift(rayoptics::optical::OpticalModel *opticalModel,
                                                      double wavelength,
                                                      double spatialFrequency) {
    double wavelengthInSystemUnits = opticalModel->nm_to_sys_units(wavelength);
    double fNumber = std::abs(opticalModel->optical_spec->parax_data->fod.fno);
    return 2.0 * wavelengthInSystemUnits * fNumber * spatialFrequency;
}

/**
 * Scale for the entrance-pupil shift so that the traced pair actually realises the
 * requested image-space spatial frequency.
 *
 * REL_PUPIL coordinates are normalised on the <em>entrance</em> pupil, but
 * the OTF is the autocorrelation of the <em>exit</em> pupil, and
 * #normalized_entry_pupil_shift derives its displacement from an exit-pupil relation
 * (2*lambda*F#*nu). The two agree only where the pupil imaging is aberration
 * free; real pupil aberration lands a rigid entrance-pupil shift as a smaller,
 * field-dependent exit-pupil shift - around 8% low at full field on a fast lens.
 *
 * A pair of rays forms image-plane fringes of frequency nu exactly when
 * their image-space direction cosines differ by lambda*nu, and that holds
 * wherever the exit pupil happens to lie. So one probe pair per field, wavelength and
 * direction measures the frequency the shift really delivers, and its reciprocal
 * corrects it.
 *
 * This removes the field-dependent bias, which is the dominant term. The remaining
 * variation across the pupil (a few percent, from pupil spherical aberration) would
 * need per-ray aiming and is left uncorrected. Cost is four rays per field and
 * wavelength against the couple of hundred used for the samples.
 */
double ContrastAnalysis::exit_pupil_frequency_calibration(
    rayoptics::optical::OpticalModel *opticalModel, rayoptics::specs::Field &field, double wavelength,
    double shift, int axis, const ContrastOptions &options) {
    if (!options.calibrateFrequency || !(shift > 0.0))
        return 1.0;
    double required = opticalModel->nm_to_sys_units(wavelength) * options.spatialFrequency;
    if (!(required > 0.0))
        return 1.0;
    auto traceOptions = options.traceOptions.copy();
    traceOptions.pupil_type = PupilType::REL_PUPIL;
    traceOptions.apply_vignetting = false;
    // Set the pupil up exactly as trace_by_wavelength does for the sampling
    // rays: take the reference image point from the central wavelength, then establish
    // this wavelength's chief ray and reference sphere against that point. Probing
    // through a differently configured pupil would measure a mapping the samples never
    // see, which is the one way this correction could do harm rather than nothing.
    double focus = opticalModel->optical_spec->defocus()->get_focus();
    auto reference = Trace::setup_pupil_coords(
        opticalModel, field, opticalModel->seq_model->central_wavelength(), focus,
        std::nullopt, std::nullopt);
    auto coordinates =
        Trace::setup_pupil_coords(opticalModel, field, wavelength, focus,
                                  reference.ref_sphere->image_pt.project_xy(),
                                  std::nullopt);
    field.chief_ray =
        std::const_pointer_cast<rayoptics::raytr::ChiefRayPkg>(coordinates.chief_ray_pkg);
    field.ref_sphere =
        std::const_pointer_cast<rayoptics::raytr::ReferenceSphere>(coordinates.ref_sphere);
    double half = 0.5 * shift;
    std::optional<double> low = image_direction(
        opticalModel, field, wavelength, traceOptions, axis,
        axis == Orientation::X ? Vector2(-half, 0.0) : Vector2(0.0, -half));
    std::optional<double> high = image_direction(
        opticalModel, field, wavelength, traceOptions, axis,
        axis == Orientation::X ? Vector2(half, 0.0) : Vector2(0.0, half));
    if (!low.has_value() || !high.has_value())
        return 1.0;
    double realized = std::abs(*high - *low);
    if (!std::isfinite(realized) || realized < 1.0e-12)
        return 1.0;
    double scale = required / realized;
    // A correction this far from unity says the probe failed rather than that the
    // pupil mapping is unusual; leave the shift alone rather than destabilise the merit.
    return scale > 0.5 && scale < 2.0 ? scale : 1.0;
}

std::optional<double> ContrastAnalysis::image_direction(
    rayoptics::optical::OpticalModel *opticalModel, rayoptics::specs::Field &field, double wavelength,
    const rayoptics::raytr::TraceOptions &traceOptions, int axis, const Vector2 &pupil) {
    auto pkg = Trace::trace_safe(opticalModel, pupil, field, wavelength, traceOptions).pkg;
    if (pkg == nullptr || pkg->ray.size() < 2)
        return std::nullopt;
    const auto &segment = rayoptics::util::Lists::get(pkg->ray, -2);
    return axis == Orientation::X ? segment.d.x : segment.d.y;
}

ContrastAnalysisResult::Sample ContrastAnalysis::sample(
    rayoptics::optical::OpticalModel *opticalModel, const ContrastRayTriplet &rays,
    rayoptics::specs::Field &field, double wavelength, double focus,
    const ContrastOptions &options) {
    (void)options;
    if (rays.reference == nullptr || rays.sagittal == nullptr ||
        rays.tangential == nullptr) {
        return ContrastAnalysisResult::Sample(rays.pupil, 0.0, 0.0, rays.weight, false,
                                              failure(rays));
    }
    // Java unboxes the Double straight away, so a null would NPE; opd only
    // answers null when the ray is missing, which the guard above rules out.
    double reference = rayoptics::analysis::WavefrontAberrationAnalysis::opd(opticalModel, rays.pupil, 0,
                                                        rays.reference, field, wavelength,
                                                        focus)
                           .value();
    double sagittal = opd(opticalModel, rays.sagittal, field, wavelength, focus,
                          rays.sagittal->input_pupil) -
                      reference;
    double tangential = opd(opticalModel, rays.tangential, field, wavelength, focus,
                            rays.tangential->input_pupil) -
                        reference;
    bool valid = std::isfinite(reference) && std::isfinite(sagittal) &&
                 std::isfinite(tangential);
    return ContrastAnalysisResult::Sample(
        rays.pupil, sagittal, tangential, rays.weight, valid,
        valid ? std::nullopt
              : std::optional<ContrastAnalysisResult::Failure>(
                    ContrastAnalysisResult::Failure("OPD",
                                                    "NonFiniteWavefrontDifference", -1)),
        reference);
}

ContrastAnalysisResult::Failure ContrastAnalysis::failure(
    const ContrastRayTriplet &rays) {
    if (rays.referenceError != nullptr)
        return failure("reference", *rays.referenceError);
    if (rays.sagittalError != nullptr)
        return failure("sagittal", *rays.sagittalError);
    if (rays.tangentialError != nullptr)
        return failure("tangential", *rays.tangentialError);
    return ContrastAnalysisResult::Failure("unknown", "MissingRay", -1);
}

ContrastAnalysisResult::Failure ContrastAnalysis::failure(
    const std::string &ray, const rayoptics::exceptions::TraceException &error) {
    return ContrastAnalysisResult::Failure(ray, error.simple_name(), error.surf);
}

double ContrastAnalysis::opd(rayoptics::optical::OpticalModel *opticalModel,
                             const std::shared_ptr<const RayPkg> &ray,
                             rayoptics::specs::Field &field, double wavelength, double focus,
                             const std::optional<Vector2> &pupil) {
    // below xy is set to 0 as it is not used by the opd calculation
    // opd ignores its pupil argument, so an absent input_pupil (Java: null)
    // changes nothing; the zero stands in for the null reference.
    return rayoptics::analysis::WavefrontAberrationAnalysis::opd(opticalModel,
                                            pupil.has_value() ? *pupil
                                                              : Vector2::vector2_0,
                                            0, ray, field, wavelength, focus)
        .value();
}

namespace {

Vector2 apply_vignetting_pt(const Vector2 &pupil, const rayoptics::specs::Field &fld) {
    auto arr = pupil.as_array();
    auto vignetted = fld.apply_vignetting(std::vector<double>{arr[0], arr[1]});
    return Vector2(vignetted[0], vignetted[1]);
}

/** Area element of Field#apply_vignetting at a nominal pupil point. */
double vignetting_jacobian(const Vector2 &pupil, const rayoptics::specs::Field &fld) {
    return fld.vignetting_scale_x(pupil.x) * fld.vignetting_scale_y(pupil.y);
}

/**
 * Smallest usable contraction of the contrast quadrature pattern. Below this
 * the samples crowd around the overlap centre, all three rays of a triplet
 * report nearly the same wavefront, and the residuals collapse towards zero -
 * which an optimizer would read as perfect contrast rather than as an
 * unusable request. The pattern reaches this limit a little below the
 * 0.707/(lambda*F#) frequency at which the overlap centre itself leaves the
 * pupil, so both bounds are reported as errors rather than silently honoured.
 */
const double MIN_CONTRACTION = 0.1;

bool valid_point(const Vector2 &pupil, const Vector2 &sagittalShift,
                          const Vector2 &tangentialShift, const rayoptics::specs::Field &fld) {
    return ContrastAnalysis::inside_vignetted_pupil(pupil, fld) &&
           ContrastAnalysis::inside_vignetted_pupil(pupil.plus(sagittalShift), fld) &&
           ContrastAnalysis::inside_vignetted_pupil(pupil.plus(tangentialShift), fld);
}

bool valid_pattern(const std::vector<Vector2> &offsets,
                            const Vector2 &overlapCenter, double contraction,
                            const Vector2 &sagittalShift, const Vector2 &tangentialShift,
                            const rayoptics::specs::Field &fld) {
    for (const auto &offset : offsets) {
        auto pupil = overlapCenter.plus(offset.times(contraction));
        if (!valid_point(pupil, sagittalShift, tangentialShift, fld)) {
            return false;
        }
    }
    return true;
}

} // namespace

bool ContrastAnalysis::inside_vignetted_pupil(const Vector2 &pupil, const rayoptics::specs::Field &fld) {
    double xScale = fld.vignetting_scale_x(pupil.x);
    double yScale = fld.vignetting_scale_y(pupil.y);
    if (!(xScale > 0.0) || !(yScale > 0.0)) {
        return false;
    }
    double x = pupil.x / xScale;
    double y = pupil.y / yScale;
    return x * x + y * y <= 1.0 + 1.0e-14;
}

std::vector<rayoptics::raytr::GaussianQuadraturePoint> ContrastAnalysis::generate_quadrature(
    const rayoptics::raytr::TraceRingsDef &grid_rng, std::optional<int> num_spokes,
    const Vector2 &sagittal_shift, const Vector2 &tangential_shift, rayoptics::specs::Field &fld) {
    auto nominal = Trace::generate_gaussian_quadrature(grid_rng, grid_rng.num_rings, num_spokes);
    Vector2 nominalCenter(grid_rng.cx, grid_rng.cy);
    auto physicalCenter = apply_vignetting_pt(nominalCenter, fld);
    auto overlapCenter =
        physicalCenter.minus(sagittal_shift.plus(tangential_shift).times(0.5));
    if (!valid_point(overlapCenter, sagittal_shift, tangential_shift, fld)) {
        throw IllegalArgumentException(
            "The requested contrast shear has no common vignetted pupil overlap");
    }
    std::vector<Vector2> offsets;
    // Map the nominal pattern into the physical pupil once. The contraction
    // search only rescales these offsets, so the map must not be repeated
    // inside the bisection loop.
    offsets.reserve(nominal.size());
    for (const auto &point : nominal)
        offsets.push_back(apply_vignetting_pt(point.pupil, fld).minus(physicalCenter));
    double contraction = 1.0;
    if (!valid_pattern(offsets, overlapCenter, contraction, sagittal_shift,
                                tangential_shift, fld)) {
        double low = 0.0;
        double high = 1.0;
        // Converge to machine precision: the contraction must be a
        // repeatable function of the design so that the finite-difference
        // Jacobian does not see quantisation steps.
        for (int iteration = 0; iteration < 60; iteration++) {
            double trial = 0.5 * (low + high);
            if (valid_pattern(offsets, overlapCenter, trial, sagittal_shift,
                                       tangential_shift, fld)) {
                low = trial;
            } else {
                high = trial;
            }
        }
        contraction = low;
    }
    if (contraction < MIN_CONTRACTION) {
        throw IllegalArgumentException(
            "The requested contrast shear leaves too little vignetted pupil overlap to "
            "sample: the pattern would contract to " +
            doubleToString(contraction) + " of its nominal extent");
    }
    std::vector<rayoptics::raytr::GaussianQuadraturePoint> points;
    points.reserve(nominal.size());
    double weightSum = 0.0;
    for (std::size_t i = 0; i < nominal.size(); i++) {
        const auto &point = nominal[i];
        auto pupil = overlapCenter.plus(offsets[i].times(contraction));
        double weight = point.weight * vignetting_jacobian(point.pupil, fld);
        points.push_back(rayoptics::raytr::GaussianQuadraturePoint(pupil, weight));
        weightSum += weight;
    }
    if (!(weightSum > 0.0)) {
        throw IllegalArgumentException("The vignetted pupil has zero area");
    }
    const double normalization = weightSum;
    std::vector<rayoptics::raytr::GaussianQuadraturePoint> normalized;
    normalized.reserve(points.size());
    for (const auto &point : points)
        normalized.push_back(
            rayoptics::raytr::GaussianQuadraturePoint(point.pupil, point.weight / normalization));
    return normalized;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Contrast trace is a Beam42 enhancement.
/// Its main use is in optimization
/// See Documentation/OPTIMIZER_NOTES.md
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<ContrastRayTriplet> ContrastAnalysis::trace_triplets(
    rayoptics::optical::OpticalModel *opt_model, const rayoptics::raytr::TraceRingsDef &grid_rng,
    std::optional<int> num_spokes, const Vector2 &sagittal_shift,
    const Vector2 &tangential_shift, rayoptics::specs::Field &fld, double wvl,
    const rayoptics::raytr::TraceOptions &trace_options) {
    return trace_triplets(opt_model, grid_rng, num_spokes, sagittal_shift,
                          tangential_shift, std::nullopt, std::nullopt, fld, wvl,
                          trace_options, false);
}

std::vector<ContrastRayTriplet> ContrastAnalysis::trace_triplets(
    rayoptics::optical::OpticalModel *opt_model, const rayoptics::raytr::TraceRingsDef &grid_rng,
    std::optional<int> num_spokes, const Vector2 &sagittal_shift,
    const Vector2 &tangential_shift, std::optional<Vector2> sagittal_exit_shift,
    std::optional<Vector2> tangential_exit_shift, rayoptics::specs::Field &fld, double wvl,
    const rayoptics::raytr::TraceOptions &trace_options_in, bool aim_exit_pupil) {
    if (aim_exit_pupil &&
        (!sagittal_exit_shift.has_value() || !tangential_exit_shift.has_value())) {
        throw IllegalArgumentException(
            "Physical exit-pupil shifts are required when aiming is enabled");
    }
    rayoptics::raytr::TraceOptions trace_options = trace_options_in.copy();
    // The quadrature generator explicitly maps samples into the physical
    // vignetted pupil. Applying Field vignetting in trace_base as well
    // would independently scale the displaced rays and change their MTF
    // shear, especially when a pair straddles a pupil axis.
    // ContrastOptions defaults aperture checking off so temporary clipping does
    // not become a discontinuous optimizer failure, but honour an explicit caller
    // choice to validate against the physical surface apertures.
    trace_options.pupil_type = PupilType::REL_PUPIL;
    trace_options.apply_vignetting = false;
    trace_options.rayerr_filter = std::string("summary");
    std::vector<ContrastRayTriplet> samples;
    auto points = generate_quadrature(grid_rng, num_spokes, sagittal_shift,
                                               tangential_shift, fld);
    for (const auto &point : points) {
        // Each sample traces a reference ray and sagittal/tangential partners.
        // The supplied shifts are normalized entrance-pupil offsets.
        // When aim_exit_pupil is enabled, these provide initial guesses;
        // each partner is inverse-aimed to a physical reference-sphere separation
        // Delta = |R| lambda nu / |n_img| from the reference ray.
        // Here lambda is in system length units and nu is image-space frequency.
        auto pupil = point.pupil;
        auto sagittalEntrancePupil = pupil.plus(sagittal_shift);
        auto tangentialEntrancePupil = pupil.plus(tangential_shift);
        auto reference = Trace::trace_safe(opt_model, pupil, fld, wvl, trace_options);
        rayoptics::raytr::RayResult sagittal;
        rayoptics::raytr::RayResult tangential;
        if (aim_exit_pupil && reference.pkg != nullptr) {
            auto referenceCoordinate = ExitPupilAiming::sphere_coord(
                reference.pkg, fld.chief_ray, fld.ref_sphere);
            if (!referenceCoordinate.has_value()) {
                auto error = std::make_shared<ExitPupilAiming::ExitPupilAimException>(
                    "Reference ray has no finite exit-pupil coordinate");
                error->surf = -1;
                sagittal = rayoptics::raytr::RayResult(nullptr, error);
                tangential = rayoptics::raytr::RayResult(nullptr, error);
            } else {
                Vector2 sagittalExitPupilTarget(
                    referenceCoordinate->x + sagittal_exit_shift->x,
                    referenceCoordinate->y + sagittal_exit_shift->y);
                Vector2 tangentialExitPupilTarget(
                    referenceCoordinate->x + tangential_exit_shift->x,
                    referenceCoordinate->y + tangential_exit_shift->y);
                sagittal = ExitPupilAiming::aim(opt_model, sagittalEntrancePupil,
                                                sagittalExitPupilTarget, fld, wvl,
                                                trace_options)
                               .ray;
                tangential = ExitPupilAiming::aim(opt_model, tangentialEntrancePupil,
                                                  tangentialExitPupilTarget, fld, wvl,
                                                  trace_options)
                                 .ray;
            }
        } else {
            sagittal =
                Trace::trace_safe(opt_model, sagittalEntrancePupil, fld, wvl, trace_options);
            tangential =
                Trace::trace_safe(opt_model, tangentialEntrancePupil, fld, wvl, trace_options);
        }
        samples.push_back(ContrastRayTriplet(pupil, reference.pkg, sagittal.pkg,
                                             tangential.pkg, reference.err,
                                             sagittal.err, tangential.err,
                                             point.weight));
    }
    return samples;
}

} // namespace redukti::analysis
