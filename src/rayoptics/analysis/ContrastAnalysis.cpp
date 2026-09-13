// C++ port of ContrastOptions, ContrastAnalysisResult and ContrastAnalysis.
#include "redukti/rayoptics/analysis/ContrastAnalysis.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/analysis/RayAberrations.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Lists.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>

namespace redukti::rayoptics::analysis {

namespace Orientation = util::Orientation;

using mathlib::Vector2;
using raytr::ContrastRayTriplet;
using raytr::PupilType;
using raytr::RayPkg;
using raytr::Trace;

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

ContrastOptions &ContrastOptions::trace_options(const raytr::TraceOptions &value) {
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
    return std::sqrt(weight) * sagittalDifference;
}

double ContrastAnalysisResult::Sample::tangentialResidual() const {
    return std::sqrt(weight) * tangentialDifference;
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
    return std::sqrt(sample.weight) * (sample.sagittalDifference - sagittalOffset);
}

double ContrastAnalysisResult::WavelengthResult::tangentialResidual(int index) const {
    const auto &sample = samples[static_cast<std::size_t>(index)];
    return std::sqrt(sample.weight) * (sample.tangentialDifference - tangentialOffset);
}

// ---------------------------------------------------------------------------
// ContrastAnalysis
// ---------------------------------------------------------------------------

ContrastAnalysisResult ContrastAnalysis::eval(optical::OpticalModel *opticalModel,
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
            double shift = normalized_entry_pupil_shift(opticalModel, wavelength,
                                                        options.spatialFrequency);
            double sagittalShift =
                shift * exit_pupil_frequency_calibration(opticalModel, *fields[fieldIndex],
                                                         wavelength, shift,
                                                         Orientation::X, options);
            double tangentialShift =
                shift * exit_pupil_frequency_calibration(opticalModel, *fields[fieldIndex],
                                                         wavelength, shift,
                                                         Orientation::Y, options);
            raytr::ContrastTraceCallback<ContrastAnalysisResult::Sample> callback =
                [opticalModel, &options](const ContrastRayTriplet &rays,
                                         specs::Field &field, double tracedWavelength,
                                         double focus) {
                    return sample(opticalModel, rays, field, tracedWavelength, focus,
                                  options);
                };
            auto traced =
                opticalModel->seq_model->trace_contrast<ContrastAnalysisResult::Sample>(
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
                                ? sample.sagittalDifference
                                : sample.tangentialDifference;
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
double ContrastAnalysis::normalized_entry_pupil_shift(optical::OpticalModel *opticalModel,
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
    optical::OpticalModel *opticalModel, specs::Field &field, double wavelength,
    double shift, int axis, const ContrastOptions &options) {
    if (!options.calibrateFrequency || !(shift > 0.0))
        return 1.0;
    double required = opticalModel->nm_to_sys_units(wavelength) * options.spatialFrequency;
    if (!(required > 0.0))
        return 1.0;
    auto traceOptions = options.traceOptions.copy();
    traceOptions.pupil_type = PupilType::REL_PUPIL;
    traceOptions.apply_vignetting = false;
    // Set the pupil up exactly as SequentialModel.trace_contrast does for the sampling
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
        std::const_pointer_cast<raytr::ChiefRayPkg>(coordinates.chief_ray_pkg);
    field.ref_sphere =
        std::const_pointer_cast<raytr::ReferenceSphere>(coordinates.ref_sphere);
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
    optical::OpticalModel *opticalModel, specs::Field &field, double wavelength,
    const raytr::TraceOptions &traceOptions, int axis, const Vector2 &pupil) {
    auto pkg = Trace::trace_safe(opticalModel, pupil, field, wavelength, traceOptions).pkg;
    if (pkg == nullptr || pkg->ray.size() < 2)
        return std::nullopt;
    const auto &segment = util::Lists::get(pkg->ray, -2);
    return axis == Orientation::X ? segment.d.x : segment.d.y;
}

ContrastAnalysisResult::Sample ContrastAnalysis::sample(
    optical::OpticalModel *opticalModel, const ContrastRayTriplet &rays,
    specs::Field &field, double wavelength, double focus,
    const ContrastOptions &options) {
    (void)options;
    if (rays.reference == nullptr || rays.sagittal == nullptr ||
        rays.tangential == nullptr) {
        return ContrastAnalysisResult::Sample(rays.pupil, 0.0, 0.0, rays.weight, false,
                                              failure(rays));
    }
    // Java unboxes the Double straight away, so a null would NPE; opd only
    // answers null when the ray is missing, which the guard above rules out.
    double reference = WavefrontAberrationAnalysis::opd(opticalModel, rays.pupil, 0,
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
                                                    "NonFiniteWavefrontDifference", -1)));
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
    const std::string &ray, const exceptions::TraceException &error) {
    return ContrastAnalysisResult::Failure(ray, error.simple_name(), error.surf);
}

double ContrastAnalysis::opd(optical::OpticalModel *opticalModel,
                             const std::shared_ptr<const RayPkg> &ray,
                             specs::Field &field, double wavelength, double focus,
                             const std::optional<Vector2> &pupil) {
    // below xy is set to 0 as it is not used by the opd calculation
    // opd ignores its pupil argument, so an absent input_pupil (Java: null)
    // changes nothing; the zero stands in for the null reference.
    return WavefrontAberrationAnalysis::opd(opticalModel,
                                            pupil.has_value() ? *pupil
                                                              : Vector2::vector2_0,
                                            0, ray, field, wavelength, focus)
        .value();
}

} // namespace redukti::rayoptics::analysis
