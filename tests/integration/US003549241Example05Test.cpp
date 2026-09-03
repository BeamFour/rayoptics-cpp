// C++ port of org.redukti.rayoptics.integration.US003549241Example05Test.
//
// Both cases are about exit-pupil aiming on a wide-angle lens: the first that
// every field and wavelength of a dense contrast pattern traces without error,
// the second that aiming beats the paraxial block calibration by orders of
// magnitude. Neither asserts against captured numbers -- they assert relations
// between measured errors, so they port across unchanged.
#include "../TestHarness.h"
#include "IntegrationModels.h"

#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/VigCalc.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace {

using namespace redukti::rayoptics;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector3;
namespace Orientation = util::Orientation;

/** The four transverse separations Java returns as a double[4]. */
using Separations = std::array<double, 4>;

Separations exitPupilSeparations(const raytr::ContrastRayTriplet &rays,
                                 specs::Field &field, double wavelength, double focus) {
    (void)wavelength;
    (void)focus;
    CHECK(rays.referenceError == nullptr);
    CHECK(rays.sagittalError == nullptr);
    CHECK(rays.tangentialError == nullptr);
    auto reference = raytr::ExitPupilAiming::sphere_coord(rays.reference, field.chief_ray,
                                                          field.ref_sphere);
    auto sagittal = raytr::ExitPupilAiming::sphere_coord(rays.sagittal, field.chief_ray,
                                                         field.ref_sphere);
    auto tangential = raytr::ExitPupilAiming::sphere_coord(
        rays.tangential, field.chief_ray, field.ref_sphere);
    CHECK(reference.has_value());
    CHECK(sagittal.has_value());
    CHECK(tangential.has_value());
    return Separations{sagittal->x - reference->x, sagittal->y - reference->y,
                       tangential->x - reference->x, tangential->y - reference->y};
}

double maxSeparationError(const std::vector<Separations> &separations, double requested) {
    double maximum = 0.0;
    for (const auto &s : separations) {
        maximum = std::fmax(maximum, std::hypot(s[0] - requested, s[1]));
        maximum = std::fmax(maximum, std::hypot(s[2], s[3] - requested));
    }
    return maximum;
}

} // namespace

TEST(integration_dense_optimizer_pattern_aims_every_field_and_wavelength) {
    auto model = integration::build_US003549241Example05();
    std::vector<double> fields{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    model->optical_spec->fov = std::make_unique<specs::FieldSpec>(
        model->optical_spec.get(),
        util::Pair<specs::ImageKey, specs::ValueKey>(specs::ImageKey::Object,
                                                     specs::ValueKey::Angle),
        45.0, fields, true, true);
    model->update_model();
    raytr::VigCalc::set_pupil(model.get());
    model->update_model();

    for (std::size_t fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++) {
        const auto &wavelengths = model->optical_spec->wvls->wavelengths;
        for (std::size_t wi = 0; wi < wavelengths.size(); wi++) {
            double wavelength = wavelengths[wi];
            double shift = analysis::ContrastAnalysis::normalized_entry_pupil_shift(
                model.get(), wavelength, 20.0);
            raytr::ContrastTraceCallback<raytr::ContrastRayTriplet> keep =
                [](const raytr::ContrastRayTriplet &rays, specs::Field &, double,
                   double) { return rays; };
            raytr::TraceOptions options;
            auto traced =
                model->seq_model->trace_contrast<raytr::ContrastRayTriplet>(
                    keep, static_cast<int>(fieldIndex), static_cast<int>(wi), 6, 12,
                    Vector2(shift, 0.0), Vector2(0.0, shift), 20.0, options, true);
            for (const auto &rays : traced[0].samples) {
                CHECK(rays.referenceError == nullptr);
                CHECK(rays.sagittalError == nullptr);
                CHECK(rays.tangentialError == nullptr);
            }
        }
    }
}

TEST(integration_exit_pupil_aiming_corrects_full_field_shear) {
    auto model = integration::build_US003549241Example05();
    int fieldIndex = 2;
    int wavelengthIndex = model->optical_spec->wvls->reference_wvl;
    double wavelength =
        model->optical_spec->wvls->wavelengths[static_cast<std::size_t>(wavelengthIndex)];
    double normalizedShift = analysis::ContrastAnalysis::normalized_entry_pupil_shift(
        model.get(), wavelength, 40.0);
    Vector2 sagittalShift(normalizedShift, 0.0);
    Vector2 tangentialShift(0.0, normalizedShift);
    specs::Field &field =
        *model->optical_spec->fov->fields[static_cast<std::size_t>(fieldIndex)];

    raytr::ContrastTraceCallback<Separations> sep = exitPupilSeparations;
    raytr::TraceOptions options;
    auto unaimed = model->seq_model->trace_contrast<Separations>(
        sep, fieldIndex, wavelengthIndex, 1, 6, sagittalShift, tangentialShift, 40.0,
        options, false);
    auto aimed = model->seq_model->trace_contrast<Separations>(
        sep, fieldIndex, wavelengthIndex, 1, 6, sagittalShift, tangentialShift, 40.0,
        options, true);

    auto calibrationOptions = analysis::ContrastOptions(40.0);
    calibrationOptions.calibrate_frequency(true);
    double sagittalScale = analysis::ContrastAnalysis::exit_pupil_frequency_calibration(
        model.get(), field, wavelength, normalizedShift, Orientation::X,
        calibrationOptions);
    double tangentialScale = analysis::ContrastAnalysis::exit_pupil_frequency_calibration(
        model.get(), field, wavelength, normalizedShift, Orientation::Y,
        calibrationOptions);
    auto calibrated = model->seq_model->trace_contrast<Separations>(
        sep, fieldIndex, wavelengthIndex, 1, 6, sagittalShift.times(sagittalScale),
        tangentialShift.times(tangentialScale), 40.0, options, false);

    double requestedShift = raytr::ExitPupilAiming::referenceSphereShift(
        model.get(), field, wavelength, 40.0);
    double legacyParaxialShift =
        normalizedShift * std::abs(model->optical_spec->parax_data->fod.exp_radius);
    // wide-angle reference-sphere scale should differ from the paraxial pupil scale
    CHECK(std::abs(requestedShift - legacyParaxialShift) > 0.01);

    double unaimedError = maxSeparationError(unaimed[0].samples, requestedShift);
    double calibratedError = maxSeparationError(calibrated[0].samples, requestedShift);
    double aimedError = maxSeparationError(aimed[0].samples, requestedShift);

    // unaimed shear should expose pupil mapping error
    CHECK(unaimedError > 0.5);
    // block calibration should materially improve the shear
    CHECK(calibratedError < unaimedError * 0.2);
    // but should retain its across-pupil approximation error
    CHECK(calibratedError > 0.02);
    // aimed shear should reach the reference-sphere target
    CHECK(aimedError < 2.0e-6);
    // aiming should improve full-field shear by at least five orders of magnitude
    CHECK(aimedError * 100000.0 < unaimedError);
    // and be substantially more accurate than block calibration
    CHECK(aimedError * 10000.0 < calibratedError);
}
