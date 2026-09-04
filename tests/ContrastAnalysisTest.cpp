// C++ port of org.redukti.rayoptics.analysis.ContrastAnalysisTest.
//
// Covers the contrast merit's two pupil-shift corrections -- exit-pupil aiming
// and frequency calibration -- on the shared Nikkor Z 58mm f/0.95 S model,
// which is the model the Java uses here too.
#include "TestHarness.h"

#include "NikkorZ58Model.h"

#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/analysis/ContrastAnalysis.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

namespace {

using redukti::mathlib::Vector2;
using redukti::rayoptics::analysis::ContrastAnalysis;
using redukti::rayoptics::analysis::ContrastOptions;
using redukti::rayoptics::raytr::ContrastRayTriplet;
using redukti::rayoptics::raytr::ExitPupilAiming;
using redukti::rayoptics::raytr::Trace;
using redukti::rayoptics::raytr::TraceOptions;
using redukti::rayoptics::specs::Field;
namespace Orientation = redukti::rayoptics::util::Orientation;

/** The four sphere-coordinate separations the aiming test measures per sample. */
using Separation = std::array<double, 4>;

} // namespace

TEST(contrast_exit_pupil_aiming_realises_both_requested_shear_vectors) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    int fieldIndex = 1;
    int wavelengthIndex = 1;
    double wavelength =
        model->optical_spec->wvls->wavelengths[static_cast<std::size_t>(wavelengthIndex)];
    double normalizedShift =
        ContrastAnalysis::normalized_entry_pupil_shift(model.get(), wavelength, 40.0);

    auto traced = model->seq_model->trace_contrast<Separation>(
        [](const ContrastRayTriplet &rays, Field &field, double, double) -> Separation {
            CHECK(rays.referenceError == nullptr);
            CHECK(rays.sagittalError == nullptr);
            CHECK(rays.tangentialError == nullptr);
            auto reference =
                ExitPupilAiming::sphere_coord(rays.reference, field.chief_ray,
                                              field.ref_sphere);
            auto sagittal = ExitPupilAiming::sphere_coord(rays.sagittal, field.chief_ray,
                                                          field.ref_sphere);
            auto tangential = ExitPupilAiming::sphere_coord(
                rays.tangential, field.chief_ray, field.ref_sphere);
            CHECK(reference.has_value());
            CHECK(sagittal.has_value());
            CHECK(tangential.has_value());
            if (!reference.has_value() || !sagittal.has_value() || !tangential.has_value())
                return Separation{0.0, 0.0, 0.0, 0.0};
            return Separation{sagittal->x - reference->x, sagittal->y - reference->y,
                              tangential->x - reference->x, tangential->y - reference->y};
        },
        fieldIndex, wavelengthIndex, 1, 6, Vector2(normalizedShift, 0.0),
        Vector2(0.0, normalizedShift), 40.0, TraceOptions(), true);

    double physicalShift = ExitPupilAiming::referenceSphereShift(
        model.get(), *model->optical_spec->fov->fields[static_cast<std::size_t>(fieldIndex)],
        wavelength, 40.0);

    double tolerance = std::max(
        1.0e-9, std::abs(model->optical_spec->parax_data->fod.exp_radius) * 2.0e-7);
    CHECK(!traced.empty());
    if (traced.empty())
        return;
    for (const auto &separation : traced[0].samples) {
        CHECK_CLOSE(separation[0], physicalShift, tolerance);
        CHECK_CLOSE(separation[1], 0.0, tolerance);
        CHECK_CLOSE(separation[2], 0.0, tolerance);
        CHECK_CLOSE(separation[3], physicalShift, tolerance);
    }
}

TEST(contrast_exit_pupil_aiming_cannot_be_combined_with_block_calibration) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    ContrastOptions options(40.0);
    options.aim_exit_pupil(true).calibrate_frequency(true);

    CHECK_THROWS(ContrastAnalysis::eval(model.get(), options),
                 redukti::IllegalArgumentException);
}

TEST(contrast_converts_image_frequency_to_normalized_pupil_shift) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    double wavelength = 550.0;
    double frequency = 50.0;

    double shift = ContrastAnalysis::normalized_entry_pupil_shift(model.get(), wavelength,
                                                                  frequency);

    CHECK_CLOSE(shift, 0.05413662972175293, 1.0e-12);
}

TEST(contrast_evaluates_weighted_wavefront_differences_for_every_field_and_wavelength) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    ContrastOptions options(40.0);
    options.num_rings(2).num_spokes(6);

    auto result = ContrastAnalysis::eval(model.get(), options);

    CHECK_EQ(static_cast<int>(result.fields.size()), 3);
    for (const auto &field : result.fields) {
        CHECK_EQ(static_cast<int>(field.wavelengths.size()), 3);
        for (const auto &wavelength : field.wavelengths) {
            CHECK_EQ(static_cast<int>(wavelength.samples.size()), 12);
            double weightSum = 0.0;
            for (const auto &sample : wavelength.samples)
                weightSum += sample.weight;
            CHECK_CLOSE(weightSum, 1.0, 1.0e-14);
            for (std::size_t i = 0; i < wavelength.samples.size(); i++) {
                const auto &sample = wavelength.samples[i];
                CHECK(std::isfinite(sample.sagittalDifference));
                CHECK(std::isfinite(sample.tangentialDifference));
                CHECK(std::isfinite(sample.sagittalResidual()));
                CHECK(std::isfinite(sample.tangentialResidual()));
            }
        }
    }
}

TEST(contrast_frequency_calibration_uses_the_same_central_reference_setup) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto &field = *model->optical_spec->fov->fields[2];
    const auto &wavelengths = model->optical_spec->wvls->wavelengths;
    double focus = model->optical_spec->defocus()->get_focus();

    // Seed a different wavelength to ensure calibration does not merely inherit
    // the wavelength-specific chief ray and reference sphere left by an earlier
    // analysis.
    auto stale = Trace::setup_pupil_coords(model.get(), field, wavelengths[2], focus,
                                           std::nullopt, std::nullopt);
    field.chief_ray =
        std::const_pointer_cast<redukti::rayoptics::raytr::ChiefRayPkg>(
            stale.chief_ray_pkg);
    field.ref_sphere =
        std::const_pointer_cast<redukti::rayoptics::raytr::ReferenceSphere>(
            stale.ref_sphere);

    double wavelength = wavelengths[1];
    auto central = Trace::setup_pupil_coords(model.get(), field,
                                             model->seq_model->central_wavelength(), focus,
                                             std::nullopt, std::nullopt);
    auto expected = Trace::setup_pupil_coords(model.get(), field, wavelength, focus,
                                              central.ref_sphere->image_pt.project_xy(),
                                              std::nullopt);

    ContrastOptions options(40.0);
    options.calibrate_frequency(true);
    double shift =
        ContrastAnalysis::normalized_entry_pupil_shift(model.get(), wavelength, 40.0);
    double scale = ContrastAnalysis::exit_pupil_frequency_calibration(
        model.get(), field, wavelength, shift, Orientation::X, options);

    CHECK(std::isfinite(scale));
    CHECK_EQ(field.chief_ray->chief_ray->wvl, wavelength);
    CHECK_CLOSE(field.ref_sphere->image_pt.x, expected.ref_sphere->image_pt.x, 1.0e-12);
    CHECK_CLOSE(field.ref_sphere->image_pt.y, expected.ref_sphere->image_pt.y, 1.0e-12);
    CHECK_CLOSE(field.ref_sphere->image_pt.z, expected.ref_sphere->image_pt.z, 1.0e-12);
}

TEST(contrast_frequency_calibration_honours_aperture_checking) {
    auto model = redukti::test::buildNikkorZ58TestModel();
    auto &field = *model->optical_spec->fov->fields[1];
    double wavelength = model->optical_spec->wvls->wavelengths[1];
    double shift =
        ContrastAnalysis::normalized_entry_pupil_shift(model.get(), wavelength, 40.0);

    // Make the first real surface reject both calibration probes when physical
    // aperture checking is requested. With checking disabled the same rays trace.
    model->seq_model->ifcs[1]->max_aperture = 1.0e-6;
    ContrastOptions checkedOptions(40.0);
    checkedOptions.calibrate_frequency(true).check_apertures(true);
    double checked = ContrastAnalysis::exit_pupil_frequency_calibration(
        model.get(), field, wavelength, shift, Orientation::X, checkedOptions);
    ContrastOptions uncheckedOptions(40.0);
    uncheckedOptions.calibrate_frequency(true).check_apertures(false);
    double unchecked = ContrastAnalysis::exit_pupil_frequency_calibration(
        model.get(), field, wavelength, shift, Orientation::X, uncheckedOptions);

    // a clipped calibration probe should fall back to no correction
    CHECK_EQ(checked, 1.0);
    // disabling aperture checks should allow the calibration probes to trace
    CHECK(std::abs(unchecked - 1.0) > 1.0e-6);
}
