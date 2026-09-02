// C++ port of org.redukti.rayoptics.raytr.ExitPupilAiming
#include "redukti/rayoptics/raytr/ExitPupilAiming.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/raytr/WaveAbr.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace redukti::rayoptics::raytr {

namespace M = mathlib::M;
using mathlib::Vector2;
using mathlib::Vector3;

const double ExitPupilAiming::ENTRANCE_STEP = 1.0e-4;

double ExitPupilAiming::referenceSphereShift(optical::OpticalModel *opticalModel,
                                             specs::Field &field, double wavelength,
                                             double spatialFrequency) {
    if (field.ref_sphere == nullptr ||
        M::is_kinda_big(field.ref_sphere->ref_sphere_radius)) {
        throw IllegalArgumentException("Finite reference sphere required");
    }
    double wavelengthInSystemUnits = opticalModel->nm_to_sys_units(wavelength);
    double imageIndex = std::abs(opticalModel->optical_spec->parax_data->fod.n_img);
    if (!(imageIndex > 0.0)) {
        throw IllegalArgumentException("Positive image-space index required");
    }
    return wavelengthInSystemUnits * std::abs(spatialFrequency) *
           std::abs(field.ref_sphere->ref_sphere_radius) / imageIndex;
}

ExitPupilAiming::Result ExitPupilAiming::failed(
    const Vector2 &pupil, std::shared_ptr<exceptions::TraceException> error) {
    error->surf = -1;
    return Result(pupil, RayResult(nullptr, error), std::nullopt, MAX_ITERATIONS,
                  std::numeric_limits<double>::quiet_NaN());
}

ExitPupilAiming::Evaluation ExitPupilAiming::evaluate(
    optical::OpticalModel *opticalModel, const Vector2 &pupil, const Vector2 &target,
    specs::Field &field, double wavelength, const TraceOptions &traceOptions) {
    Evaluation e;
    e.ray = Trace::trace_safe(opticalModel, pupil, field, wavelength, traceOptions);
    if (e.ray.pkg == nullptr)
        return e;
    auto coordinate = sphere_coord(e.ray.pkg, field.chief_ray, field.ref_sphere);
    if (!coordinate.has_value())
        return e;
    e.coordinate = coordinate;
    e.residual = Vector2(coordinate->x - target.x, coordinate->y - target.y);
    return e;
}

ExitPupilAiming::Result ExitPupilAiming::aim(optical::OpticalModel *opticalModel,
                                             const Vector2 &initialPupil,
                                             const Vector2 &target, specs::Field &field,
                                             double wavelength,
                                             const TraceOptions &traceOptions) {
    if (field.chief_ray == nullptr || field.ref_sphere == nullptr ||
        M::is_kinda_big(field.ref_sphere->ref_sphere_radius)) {
        return failed(initialPupil,
                      std::make_shared<ExitPupilAimException>(
                          "Finite exit-pupil reference sphere required"));
    }
    double pupilRadius = std::abs(opticalModel->optical_spec->parax_data->fod.exp_radius);
    // Sub-micrometre accuracy on a typical photographic exit pupil is already
    // far below the frequency resolution useful to the contrast merit function.
    // A tighter threshold stalls on trace/intersection roundoff for off-axis rays.
    double targetTolerance = std::max(1.0e-10, pupilRadius * 2.0e-7);
    double acceptableTolerance = std::max(targetTolerance, pupilRadius * 1.0e-6);
    Vector2 pupil = initialPupil;
    Evaluation current =
        evaluate(opticalModel, pupil, target, field, wavelength, traceOptions);
    if (current.ray.err != nullptr) {
        return Result(pupil, current.ray, current.coordinate, 0,
                      std::numeric_limits<double>::quiet_NaN());
    }
    if (!current.coordinate.has_value())
        return failed(pupil, std::make_shared<ExitPupilAimException>(
                                 "Exit-pupil coordinate could not be computed"));
    for (int iteration = 0; iteration <= MAX_ITERATIONS; iteration++) {
        double error = current.residual->len();
        if (error <= targetTolerance) {
            return Result(pupil, current.ray, current.coordinate, iteration, error);
        }
        if (iteration == MAX_ITERATIONS)
            break;
        Evaluation xProbe = evaluate(opticalModel,
                                     pupil.plus(Vector2(ENTRANCE_STEP, 0.0)), target,
                                     field, wavelength, traceOptions);
        Evaluation yProbe = evaluate(opticalModel,
                                     pupil.plus(Vector2(0.0, ENTRANCE_STEP)), target,
                                     field, wavelength, traceOptions);
        if (!xProbe.coordinate.has_value() || !yProbe.coordinate.has_value())
            break;
        double j00 = (xProbe.coordinate->x - current.coordinate->x) / ENTRANCE_STEP;
        double j10 = (xProbe.coordinate->y - current.coordinate->y) / ENTRANCE_STEP;
        double j01 = (yProbe.coordinate->x - current.coordinate->x) / ENTRANCE_STEP;
        double j11 = (yProbe.coordinate->y - current.coordinate->y) / ENTRANCE_STEP;
        double determinant = j00 * j11 - j01 * j10;
        if (!std::isfinite(determinant) || std::abs(determinant) < 1.0e-14)
            break;
        // J.dp = -residual
        double dx =
            (-j11 * current.residual->x + j01 * current.residual->y) / determinant;
        double dy =
            (j10 * current.residual->x - j00 * current.residual->y) / determinant;
        if (!std::isfinite(dx) || !std::isfinite(dy))
            break;
        // Backtracking keeps the inverse map stable near a strongly aberrated edge.
        bool haveAccepted = false;
        Evaluation accepted;
        Vector2 acceptedPupil(0.0, 0.0);
        for (int lineSearch = 0; lineSearch < 8; lineSearch++) {
            double scale = std::scalbn(1.0, -lineSearch);
            Vector2 trialPupil = pupil.plus(Vector2(dx * scale, dy * scale));
            Evaluation trial = evaluate(opticalModel, trialPupil, target, field,
                                        wavelength, traceOptions);
            if (trial.coordinate.has_value() && trial.residual->len() < error) {
                accepted = trial;
                acceptedPupil = trialPupil;
                haveAccepted = true;
                break;
            }
        }
        if (!haveAccepted)
            break;
        pupil = acceptedPupil;
        current = accepted;
    }
    double finalError = current.residual->len();
    if (finalError <= acceptableTolerance)
        return Result(pupil, current.ray, current.coordinate, MAX_ITERATIONS,
                      finalError);
    return failed(pupil, std::make_shared<ExitPupilAimException>(
                             "Exit-pupil aiming did not converge; residual=" +
                             doubleToString(current.residual->len()) + ", pupil=" +
                             pupil.toString()));
}

std::optional<Vector3> ExitPupilAiming::sphere_coord(
    const std::shared_ptr<const RayPkg> &rayPkg,
    const std::shared_ptr<const ChiefRayPkg> &chiefRayPkg,
    const std::shared_ptr<const ReferenceSphere> &referenceSphere) {
    if (rayPkg == nullptr)
        return std::nullopt;
    if (chiefRayPkg == nullptr || referenceSphere == nullptr)
        return std::nullopt;
    if (M::is_kinda_big(referenceSphere->ref_sphere_radius))
        return std::nullopt;
    auto &chiefRay = chiefRayPkg->chief_ray->ray;
    if (chiefRay.size() < 2)
        return std::nullopt;
    auto &testRay = rayPkg->ray;
    if (testRay.size() != chiefRay.size())
        return std::nullopt;
    return WaveAbr::wave_abr_calc_finite_pupil(rayPkg, chiefRayPkg, referenceSphere)
        .ray_exit_pupil_coord;
}

} // namespace redukti::rayoptics::raytr
