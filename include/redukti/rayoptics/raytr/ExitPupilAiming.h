// C++ port of org.redukti.rayoptics.raytr.ExitPupilAiming
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_EXITPUPILAIMING_H
#define REDUKTI_RAYOPTICS_RAYTR_EXITPUPILAIMING_H

#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <optional>
#include <string>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::raytr {

class ExitPupilAiming {
public:
    /** Java's nested `ExitPupilAimException extends TraceException`. */
    class ExitPupilAimException : public exceptions::TraceException {
    public:
        explicit ExitPupilAimException(std::string message)
            : exceptions::TraceException(std::move(message)) {}

        std::string simple_name() const override { return "ExitPupilAimException"; }

        std::shared_ptr<exceptions::TraceException> clone() const override {
            return std::make_shared<ExitPupilAimException>(*this);
        }
    };

    /** Result of aiming one ray at a transverse coordinate on the reference sphere. */
    class Result {
    public:
        mathlib::Vector2 pupil;
        RayResult ray;
        /** Null when the exit-pupil coordinate could not be computed. */
        std::optional<mathlib::Vector3> exitCoordinate;
        int iterations;
        double error;

        Result(const mathlib::Vector2 &pupil_, RayResult ray_,
               std::optional<mathlib::Vector3> exitCoordinate_, int iterations_,
               double error_)
            : pupil(pupil_), ray(std::move(ray_)),
              exitCoordinate(exitCoordinate_), iterations(iterations_),
              error(error_) {}
    };

    /** Defocus shift of the reference sphere for a given spatial frequency. */
    static double referenceSphereShift(optical::OpticalModel *opticalModel,
                                       specs::Field &field, double wavelength,
                                       double spatialFrequency);

    /** Iterate the entrance pupil coordinate until the ray lands on `target`. */
    static Result aim(optical::OpticalModel *opticalModel,
                      const mathlib::Vector2 &initialPupil,
                      const mathlib::Vector2 &target, specs::Field &field,
                      double wavelength, const TraceOptions &traceOptions);

    /** Exit-pupil coordinate of a traced ray, or null if not computable. */
    static std::optional<mathlib::Vector3> sphere_coord(
        const std::shared_ptr<const RayPkg> &rayPkg,
        const std::shared_ptr<const ChiefRayPkg> &chiefRayPkg,
        const std::shared_ptr<const ReferenceSphere> &referenceSphere);

private:
    static const int MAX_ITERATIONS = 10;
    static const double ENTRANCE_STEP;

    class Evaluation {
    public:
        RayResult ray;
        std::optional<mathlib::Vector3> coordinate;
        std::optional<mathlib::Vector2> residual;
    };

    static Result failed(const mathlib::Vector2 &pupil,
                         std::shared_ptr<exceptions::TraceException> error);

    static Evaluation evaluate(optical::OpticalModel *opticalModel,
                               const mathlib::Vector2 &pupil,
                               const mathlib::Vector2 &target, specs::Field &field,
                               double wavelength, const TraceOptions &traceOptions);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_EXITPUPILAIMING_H
