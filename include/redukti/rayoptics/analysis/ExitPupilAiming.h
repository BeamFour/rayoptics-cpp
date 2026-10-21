// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.rayoptics.analysis.ExitPupilAiming
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_EXITPUPILAIMING_H
#define REDUKTI_RAYOPTICS_ANALYSIS_EXITPUPILAIMING_H

#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <optional>
#include <string>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::analysis {

/** Inverse entrance-to-exit-pupil mapping used by contrast tracing. */
class ExitPupilAiming {
public:
    /** A ray traced successfully, but its requested exit-pupil coordinate was not found. */
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
        raytr::RayResult ray;
        /** Null when the exit-pupil coordinate could not be computed. */
        std::optional<mathlib::Vector3> exitCoordinate;
        int iterations;
        double error;

        Result(const mathlib::Vector2 &pupil_, raytr::RayResult ray_,
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

    /**
     * Iterate the entrance pupil coordinate until the ray lands on
     * `targetExitPupilCoord`.
     */
    static Result aim(optical::OpticalModel *opticalModel,
                      const mathlib::Vector2 &initialEntrancePupilCoord,
                      const mathlib::Vector2 &targetExitPupilCoord, specs::Field &field,
                      double wavelength, const raytr::TraceOptions &traceOptions);

    /** Exit-pupil coordinate of a traced ray, or null if not computable. */
    static std::optional<mathlib::Vector3> sphere_coord(
        const std::shared_ptr<const raytr::RayPkg> &rayPkg,
        const std::shared_ptr<const raytr::ChiefRayPkg> &chiefRayPkg,
        const std::shared_ptr<const raytr::ReferenceSphere> &referenceSphere);

private:
    static const int MAX_ITERATIONS = 10;
    static const double ENTRANCE_STEP;

    class Evaluation {
    public:
        raytr::RayResult ray;
        std::optional<mathlib::Vector3> coordinate;
        std::optional<mathlib::Vector2> residual;
    };

    static Result failed(const mathlib::Vector2 &pupil,
                         std::shared_ptr<exceptions::TraceException> error);

    /**
     * Aim a ray so its transverse coordinate on the exit-pupil reference sphere equals
     * target. The initial entrance-pupil coordinate is normally the traditional
     * entrance-pupil displacement and is therefore already a close guess.
     */
    static Evaluation evaluate(optical::OpticalModel *opticalModel,
                               const mathlib::Vector2 &pupil,
                               const mathlib::Vector2 &target, specs::Field &field,
                               double wavelength, const raytr::TraceOptions &traceOptions);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_EXITPUPILAIMING_H
