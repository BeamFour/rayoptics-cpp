// C++ port of org.redukti.rayoptics.exceptions.TraceException and subclasses
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_EXCEPTIONS_TRACEEXCEPTION_H
#define REDUKTI_RAYOPTICS_EXCEPTIONS_TRACEEXCEPTION_H

#include "redukti/Exceptions.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/math/Tfm3d.h"

#include <memory>
#include <optional>
#include <string>

namespace redukti::rayoptics::raytr {
class RayPkg;
}

namespace redukti::rayoptics::exceptions {

/**
 * These are control flow, not just errors. A ray that misses a surface, is
 * blocked by an aperture, or totally internally reflects carries a payload
 * describing how far it got, and callers routinely catch one and go on to use
 * that payload -- see Trace, VigCalc and Wideangle, which read `ray_pkg` and
 * `surf` off a caught exception and continue.
 *
 * Deriving from RuntimeException matters: the optimizer catches broadly
 * (LMDerMeritFunction wraps analysis.compute() and maps any failure onto a
 * penalty value), and ConstraintEdgeThickness catches RuntimeException around
 * a sag computation. Those catch clauses port verbatim only because the whole
 * Java hierarchy is reproduced in redukti/Exceptions.h.
 *
 * Two deliberate differences from the Java:
 *
 *  - The `ifc` field is dropped. It is assigned once, in RayTrace, and read
 *    nowhere in the entire codebase. Keeping it would put a pointer to an
 *    Interface owned by the SequentialModel inside an object that outlives the
 *    trace and gets stored in RayResult, for no benefit.
 *
 *  - `ray_pkg` is a shared_ptr to a forward-declared RayPkg. RayPkg lives in
 *    raytr and participates in a dependency cycle with specs::Field, so that
 *    cluster is ported as a unit later; the pointer lets this header be
 *    complete now. RayPkg is immutable once constructed, so sharing is safe.
 */
class TraceException : public RuntimeException {
public:
    int surf = 0;
    std::optional<mathlib::Vector3> int_pt;
    std::optional<math::Tfm3d> prev_tfrm;
    std::shared_ptr<const raytr::RayPkg> ray_pkg;
    std::optional<mathlib::Vector2> rel_p1;

    TraceException() = default;
    explicit TraceException(std::string message) : RuntimeException(std::move(message)) {}
};

class TraceMissedSurfaceException : public TraceException {
public:
    TraceMissedSurfaceException() = default;
    explicit TraceMissedSurfaceException(std::string message)
        : TraceException(std::move(message)) {}
};

class TraceRayBlockedException : public TraceException {
public:
    TraceRayBlockedException(const mathlib::Vector3 &int_pt_) { int_pt = int_pt_; }
};

class TraceTIRException : public TraceException {
public:
    /**
     * The Java constructor takes (d_in, normal, n_in, n_out) and discards all
     * four -- it only calls super(). Kept as a default constructor; the caller
     * in RayTrace fills the payload fields in directly afterwards.
     */
    TraceTIRException() = default;
};

} // namespace redukti::rayoptics::exceptions

#endif // REDUKTI_RAYOPTICS_EXCEPTIONS_TRACEEXCEPTION_H
