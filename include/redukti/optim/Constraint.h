// C++ port of org.redukti.optim.Constraint and its subclasses
// (ConstraintCurvature, ConstraintThickness, ConstraintEdgeThickness).
#ifndef REDUKTI_OPTIM_CONSTRAINT_H
#define REDUKTI_OPTIM_CONSTRAINT_H

#include "redukti/optim/Goal.h"

namespace redukti::optim {

/**
 * Anchors a prescription parameter to the value it started at.
 *
 * Contrast optimization sees the design clearly enough to rearrange it
 * wholesale, and left to itself it will: on the Leica 75/2 with every air space
 * free it drove three gaps negative, passing elements through each other and
 * through the stop. Nothing in an optical merit function has an opinion about
 * mechanical layout, so something has to.
 *
 * These goals supply that opinion the same way GoalParax holds focal length and
 * f-number: a target, a deviation, and a weight. There is no dead band and no
 * bound. The parameter is free to move, it simply costs merit to do so, and the
 * weight decides how much. Raise it to hold the design close, lower it to let
 * the optimizer explore.
 *
 * The residual is a FRACTION of the starting value rather than an absolute
 * deviation, which is what makes one weight sensible across a whole
 * prescription: a lens has 0.1mm air gaps beside 39mm ones, and surfaces at
 * r=14 beside r=2009. An absolute residual would effectively freeze the small
 * ones and ignore the large.
 *
 * Deliberately never reports LMLSolver::BIGVAL. A goal like this exists to steer
 * the solver, not to end the run, and a single BIGVAL raised during a Jacobian
 * probe step aborts the whole solve.
 */
class Constraint : public Goal {
public:
    const int _surface_id;

    /** The parameter's present value, in the same units as _target. */
    virtual double current_value() const = 0;

    double value() override;

    /** Signed change from the starting value, as a fraction of it. */
    double fractional_deviation() { return value() / _target - 1.0; }

    std::string toString() override;

protected:
    /**
     * @param base   the parameter's starting value, which becomes _target
     * @param weight the caller's weight, before fractional normalization
     */
    Constraint(Analysis *analysis, int surfaceId, double base, double weight);

    /** Class name, for toString; the Java uses getClass().getSimpleName(). */
    virtual const char *kind() const = 0;

private:
    /** How far out of range an undefined parameter reads: repellent, but finite. */
    static constexpr double UNDEFINED_DEVIATION = 1.0e6;

    /**
     * Fold the fractional normalization into the weight so the starting value
     * can serve as the target.
     *
     * The solver forms (value - target) * sqrt(weight). Holding a parameter to a
     * FRACTION of where it started means
     *
     *   (v/base - 1) * sqrt(w) = (v - base) * sqrt(w / base^2)
     *
     * so scaling the weight by 1/base^2 is exactly equivalent, and lets value()
     * report the parameter itself against a target of its starting value - the
     * same shape as GoalParax - instead of carrying a separate base.
     *
     * The consequence is that _weight is not the number the caller passed. It is
     * larger for small parameters and smaller for large ones, which is the
     * normalization doing its job: a 0.1mm air gap and a 39mm back focus then
     * resist a given PROPORTIONAL change equally.
     */
    static double normalized_weight(double base, double weight);
};

/**
 * Holds a surface near its starting CURVATURE.
 *
 * Curvature rather than radius, deliberately. The optimizer varies radius
 * through VarRadius, but radius is a poor measure of how much a surface has
 * really changed: on a near-flat surface it runs away towards infinity for a
 * negligible optical change, so constraining it fractionally would barely
 * restrain that surface while over-restraining a strongly curved one. Curvature
 * c = 1/r tracks optical effect.
 *
 * Target and value are therefore both curvatures, and the fractional change the
 * constraint resists is
 *
 *   c/c0 - 1 = (1/r)/(1/r0) - 1 = r0/r - 1
 *
 * It behaves sensibly at both extremes. As the surface flattens, r -> infinity
 * and the change tends to -1, so fully flattening a surface reads as 100%. As it
 * curves up, r -> 0 and the change grows without bound, which is exactly where
 * strong resistance is wanted. A surface that starts flat has no curvature
 * variable in the first place, so a zero starting radius cannot arise here.
 */
class ConstraintCurvature : public Constraint {
public:
    ConstraintCurvature(Analysis *analysis, int surfaceId, double weight)
        : Constraint(analysis, surfaceId, curvature(analysis, surfaceId), weight) {}

    double current_value() const override {
        // A radius driven to zero gives a non-finite curvature, which the base
        // class turns into a large finite miss rather than a solve-killing BIGVAL.
        return curvature(_analysis, _surface_id);
    }

protected:
    const char *kind() const override { return "ConstraintCurvature"; }

private:
    /** Curvature, so that target, value and residual are all in the same quantity. */
    static double curvature(Analysis *analysis, int surfaceId);
};

/**
 * Holds an air space or element thickness near its starting value.
 *
 * This is what stops the solver collapsing a gap or driving elements through one
 * another. Note it constrains the axial (centre) thickness only, so it does not
 * by itself guarantee positive EDGE separation, which also depends on the sag of
 * the two bounding surfaces. Keeping the layout recognisable is what makes it
 * effective in practice rather than any guarantee.
 */
class ConstraintThickness : public Constraint {
public:
    ConstraintThickness(Analysis *analysis, int surfaceId, double weight)
        : Constraint(analysis, surfaceId, thickness(analysis, surfaceId), weight) {}

    double current_value() const override { return thickness(_analysis, _surface_id); }

protected:
    const char *kind() const override { return "ConstraintThickness"; }

private:
    static double thickness(Analysis *analysis, int surfaceId);
};

/**
 * Holds the EDGE separation of a gap near its starting value.
 *
 * ConstraintThickness holds axial centre thickness, which is not the same thing:
 * two surfaces can keep their axial gap and still pass through one another away
 * from the axis, because the separation at height h is
 *
 *   gap(h) = t + sag_next(h) - sag_this(h)
 *
 * and curvature is free to move. That is how a solve with thickness constraints
 * in place still produced overlapping first and second surfaces on the Leica
 * 75/2. This constraint watches the quantity that actually goes negative.
 *
 * Measured by default at the smaller of the two bounding semi-diameters, which
 * is the outermost height at which both surfaces exist. Pass an explicit height
 * to check somewhere else - a mount land outside the clear aperture, say.
 *
 * Like every Constraint this is a penalty, not a bound. It makes crossing
 * expensive rather than impossible, and it anchors to the starting separation
 * rather than to zero, so it resists CHANGE in either direction. A design that
 * needs its edges opened up rather than preserved wants a different goal.
 */
class ConstraintEdgeThickness : public Constraint {
public:
    /** Height at which the separation is measured, in system units. */
    const double _height;

    /**
     * Constrain the gap following surfaceId, measured at the smaller of the two
     * bounding semi-diameters.
     */
    ConstraintEdgeThickness(Analysis *analysis, int surfaceId, double weight)
        : ConstraintEdgeThickness(analysis, surfaceId, default_height(analysis, surfaceId),
                                  weight) {}

    ConstraintEdgeThickness(Analysis *analysis, int surfaceId, double height,
                            double weight)
        : Constraint(analysis, surfaceId,
                     edge_gap(analysis, surfaceId, checked_height(height)), weight),
          _height(height) {}

    double current_value() const override {
        return edge_gap(_analysis, _surface_id, _height);
    }

    /**
     * True when the gap after surfaceId can be constrained: there is a following
     * surface, and the starting separation is positive and finite.
     *
     * A design that already starts with coincident or crossed surfaces has
     * nothing useful to anchor to, and a fractional constraint cannot be formed
     * around zero.
     */
    static bool is_constrainable(Analysis *analysis, int surfaceId);

    std::string toString() override;

protected:
    const char *kind() const override { return "ConstraintEdgeThickness"; }

private:
    /** The smaller of the two bounding semi-diameters. */
    static double default_height(Analysis *analysis, int surfaceId);
    static double checked_height(double height);

    /**
     * Separation between this surface and the next at `height`. Returns NaN when
     * either surface is undefined there - a radius driven below the
     * semi-diameter, say - which Constraint::value() turns into a large finite
     * miss rather than a solve-ending BIGVAL.
     */
    static double edge_gap(Analysis *analysis, int surfaceId, double height);
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_CONSTRAINT_H
