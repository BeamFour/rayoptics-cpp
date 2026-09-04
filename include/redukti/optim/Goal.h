// C++ port of org.redukti.optim.Goal
#ifndef REDUKTI_OPTIM_GOAL_H
#define REDUKTI_OPTIM_GOAL_H

#include "redukti/optim/Analysis.h"

#include <string>

namespace redukti::optim {

/**
 * Goal represents a target we would like to achieve.
 */
class Goal {
public:
    /**
     * Borrowed. The Java holds a reference to the Analysis the optimization
     * setup owns; that Analysis outlives every goal built against it.
     */
    Analysis *const _analysis;
    const double _weight;
    const double _target;

    Goal(Analysis *analysis, double target, double weight)
        : _analysis(analysis), _weight(weight), _target(target) {}
    virtual ~Goal() = default;

    virtual double value() = 0;

    /** As for Var::toString: the Java inherits Object's when not overridden. */
    virtual std::string toString();
};

} // namespace redukti::optim

#endif // REDUKTI_OPTIM_GOAL_H
