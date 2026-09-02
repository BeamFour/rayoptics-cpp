// C++ port of org.redukti.mathlib.Line3
#ifndef REDUKTI_MATHLIB_LINE3_H
#define REDUKTI_MATHLIB_LINE3_H

#include "redukti/mathlib/Vector3.h"

namespace redukti::mathlib {

class Line3 {
public:
    /** intersection point with interface */
    Vector3 origin;
    /** direction cosine exiting the interface */
    Vector3 direction;

    Line3(const Vector3 &origin_, const Vector3 &direction_)
        : origin(origin_), direction(direction_) {}
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_LINE3_H
