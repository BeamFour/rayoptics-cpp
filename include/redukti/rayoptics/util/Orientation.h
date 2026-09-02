// C++ port of org.redukti.rayoptics.util.Orientation
#ifndef REDUKTI_RAYOPTICS_UTIL_ORIENTATION_H
#define REDUKTI_RAYOPTICS_UTIL_ORIENTATION_H

#include "redukti/Exceptions.h"

#include <string>

namespace redukti::rayoptics::util {

/** Java's final class with a private constructor becomes a namespace. */
namespace Orientation {

/** Sagittal: the x meridian, perpendicular to the field direction. */
inline constexpr int SAGITTAL = 0;
/** Tangential: the y meridian, along the field direction. */
inline constexpr int TANGENTIAL = 1;
/** The sagittal meridian, named as a ray coordinate. */
inline constexpr int X = SAGITTAL;
/** The tangential meridian, named as a ray coordinate. */
inline constexpr int Y = TANGENTIAL;
/** Number of meridians, i.e. the exclusive bound for a loop over both. */
inline constexpr int COUNT = 2;

/** Returns the orientation unchanged, rejecting anything that is not one of the two meridians. */
inline int checked(int orientation) {
    if (orientation != SAGITTAL && orientation != TANGENTIAL)
        throw IllegalArgumentException(
            "orientation must be SAGITTAL (" + std::to_string(SAGITTAL) + ") or TANGENTIAL (" +
            std::to_string(TANGENTIAL) + ") but was " + std::to_string(orientation));
    return orientation;
}

/** "sag" or "tan", for labels and descriptions. */
inline std::string name(int orientation) {
    return orientation == TANGENTIAL ? "tan" : "sag";
}

} // namespace Orientation
} // namespace redukti::rayoptics::util

#endif // REDUKTI_RAYOPTICS_UTIL_ORIENTATION_H
