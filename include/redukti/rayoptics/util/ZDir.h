// C++ port of org.redukti.rayoptics.util.ZDir
#ifndef REDUKTI_RAYOPTICS_UTIL_ZDIR_H
#define REDUKTI_RAYOPTICS_UTIL_ZDIR_H

namespace redukti::rayoptics::util {

/**
 * Java enums carry fields and methods; a C++ enum class cannot, so the value
 * and the two helpers become free functions in the same namespace. Call sites
 * that read `z_dir.value` in Java read `value(z_dir)` here.
 */
enum class ZDir {
    PROPAGATE_RIGHT,
    PROPAGATE_LEFT,
};

inline double value(ZDir d) { return d == ZDir::PROPAGATE_LEFT ? -1.0 : 1.0; }

inline ZDir opposite(ZDir d) {
    if (d == ZDir::PROPAGATE_LEFT)
        return ZDir::PROPAGATE_RIGHT;
    else
        return ZDir::PROPAGATE_LEFT;
}

inline ZDir zdir_from(double v) {
    if (v >= 0.0)
        return ZDir::PROPAGATE_RIGHT;
    return ZDir::PROPAGATE_LEFT;
}

} // namespace redukti::rayoptics::util

#endif // REDUKTI_RAYOPTICS_UTIL_ZDIR_H
