// C++ port of org.redukti.mathlib.M
#ifndef REDUKTI_MATHLIB_M_H
#define REDUKTI_MATHLIB_M_H

#include <cmath>

namespace redukti::mathlib {

/** Free-function equivalents of the static methods on Java's `M` class. */
namespace M {

inline constexpr double EPSILON = 2.2204460492503131e-016;

inline bool isZero(double d) { return std::abs(d) <= EPSILON; }

/**
 * Test for |x| < fuzz.
 *
 * The choice of fuzz depends on the context of the test. The default is
 * appropriate for ignoring small double precision rounding errors. Chunkier
 * values are appropriate for different calculations: ray trace values are
 * typically good to 1e-10 to 1e-12, geometric modelling or rendering might
 * be as loose as 1e-8 or 1e-6.
 */
inline bool is_fuzzy_zero(double x, double fuzz) { return std::abs(x) < fuzz; }
inline bool is_fuzzy_zero(double x) { return is_fuzzy_zero(x, 1e-14); }

inline bool is_kinda_big(double x, double kinda_big) {
    if (std::isinf(x))
        return true;
    if (std::abs(x) > kinda_big)
        return true;
    return false;
}
inline bool is_kinda_big(double x) { return is_kinda_big(x, 1e8); }

/** Replace IEEE inf with a signed big number. */
inline double infinity_guard(double x, double big) {
    if (std::isinf(x))
        return x < 0.0 ? -big : big;
    return x;
}
inline double infinity_guard(double x) { return infinity_guard(x, 1e12); }

inline double square(double x) { return x * x; }

inline int trunc(double value) {
    return static_cast<int>(value < 0 ? std::ceil(value) : std::floor(value));
}

/** Java's Math.toRadians / Math.toDegrees. */
inline double toRadians(double deg) { return deg * (3.14159265358979323846 / 180.0); }
inline double toDegrees(double rad) { return rad * (180.0 / 3.14159265358979323846); }

double cosd(double deg);
double sind(double deg);
double tand(double deg);

// NYI: decimal_format / decimal_format_scientific -- java.text.DecimalFormat.
// Only used by the exporters and renderers, which are ported later.

} // namespace M
} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_M_H
