// Minimal logging for the ported ray-optics and optimizer messages.
//
// The Java uses java.util.logging, one logger per class, and LensTool2 sets the
// thresholds from --verbose and --debug. This keeps the two properties that
// matter and nothing else: a threshold per area, and a message that is only
// built when its level is enabled, so the vignetting and pupil searches pay a
// single comparison when logging is off. There is no third-party dependency.
#ifndef REDUKTI_UTIL_LOG_H
#define REDUKTI_UTIL_LOG_H

#include "redukti/Text.h"

#include <atomic>
#include <cstdio>
#include <optional>
#include <string>

namespace redukti::util::log {

/**
 * The three levels the port uses, numbered as java.util.logging numbers WARNING,
 * CONFIG and FINE. Upstream's Python warning, info and debug map onto them.
 */
enum class Level : int { Warning = 900, Info = 700, Debug = 500 };

/**
 * Loggers grouped as the Java packages group them: org.redukti.optim and
 * org.redukti.rayoptics. --verbose opens Optim alone, --debug opens both.
 */
enum class Area : int { Optim = 0, Rayoptics = 1 };

inline constexpr int AREA_COUNT = 2;

namespace detail {
// Warnings only by default, as with Python's logging defaults upstream.
inline std::atomic<int> thresholds[AREA_COUNT] = {static_cast<int>(Level::Warning),
                                                 static_cast<int>(Level::Warning)};
} // namespace detail

inline void set_level(Area area, Level level) {
    detail::thresholds[static_cast<int>(area)].store(static_cast<int>(level),
                                                     std::memory_order_relaxed);
}

inline void set_all_levels(Level level) {
    for (auto &threshold : detail::thresholds)
        threshold.store(static_cast<int>(level), std::memory_order_relaxed);
}

inline bool enabled(Area area, Level level) {
    return static_cast<int>(level) >=
           detail::thresholds[static_cast<int>(area)].load(std::memory_order_relaxed);
}

/** One bare line on stderr, which is how LensTool2's Java console handler prints. */
inline void write(const std::string &message) {
    std::fputs(message.c_str(), stderr);
    std::fputc('\n', stderr);
}

// Pieces for building a message the way the Java builds it with
// String.format(Locale.ROOT, ...). Each is named after its conversion.

/** Right-aligns in `width` columns, as a Java format width does. */
inline std::string pad(std::string s, int width) {
    if (static_cast<int>(s.size()) < width)
        s.insert(s.begin(), static_cast<std::size_t>(width) - s.size(), ' ');
    return s;
}

/** `%<width>.<precision>f` */
inline std::string f(double value, int width, int precision) {
    return pad(formatF(value, precision), width);
}

/** `%<width>.<precision>e` */
inline std::string e(double value, int width, int precision) {
    return formatE(value, width, precision);
}

/** `%<width>.<precision>g` */
inline std::string g(double value, int width, int precision) {
    return formatG(value, width, precision);
}

/** `%d` of a Java Integer, which prints "null" when there is none. */
inline std::string d(std::optional<int> value) {
    return value.has_value() ? std::to_string(*value) : std::string("null");
}

/** `%s` of a Java boolean. */
inline std::string b(bool value) { return value ? "true" : "false"; }

} // namespace redukti::util::log

/**
 * Logs `message` at `level` in `area` (the enumerator names, unqualified). The
 * message expression is evaluated only when the level is enabled.
 */
#define REDUKTI_LOG(area, level, message)                                                  \
    do {                                                                                   \
        if (::redukti::util::log::enabled(::redukti::util::log::Area::area,                \
                                          ::redukti::util::log::Level::level))             \
            ::redukti::util::log::write(message);                                          \
    } while (false)

#define REDUKTI_LOG_WARNING(area, message) REDUKTI_LOG(area, Warning, message)
#define REDUKTI_LOG_INFO(area, message) REDUKTI_LOG(area, Info, message)
#define REDUKTI_LOG_DEBUG(area, message) REDUKTI_LOG(area, Debug, message)

#endif // REDUKTI_UTIL_LOG_H
