// Minimal zero-dependency test harness.
//
//   TEST(name) { CHECK(cond); CHECK_CLOSE(a, b, tol); }
//
// Each TEST self-registers, so adding a test file to CMakeLists is enough.
#ifndef REDUKTI_TESTHARNESS_H
#define REDUKTI_TESTHARNESS_H

#include "redukti/Text.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace redukti::test {

struct TestCase {
    const char *name;
    void (*fn)();
};

std::vector<TestCase> &registry();
/** Runs every registered test, or only those whose name contains `filter`. */
int runAll(const char *filter = nullptr);

/** Failures within the currently running test. */
extern int currentFailures;

inline void reportFailure(const char *file, int line, const std::string &what) {
    std::printf("    FAIL %s:%d: %s\n", file, line, what.c_str());
    ++currentFailures;
}

struct Registrar {
    Registrar(const char *name, void (*fn)()) { registry().push_back({name, fn}); }
};

} // namespace redukti::test

#define TEST(name)                                                                       \
    static void name();                                                                  \
    static ::redukti::test::Registrar name##_registrar(#name, &name);                    \
    static void name()

#define CHECK(cond)                                                                      \
    do {                                                                                 \
        if (!(cond))                                                                     \
            ::redukti::test::reportFailure(__FILE__, __LINE__, "CHECK(" #cond ")");      \
    } while (0)

#define CHECK_EQ(actual, expected)                                                       \
    do {                                                                                 \
        auto _a = (actual);                                                              \
        auto _e = (expected);                                                            \
        if (!(_a == _e))                                                                 \
            ::redukti::test::reportFailure(__FILE__, __LINE__,                           \
                                           std::string(#actual) + " != " + #expected);   \
    } while (0)

#define CHECK_STR_EQ(actual, expected)                                                   \
    do {                                                                                 \
        std::string _a = (actual);                                                       \
        std::string _e = (expected);                                                     \
        if (_a != _e)                                                                    \
            ::redukti::test::reportFailure(__FILE__, __LINE__,                           \
                                           "got \"" + _a + "\" want \"" + _e + "\"");    \
    } while (0)

#define CHECK_CLOSE(actual, expected, tol)                                               \
    do {                                                                                 \
        double _a = (actual);                                                            \
        double _e = (expected);                                                          \
        if (!(std::abs(_a - _e) <= (tol)))                                               \
            ::redukti::test::reportFailure(                                              \
                __FILE__, __LINE__,                                                      \
                "got " + ::redukti::doubleToString(_a) + " want " +                      \
                    ::redukti::doubleToString(_e) + " (delta " +                         \
                    ::redukti::doubleToString(_a - _e) + ")");                           \
    } while (0)

#define CHECK_THROWS(expr, exceptionType)                                                \
    do {                                                                                 \
        bool _caught = false;                                                            \
        try {                                                                            \
            (void)(expr);                                                                \
        } catch (const exceptionType &) {                                                \
            _caught = true;                                                              \
        } catch (...) {                                                                  \
        }                                                                                \
        if (!_caught)                                                                    \
            ::redukti::test::reportFailure(__FILE__, __LINE__,                           \
                                           #expr " did not throw " #exceptionType);      \
    } while (0)

#endif // REDUKTI_TESTHARNESS_H
