// Checks redukti::DecimalFormat against java.text.DecimalFormat as configured
// by M.decimal_format, over 3207 cases dumped from JDK 25.
//
// The sample deliberately includes exact halves and near-halves at several
// magnitudes: that is where HALF_EVEN, and Java rounding the shortest
// round-tripping decimal rather than the exact binary value, both show up.
#include "DecimalFormatExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"
#include "redukti/mathlib/M.h"

#include <cstdlib>
#include <string>

TEST(decimal_format_matches_jvm) {
    const std::size_t n = sizeof(DECIMAL_FORMAT_CASES) / sizeof(DECIMAL_FORMAT_CASES[0]);
    for (std::size_t i = 0; i < n; i++) {
        const auto &c = DECIMAL_FORMAT_CASES[i];
        double value = std::strtod(c.value, nullptr);
        auto fmt = redukti::mathlib::M::decimal_format(c.max_fraction_digits);
        std::string got = fmt.format(value);
        if (got != c.expected)
            ::redukti::test::reportFailure(
                __FILE__, __LINE__,
                std::string("decimal_format(") + std::to_string(c.max_fraction_digits) +
                    ").format(" + c.value + "): got \"" + got + "\" want \"" +
                    c.expected + "\"");
    }
}
