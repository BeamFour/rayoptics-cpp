// Checks redukti::formatF against Java's String.format("%.Nf", value).
//
// Java rounds HALF_UP where C printf rounds half to even, so this is not a
// wrapper test: the sample is deliberately loaded with dyadic rationals like
// 1/128 = 0.0078125, whose exact expansion ends in a 5 at the rounding
// position. Java answers 0.007813 there and printf answers 0.007812.
#include "FormatFExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"

#include <cstdlib>
#include <string>

TEST(format_f_matches_jvm) {
    const std::size_t n = sizeof(FORMATF_CASES) / sizeof(FORMATF_CASES[0]);
    for (std::size_t i = 0; i < n; i++) {
        const auto &c = FORMATF_CASES[i];
        double value = std::strtod(c.value, nullptr);
        std::string got = redukti::formatF(value, c.precision);
        if (got != c.expected)
            ::redukti::test::reportFailure(
                __FILE__, __LINE__,
                std::string("formatF(") + c.value + ", " +
                    std::to_string(c.precision) + "): got \"" + got + "\" want \"" +
                    c.expected + "\"");
    }
}
