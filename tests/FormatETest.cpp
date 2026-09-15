// Checks redukti::formatE against Java's String.format("%<width>.<precision>e", value).
//
// Java rounds the shortest round-tripping decimal HALF_UP where C printf rounds
// the exact binary value half to even, so the sample is loaded with dyadic
// rationals whose exact expansion ends in a 5, and with values like 9.995 that
// carry into the next power of ten.
#include "FormatEExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"

#include <cstdlib>
#include <string>

TEST(format_e_matches_jvm) {
    const std::size_t n = sizeof(FORMATE_CASES) / sizeof(FORMATE_CASES[0]);
    for (std::size_t i = 0; i < n; i++) {
        const auto &c = FORMATE_CASES[i];
        double value = std::strtod(c.value, nullptr);
        std::string got = redukti::formatE(value, c.width, c.precision);
        if (got != c.expected)
            ::redukti::test::reportFailure(
                __FILE__, __LINE__,
                std::string("formatE(") + c.value + ", " + std::to_string(c.width) + ", " +
                    std::to_string(c.precision) + "): got \"" + got + "\" want \"" +
                    c.expected + "\"");
    }
}
