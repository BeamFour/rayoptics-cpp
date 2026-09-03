// End-to-end check of the org.redukti.data port against the JVM.
//
// Builds the same irregularly spaced DiscreteSet as scratchpad/DumpData.java,
// then renders count, ranges, intervals, the two search functions and every
// interpolation method at four derivative orders, and compares the whole thing
// with output dumped from JDK 25.
//
// This layer is pure arithmetic on stored samples -- no trig anywhere -- so the
// comparison is exact, character for character.
#include "DataExpected.h"
#include "TestHarness.h"

#include "redukti/Text.h"
#include "redukti/data/DataSet.h"

#include <string>
#include <vector>

using redukti::doubleToString;
using redukti::data::DiscreteSet;
using redukti::data::Interpolation;

namespace {

/** The Java writes String.valueOf(double); this is the same text. */
std::string d(double v) { return doubleToString(v); }

DiscreteSet build() {
    DiscreteSet s;
    const double pts[][3] = {
        {0.0, 1.0, 0.5},    {0.7, 2.3, -0.25}, {1.9, 0.4, 1.75}, {3.1, -1.2, 0.0},
        {4.05, 2.9, -1.5},  {6.3, 3.4, 0.125}, {7.0, -0.75, 2.0},
    };
    for (const auto &p : pts)
        s.add_data(p[0], p[1], p[2]);
    return s;
}

/** Java prints an enum by its name; these must match those names exactly. */
const char *name(Interpolation m) {
    switch (m) {
    case Interpolation::Nearest: return "Nearest";
    case Interpolation::Linear: return "Linear";
    case Interpolation::Quadratic: return "Quadratic";
    case Interpolation::CubicSimple: return "CubicSimple";
    case Interpolation::Cubic: return "Cubic";
    case Interpolation::Cubic2: return "Cubic2";
    case Interpolation::CubicDerivInit: return "CubicDerivInit";
    case Interpolation::Cubic2DerivInit: return "Cubic2DerivInit";
    case Interpolation::CubicDeriv: return "CubicDeriv";
    case Interpolation::Cubic2Deriv: return "Cubic2Deriv";
    default: return "?";
    }
}

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        if (nl == std::string::npos) {
            out.push_back(text.substr(start));
            break;
        }
        out.push_back(text.substr(start, nl - start));
        start = nl + 1;
    }
    return out;
}

} // namespace

TEST(data_discrete_set_matches_jvm) {
    const Interpolation methods[] = {
        Interpolation::Nearest,        Interpolation::Linear,
        Interpolation::Quadratic,      Interpolation::CubicSimple,
        Interpolation::Cubic,          Interpolation::Cubic2,
        Interpolation::CubicDerivInit, Interpolation::Cubic2DerivInit,
        Interpolation::CubicDeriv,     Interpolation::Cubic2Deriv,
    };
    DiscreteSet s = build();
    std::string sb;
    sb += "count=" + std::to_string(s.get_count()) + "\n";
    sb += "version=" + std::to_string(s.get_version()) + "\n";
    auto xr = s.get_x_range();
    sb += "x_range=" + d(xr.first) + "," + d(xr.second) + "\n";
    auto yr = s.get_y_range();
    sb += "y_range=" + d(yr.first) + "," + d(yr.second) + "\n";
    for (int i = 0; i < s.get_count(); i++)
        sb += "pt." + std::to_string(i) + " x=" + d(s.get_x_value(i)) +
              " y=" + d(s.get_y_value(i)) + " d=" + d(s.get_d_value(i)) + "\n";
    for (int i = 0; i < s.get_count() - 1; i++)
        sb += "xi." + std::to_string(i) + " " + d(s.get_x_interval(i)) + "\n";
    for (double x = -0.5; x <= 7.6; x += 0.35)
        sb += "interval(" + d(x) + ")=" + std::to_string(s.get_interval(x)) +
              " nearest=" + std::to_string(s.get_nearest(x)) + "\n";
    for (auto m : methods) {
        s.set_interpolation(m);
        for (int deriv = 0; deriv <= 3; deriv++) {
            sb += std::string(name(m)) + ".d" + std::to_string(deriv) + ":";
            for (double x = 0.0; x <= 7.0001; x += 0.25)
                sb += " " + d(s.interpolate(x, deriv));
            sb += "\n";
        }
    }
    s.set_interpolation(Interpolation::Cubic);
    double before = s.interpolate(3.5);
    s.add_data(3.5, 9.9, 0.0);
    double after = s.interpolate(3.5);
    sb += "invalidate before=" + d(before) + " after=" + d(after) + "\n";
    sb += "count_after=" + std::to_string(s.get_count()) +
          " version=" + std::to_string(s.get_version()) + "\n";
    s.add_data(3.5, -9.9, 0.0);
    sb += "count_replaced=" + std::to_string(s.get_count()) +
          " y=" + d(s.get_y_value(s.get_interval(3.5) - 1)) + "\n";

    auto actual = lines(sb);
    const std::size_t expected_n =
        sizeof(EXPECTED_DATA_LINES) / sizeof(EXPECTED_DATA_LINES[0]);
    CHECK_EQ(actual.size(), expected_n);
    for (std::size_t i = 0; i < actual.size() && i < expected_n; i++)
        CHECK_STR_EQ(actual[i], std::string(EXPECTED_DATA_LINES[i]));
}
