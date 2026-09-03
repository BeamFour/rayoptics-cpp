#include "redukti/Text.h"

#include "ryu/ryu.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace redukti {

std::string intToString(int i) { return std::to_string(i); }

std::string doubleToString(double d) {
    if (std::isnan(d))
        return "NaN";
    if (std::isinf(d))
        return d < 0 ? "-Infinity" : "Infinity";
    if (d == 0.0)
        return std::signbit(d) ? "-0.0" : "0.0";

    const bool negative = std::signbit(d);
    const double m = std::abs(d);

    // Ryu emits the shortest round-tripping decimal in scientific form with the
    // exponent written bare, e.g. "1E0", "1.23456E2", "5E-324".
    std::array<char, 32> buf{};
    const int len = d2s_buffered_n(m, buf.data());
    const std::string sci(buf.data(), static_cast<std::size_t>(len));

    const std::size_t epos = sci.find('E');
    std::string digits;
    for (char c : sci.substr(0, epos)) {
        if (c != '.')
            digits.push_back(c);
    }
    while (digits.size() > 1 && digits.back() == '0')
        digits.pop_back();
    const int exponent = std::stoi(sci.substr(epos + 1));

    std::string out;
    if (negative)
        out.push_back('-');

    if (exponent >= -3 && exponent < 7) {
        // Plain decimal notation.
        if (exponent >= 0) {
            const std::size_t intDigits = static_cast<std::size_t>(exponent) + 1;
            if (digits.size() <= intDigits) {
                out += digits;
                out.append(intDigits - digits.size(), '0');
                out += ".0";
            } else {
                out += digits.substr(0, intDigits);
                out.push_back('.');
                out += digits.substr(intDigits);
            }
        } else {
            out += "0.";
            out.append(static_cast<std::size_t>(-exponent) - 1, '0');
            out += digits;
        }
    } else {
        // Computerized scientific notation: one digit, point, rest, 'E', exponent.
        out.push_back(digits[0]);
        out.push_back('.');
        if (digits.size() > 1)
            out += digits.substr(1);
        else
            out.push_back('0');
        out.push_back('E');
        out += std::to_string(exponent);
    }
    return out;
}


std::string formatG(double value, int width, int precision) {
    std::string out;
    if (std::isnan(value)) {
        out = "NaN";
    } else if (std::isinf(value)) {
        out = value < 0 ? "-Infinity" : "Infinity";
    } else if (value == 0.0) {
        // Java renders zero in decimal form with precision-1 fraction digits.
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.*f", precision - 1, value);
        out = buf;
    } else {
        // Decide the form from the magnitude AFTER rounding to `precision`
        // significant digits: 9.99999e-5 rounds up to 1.000e-4 and is then
        // rendered in decimal form.
        char sci[64];
        std::snprintf(sci, sizeof(sci), "%.*e", precision - 1, value);
        const double rounded = std::strtod(sci, nullptr);
        const double m = std::abs(rounded);
        const double upper = std::pow(10.0, static_cast<double>(precision));
        if (m >= 1e-4 && m < upper) {
            const int exp10 = static_cast<int>(std::floor(std::log10(m)));
            int frac = precision - 1 - exp10;
            if (frac < 0)
                frac = 0;
            char buf[512];
            std::snprintf(buf, sizeof(buf), "%.*f", frac, value);
            out = buf;
        } else {
            out = sci;
        }
    }
    if (static_cast<int>(out.size()) < width)
        out.insert(out.begin(), static_cast<std::size_t>(width) - out.size(), ' ');
    return out;
}

namespace {

/**
 * The shortest round-tripping decimal of |m|: the digit string with no
 * trailing zeros, plus how many of those digits sit before the decimal point,
 * so m = 0.DIGITS * 10^decimalAt. This is the same form Java DigitList holds.
 */
struct ShortestDecimal {
    std::string digits;
    int decimalAt;
};

ShortestDecimal shortestDecimal(double m) {
    std::array<char, 32> buf{};
    const int len = d2s_buffered_n(m, buf.data());
    const std::string sci(buf.data(), static_cast<std::size_t>(len));
    const std::size_t epos = sci.find('E');
    ShortestDecimal out;
    for (char c : sci.substr(0, epos)) {
        if (c != '.')
            out.digits.push_back(c);
    }
    while (out.digits.size() > 1 && out.digits.back() == '0')
        out.digits.pop_back();
    out.decimalAt = std::stoi(sci.substr(epos + 1)) + 1;
    return out;
}

} // namespace

std::string DecimalFormat::format(double value) const {
    if (std::isnan(value))
        return "NaN";
    // DecimalFormat writes the infinity sign, not the word.
    if (std::isinf(value))
        return value < 0 ? "-∞" : "∞";
    if (value == 0.0)
        return std::signbit(value) ? "-0" : "0";

    // Java feeds DecimalFormat from the shortest round-tripping decimal, so
    // anything past those digits reads as zero: 1e100 formats as a 1 and a
    // hundred zeros, not the exact binary expansion. When the rounding
    // position falls at or past the end of those digits nothing is rounded at
    // all, and the answer is just the digits placed against the decimal point.
    const ShortestDecimal sd = shortestDecimal(std::abs(value));
    const int len = static_cast<int>(sd.digits.size());
    if (sd.decimalAt + max_fraction_digits_ >= len) {
        std::string out;
        if (std::signbit(value))
            out.push_back('-');
        if (sd.decimalAt <= 0) {
            out.push_back('0');
            out.push_back('.');
            out.append(static_cast<std::size_t>(-sd.decimalAt), '0');
            out += sd.digits;
        } else if (sd.decimalAt >= len) {
            out += sd.digits;
            out.append(static_cast<std::size_t>(sd.decimalAt - len), '0');
        } else {
            out += sd.digits.substr(0, static_cast<std::size_t>(sd.decimalAt));
            out.push_back('.');
            out += sd.digits.substr(static_cast<std::size_t>(sd.decimalAt));
        }
        return out;
    }

    // Rounding does happen. Java rounds the exact binary value with HALF_EVEN
    // -- its tie-break consults whether the shortest decimal sits above or
    // below the true value, which comes to the same thing -- and that is
    // exactly what a correctly rounded printf gives. 0.15 is really
    // 0.1499...94, so one fraction digit yields "0.1", while an exact tie like
    // 2.5 goes to the even neighbour.
    const int digits = max_fraction_digits_;
    const int needed = std::snprintf(nullptr, 0, "%.*f", digits, value);
    std::string out(static_cast<std::size_t>(needed) + 1, ' ');
    std::snprintf(out.data(), out.size(), "%.*f", digits, value);
    out.resize(static_cast<std::size_t>(needed));

    // minimumFractionDigits is 0 and the separator is only shown when needed.
    if (out.find('.') != std::string::npos) {
        while (!out.empty() && out.back() == '0')
            out.pop_back();
        if (!out.empty() && out.back() == '.')
            out.pop_back();
    }
    return out;
}

} // namespace redukti
