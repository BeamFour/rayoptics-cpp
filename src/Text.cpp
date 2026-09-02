#include "redukti/Text.h"

#include "ryu/ryu.h"

#include <array>
#include <cmath>
#include <cstddef>
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

} // namespace redukti
