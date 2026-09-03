// Numeric formatting for the ported reports and dumps.
#ifndef REDUKTI_TEXT_H
#define REDUKTI_TEXT_H

#include <string>

namespace redukti {

/**
 * Formats a double the way the Java version presents one.
 *
 * Presentation follows java.lang.Double.toString: always at least one
 * fractional digit ("1.0", never "1"), plain notation within [1e-3, 1e7),
 * scientific outside it written as `1.0E7` / `1.0E-4` -- no `+`, no zero
 * padding on the exponent.
 *
 * Digits come from Ryu (third_party/ryu), which yields the shortest decimal
 * that round-trips. That agrees with Java for every value the models actually
 * produce, but it is deliberately NOT bit-for-bit identical to Java in every
 * corner: JDK 19+ widens its candidate set when the shortest form has a single
 * significant digit, so e.g. Double.MIN_VALUE prints as 4.9E-324 there and
 * 5.0E-324 here. Chasing that is pointless -- JVM and libm differ in the last
 * bits of sin/cos/pow anyway, so ported results are compared numerically with a
 * tolerance rather than by diffing text.
 */
std::string doubleToString(double d);

/** Reproduces java.lang.Integer.toString(int) -- provided for symmetry. */
std::string intToString(int i);

/**
 * Reproduces Java's String.format("%<width>.<precision>g", value).
 *
 * This is NOT C's %g. Java keeps trailing zeros out to the stated number of
 * significant digits, where C strips them: at precision 4, Java renders 100.0
 * as "100.0" and 0.001 as "0.001000", while C gives "100" and "0.001". The
 * first-order and third-order reports are formatted entirely with %12.4g, so
 * using C's %g would quietly change every number in them.
 *
 * Java switches to scientific notation when the value, after rounding to
 * `precision` significant digits, falls outside [1e-4, 10^precision).
 * `width` right-aligns with spaces; pass 0 for no padding.
 */
std::string formatG(double value, int width, int precision);

/**
 * Java's `String.format("%.<precision>f", value)`.
 *
 * Java rounds HALF_UP here, where C printf rounds half to even, so this is not
 * a thin wrapper around snprintf. See the note in the implementation.
 */
std::string formatF(double value, int precision);

/**
 * The subset of java.text.DecimalFormat that M::decimal_format configures:
 * minimum 1 integer digit, at most `maxFractionDigits` fraction digits,
 * minimum 0 fraction digits, no grouping, decimal separator only when needed.
 *
 * Java rounds the *shortest round-tripping decimal* of the double, not the
 * exact binary value -- DecimalFormat feeds DigitList from
 * FloatingDecimal.getBinaryToASCIIConverter -- and then applies HALF_EVEN. The
 * two differ, so this rounds the Ryu digits the same way.
 */
class DecimalFormat {
public:
    explicit DecimalFormat(int maxFractionDigits)
        : max_fraction_digits_(maxFractionDigits) {}

    std::string format(double value) const;

private:
    int max_fraction_digits_;
};

} // namespace redukti

#endif // REDUKTI_TEXT_H
