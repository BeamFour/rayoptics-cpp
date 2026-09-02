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

} // namespace redukti

#endif // REDUKTI_TEXT_H
