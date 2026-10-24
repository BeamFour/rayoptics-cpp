// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.mathlib.FFTUtil.
#ifndef REDUKTI_MATHLIB_FFTUTIL_H
#define REDUKTI_MATHLIB_FFTUTIL_H

#include <vector>

namespace redukti::mathlib {

/**
 * Replacements for the numpy.fft functions used by the ray-optics port.
 *
 * Complex arrays are square, with rows interleaved re, im: element [i][j] is
 * x[i][2*j] + 1j*x[i][2*j+1].
 */
class FFTUtil {
public:
    /** Java's double[][]: a square complex array, rows interleaved re, im. */
    using ComplexArray = std::vector<std::vector<double>>;

    /** numpy.fft.fftshift over both axes of a square complex array, rows interleaved re, im. */
    static ComplexArray fftshift(const ComplexArray &x);

    /** numpy.fft.fft2 of a square complex array, rows interleaved re, im. */
    static ComplexArray fft2(const ComplexArray &x);
};

} // namespace redukti::mathlib

#endif // REDUKTI_MATHLIB_FFTUTIL_H
