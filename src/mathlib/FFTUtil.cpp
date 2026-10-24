// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.mathlib.FFTUtil.
#include "redukti/mathlib/FFTUtil.h"

#include "redukti/mathlib/fftpack/ComplexDoubleFFT.h"

namespace redukti::mathlib {

FFTUtil::ComplexArray FFTUtil::fftshift(const ComplexArray &x) {
    int n = static_cast<int>(x.size());
    int shift = n / 2;
    ComplexArray y(static_cast<std::size_t>(n), std::vector<double>(static_cast<std::size_t>(2 * n)));
    for (int i = 0; i < n; i++) {
        auto si = static_cast<std::size_t>((i + shift) % n);
        const auto &row = x[static_cast<std::size_t>(i)];
        for (int j = 0; j < n; j++) {
            auto sj = static_cast<std::size_t>((j + shift) % n);
            y[si][2 * sj] = row[static_cast<std::size_t>(2 * j)];
            y[si][2 * sj + 1] = row[static_cast<std::size_t>(2 * j + 1)];
        }
    }
    return y;
}

FFTUtil::ComplexArray FFTUtil::fft2(const ComplexArray &x) {
    int n = static_cast<int>(x.size());
    fftpack::ComplexDoubleFFT fft(n);
    ComplexArray y(x);
    for (auto &row : y)
        fft.ft(row);
    std::vector<double> column(static_cast<std::size_t>(2 * n));
    for (std::size_t j = 0; j < static_cast<std::size_t>(n); j++) {
        for (std::size_t i = 0; i < static_cast<std::size_t>(n); i++) {
            column[2 * i] = y[i][2 * j];
            column[2 * i + 1] = y[i][2 * j + 1];
        }
        fft.ft(column);
        for (std::size_t i = 0; i < static_cast<std::size_t>(n); i++) {
            y[i][2 * j] = column[2 * i];
            y[i][2 * j + 1] = column[2 * i + 1];
        }
    }
    return y;
}

} // namespace redukti::mathlib
