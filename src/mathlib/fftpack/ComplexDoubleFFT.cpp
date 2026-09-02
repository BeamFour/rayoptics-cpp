// C++ port of org.redukti.mathlib.fftpack.ComplexDoubleFFT
#include "redukti/mathlib/fftpack/ComplexDoubleFFT.h"

#include "redukti/Exceptions.h"

namespace redukti::mathlib::fftpack {

ComplexDoubleFFT::ComplexDoubleFFT(int n) {
    ndim = n;
    norm_factor = n;
    if (wavetable.size() != static_cast<std::size_t>(4 * ndim + 15)) {
        wavetable.assign(static_cast<std::size_t>(4 * ndim + 15), 0.0);
    }
    cffti(ndim, wavetable);
}

void ComplexDoubleFFT::ft(std::vector<double> &x) {
    if (x.size() != static_cast<std::size_t>(2 * ndim))
        throw IllegalArgumentException(
            "The length of data can not match that of the wavetable");
    cfftf(ndim, x, wavetable);
}

void ComplexDoubleFFT::ft(Complex1D &x) {
    if (x.x.size() != static_cast<std::size_t>(ndim))
        throw IllegalArgumentException(
            "The length of data can not match that of the wavetable");
    std::vector<double> y(static_cast<std::size_t>(2 * ndim), 0.0);
    for (int i = 0; i < ndim; i++) {
        y[2 * i] = x.x[i];
        y[2 * i + 1] = x.y[i];
    }
    cfftf(ndim, y, wavetable);
    for (int i = 0; i < ndim; i++) {
        x.x[i] = y[2 * i];
        x.y[i] = y[2 * i + 1];
    }
}

void ComplexDoubleFFT::bt(std::vector<double> &x) {
    if (x.size() != static_cast<std::size_t>(2 * ndim))
        throw IllegalArgumentException(
            "The length of data can not match that of the wavetable");
    cfftb(ndim, x, wavetable);
}

void ComplexDoubleFFT::bt(Complex1D &x) {
    if (x.x.size() != static_cast<std::size_t>(ndim))
        throw IllegalArgumentException(
            "The length of data can not match that of the wavetable");
    std::vector<double> y(static_cast<std::size_t>(2 * ndim), 0.0);
    for (int i = 0; i < ndim; i++) {
        y[2 * i] = x.x[i];
        y[2 * i + 1] = x.y[i];
    }
    cfftb(ndim, y, wavetable);
    for (int i = 0; i < ndim; i++) {
        x.x[i] = y[2 * i];
        x.y[i] = y[2 * i + 1];
    }
}

} // namespace redukti::mathlib::fftpack
