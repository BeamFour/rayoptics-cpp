// C++ port of org.redukti.mathlib.fftpack.ComplexDoubleFFT
//
// FFT transform of a complex periodic sequence.
// @author Baoshe Zhang
// @author Astronomical Instrument Group of University of Lethbridge.
#ifndef REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_H
#define REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_H

#include "redukti/mathlib/fftpack/Complex1D.h"
#include "redukti/mathlib/fftpack/ComplexDoubleFFT_Mixed.h"

#include <vector>

namespace redukti::mathlib::fftpack {

/**
 * This is the only part of the Java fftpack package that anything uses:
 * BaseMTF::compute_fft does `new ComplexDoubleFFT(fft_size).ft(fft)` and that
 * is the whole of it. The RealDoubleFFT family in the Java package has no
 * callers and is not ported.
 *
 * fft_size is `num_bins * 2` and so is not necessarily a power of two, which
 * is why the mixed-radix implementation is required rather than a radix-2 one.
 *
 * The Java pom also declares a JTransforms dependency, but the only reference
 * to it in BaseMTF is commented out; the C++ port needs no FFT library.
 */
class ComplexDoubleFFT : public ComplexDoubleFFT_Mixed {
public:
    /**
     * norm_factor can be used to normalize this FFT transform: a call of the
     * forward transform (ft) followed by a call of the backward transform (bt)
     * multiplies the input sequence by norm_factor.
     */
    double norm_factor;

    /**
     * Construct a wavenumber table of size n. Sequences of the same size can
     * share a wavenumber table. The prime factorization of n together with a
     * tabulation of the trigonometric functions are computed and stored.
     *
     * @param n the size of a complex data sequence. When n is a product of
     *          small numbers (4, 2, 3, 5) this transform is very efficient.
     */
    explicit ComplexDoubleFFT(int n);

    /**
     * Forward complex FFT transform, in place.
     *
     * @param x 2*n reals representing n complex data: x[2*i] is the real part
     *          and x[2*i+1] the imaginary part of the i-th datum.
     */
    void ft(std::vector<double> &x);

    /** Forward complex FFT transform of n Complex data. */
    void ft(Complex1D &x);

    /** Backward complex FFT transform; the unnormalized inverse of ft. */
    void bt(std::vector<double> &x);

    /** Backward complex FFT transform of n Complex data. */
    void bt(Complex1D &x);

private:
    std::vector<double> wavetable;
    int ndim;
};

} // namespace redukti::mathlib::fftpack

#endif // REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_H
