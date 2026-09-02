// C++ port of org.redukti.mathlib.fftpack.ComplexDoubleFFT_Mixed
//
// Mixed-radix complex FFT, after FFTPACK.
// @author Baoshe Zhang
// @author Astronomical Instrument Group of University of Lethbridge.
//
// Only the complex transform is ported. The Java package also carries a
// RealDoubleFFT family (RealDoubleFFT, _Mixed, _Even, _Odd, _Even_Odd,
// _Odd_Odd -- about 1800 lines); nothing in the codebase references it, so it
// is deliberately left unported. See the note in ComplexDoubleFFT.h.
#ifndef REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_MIXED_H
#define REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_MIXED_H

#include <vector>

namespace redukti::mathlib::fftpack {

/**
 * The forward and backward transforms share the same butterflies, selected by
 * the `isign` argument, so both directions come essentially for free once
 * cfftf1 is ported.
 *
 * Arrays carry explicit `offset` arguments in place of the pointer arithmetic
 * of the original Fortran, exactly as the Java does.
 */
class ComplexDoubleFFT_Mixed {
public:
    virtual ~ComplexDoubleFFT_Mixed() = default;

protected:
    void passf2(int ido, int l1, std::vector<double> &cc, std::vector<double> &ch,
                std::vector<double> &wtable, int offset, int isign);
    void passf3(int ido, int l1, std::vector<double> &cc, std::vector<double> &ch,
                std::vector<double> &wtable, int offset, int isign);
    void passf4(int ido, int l1, std::vector<double> &cc, std::vector<double> &ch,
                std::vector<double> &wtable, int offset, int isign);
    void passf5(int ido, int l1, std::vector<double> &cc, std::vector<double> &ch,
                std::vector<double> &wtable, int offset, int isign);
    void passfg(std::vector<int> &nac, int ido, int ip, int l1, int idl1,
                std::vector<double> &cc, std::vector<double> &c1,
                std::vector<double> &c2, std::vector<double> &ch,
                std::vector<double> &ch2, std::vector<double> &wtable, int offset,
                int isign);

    void cfftf1(int n, std::vector<double> &c, std::vector<double> &wtable, int isign);

    /** Forward complex transform. */
    void cfftf(int n, std::vector<double> &c, std::vector<double> &wtable);

    /** Backward (unnormalized inverse) complex transform. */
    void cfftb(int n, std::vector<double> &c, std::vector<double> &wtable);

    void cffti1(int n, std::vector<double> &wtable);

    /** Builds the wavenumber table for size n. */
    void cffti(int n, std::vector<double> &wtable);
};

} // namespace redukti::mathlib::fftpack

#endif // REDUKTI_MATHLIB_FFTPACK_COMPLEXDOUBLEFFT_MIXED_H
