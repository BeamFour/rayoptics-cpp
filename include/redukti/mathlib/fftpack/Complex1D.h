// C++ port of org.redukti.mathlib.fftpack.Complex1D
#ifndef REDUKTI_MATHLIB_FFTPACK_COMPLEX1D_H
#define REDUKTI_MATHLIB_FFTPACK_COMPLEX1D_H

#include <vector>

namespace redukti::mathlib::fftpack {

/**
  * Construct a 1-D complex data sequence.
*/
/** A 1-D complex data sequence held as separate real and imaginary arrays. */
class Complex1D {
public:
    /**
      * <em>x</em>[<em>i</em>] is the real part of <em>i</em>-th complex data.
    */
    /** x[i] is the real part of the i-th complex datum. */
    std::vector<double> x;
    /**
      * <em>y</em>[<em>i</em>] is the imaginary part of <em>i</em>-th complex data.
    */
    /** y[i] is the imaginary part of the i-th complex datum. */
    std::vector<double> y;
};

} // namespace redukti::mathlib::fftpack

#endif // REDUKTI_MATHLIB_FFTPACK_COMPLEX1D_H
