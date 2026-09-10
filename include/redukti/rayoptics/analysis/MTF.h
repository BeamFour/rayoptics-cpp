// C++ port of the analysis MTF classes:
//   BaseMTF, MTF, PolyMTF, MonochromaticGeometricMTF,
//   PolyChromaticGeometricMTF and MTFResultByFreq.
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_MTF_H
#define REDUKTI_RAYOPTICS_ANALYSIS_MTF_H

#include "redukti/rayoptics/analysis/Histogram.h"
#include "redukti/rayoptics/analysis/SpotIntercepts.h"

#include <memory>
#include <vector>

namespace redukti::rayoptics::analysis {

/** Frequency axis, complex spectra and normalized magnitudes for both meridians. */
class BaseMTF {
public:
    double pixel_size;
    /**
     * The FFT is calculated on LSF padded to twice the size
     * of the measurements - this is apparently required for FFT correctness.
     */
    int fft_size;
    /**
     * The MTF uses only the non-negative bins
     */
    int mtf_size;
    // frequencies
    std::vector<double> freq;
    // fourier transforms - complex numbers so has real and imaginary pairs
    /** Interleaved re/im, which is what ComplexDoubleFFT::ft expects. */
    std::vector<double> fft_x;
    std::vector<double> fft_y;
    // magnitudes - computed for the positive half
    std::vector<double> mag_x;
    std::vector<double> mag_y;

    virtual ~BaseMTF() = default;

protected:
    BaseMTF(int fft_size_, double pixel_size_);

    void compute_freq();
    void compute_fft(int xy);
    static void compute_magnitude(std::vector<double> &mag,
                                  const std::vector<double> &fft);
    void compute_magnitude(int xy);
};

/** The geometric MTF of one histogram: pad each LSF, transform, take magnitudes. */
class MTF : public BaseMTF {
public:
    // line spreads get padded to twice the size to improve FFT correctness
    std::vector<double> padded_lsf_x;
    std::vector<double> padded_lsf_y;

    explicit MTF(const Histogram &h2d_);

private:
    void pad_lfs(const std::vector<double> &lsf, std::vector<double> &padded_lsf);
    void compute_mtf(int xy);
    void compute_mtfs(const Histogram &histogram);
};

/**
 * Combines monochromatic MTFs for a field
 */
/** A weighted sum of monochromatic spectra, magnitudes taken once at the end. */
class PolyMTF : public BaseMTF {
public:
    PolyMTF(int fft_size_, double pixel_size_) : BaseMTF(fft_size_, pixel_size_) {}

    void add(const MTF &mono_mtf, double wt);

    void compute();

private:
    void add(const MTF &mono_mtf, int xy, double wt);
};

/**
 * Input : Traced rays through the system for a given field point and wavelength. x and y intersections with the image plane.
 *
 * 1. Bin these hits into a 2D intensity histogram to build the geometric PSF.
 * 2. Integrate the PSF along the perpendicular axis -> generate LSF for each axis.
 * 3. Compute the 1D Fourier transform for each LSF, then compute magnitude and normalize.
 * 4. Compute the related frequencies.
 */
class MonochromaticGeometricMTF {
public:
    double wvl;
    // 2d histogram
    Histogram h2d;
    MTF mtf;

    /**
     * Uses the default fixed grid. Prefer #MonochromaticGeometricMTF(SpotIntercepts, Histogram.Config)
     * with a field-level Histogram.Config so the window adapts to the spot
     * size and so all wavelengths of a field share one grid (required when the
     * results are combined into a PolyMTF).
     */
    explicit MonochromaticGeometricMTF(const SpotIntercepts &intercepts)
        : MonochromaticGeometricMTF(
              intercepts,
              Histogram::Config(Histogram::DEFAULT_NUM_BINS,
                                Histogram::DEFAULT_PIXEL_SIZE)) {}

    MonochromaticGeometricMTF(const SpotIntercepts &intercepts,
                              const Histogram::Config &cfg);
};

/**
 * Combines mono chromatic MTFs for a field
 */
class PolyChromaticGeometricMTF {
public:
    Histogram h2d;
    /** Null until compute(); MTF has no meaningful empty state. */
    std::unique_ptr<MTF> mtf;

    PolyChromaticGeometricMTF()
        : PolyChromaticGeometricMTF(Histogram::Config(Histogram::DEFAULT_NUM_BINS,
                                                      Histogram::DEFAULT_PIXEL_SIZE)) {}

    explicit PolyChromaticGeometricMTF(const Histogram::Config &cfg) : h2d(cfg) {}

    void add(const SpotIntercepts &intercepts, double wt) {
        h2d.accumulate(intercepts, wt);
    }

    void compute();
};

/** One spatial frequency sampled across every field, sagittal and tangential. */
class MTFResultByFreq {
public:
    int freq;
    // data for above freq across fields
    std::vector<double> sag_mtf_by_field;
    std::vector<double> tan_mtf_by_field;

    MTFResultByFreq(const std::vector<PolyMTF> &mtfs_by_field, int freq_);

    /**
     * Linearly interpolate an MTF magnitude at the requested spatial frequency.
     * The frequency axis is monotonically increasing and uniformly spaced at
     * 1/(fft_size*pixel_size), so the requested frequency almost never lands
     * exactly on a bin; interpolating avoids the ~0.5-1 cycle/mm error of picking
     * the nearest bin. Frequencies outside the sampled range are clamped to the
     * endpoints.
     */
    static double interpolate(const std::vector<double> &freq,
                              const std::vector<double> &mag, double f);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_MTF_H
