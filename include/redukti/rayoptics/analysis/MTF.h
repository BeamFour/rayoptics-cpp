// C++ port of the analysis MTF classes:
//   BaseMTF, MTF, PolyMTF, MonochromaticGeometricMTF,
//   PolyChromaticGeometricMTF and MTFResultByFreq.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
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
    int fft_size;
    int mtf_size;
    std::vector<double> freq;
    /** Interleaved re/im, which is what ComplexDoubleFFT::ft expects. */
    std::vector<double> fft_x;
    std::vector<double> fft_y;
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
    /** Borrowed; the owner of the histogram outlives this. */
    const Histogram *h2d;
    std::vector<double> padded_lsf_x;
    std::vector<double> padded_lsf_y;

    explicit MTF(const Histogram &h2d_);

private:
    void pad_lfs(const std::vector<double> &lsf, std::vector<double> &padded_lsf);
    void compute_mtf(int xy);
    void compute_mtfs();
};

/** A weighted sum of monochromatic spectra, magnitudes taken once at the end. */
class PolyMTF : public BaseMTF {
public:
    PolyMTF(int fft_size_, double pixel_size_) : BaseMTF(fft_size_, pixel_size_) {}

    void add(const MTF &mono_mtf, double wt);

    void compute();

private:
    void add(const MTF &mono_mtf, int xy, double wt);
};

class MonochromaticGeometricMTF {
public:
    double wvl;
    Histogram h2d;
    MTF mtf;

    explicit MonochromaticGeometricMTF(const SpotIntercepts &intercepts)
        : MonochromaticGeometricMTF(
              intercepts,
              Histogram::Config(Histogram::DEFAULT_NUM_BINS,
                                Histogram::DEFAULT_PIXEL_SIZE)) {}

    MonochromaticGeometricMTF(const SpotIntercepts &intercepts,
                              const Histogram::Config &cfg);
};

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
    std::vector<double> sag_mtf_by_field;
    std::vector<double> tan_mtf_by_field;

    MTFResultByFreq(const std::vector<PolyMTF> &mtfs_by_field, int freq_);

    static double interpolate(const std::vector<double> &freq,
                              const std::vector<double> &mag, double f);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_MTF_H
