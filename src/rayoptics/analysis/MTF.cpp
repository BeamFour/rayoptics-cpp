// C++ port of the analysis MTF classes.
#include "redukti/rayoptics/analysis/MTF.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/fftpack/ComplexDoubleFFT.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>

namespace redukti::rayoptics::analysis {

namespace Orientation = util::Orientation;

// ---------------------------------------------------------------------------
// BaseMTF
// ---------------------------------------------------------------------------

BaseMTF::BaseMTF(int fft_size_, double pixel_size_)
    : pixel_size(pixel_size_), fft_size(fft_size_) {
    mtf_size = fft_size / 2 + 1;
    freq.assign(static_cast<std::size_t>(mtf_size), 0.0);
    fft_x.assign(static_cast<std::size_t>(fft_size) * 2, 0.0);
    fft_y.assign(static_cast<std::size_t>(fft_size) * 2, 0.0);
    mag_x.assign(static_cast<std::size_t>(mtf_size), 0.0);
    mag_y.assign(static_cast<std::size_t>(mtf_size), 0.0);
    compute_freq();
}

void BaseMTF::compute_freq() {
    for (std::size_t i = 0; i < freq.size(); i++)
        freq[i] = static_cast<double>(i) / (fft_size * pixel_size);
}

void BaseMTF::compute_fft(int xy) {
    std::vector<double> &fft = xy == Orientation::SAGITTAL ? fft_x : fft_y;
    mathlib::fftpack::ComplexDoubleFFT fft2d(fft_size);
    fft2d.ft(fft);
}

void BaseMTF::compute_magnitude(std::vector<double> &mag, const std::vector<double> &fft) {
    for (std::size_t i = 0; i < mag.size(); i++)
        mag[i] = std::hypot(fft[2 * i], fft[2 * i + 1]);
    double dc = mag[0];
    for (std::size_t k = 0; k < mag.size(); k++) {
        mag[k] = mag[k] / dc;
    }
}

void BaseMTF::compute_magnitude(int xy) {
    std::vector<double> &fft = xy == Orientation::SAGITTAL ? fft_x : fft_y;
    std::vector<double> &mag = xy == Orientation::SAGITTAL ? mag_x : mag_y;
    compute_magnitude(mag, fft);
}

// ---------------------------------------------------------------------------
// MTF
// ---------------------------------------------------------------------------

MTF::MTF(const Histogram &h2d_)
    : BaseMTF(h2d_.num_bins * 2, h2d_.pixel_size), h2d(&h2d_) {
    padded_lsf_x.assign(static_cast<std::size_t>(fft_size), 0.0);
    padded_lsf_y.assign(static_cast<std::size_t>(fft_size), 0.0);
    compute_mtfs();
}

void MTF::pad_lfs(const std::vector<double> &lsf, std::vector<double> &padded_lsf) {
    for (int i = 0; i < h2d->num_bins; i++)
        padded_lsf[static_cast<std::size_t>(i + h2d->num_bins / 2)] =
            lsf[static_cast<std::size_t>(i)];
}

void MTF::compute_mtf(int xy) {
    std::vector<double> &lsf = xy == Orientation::SAGITTAL ? padded_lsf_x : padded_lsf_y;
    std::vector<double> &fft = xy == Orientation::SAGITTAL ? fft_x : fft_y;
    for (std::size_t i = 0; i < static_cast<std::size_t>(fft_size); i++) {
        fft[2 * i] = lsf[i];
        fft[2 * i + 1] = 0.0;
    }
    compute_fft(xy);
}

void MTF::compute_mtfs() {
    pad_lfs(h2d->lsf_x, padded_lsf_x);
    pad_lfs(h2d->lsf_y, padded_lsf_y);
    compute_mtf(0);
    compute_mtf(1);
    compute_magnitude(0);
    compute_magnitude(1);
}

// ---------------------------------------------------------------------------
// PolyMTF
// ---------------------------------------------------------------------------

void PolyMTF::add(const MTF &mono_mtf, int xy, double wt) {
    std::vector<double> &fft = xy == Orientation::SAGITTAL ? fft_x : fft_y;
    const std::vector<double> &mono_fft =
        xy == Orientation::SAGITTAL ? mono_mtf.fft_x : mono_mtf.fft_y;
    for (std::size_t i = 0; i < static_cast<std::size_t>(fft_size); i++) {
        fft[2 * i] += wt * mono_fft[2 * i];
        fft[2 * i + 1] += wt * mono_fft[2 * i + 1];
    }
}

void PolyMTF::add(const MTF &mono_mtf, double wt) {
    add(mono_mtf, 0, wt);
    add(mono_mtf, 1, wt);
}

void PolyMTF::compute() {
    compute_magnitude(0);
    compute_magnitude(1);
}

// ---------------------------------------------------------------------------
// MonochromaticGeometricMTF / PolyChromaticGeometricMTF
// ---------------------------------------------------------------------------

namespace {

/**
 * Java builds the histogram in the constructor body before handing it to the
 * MTF; C++ initialises members in declaration order, so the accumulate/compute
 * has to happen inside the initialiser for `h2d` that `mtf` then reads.
 */
Histogram build_histogram(const SpotIntercepts &intercepts, const Histogram::Config &cfg) {
    Histogram h(cfg);
    h.accumulate(intercepts, 1.0);
    h.compute();
    return h;
}

} // namespace

MonochromaticGeometricMTF::MonochromaticGeometricMTF(const SpotIntercepts &intercepts,
                                                     const Histogram::Config &cfg)
    : wvl(intercepts.wvl), h2d(build_histogram(intercepts, cfg)), mtf(h2d) {}

void PolyChromaticGeometricMTF::compute() {
    h2d.compute();
    mtf = std::make_unique<MTF>(h2d);
}

// ---------------------------------------------------------------------------
// MTFResultByFreq
// ---------------------------------------------------------------------------

MTFResultByFreq::MTFResultByFreq(const std::vector<PolyMTF> &mtfs_by_field, int freq_)
    : freq(freq_) {
    sag_mtf_by_field.assign(mtfs_by_field.size(), 0.0);
    tan_mtf_by_field.assign(mtfs_by_field.size(), 0.0);
    for (std::size_t fi = 0; fi < mtfs_by_field.size(); fi++) {
        const auto &mtf = mtfs_by_field[fi];
        sag_mtf_by_field[fi] = interpolate(mtf.freq, mtf.mag_x, freq);
        tan_mtf_by_field[fi] = interpolate(mtf.freq, mtf.mag_y, freq);
    }
}

double MTFResultByFreq::interpolate(const std::vector<double> &freq,
                                    const std::vector<double> &mag, double f) {
    int n = static_cast<int>(freq.size());
    if (n == 0)
        throw IllegalArgumentException("empty MTF");
    if (f <= freq[0])
        return mag[0];
    if (f >= freq[static_cast<std::size_t>(n - 1)])
        return mag[static_cast<std::size_t>(n - 1)];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        // Java uses >>> here; both operands are non-negative and small, so an
        // unsigned shift is the same value, but it is kept unsigned anyway.
        int mid = static_cast<int>(static_cast<unsigned int>(lo + hi) >> 1);
        if (freq[static_cast<std::size_t>(mid)] <= f)
            lo = mid;
        else
            hi = mid;
    }
    double t = (f - freq[static_cast<std::size_t>(lo)]) /
               (freq[static_cast<std::size_t>(lo + 1)] - freq[static_cast<std::size_t>(lo)]);
    return mag[static_cast<std::size_t>(lo)] +
           t * (mag[static_cast<std::size_t>(lo + 1)] - mag[static_cast<std::size_t>(lo)]);
}

} // namespace redukti::rayoptics::analysis
