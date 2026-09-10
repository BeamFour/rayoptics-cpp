// C++ port of org.redukti.rayoptics.analysis.Histogram
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_HISTOGRAM_H
#define REDUKTI_RAYOPTICS_ANALYSIS_HISTOGRAM_H

#include "redukti/rayoptics/analysis/SpotIntercepts.h"

#include <vector>

namespace redukti::rayoptics::analysis {

/**
 * A square 2D histogram of spot intercepts, plus the two line spread functions
 * obtained by summing it along each axis. This is what feeds the geometric MTF.
 */
class Histogram {
public:
    /**
     * Defines the size of the grid used as histogram
     */
    int num_bins;
    double pixel_size;
    double hmin, hmax;
    // 2d histogram
    std::vector<std::vector<double>> h2d;
    // line spreads
    std::vector<double> lsf_x;
    std::vector<double> lsf_y;

    // Default grid used when no spot extent is available: a 0.512 mm window
    // (+/-0.256 mm) sampled at 1 micron. This reproduces the historical fixed grid.
    static constexpr double DEFAULT_PIXEL_SIZE = 0.001;
    static constexpr int DEFAULT_NUM_BINS = 512;

    /**
     * Bin count and bin size for the spot histogram. All monochromatic MTFs that
     * are combined into a single polychromatic PolyMTF for one field MUST
     * share the same Config, otherwise their FFT sizes and frequency axes differ
     * and the complex-OTF summation in PolyMTF#add is invalid.
     */
    class Config {
    public:
        int num_bins;
        double pixel_size;

        Config(int num_bins_, double pixel_size_)
            : num_bins(num_bins_), pixel_size(pixel_size_) {}
    };

    Histogram(int num_bins_, double pixel_size_);

    explicit Histogram(const Config &cfg) : Histogram(cfg.num_bins, cfg.pixel_size) {}

    /**
     * Adaptive grid sized to contain the geometric spot with margin, using the
     * default bin size and limits. See the overload for the meaning of the tuning
     * parameters.
     *
     * @param max_radius maximum spot radius (in lens units, i.e. mm) for the field,
     *                   taken across all wavelengths so every wavelength shares one grid
     */
    static Config adaptiveConfig(double max_radius);

    /**
     * Choose a histogram grid that:
     * <ul>
     *   <li>keeps pixel_size fine enough to reach the desired maximum
     *       spatial frequency (Nyquist = 1/(2*pixel_size)),</li>
     *   <li>makes the window wide enough to contain the spot (radius * margin), so
     *       rays are not clipped by #accumulate and the LSF tapers to zero,</li>
     *   <li>never shrinks the window below min_window, so the MTF frequency
     *       step (1/(2*window)) stays small enough to sample low frequencies, and</li>
     *   <li>caps the bin count at max_bins to bound the O(num_bins^2) memory,
     *       coarsening the bin size instead if the spot is very large.</li>
     * </ul>
     * When the spot is small this returns the historical default grid unchanged.
     *
     * @param max_radius maximum spot radius (mm) across all wavelengths of the field
     * @param pixel_size preferred bin size (mm)
     * @param margin     window half-width as a multiple of max_radius (>= 1)
     * @param min_window minimum full window width (mm)
     * @param max_bins   upper bound on bin count (power of two recommended)
     */
    static Config adaptiveConfig(double max_radius, double pixel_size, double margin,
                                 double min_window, int max_bins);

    // keep the whole spot in the window by coarsening the bin size
    void accumulate(const SpotIntercepts &intercepts, double wt);

    void compute();

private:
    static int nextPow2(int n);

    void normalize_histogram();
    void build_lsf(int xy);
    void build_lsfs();
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_HISTOGRAM_H
