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
    int num_bins;
    double pixel_size;
    double hmin, hmax;
    std::vector<std::vector<double>> h2d;
    std::vector<double> lsf_x;
    std::vector<double> lsf_y;

    static constexpr double DEFAULT_PIXEL_SIZE = 0.001;
    static constexpr int DEFAULT_NUM_BINS = 512;

    class Config {
    public:
        int num_bins;
        double pixel_size;

        Config(int num_bins_, double pixel_size_)
            : num_bins(num_bins_), pixel_size(pixel_size_) {}
    };

    Histogram(int num_bins_, double pixel_size_);

    explicit Histogram(const Config &cfg) : Histogram(cfg.num_bins, cfg.pixel_size) {}

    static Config adaptiveConfig(double max_radius);

    static Config adaptiveConfig(double max_radius, double pixel_size, double margin,
                                 double min_window, int max_bins);

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
