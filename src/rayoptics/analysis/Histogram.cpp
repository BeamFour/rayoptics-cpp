// C++ port of org.redukti.rayoptics.analysis.Histogram
#include "redukti/rayoptics/analysis/Histogram.h"

#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>

namespace redukti::rayoptics::analysis {

namespace Orientation = util::Orientation;

Histogram::Histogram(int num_bins_, double pixel_size_)
    : num_bins(num_bins_), pixel_size(pixel_size_) {
    auto width = pixel_size * num_bins;
    hmin = -width / 2;
    hmax = width / 2;
    h2d.assign(static_cast<std::size_t>(num_bins),
               std::vector<double>(static_cast<std::size_t>(num_bins), 0.0));
    lsf_x.assign(static_cast<std::size_t>(num_bins), 0.0);
    lsf_y.assign(static_cast<std::size_t>(num_bins), 0.0);
}

Histogram::Config Histogram::adaptiveConfig(double max_radius) {
    return adaptiveConfig(max_radius, DEFAULT_PIXEL_SIZE, 2.0,
                          DEFAULT_PIXEL_SIZE * DEFAULT_NUM_BINS, 2048);
}

Histogram::Config Histogram::adaptiveConfig(double max_radius, double pixel_size,
                                            double margin, double min_window,
                                            int max_bins) {
    if (!std::isfinite(max_radius) || max_radius < 0)
        max_radius = 0;
    double half_width = std::fmax(max_radius * margin, min_window / 2.0);
    double window = 2.0 * half_width;
    int bins = nextPow2(static_cast<int>(std::ceil(window / pixel_size)));
    if (bins > max_bins) {
        bins = max_bins;
        pixel_size = window / bins;
    }
    return Config(bins, pixel_size);
}

int Histogram::nextPow2(int n) {
    int p = 1;
    while (p < n)
        p <<= 1;
    return p;
}

void Histogram::accumulate(const SpotIntercepts &intercepts, double wt) {
    for (std::size_t i = 0; i < intercepts.x.size(); i++) {
        auto x = intercepts.x[i];
        auto y = intercepts.y[i];
        if (x < hmin || x > hmax || y < hmin || y > hmax)
            continue;
        int ix = static_cast<int>(std::floor(num_bins * (x - hmin) / (hmax - hmin)));
        int iy = static_cast<int>(std::floor(num_bins * (y - hmin) / (hmax - hmin)));
        if (ix < 0 || ix >= num_bins || iy < 0 || iy >= num_bins)
            continue;
        h2d[static_cast<std::size_t>(ix)][static_cast<std::size_t>(iy)] +=
            wt * intercepts.weights[i];
    }
}

void Histogram::normalize_histogram() {
    double sum = 0.0;
    for (int i = 0; i < num_bins; i++)
        for (int j = 0; j < num_bins; j++)
            sum += h2d[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    for (int i = 0; i < num_bins; i++)
        for (int j = 0; j < num_bins; j++)
            h2d[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] /= sum;
}

void Histogram::build_lsf(int xy) {
    std::vector<double> &lsf = xy == Orientation::SAGITTAL ? lsf_x : lsf_y;
    if (xy == Orientation::SAGITTAL) {
        for (int i = 0; i < num_bins; i++) {
            double s = 0;
            for (int j = 0; j < num_bins; j++) {
                s += h2d[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
            lsf[static_cast<std::size_t>(i)] = s;
        }
    } else {
        for (int j = 0; j < num_bins; j++) {
            double s = 0;
            for (int i = 0; i < num_bins; i++) {
                s += h2d[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
            }
            lsf[static_cast<std::size_t>(j)] = s;
        }
    }
    double lsfSum = 0;
    for (double v : lsf)
        lsfSum += v;
    for (int i = 0; i < num_bins; i++)
        lsf[static_cast<std::size_t>(i)] /= lsfSum;
}

void Histogram::build_lsfs() {
    build_lsf(0);
    build_lsf(1);
}

void Histogram::compute() {
    normalize_histogram();
    build_lsfs();
}

} // namespace redukti::rayoptics::analysis
