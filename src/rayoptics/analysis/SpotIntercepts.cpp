// C++ port of org.redukti.rayoptics.analysis.SpotIntercepts
#include "redukti/rayoptics/analysis/SpotIntercepts.h"

#include <cmath>
#include <limits>

namespace redukti::rayoptics::analysis {

using mathlib::Vector2;

SpotIntercepts::SpotIntercepts(const raytr::TraceGridByWvl &trace_data_)
    : trace_data(&trace_data_) {
    const auto n = trace_data_.grid.size();
    wvl = trace_data_.wvl;
    x.resize(n);
    y.resize(n);
    weights.resize(n);
    valid.resize(n);
    for (std::size_t i = 0; i < n; i++) {
        const auto &item = trace_data_.grid[i];
        valid[i] = item.valid ? 1 : 0;
        x[i] = valid[i] ? item.pupil.x : std::numeric_limits<double>::quiet_NaN();
        y[i] = valid[i] ? item.pupil.y : std::numeric_limits<double>::quiet_NaN();
        weights[i] = item.weight;
    }
}

Vector2 SpotIntercepts::compute_centroid() const {
    double cx = 0, cy = 0;
    double totalWeight = 0.0;
    for (std::size_t i = 0; i < trace_data->grid.size(); i++) {
        if (!valid[i])
            continue;
        cx += this->weights[i] * this->x[i];
        cy += this->weights[i] * this->y[i];
        totalWeight += this->weights[i];
    }
    // Written as !(w > 0) so a NaN total also takes this branch.
    if (!(totalWeight > 0.0))
        return Vector2(std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN());
    cx = cx / totalWeight;
    cy = cy / totalWeight;
    return Vector2(cx, cy);
}

void SpotIntercepts::adjust_to_centroid(const Vector2 &centroid) {
    for (std::size_t i = 0; i < trace_data->grid.size(); i++) {
        if (!valid[i])
            continue;
        this->x[i] -= centroid.x;
        this->y[i] -= centroid.y;
    }
}

} // namespace redukti::rayoptics::analysis
