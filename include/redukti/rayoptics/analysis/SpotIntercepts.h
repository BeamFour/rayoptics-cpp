// C++ port of org.redukti.rayoptics.analysis.SpotIntercepts
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_ANALYSIS_SPOTINTERCEPTS_H
#define REDUKTI_RAYOPTICS_ANALYSIS_SPOTINTERCEPTS_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <vector>

namespace redukti::rayoptics::analysis {

/**
 * The image-plane intercepts of one traced grid, unpacked into parallel arrays.
 *
 * Note that the GridItem `pupil` this reads is not a pupil coordinate: the spot
 * callback stores the transverse aberration there, so x and y are image-plane
 * offsets. Invalid entries are NaN, which is why every consumer tests `valid`
 * before touching them.
 */
class SpotIntercepts {
public:
    double wvl;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> weights;
    /** Not vector<bool>: that is a bitfield and hands out proxies, not refs. */
    std::vector<char> valid;
    /**
     * Borrowed. The SpotAnalysisResult that builds this owns the traced grid
     * and outlives every SpotIntercepts made from it.
     */
    const raytr::TraceGridByWvl *trace_data;

    explicit SpotIntercepts(const raytr::TraceGridByWvl &trace_data_);

    mathlib::Vector2 compute_centroid() const;

    void adjust_to_centroid(const mathlib::Vector2 &centroid);
};

} // namespace redukti::rayoptics::analysis

#endif // REDUKTI_RAYOPTICS_ANALYSIS_SPOTINTERCEPTS_H
