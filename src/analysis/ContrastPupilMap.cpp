// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.ContrastPupilMap.
#include "redukti/analysis/ContrastPupilMap.h"

#include "redukti/Exceptions.h"
#include "redukti/analysis/ContrastReport.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <cmath>

namespace redukti::analysis {

namespace Orientation = rayoptics::util::Orientation;
using mathlib::M::PI;

ContrastPupilMap::ContrastPupilMap(const ContrastAnalysisResult::WavelengthResult &block,
                                   int orientation_) {
    Orientation::checked(orientation_);
    this->wavelength = block.wavelength;
    this->orientation = orientation_;
    bool sagittal = orientation_ == Orientation::SAGITTAL;

    double re = 0, im = 0, w = 0;
    for (const auto &s : block.samples) {
        double d = difference(s, sagittal);
        if (!usable(s, d)) continue;
        re += s.weight * std::cos(2 * PI * d);
        im += s.weight * std::sin(2 * PI * d);
        w += s.weight;
    }
    if (!(w > 0))
        throw IllegalArgumentException("the contrast block has no valid samples");
    mean_phase = std::atan2(im, re);
    phasor = std::hypot(re, im) / w;
    rms = ContrastReport::rms(block, orientation_);

    for (const auto &s : block.samples) {
        double d = difference(s, sagittal);
        if (!usable(s, d)) continue;
        double loss = 0.5 * (1 - std::cos(2 * PI * d - mean_phase));
        double average = PI * (2 * s.referenceOpd + d);
        points.emplace_back(s.pupil, loss, average, s.weight / w);
    }
}

double ContrastPupilMap::difference(const ContrastAnalysisResult::Sample &s, bool sagittal) {
    return sagittal ? s.sagittalOpdDifference : s.tangentialOpdDifference;
}

bool ContrastPupilMap::usable(const ContrastAnalysisResult::Sample &s, double d) {
    return s.valid && std::isfinite(d) && s.weight > 0;
}

} // namespace redukti::analysis
