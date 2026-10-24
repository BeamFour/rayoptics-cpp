// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.ContrastPupilMap.
#ifndef REDUKTI_ANALYSIS_CONTRASTPUPILMAP_H
#define REDUKTI_ANALYSIS_CONTRASTPUPILMAP_H

#include "redukti/analysis/ContrastAnalysis.h"
#include "redukti/mathlib/Vector2.h"

#include <vector>

namespace redukti::analysis {

/**
 * Where in the pupil one contrast block loses contrast, sample by sample.
 *
 * Each valid sample gets a loss m = (1 - cos(dphi - phi_mean)) / 2, where
 * dphi = 2 pi dW is its phase difference across the shear and phi_mean the
 * phase of the block's weighted mean phasor. m is 0 where the pair is in phase with
 * the block, 1/2 at a quarter wave and 1 at half a wave. The weighted losses add up to the
 * phasor modulus: |sum w exp(i dphi)| / sum w = 1 - 2 sum w m / sum w, the MTF
 * relative to an aberration-free pupil over the sampled region.
 *
 * Zemax's contrast pupil map uses dphi itself, not dphi - phi_mean. The
 * two agree when the mean phase is zero, as it is sagittally on a rotationally symmetric
 * lens. Referencing the mean phase keeps a pure image shift, which costs no contrast, from
 * showing as loss.
 *
 * Each sample also carries the average phase of its pair, pi (W(p) + W(p+s)), so a
 * plot can show the wavefront shape. It needs the reference OPD, which
 * ContrastAnalysisResult::Sample::referenceOpd records; it is NaN otherwise.
 *
 * The samples are those of the contrast analysis, so the map covers the region it
 * samples, not the full overlap of the pupil with its shifted copy.
 */
class ContrastPupilMap {
public:
    /**
     * @param loss          (1 - cos(dphi - phi_mean)) / 2, between 0 and 1
     * @param average_phase pi (W(p) + W(p+s)) in radians, or NaN
     * @param weight        the sample's quadrature weight
     */
    /** Java record Point(Vector2 pupil, double loss, double average_phase, double weight). */
    class Point {
    public:
        mathlib::Vector2 pupil;
        double loss;
        double average_phase;
        double weight;

        Point(const mathlib::Vector2 &pupil_, double loss_, double average_phase_, double weight_)
            : pupil(pupil_), loss(loss_), average_phase(average_phase_), weight(weight_) {}
    };

    double wavelength;
    int orientation;
    std::vector<Point> points;
    /** phase of the weighted mean phasor, radians */
    double mean_phase;
    /** |sum w exp(i dphi)| / sum w, equal to 1 - 2 * weighted mean loss */
    double phasor;
    /** weighted RMS of the OPD difference about the block's offset, waves */
    double rms;

    ContrastPupilMap(const ContrastAnalysisResult::WavelengthResult &block, int orientation);

private:
    static double difference(const ContrastAnalysisResult::Sample &s, bool sagittal);

    static bool usable(const ContrastAnalysisResult::Sample &s, double d);
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_CONTRASTPUPILMAP_H
