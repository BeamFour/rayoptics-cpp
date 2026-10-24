// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.PolyDiffractionMTF.
#ifndef REDUKTI_ANALYSIS_POLYDIFFRACTIONMTF_H
#define REDUKTI_ANALYSIS_POLYDIFFRACTIONMTF_H

#include "redukti/analysis/DiffractionMTF.h"
#include "redukti/mathlib/Vector2.h"

#include <vector>

namespace redukti::analysis {

/**
 * Polychromatic diffraction MTF of one field: the modulus of the weighted sum of the
 * monochromatic OTFs, divided by the sum of the weights.
 *
 * Each monochromatic PSF is computed about its own chief-ray image point. Its complex
 * OTF is translated analytically to the reference wavelength's chief-ray image point,
 * so lateral colour appears as a phase difference without shifting the PSF through the
 * finite FFT window. The OTFs are evaluated at the requested frequency, since the FFT bins of each
 * wavelength fall at different frequencies. A wavelength contributes nothing beyond its
 * own cutoff.
 */
class PolyDiffractionMTF {
public:
    /** Borrowed. */
    rayoptics::specs::Field *fld;
    /** the point every wavelength is referenced to */
    mathlib::Vector2 image_pt_2d;
    /** one per wavelength with a non-zero weight, in the order given */
    std::vector<DiffractionMTF> mono;
    std::vector<double> weights;

    /**
     * @param wvls    wavelengths, nm; those with zero weight are skipped
     * @param weights one per wavelength
     * @param n_pupil rays across the pupil; even and at least 4
     * @param maxdim  FFT size; at least 2 * n_pupil
     */
    PolyDiffractionMTF(rayoptics::optical::OpticalModel *opm, rayoptics::specs::Field &fld,
                       const std::vector<double> &wvls, const std::vector<double> &weights,
                       double foc, int n_pupil, int maxdim);

    /** Sagittal MTF at frequency nu in cycles per system length unit. */
    double sagittal(double nu) const;

    /** Tangential MTF at frequency nu in cycles per system length unit. */
    double tangential(double nu) const;

private:
    double along(double nu, bool x) const;
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_POLYDIFFRACTIONMTF_H
