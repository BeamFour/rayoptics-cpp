// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.ContrastReport.
#ifndef REDUKTI_ANALYSIS_CONTRASTREPORT_H
#define REDUKTI_ANALYSIS_CONTRASTREPORT_H

#include "redukti/analysis/ContrastAnalysis.h"

#include <string>
#include <vector>

namespace redukti::analysis {

/**
 * Tabulates the wavefront differences the contrast optimizer works on, one row per
 * frequency, wavelength and orientation, one column per field.
 *
 * Each value is the quadrature-weighted RMS of the OPD difference over the valid
 * samples of that block, in waves: sqrt(sum w dW^2 / sum w). Its square is the
 * block's contribution to the contrast merit before field, orientation and wavelength
 * weights, normalized by the weight of the samples that traced.
 */
class ContrastReport {
public:
    /** Weighted RMS of one orientation's OPD differences, or NaN if no sample is valid. */
    static double rms(const ContrastAnalysisResult::WavelengthResult &block, int orientation);

    /**
     * Evaluates the contrast analysis at each frequency with options and formats
     * the table. A frequency the analysis cannot sample, for instance one with no common
     * pupil overlap, is reported as a line of its own rather than failing the report.
     *
     * @param options supplies everything but the frequency: sampling, calibration,
     *                exit-pupil aiming and centering
     */
    static std::string report(rayoptics::optical::OpticalModel *opticalModel,
                              const std::vector<int> &frequencies,
                              const ContrastOptions &options);

private:
    static constexpr int COLUMN_WIDTH = 9;

    /** True if any sample of the block failed to trace. */
    static bool incomplete(const ContrastAnalysisResult::WavelengthResult &block);

    static void cell(std::string &sb, const std::string &value);

    static std::string yesNo(bool value);
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_CONTRASTREPORT_H
