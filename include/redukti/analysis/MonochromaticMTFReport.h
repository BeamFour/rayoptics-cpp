// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.MonochromaticMTFReport.
#ifndef REDUKTI_ANALYSIS_MONOCHROMATICMTFREPORT_H
#define REDUKTI_ANALYSIS_MONOCHROMATICMTFREPORT_H

#include "redukti/analysis/MTF.h"

#include <string>
#include <vector>

namespace redukti::analysis {

/**
 * Tabulates the monochromatic geometric MTF, one row per frequency, wavelength and
 * orientation, one column per field. The layout is that of ContrastReport, so
 * the two can be read side by side.
 */
class MonochromaticMTFReport {
public:
    /**
     * @param fields      field labels, one per entry of mtfs
     * @param mtfs        per field, the monochromatic MTF of each wavelength, in the same
     *                    wavelength order for every field
     * @param frequencies cycles/mm
     */
    static std::string report(const std::vector<double> &fields,
                              const std::vector<std::vector<MonochromaticGeometricMTF>> &mtfs,
                              const std::vector<int> &frequencies);

private:
    static constexpr int COLUMN_WIDTH = 9;

    /** Null when the field has no MTF at wvl. */
    static const MonochromaticGeometricMTF *find(const std::vector<MonochromaticGeometricMTF> &field,
                                                 double wvl);

    static void cell(std::string &sb, const std::string &value);
};

} // namespace redukti::analysis

#endif // REDUKTI_ANALYSIS_MONOCHROMATICMTFREPORT_H
