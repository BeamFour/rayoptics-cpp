// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.MonochromaticMTFReport.
#include "redukti/analysis/MonochromaticMTFReport.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <algorithm>

namespace redukti::analysis {

namespace Orientation = rayoptics::util::Orientation;

std::string MonochromaticMTFReport::report(
    const std::vector<double> &fields,
    const std::vector<std::vector<MonochromaticGeometricMTF>> &mtfs,
    const std::vector<int> &frequencies) {
    const auto VALUE_FORMAT = mathlib::M::decimal_format(3);
    const auto LABEL_FORMAT = mathlib::M::decimal_format(4);
    if (fields.size() != mtfs.size())
        throw IllegalArgumentException("one field label is needed per field");
    std::string sb;
    sb += "Monochromatic geometric MTF\n\n";
    cell(sb, "freq");
    cell(sb, "wvl");
    cell(sb, "dir");
    for (double field : fields)
        cell(sb, LABEL_FORMAT.format(field));
    sb += '\n';
    if (mtfs.empty())
        return sb;

    // Rows follow the wavelengths of the first field; a field without one shows "-".
    const auto &wavelengths = mtfs[0];
    for (int frequency : frequencies) {
        for (const auto &wavelength : wavelengths) {
            for (int orientation : {Orientation::TANGENTIAL, Orientation::SAGITTAL}) {
                cell(sb, std::to_string(frequency));
                cell(sb, LABEL_FORMAT.format(wavelength.wvl));
                cell(sb, Orientation::name(orientation));
                for (const auto &field : mtfs) {
                    const auto *mono = find(field, wavelength.wvl);
                    if (mono == nullptr) {
                        cell(sb, "-");
                        continue;
                    }
                    const auto &mtf = mono->mtf;
                    double value = MTFResultByFreq::interpolate(
                        mtf.freq, orientation == Orientation::SAGITTAL ? mtf.mag_x : mtf.mag_y,
                        frequency);
                    cell(sb, VALUE_FORMAT.format(value));
                }
                sb += '\n';
            }
        }
    }
    return sb;
}

const MonochromaticGeometricMTF *MonochromaticMTFReport::find(
    const std::vector<MonochromaticGeometricMTF> &field, double wvl) {
    for (const auto &mono : field)
        if (mono.wvl == wvl) return &mono;
    return nullptr;
}

void MonochromaticMTFReport::cell(std::string &sb, const std::string &value) {
    sb.append(static_cast<std::size_t>(
                  std::max(1, COLUMN_WIDTH - static_cast<int>(value.size()))),
              ' ');
    sb += value;
}

} // namespace redukti::analysis
