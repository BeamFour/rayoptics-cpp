// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.analysis.ContrastReport.
#include "redukti/analysis/ContrastReport.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/rayoptics/util/Orientation.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace redukti::analysis {

namespace Orientation = rayoptics::util::Orientation;

double ContrastReport::rms(const ContrastAnalysisResult::WavelengthResult &block, int orientation) {
    double sum = 0.0;
    double weight = 0.0;
    for (std::size_t i = 0; i < block.samples.size(); i++) {
        const auto &sample = block.samples[i];
        if (!sample.valid) continue;
        double residual = orientation == Orientation::SAGITTAL
                              ? block.sagittalResidual(static_cast<int>(i))
                              : block.tangentialResidual(static_cast<int>(i));
        if (!std::isfinite(residual)) continue;
        sum += residual * residual;
        weight += sample.weight;
    }
    return weight > 0.0 ? std::sqrt(sum / weight) : std::numeric_limits<double>::quiet_NaN();
}

bool ContrastReport::incomplete(const ContrastAnalysisResult::WavelengthResult &block) {
    for (const auto &sample : block.samples)
        if (!sample.valid) return true;
    return false;
}

std::string ContrastReport::report(rayoptics::optical::OpticalModel *opticalModel,
                                   const std::vector<int> &frequencies,
                                   const ContrastOptions &options) {
    const auto VALUE_FORMAT = mathlib::M::decimal_format(4);
    const auto LABEL_FORMAT = mathlib::M::decimal_format(4);
    std::string sb;
    // Java appends the Integer numSpokes, which prints "null" when unset.
    sb += "RMS OPD difference (waves), sampling " + std::to_string(options.numRings) +
          " rings x " +
          (options.numSpokes.has_value() ? std::to_string(*options.numSpokes) : "null") +
          " spokes" + ", calibrate " + yesNo(options.calibrateFrequency) +
          ", exit-pupil-aiming " + yesNo(options.aimExitPupil) + ", centering " +
          yesNo(options.centerResiduals) + "\n\n";

    const auto &fields = opticalModel->optical_spec->fov->fields;
    cell(sb, "freq");
    cell(sb, "wvl");
    cell(sb, "dir");
    for (const auto &field : fields)
        cell(sb, LABEL_FORMAT.format(field->y));
    sb += '\n';

    std::string notes;
    bool anyIncomplete = false;
    for (int frequency : frequencies) {
        std::optional<ContrastAnalysisResult> result;
        try {
            result.emplace(ContrastAnalysis::eval(opticalModel, options.with_frequency(frequency)));
        } catch (const IllegalArgumentException &e) {
            notes += "frequency " + std::to_string(frequency) + ": " + e.what() + '\n';
            continue;
        }
        const auto &wavelengths = opticalModel->optical_spec->wvls->wavelengths;
        for (std::size_t wi = 0; wi < wavelengths.size(); wi++) {
            for (int orientation : {Orientation::TANGENTIAL, Orientation::SAGITTAL}) {
                cell(sb, std::to_string(frequency));
                cell(sb, LABEL_FORMAT.format(wavelengths[wi]));
                cell(sb, Orientation::name(orientation));
                for (const auto &field : result->fields) {
                    const auto &block = field.wavelengths[wi];
                    double value = rms(block, orientation);
                    bool partial = incomplete(block);
                    anyIncomplete |= partial;
                    cell(sb, std::isnan(value) ? std::string("-")
                                               : VALUE_FORMAT.format(value) + (partial ? "*" : ""));
                }
                sb += '\n';
            }
        }
    }
    if (anyIncomplete)
        sb += "\n* some samples failed to trace; the value is over the samples that did\n";
    if (!notes.empty())
        sb += '\n' + notes;
    return sb;
}

void ContrastReport::cell(std::string &sb, const std::string &value) {
    sb.append(static_cast<std::size_t>(
                  std::max(1, COLUMN_WIDTH - static_cast<int>(value.size()))),
              ' ');
    sb += value;
}

std::string ContrastReport::yesNo(bool value) {
    return value ? "yes" : "no";
}

} // namespace redukti::analysis
