// C++ port of org.redukti.rayoptics.analysis.PupilMapAnalysis
#include "redukti/rayoptics/analysis/PupilMapAnalysis.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/exceptions/TraceException.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"

#include <algorithm>
#include <cmath>
#include <memory>

namespace redukti::rayoptics::analysis {

namespace {

/** Java's `%<width>.<precision>f`: formatF, right-aligned in the width. */
std::string padded(double value, int width, int precision) {
    std::string text = formatF(value, precision);
    return text.size() < static_cast<std::size_t>(width)
               ? std::string(static_cast<std::size_t>(width) - text.size(), ' ') + text
               : text;
}

} // namespace

std::string PupilMapAnalysis::MappingQuality::toString() const {
    return "MappingQuality[sampled=" + doubleToString(sampled) +
           ", covered=" + doubleToString(covered) + "]";
}

PupilMapAnalysis::PupilMapForField::PupilMapForField(int fi_, specs::Field *fld_,
                                                     double reach_, int num_samples_,
                                                     std::vector<Sample> samples_)
    : fi(fi_), fld(fld_), reach(reach_), num_samples(num_samples_),
      samples(std::move(samples_)) {
    double maxX = 0.0, maxY = 0.0;
    int nominal = 0, nominalPassed = 0;
    int passed = 0;
    int inPiecewise = 0, inPiecewisePassed = 0;
    int inEllipse = 0, inEllipsePassed = 0;
    for (const Sample &s : samples) {
        bool piecewise = inside_piecewise(s.x, s.y, *fld);
        bool ellipse = inside_ellipse(s.x, s.y, *fld);
        if (piecewise)
            inPiecewise++;
        if (ellipse)
            inEllipse++;
        if (s.x * s.x + s.y * s.y <= 1.0)
            nominal++;
        if (!s.passed)
            continue;
        maxX = std::max(maxX, std::abs(s.x));
        maxY = std::max(maxY, std::abs(s.y));
        passed++;
        if (piecewise)
            inPiecewisePassed++;
        if (ellipse)
            inEllipsePassed++;
        if (s.x * s.x + s.y * s.y <= 1.0)
            nominalPassed++;
    }
    passed_x = maxX;
    passed_y = maxY;
    nominal_passed = nominal == 0 ? 0.0 : static_cast<double>(nominalPassed) / nominal;
    piecewise_quality = MappingQuality{
        inPiecewise == 0 ? 0.0 : static_cast<double>(inPiecewisePassed) / inPiecewise,
        passed == 0 ? 0.0 : static_cast<double>(inPiecewisePassed) / passed};
    ellipse_quality = MappingQuality{
        inEllipse == 0 ? 0.0 : static_cast<double>(inEllipsePassed) / inEllipse,
        passed == 0 ? 0.0 : static_cast<double>(inEllipsePassed) / passed};
}

std::string PupilMapAnalysis::PupilMapResult::toString() const {
    std::string sb = "field  vig scales x-,x+,y-,y+       passed |x|,|y|   nominal pupil"
                     "   piecewise sampled/covered   ellipse sampled/covered\n";
    for (const PupilMapForField &map : maps) {
        // Java's %n is the platform line separator; the report keeps the bare newline
        // every other output here uses.
        sb += padded(map.fld->yv(), 5, 2) + "  " + padded(scale(map.fld->vlx), 6, 3) + " " +
              padded(scale(map.fld->vux), 6, 3) + " " + padded(scale(map.fld->vly), 6, 3) +
              " " + padded(scale(map.fld->vuy), 6, 3) + "   " + padded(map.passed_x, 6, 3) +
              " " + padded(map.passed_y, 6, 3) + "   " +
              padded(100 * map.nominal_passed, 11, 0) + "%   " +
              padded(100 * map.piecewise_quality.sampled, 14, 0) + "% " +
              formatF(100 * map.piecewise_quality.covered, 0) + "%   " +
              padded(100 * map.ellipse_quality.sampled, 12, 0) + "% " +
              formatF(100 * map.ellipse_quality.covered, 0) + "%\n";
    }
    return sb;
}

PupilMapAnalysis::PupilMapResult PupilMapAnalysis::eval(optical::OpticalModel *opm) {
    return eval(opm, DEFAULT_NUM_SAMPLES, std::nullopt);
}

PupilMapAnalysis::PupilMapResult PupilMapAnalysis::eval(
    optical::OpticalModel *opm, int num_samples, const std::optional<std::vector<int>> &fields) {
    if (num_samples < 3)
        throw IllegalArgumentException("a pupil map needs at least 3 samples per axis");
    auto *osp = opm->optical_spec.get();
    double wvl = opm->seq_model->central_wavelength();
    raytr::TraceOptions options;
    options.check_apertures = true;
    // Raw pupil coordinates: the map measures the bundle, so it must not be told
    // where the bundle is supposed to be.
    options.apply_vignetting = false;
    // trace_safe only reports the blocking surface when asked for the error.
    options.rayerr_filter = "summary";

    PupilMapResult result;
    auto &fieldList = osp->fov->fields;
    for (int fi = 0; fi < static_cast<int>(fieldList.size()); fi++) {
        if (fields.has_value() &&
            std::find(fields->begin(), fields->end(), fi) == fields->end())
            continue;
        specs::Field &fld = *fieldList[static_cast<std::size_t>(fi)];
        double reach = reach_for(fld);
        // The factors are the hypothesis being checked, so cleared or underestimated
        // factors can make this initial square smaller than the transmitted bundle.
        // Passing boundary rays signal that we must expand: otherwise the map would
        // miss usable pupil area and overstate the mapping's coverage of the bundle.
        for (int expansion = 0;; expansion++) {
            auto samples = sample_grid(opm, fld, wvl, options, reach, num_samples);
            bool boundaryPassed = false;
            auto at = [&](int index) { return samples[static_cast<std::size_t>(index)].passed; };
            for (int k = 0; k < num_samples; k++) {
                if (at(k) || at((num_samples - 1) * num_samples + k) || at(k * num_samples) ||
                    at(k * num_samples + num_samples - 1)) {
                    boundaryPassed = true;
                    break;
                }
            }
            if (!boundaryPassed) {
                result.maps.emplace_back(fi, &fld, reach, num_samples, std::move(samples));
                break;
            }
            if (expansion == MAX_REACH_EXPANSIONS)
                throw IllegalStateException("Pupil map for field " + intToString(fi) +
                                            " still has passing boundary rays at reach " +
                                            doubleToString(reach));
            reach *= 1.5;
        }
    }
    return result;
}

std::vector<PupilMapAnalysis::Sample> PupilMapAnalysis::sample_grid(
    optical::OpticalModel *opm, specs::Field &fld, double wvl,
    const raytr::TraceOptions &options, double reach, int num_samples) {
    std::vector<Sample> samples;
    samples.reserve(static_cast<std::size_t>(num_samples) * static_cast<std::size_t>(num_samples));
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_samples; j++) {
            double x = -reach + 2 * reach * i / (num_samples - 1.0);
            double y = -reach + 2 * reach * j / (num_samples - 1.0);
            auto ray = raytr::Trace::trace_safe(opm, mathlib::Vector2(x, y), fld, wvl, options);
            bool passed = ray.pkg != nullptr && ray.err == nullptr;
            int blocked_by =
                std::dynamic_pointer_cast<exceptions::TraceRayBlockedException>(ray.err)
                    ? ray.err->surf
                    : -1;
            samples.push_back(Sample{x, y, passed, blocked_by});
        }
    }
    return samples;
}

double PupilMapAnalysis::reach_for(const specs::Field &fld) {
    double widest = std::max(std::max(scale(fld.vlx), scale(fld.vux)),
                             std::max(scale(fld.vly), scale(fld.vuy)));
    return 1.15 * std::max(1.0, widest);
}

bool PupilMapAnalysis::inside_piecewise(double x, double y, const specs::Field &fld) {
    double xScale = fld.vignetting_scale_x(x);
    double yScale = fld.vignetting_scale_y(y);
    if (!(xScale > 0.0) || !(yScale > 0.0))
        return false;
    double u = x / xScale, v = y / yScale;
    return u * u + v * v <= 1.0;
}

bool PupilMapAnalysis::inside_ellipse(double x, double y, const specs::Field &fld) {
    double xScale = ellipse_scale(fld.vlx, fld.vux);
    double yScale = ellipse_scale(fld.vly, fld.vuy);
    if (!(xScale > 0.0) || !(yScale > 0.0))
        return false;
    double u = (x - ellipse_offset(fld.vlx, fld.vux)) / xScale;
    double v = (y - ellipse_offset(fld.vly, fld.vuy)) / yScale;
    return u * u + v * v <= 1.0;
}

} // namespace redukti::rayoptics::analysis
