// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.plotter.ContrastPupilMapPlot.
#ifndef REDUKTI_PLOTTER_CONTRASTPUPILMAPPLOT_H
#define REDUKTI_PLOTTER_CONTRASTPUPILMAPPLOT_H

#include "redukti/analysis/ContrastPupilMap.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/specs/Field.h"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace redukti::render {
class RendererSvg;
}

namespace redukti::plotter {

/**
 * Draws ContrastPupilMaps for one field: a row per frequency, sagittal on the
 * left and tangential on the right.
 *
 * Each contrast sample is a circle at its reference ray's pupil position. The circle's
 * radius is proportional to the sample's contrast loss, from nothing when the pair is in
 * phase with the block to the full radius at half a wave. The line through each sample is
 * a clock hand at the pair's average phase, so neighbouring hands turning shows the shape of
 * the wavefront. The outline is the vignetted pupil, scaled in each panel to fill it.
 */
class ContrastPupilMapPlot {
public:
    /** One row of ContrastPupilMap[][]: sagittal, tangential; empty where Java holds null. */
    using Row = std::array<std::optional<analysis::ContrastPupilMap>, 2>;

    std::shared_ptr<const rayoptics::specs::FieldSnapshot> fld;
    std::vector<int> frequencies;
    /** [frequency][orientation] */
    std::vector<Row> maps;

    ContrastPupilMapPlot(std::shared_ptr<const rayoptics::specs::FieldSnapshot> fld,
                         std::vector<int> frequencies, std::vector<Row> maps);

    std::string plot() const;

private:
    void panel(render::RendererSvg &r, const analysis::ContrastPupilMap &map, int frequency,
               const mathlib::Vector2 &centre) const;

    /** The vignetted pupil, as the unit circle mapped through the piecewise scales. */
    void outline(render::RendererSvg &r, const mathlib::Vector2 &centre,
                 const mathlib::Vector2 &middle, double zoom) const;

    /** Mean of the two pupil half-widths along an axis, from the vignetting factors. */
    static double scale(double lower, double upper);

    static double half_width(double factor);
};

} // namespace redukti::plotter

#endif // REDUKTI_PLOTTER_CONTRASTPUPILMAPPLOT_H
