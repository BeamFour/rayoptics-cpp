// This code is part of Beam42 project (https://github.com/BeamFour/Beam42)
// Copyright 2025-2026 by Dibyendu Majumdar
// License GPL v3
// See LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.plotter.ContrastPupilMapPlot.
#include "redukti/plotter/ContrastPupilMapPlot.h"

#include "redukti/Exceptions.h"
#include "redukti/Text.h"
#include "redukti/mathlib/M.h"
#include "redukti/mathlib/Vector2Pair.h"
#include "redukti/rayoptics/util/Orientation.h"
#include "redukti/render/RendererSvg.h"

#include <algorithm>
#include <cmath>

namespace redukti::plotter {

using mathlib::Vector2;
using mathlib::Vector2Pair;
using mathlib::M::PI;
using render::Renderer;
using render::RendererSvg;
using render::Rgb;
namespace Orientation = rayoptics::util::Orientation;

namespace {
const Rgb LOSS(0.84f, 0.15f, 0.16f, 1.0f);
const Rgb HAND(0.12f, 0.47f, 0.71f, 1.0f);
const Rgb &OUTLINE = Rgb::rgb_black;
constexpr double PANEL_WIDTH = 2.6;
constexpr double PANEL_HEIGHT = 2.9;
} // namespace

ContrastPupilMapPlot::ContrastPupilMapPlot(
    std::shared_ptr<const rayoptics::specs::FieldSnapshot> fld_, std::vector<int> frequencies_,
    std::vector<Row> maps_) {
    if (frequencies_.size() != maps_.size())
        throw IllegalArgumentException("one row of maps is needed per frequency");
    this->fld = std::move(fld_);
    this->frequencies = std::move(frequencies_);
    this->maps = std::move(maps_);
}

std::string ContrastPupilMapPlot::plot() const {
    int rows = static_cast<int>(frequencies.size());
    RendererSvg r(660, 370 * rows, Rgb::rgb_white);
    r.set_window(Vector2Pair(Vector2(-PANEL_WIDTH / 2, -(rows - 1) * PANEL_HEIGHT - 1.3),
                             Vector2(1.5 * PANEL_WIDTH, 1.6)),
                 true);
    for (int row = 0; row < rows; row++)
        for (int column = 0; column < 2; column++) {
            const auto &map = maps[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)];
            if (map.has_value())
                panel(r, *map, frequencies[static_cast<std::size_t>(row)],
                      Vector2(column * PANEL_WIDTH, -row * PANEL_HEIGHT));
        }
    std::string out;
    r.write(out);
    return out;
}

void ContrastPupilMapPlot::panel(RendererSvg &r, const analysis::ContrastPupilMap &map,
                                 int frequency, const Vector2 &centre) const {
    // Scale the vignetted pupil to fill the panel, so off-axis maps stay legible.
    double x0 = -half_width(fld->vlx), x1 = half_width(fld->vux);
    double y0 = -half_width(fld->vly), y1 = half_width(fld->vuy);
    double zoom = 1 / std::max((x1 - x0) / 2, (y1 - y0) / 2);
    Vector2 middle((x0 + x1) / 2, (y0 + y1) / 2);
    outline(r, centre, middle, zoom);
    // Scale circles to the mean sample spacing so that full-size ones only just touch.
    double radius = 0.5 * zoom *
                    std::sqrt(PI * scale(fld->vlx, fld->vux) * scale(fld->vly, fld->vuy) /
                              std::max(1, static_cast<int>(map.points.size())));
    for (const auto &p : map.points) {
        auto at = centre.plus(p.pupil.minus(middle).times(zoom));
        if (p.loss > 0)
            r.draw_circle(at, radius * p.loss, LOSS, false);
        if (std::isfinite(p.average_phase)) {
            auto hand =
                Vector2(std::cos(p.average_phase), std::sin(p.average_phase)).times(0.8 * radius);
            r.draw_segment(at, at.plus(hand), HAND);
        }
    }
    double x = centre.x - 1.2;
    r.draw_text(Vector2(x, centre.y + 1.42), Vector2::vector2_10,
                "field " + formatF(fld->y, 2) + "  " + std::to_string(frequency) + " cyc/mm  " +
                    Orientation::name(map.orientation),
                Renderer::TextAlignLeft, 13, Rgb::rgb_black);
    r.draw_text(Vector2(x, centre.y + 1.22), Vector2::vector2_10,
                "phasor " + formatF(map.phasor, 3) + "   rms " + formatF(map.rms, 3) + " waves",
                Renderer::TextAlignLeft, 11, Rgb::rgb_gray);
}

void ContrastPupilMapPlot::outline(RendererSvg &r, const Vector2 &centre, const Vector2 &middle,
                                   double zoom) const {
    const int steps = 180;
    std::optional<Vector2> previous;
    for (int i = 0; i <= steps; i++) {
        double t = 2 * PI * i / steps;
        double x = std::cos(t), y = std::sin(t);
        Vector2 pupil(x * half_width(x < 0 ? fld->vlx : fld->vux),
                      y * half_width(y < 0 ? fld->vly : fld->vuy));
        auto point = centre.plus(pupil.minus(middle).times(zoom));
        if (previous.has_value())
            r.draw_segment(*previous, point, OUTLINE);
        previous = point;
    }
}

double ContrastPupilMapPlot::scale(double lower, double upper) {
    return 0.5 * (half_width(lower) + half_width(upper));
}

double ContrastPupilMapPlot::half_width(double factor) {
    return factor == 0.0 ? 1.0 : 1.0 - factor;
}

} // namespace redukti::plotter
