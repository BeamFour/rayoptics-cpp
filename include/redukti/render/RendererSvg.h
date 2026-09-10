// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Goptical: Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
// Licensed under the GNU General Public License, version 3 or later; see LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.render.rendering.RendererSvg
#ifndef REDUKTI_RENDER_RENDERERSVG_H
#define REDUKTI_RENDER_RENDERERSVG_H

#include "redukti/Text.h"
#include "redukti/render/Renderer.h"

#include <string>
#include <vector>

namespace redukti::render {

/**
 * SVG file rendering driver
 *
 * This class implements a SVG graphic output driver.
 */
/** Draws into an SVG document held as text. */
class RendererSvg : public Renderer2d {
public:
    /**
     * Create a new svg renderer with given resolution. The
     *
     * write function must be used to write svg to output
     * stream.
     */
    RendererSvg(double width, double height, const Rgb &bg);

    RendererSvg(double width, double height) : RendererSvg(width, height, Rgb::rgb_white) {}

    RendererSvg() : RendererSvg(800, 600, Rgb::rgb_white) {}

    void clear();

    void group_begin(const std::string &name) override;
    void group_end() override;

    void draw_point(const mathlib::Vector2 &p, const Rgb &rgb, PointStyle s) override;
    void draw_segment(const mathlib::Vector2Pair &l, const Rgb &rgb) override;
    void draw_circle(const mathlib::Vector2 &c, double r, const Rgb &rgb,
                     bool filled) override;
    void draw_text(const mathlib::Vector2 &v, const mathlib::Vector2 &dir,
                   const std::string &str, int a, int size, const Rgb &rgb) override;
    void draw_polygon(const std::vector<mathlib::Vector2> &array, const Rgb &rgb,
                      bool filled, bool closed) override;

    using Renderer2d::draw_point;
    using Renderer2d::draw_segment;
    using Renderer2d::draw_text;

    double y_trans_pos(double y) const override;

    /** Java's `write(StringBuilder)`; appends and returns the same buffer. */
    std::string &write(std::string &s) const;

private:
    std::string format(double value) const { return _decimal_format.format(value); }

    mathlib::Vector2 trans_pos(const mathlib::Vector2 &v) const;

    void svg_begin_rect(double x1, double y1, double x2, double y2, bool terminate);
    void svg_begin_line(double x1, double y1, double x2, double y2, bool terminate);
    void svg_begin_ellipse(double x, double y, double rx, double ry, bool terminate);
    void svg_begin_use(const std::string &id, double x, double y, bool terminate);
    void write_srgb(const Rgb &rgb);
    void svg_add_fill(const Rgb &rgb);
    void svg_add_stroke(const Rgb &rgb);
    void svg_add_dasharray();
    void svg_add_stroke_width();
    void svg_end();

    std::string _out;
    DecimalFormat _decimal_format;
};

} // namespace redukti::render

#endif // REDUKTI_RENDER_RENDERERSVG_H
