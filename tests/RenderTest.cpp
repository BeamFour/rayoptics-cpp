// End-to-end check of the org.redukti.render port against the JVM.
//
// Three SVG documents -- a full plot with axes, tics and value labels; the
// drawing primitives; and the 2D/perspective projections across a page layout
// -- are generated exactly as scratchpad/DumpRender.java generates them, and
// compared line by line with the JDK 25 output.
//
// Comparison is exact. The plot path runs through DecimalFormat for every tic
// label and through the data package for every interpolated point, so this
// pins those two together with the renderer.
#include "RenderExpected.h"
#include "TestHarness.h"

#include "redukti/data/DataSet.h"
#include "redukti/render/Plot.h"
#include "redukti/render/RendererSvg.h"

#include <string>
#include <vector>

namespace {

using redukti::data::DiscreteSet;
using redukti::data::Interpolation;
using redukti::mathlib::Vector2;
using redukti::mathlib::Vector2Pair;
using redukti::mathlib::Vector3;
using namespace redukti::render;

DiscreteSet build() {
    DiscreteSet s;
    const double pts[][2] = {
        {0.0, 1.0},  {0.7, 2.3},  {1.9, 0.4}, {3.1, -1.2},
        {4.05, 2.9}, {6.3, 3.4},  {7.0, -0.75},
    };
    for (const auto &p : pts)
        s.add_data(p[0], p[1]);
    s.set_interpolation(Interpolation::Cubic);
    return s;
}

std::vector<std::string> lines(const std::string &text) {
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start < text.size()) {
        auto nl = text.find('\n', start);
        if (nl == std::string::npos) {
            out.push_back(text.substr(start));
            break;
        }
        out.push_back(text.substr(start, nl - start));
        start = nl + 1;
    }
    while (!out.empty() && out.back().empty())
        out.pop_back();
    return out;
}

} // namespace

#define CHECK_SVG(actual, expected)                                                      \
    do {                                                                                 \
        auto _al = lines(actual);                                                        \
        const std::size_t _en = sizeof(expected) / sizeof(expected[0]);                  \
        CHECK_EQ(_al.size(), _en);                                                       \
        for (std::size_t _i = 0; _i < _al.size() && _i < _en; _i++)                       \
            CHECK_STR_EQ(_al[_i], std::string(expected[_i]));                            \
    } while (0)

TEST(render_plot_svg_matches_jvm) {
    DiscreteSet data = build();
    Plot plot;
    plot.set_title("Interpolated data");
    PlotData *pd = plot.add_plot_data(&data, Rgb::rgb_blue, "series",
                                      static_cast<int>(InterpolatePlot) |
                                          static_cast<int>(PointPlot) |
                                          static_cast<int>(LinePlot) |
                                          static_cast<int>(ValuePlot));
    pd->set_line_style(PlotData::LineStyle::Dashed);
    PlotAxes &axes = plot.get_axes();
    axes.set_label("x axis", PlotAxes::AxisMask::X);
    axes.set_label("y axis", PlotAxes::AxisMask::Y);
    axes.set_unit("m", true, true, 0, PlotAxes::AxisMask::XY);
    axes.set_show_grid(true);
    axes.set_show_frame(true);
    axes.set_tics_count(6, PlotAxes::AxisMask::XY);
    plot.fit_axes_range();

    RendererSvg r(800, 600);
    PlotRenderer().draw_plot(r, plot);
    std::string svg;
    r.write(svg);
    CHECK_SVG(svg, EXPECTED_PLOT);
}

TEST(render_primitives_svg_matches_jvm) {
    RendererSvg r(400, 300, Rgb::rgb_black);
    r.set_window(Vector2Pair(Vector2(-10, -10), Vector2(10, 10)), false);
    r.group_begin("shapes");
    for (int i = 0; i < 5; i++)
        r.draw_point(Vector2(-8 + i * 4, 8), Rgb::rgb_red, Renderer::point_style_of(i));
    r.draw_segment(Vector2Pair(Vector2(-9, -9), Vector2(9, 9)), Rgb::rgb_green);
    r.set_stroke_width(2.5);
    r.set_stroke_dasharray("6,4");
    r.draw_segment(Vector2Pair(Vector2(-9, 9), Vector2(9, -9)), Rgb::rgb_cyan);
    r.clear_stroke_dasharray();
    r.set_stroke_width(1.0);
    r.draw_circle(Vector2(0, 0), 5.0, Rgb::rgb_yellow, false);
    r.draw_circle(Vector2(3, 3), 2.0, Rgb::rgb_magenta, true);
    std::vector<Vector2> poly{Vector2(-6, -2), Vector2(-3, -6), Vector2(1, -4),
                              Vector2(-2, 0)};
    r.draw_polygon(poly, Rgb::rgb_white, false, true);
    r.draw_polygon(poly, Rgb::rgb_gray, true, false);
    r.draw_box(Vector2Pair(Vector2(-7, -8), Vector2(-1, -7)), Rgb::rgb_blue);
    r.draw_text(Vector2(0, -8), Vector2(1, 0), "flat",
                Renderer::TextAlignCenter | Renderer::TextAlignMiddle, 14, Rgb::rgb_white);
    r.draw_text(Vector2(6, 0), Vector2(1, 1), "rotated",
                Renderer::TextAlignLeft | Renderer::TextAlignTop, 10, Rgb::rgb_red);
    r.draw_text(Vector2(-6, 0), Vector2(0, 1), "vertical",
                Renderer::TextAlignRight | Renderer::TextAlignBottom, 12, Rgb::rgb_green);
    r.group_end();
    r.draw_frame_2d();
    std::string svg;
    r.write(svg);
    CHECK_SVG(svg, EXPECTED_PRIMITIVES);
}

TEST(render_projection_svg_matches_jvm) {
    RendererSvg r(320, 240, Rgb::rgb_white);
    r.set_page_layout(2, 2);
    r.set_window(Vector2Pair(Vector2(-5, -5), Vector2(5, 5)), true);
    r.set_camera_position(Vector3(0, 0, 10));
    r.set_camera_direction(Vector3(0, 0, -1));
    for (int i = 0; i < 4; i++) {
        r.set_page(i);
        r.draw_segment(Vector3(-4, -4, 0), Vector3(4, 4, 0), Rgb::rgb_red);
        r.draw_point(Vector3(0, 0, 0), Rgb::rgb_blue,
                     Renderer::PointStyle::PointStyleRound);
    }
    r.set_perspective();
    r.set_page(0);
    r.draw_segment(Vector3(-4, -4, 2), Vector3(4, 4, -2), Rgb::rgb_green);
    std::string svg;
    r.write(svg);
    CHECK_SVG(svg, EXPECTED_PROJECTION);
}
