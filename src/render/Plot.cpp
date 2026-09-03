// C++ port of org.redukti.render.plotting
#include "redukti/render/Plot.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"

#include <cmath>
#include <cstdio>
#include <limits>

namespace redukti::render {

using data::Range;
using data::Set1d;
using mathlib::Vector2;
using mathlib::Vector2Pair;
using mathlib::Vector3;

// ---------------------------------------------------------------------------
// PlotData
// ---------------------------------------------------------------------------

const char *PlotData::dasharray(LineStyle s) {
    switch (s) {
    case LineStyle::Solid: return nullptr; // Java: null
    case LineStyle::Dashed: return "6,4";
    case LineStyle::Dotted: return "2,3";
    case LineStyle::DashDot: return "6,3,2,3";
    }
    return nullptr;
}

PlotData::PlotData(data::DataSet *s)
    : _set(s), _color(Rgb::rgb_red),
      _style(static_cast<int>(InterpolatePlot) | static_cast<int>(PointPlot)) {}

// ---------------------------------------------------------------------------
// PlotAxes
// ---------------------------------------------------------------------------

namespace {
/** Java's `_axes_bits`: only X and Y are backed by an Axis. */
const int kAxesBits[] = {static_cast<int>(PlotAxes::AxisMask::X),
                         static_cast<int>(PlotAxes::AxisMask::Y)};
constexpr std::size_t kAxesCount = sizeof(kAxesBits) / sizeof(kAxesBits[0]);
} // namespace

void PlotAxes::set_tics_step(double step, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._step_base = step;
            _axes[i]._step_mode = step_mode_e::step_interval;
        }
    }
}

void PlotAxes::set_tics_count(int count, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._count = count;
            _axes[i]._step_mode = step_mode_e::step_count;
        }
    }
}

void PlotAxes::set_tics_base(int min_count, double base, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._count = min_count;
            _axes[i]._step_base = base;
            _axes[i]._step_mode = step_mode_e::step_base;
        }
    }
}

void PlotAxes::set_show_axes(bool show, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0)
            _axes[i]._axis = show;
    }
}

void PlotAxes::set_show_tics(bool show, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._tics = show;
            _axes[i]._axis |= show;
        }
    }
}

void PlotAxes::set_show_values(bool show, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._values = show;
            _axes[i]._tics |= show;
            _axes[i]._axis |= show;
        }
    }
}

void PlotAxes::set_label(const std::string &label, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0)
            _axes[i]._label = label;
    }
}

void PlotAxes::set_unit(const std::string &unit, bool pow10_scale, bool si_prefix,
                        int pow10, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0) {
            _axes[i]._si_prefix = si_prefix;
            _axes[i]._unit = unit;
            _axes[i]._pow10_scale = pow10_scale;
            _axes[i]._pow10 = pow10;
        }
    }
}

void PlotAxes::set_range(const Range &r, AxisMask a) {
    for (std::size_t i = 0; i < kAxesCount; i++) {
        if ((static_cast<int>(a) & kAxesBits[i]) != 0)
            _axes[i]._range = r;
    }
}

double PlotAxes::get_tics_step(int index, const Range &r) const {
    const Axis &a = _axes[static_cast<std::size_t>(index)];
    double d = r.second - r.first;
    switch (a._step_mode) {
    case step_mode_e::step_interval:
        return d > 0 ? a._step_base : -a._step_base;
    case step_mode_e::step_count:
        return d / static_cast<double>(a._count);
    case step_mode_e::step_base: {
        if (d == 0.0)
            return 1;
        double da = std::abs(d);
        double p = std::floor(std::log(da) / std::log(a._step_base));
        double n = std::pow(a._step_base, p);
        int f = 1;
        while (static_cast<int>(da / n * f) < a._count) {
            if (static_cast<int>(da / n * f * 2) >= a._count) {
                f *= 2;
                break;
            } else if (static_cast<int>(da / n * f * 5) >= a._count) {
                f *= 5;
                break;
            } else {
                f *= 10;
            }
        }
        n /= f;
        return d > 0 ? n : -n;
    }
    }
    throw IllegalArgumentException("step mode");
}

// ---------------------------------------------------------------------------
// Plot
// ---------------------------------------------------------------------------

PlotData *Plot::add_plot_data(data::DataSet *data, const Rgb &color,
                              const std::string &label, int style) {
    PlotData plotdata(data);
    _plots.push_back(plotdata);
    PlotData &added = _plots.back();
    added.set_color(color);
    added.set_label(label);
    added.set_style(style);
    return &added;
}

void Plot::set_color(const Rgb &color) {
    for (auto &i : _plots)
        i.set_color(color);
}

void Plot::set_different_colors() {
    int n = 1;
    for (auto &i : _plots) {
        double r = static_cast<double>((n >> 0) & 0x01);
        double g = static_cast<double>((n >> 1) & 0x01);
        double b = static_cast<double>((n >> 2) & 0x01);
        i.set_color(Rgb(r, g, b, 1.0f));
        n++;
    }
}

void Plot::set_style(int style) {
    for (auto &i : _plots)
        i.set_style(style);
}

void Plot::fit_axes_range() {
    switch (get_dimensions()) {
    case 1:
        _axes.set_range(get_x_data_range(0), PlotAxes::AxisMask::X);
        _axes.set_range(get_y_data_range(), PlotAxes::AxisMask::Y);
        break;
    case 2:
        _axes.set_range(get_x_data_range(0), PlotAxes::AxisMask::X);
        _axes.set_range(get_x_data_range(1), PlotAxes::AxisMask::Y);
        _axes.set_range(get_y_data_range(), PlotAxes::AxisMask::Z);
        break;
    default:
        throw IllegalArgumentException("inconsistent dimensions of data sets in plot");
    }
}

int Plot::get_dimensions() const {
    int dimension = 0;
    for (const auto &i : _plots) {
        int d = i.get_set()->get_dimensions();
        if (dimension == 0)
            dimension = d;
        else if (dimension != d)
            return 0;
    }
    return dimension;
}

Range Plot::get_x_data_range(int dimension) const {
    // Double.MIN_VALUE is the smallest positive value, so the upper seed is
    // 4.9e-324; see the note on DataSet::get_y_range.
    Range r(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::denorm_min());
    for (const auto &i : _plots) {
        Range ri = i.get_set()->get_x_range(dimension);
        if (ri.first < r.first)
            r.first = ri.first;
        if (ri.second > r.second)
            r.second = ri.second;
    }
    return r;
}

Range Plot::get_y_data_range() const {
    Range r(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::denorm_min());
    for (const auto &i : _plots) {
        Range ri = i.get_set()->get_y_range();
        if (ri.first < r.first)
            r.first = ri.first;
        if (ri.second > r.second)
            r.second = ri.second;
    }
    return r;
}

// ---------------------------------------------------------------------------
// PlotRenderer
// ---------------------------------------------------------------------------

void PlotRenderer::draw_plot(RendererViewport &r, Plot &plot) {
    switch (plot.get_dimensions()) {
    case 1: {
        set_2d_plot_window(r, plot);
        draw_axes_2d(r, plot.get_axes());
        Vector2Pair _window2d = r.get_window2d();
        Vector2Pair _window2d_fit = r.get_window2d_fit();
        r.draw_text(Vector2((_window2d.v0.x + _window2d.v1.x) / 2.,
                            (_window2d_fit.v1.y + _window2d.v1.y) / 2.),
                    Vector2::vector2_10, plot.get_title(),
                    Renderer::TextAlignCenter | Renderer::TextAlignMiddle, 18,
                    r.get_style_color(Renderer::Style::StyleForeground));
        for (int i = 0; i < plot.get_plot_count(); i++) {
            PlotData &d = plot.get_plot_data(i);
            draw_plot_data_2d(r, *static_cast<Set1d *>(d.get_set()), d);
        }
        break;
    }
    default:
        throw IllegalArgumentException("Unsupported dimensions " +
                                       std::to_string(plot.get_dimensions()));
    }
}

void PlotRenderer::draw_plot_data_2d(RendererViewport &r, Set1d &data, PlotData &style) {
    Vector2Pair _window2d_fit = r.get_window2d_fit();
    Vector2Pair _window2d = r.get_window2d();
    Vector2 _2d_output_res = r.get_2d_output_res();
    const char *dash = PlotData::dasharray(style.get_line_style());
    if (dash != nullptr)
        r.set_stroke_dasharray(dash);
    else
        r.clear_stroke_dasharray();
    if ((style.get_style() & static_cast<int>(InterpolatePlot)) != 0) {
        const double x_step = (_window2d.v1.x - _window2d.v0.x) / _2d_output_res.x;
        Range xr = data.get_x_range(0);
        double x_low = std::fmax(_window2d_fit.v0.x, xr.first);
        double x_high = std::fmin(_window2d_fit.v1.x, xr.second);
        std::vector<Vector2> pts;
        pts.push_back(Vector2(x_low, data.interpolate(x_low)));
        for (double x = x_low + x_step; x < x_high + x_step / 2; x += x_step)
            pts.push_back(Vector2(x, data.interpolate(x)));
        draw_polyline(r, pts, style.get_color());
    }
    if ((style.get_style() & static_cast<int>(LinePlot)) != 0) {
        std::vector<Vector2> pts;
        for (int j = 0; j < data.get_count(); j++)
            pts.push_back(Vector2(data.get_x_value(j), data.get_y_value(j)));
        draw_polyline(r, pts, style.get_color());
    }
    r.clear_stroke_dasharray();
    if ((style.get_style() & static_cast<int>(PointPlot)) != 0) {
        for (int j = 0; j < data.get_count(); j++) {
            Vector2 p(data.get_x_value(j), data.get_y_value(j));
            r.draw_point(p, style.get_color(), Renderer::PointStyle::PointStyleCross);
        }
    }
    if ((style.get_style() & static_cast<int>(ValuePlot)) != 0) {
        for (int j = 0; j < data.get_count(); j++) {
            int a;
            Range p(data.get_x_value(j), data.get_y_value(j));
            double prev = j > 0 ? data.get_y_value(j - 1) : p.second;
            double next = j + 1 < data.get_count() ? data.get_y_value(j + 1) : p.second;
            if (p.second > prev) {
                if (p.second > next)
                    a = Renderer::TextAlignBottom | Renderer::TextAlignCenter;
                else
                    a = Renderer::TextAlignBottom | Renderer::TextAlignRight;
            } else {
                if (p.second > next)
                    a = Renderer::TextAlignTop | Renderer::TextAlignRight;
                else
                    a = Renderer::TextAlignBottom | Renderer::TextAlignLeft;
            }
            // Java: String.format("%.02f", ...), which is plain fixed-point.
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%.02f", p.second);
            r.draw_text(Vector2(p.first, p.second), Vector2::vector2_10, buf, a, 12,
                        style.get_color());
        }
    }
}

void PlotRenderer::draw_polyline(RendererViewport &r, const std::vector<Vector2> &pts,
                                 const Rgb &color) {
    if (pts.size() < 2)
        return;
    if (pts.size() >= 3)
        r.draw_polygon(pts, color, false, false);
    else
        r.draw_segment(Vector2Pair(pts[0], pts[1]), color);
}

void PlotRenderer::draw_frame_2d(RendererViewport &r) {
    std::vector<Vector2> fr(4, Vector2::vector2_0);
    Vector2Pair _window2d_fit = r.get_window2d_fit();
    fr[0] = _window2d_fit.v0;
    fr[1] = Vector2(_window2d_fit.v0.x, _window2d_fit.v1.y);
    fr[2] = _window2d_fit.v1;
    fr[3] = Vector2(_window2d_fit.v1.x, _window2d_fit.v0.y);
    r.draw_polygon(fr, r.get_style_color(Renderer::Style::StyleForeground), false, true);
}

void PlotRenderer::set_2d_plot_window(RendererViewport &r, Plot &plot) {
    Range x_range = plot.get_axes()._axes[0]._range;
    if (x_range.first == x_range.second)
        x_range = plot.get_x_data_range(0);
    Range y_range = plot.get_axes()._axes[1]._range;
    if (y_range.first == y_range.second)
        y_range = plot.get_y_data_range();
    r.set_window(Vector2Pair(Vector2(x_range.first, y_range.first),
                             Vector2(x_range.second, y_range.second)),
                 false);
}

namespace {
const char *const kSiPrefixes[] = {"y", "z", "a", "f", "p", "n", "u", "m", "",
                                   "k", "M", "G", "T", "P", "E", "Z", "Y"};
} // namespace

void PlotRenderer::draw_axes_2d(RendererViewport &renderer, PlotAxes &a) {
    const int N = 2;
    Vector2 p(a.get_position().x, a.get_position().y);
    int pow10;
    int max[N];
    int min[N];
    double step[N];
    Vector2Pair _window2d = renderer.get_window2d();
    Vector2Pair _window2d_fit = renderer.get_window2d_fit();
    if (a._frame)
        draw_frame_2d(renderer);
    for (int i = 0; i < N; i++) {
        PlotAxes::Axis &ax = a._axes[static_cast<std::size_t>(i)];
        Range r(_window2d_fit.v0.v(i), _window2d_fit.v1.v(i));
        double s = step[i] = std::abs(a.get_tics_step(i, r));
        min[i] = mathlib::M::trunc((r.first - p.v(i)) / s);
        max[i] = mathlib::M::trunc((r.second - p.v(i)) / s);
        pow10 = ax._pow10_scale ? static_cast<int>(std::floor(std::log10(s))) : 0;
        std::string si_unit;
        if (ax._si_prefix) {
            int u = (24 + pow10 + ax._pow10) / 3;
            if (u >= 0 && u < 17) {
                si_unit = std::string(kSiPrefixes[u]) + ax._unit;
                pow10 = (u - 8) * 3 - ax._pow10;
            }
        }
        Vector2 lp = Vector2::vector2_0;
        Vector2 ld = Vector2::vector2_0;
        switch (i) {
        case 0:
            lp = Vector2((_window2d.v0.x + _window2d.v1.x) / 2.,
                         (_window2d_fit.v0.y * .50 + _window2d.v0.y * 1.50) / 2.);
            ld = Vector2::vector2_10;
            break;
        case 1:
            lp = Vector2((_window2d_fit.v0.x * .50 + _window2d.v0.x * 1.50) / 2.,
                         (_window2d.v0.y + _window2d.v1.y) / 2.);
            ld = Vector2::vector2_01;
            break;
        default:
            throw IllegalArgumentException("Invalid axis " + std::to_string(i));
        }
        {
            std::string lx = ax._label;
            bool useunit = !ax._unit.empty();
            bool usep10 = pow10 != 0;
            if (!si_unit.empty()) {
                lx += " (" + si_unit + ")";
            } else if (useunit || usep10) {
                lx += " (";
                if (usep10) {
                    char buf[32];
                    std::snprintf(buf, sizeof(buf), "x10^%d", pow10);
                    lx += buf;
                }
                if (useunit && usep10)
                    lx += " ";
                if (useunit)
                    lx += ax._unit;
                lx += ")";
            }
            renderer.draw_text(lp, ld, lx,
                               Renderer::TextAlignCenter | Renderer::TextAlignMiddle, 12,
                               renderer.get_style_color(Renderer::Style::StyleForeground));
        }
        bool oor = false;
        for (int j = 0; j < N; j++)
            oor |= (j != i && ((p.v(j) <= std::fmin(_window2d_fit.v0.v(j),
                                                    _window2d_fit.v1.v(j))) ||
                               (p.v(j) >= std::fmax(_window2d_fit.v0.v(j),
                                                    _window2d_fit.v1.v(j)))));
        if (!oor && ax._axis) {
            Vector2Pair seg(p.set(i, r.first), p.set(i, r.second));
            renderer.draw_segment(seg,
                                  renderer.get_style_color(Renderer::Style::StyleForeground));
        }
        if (ax._tics && (ax._axis || a._frame)) {
            for (int j = min[i]; j <= max[i]; j++)
                draw_axes_tic2(renderer, a, i, pow10, oor, j * s);
        }
    }
    if (a._grid) {
        for (int x = min[0]; x <= max[0]; x++)
            for (int y = min[1]; y <= max[1]; y++)
                renderer.draw_point(
                    Vector2(p.v(0) + x * step[0], p.v(1) + y * step[1]),
                    renderer.get_style_color(Renderer::Style::StyleForeground),
                    Renderer::PointStyle::PointStyleDot);
    }
}

void PlotRenderer::draw_axes_tic2(RendererViewport &r, PlotAxes &a, int i, int pow10,
                                  bool oor, double x) {
    Vector2 p(a.get_position().x, a.get_position().y);
    PlotAxes::Axis &ax = a._axes[static_cast<std::size_t>(i)];
    Vector2 vtic = Vector2::vector2_0;
    Vector2Pair _window2d_fit = r.get_window2d_fit();
    if (!oor && ax._axis) {
        vtic = p;
        vtic = vtic.set(i, x + p.v(i));
        r.draw_point(vtic, r.get_style_color(Renderer::Style::StyleForeground),
                     Renderer::PointStyle::PointStyleCross);
    }
    if (a._frame) {
        vtic = _window2d_fit.v1;
        vtic = vtic.set(i, x + p.v(i));
        r.draw_point(vtic, r.get_style_color(Renderer::Style::StyleForeground),
                     Renderer::PointStyle::PointStyleCross);
        vtic = _window2d_fit.v0;
        vtic = vtic.set(i, x + p.v(i));
        r.draw_point(vtic, r.get_style_color(Renderer::Style::StyleForeground),
                     Renderer::PointStyle::PointStyleCross);
    }
    if (ax._values) {
        // Java builds three alignments but only ever indexes 0 and 1 here.
        const int alignments[] = {
            Renderer::TextAlignCenter | Renderer::TextAlignTop,
            Renderer::TextAlignRight | Renderer::TextAlignMiddle,
            Renderer::TextAlignTop | Renderer::TextAlignCenter,
        };
        std::string s =
            _decimal_format.format((x + p.v(i) - a._origin.v(i)) / std::pow(10., pow10));
        int align = alignments[static_cast<std::size_t>(i)];
        r.draw_text(vtic, Vector2::vector2_10, s, align, 12,
                    r.get_style_color(Renderer::Style::StyleForeground));
    }
}

} // namespace redukti::render
