// C++ port of org.redukti.render.rendering.RendererSvg
#include "redukti/render/RendererSvg.h"

#include "redukti/mathlib/M.h"

#include <cmath>
#include <cstdio>

namespace redukti::render {

using mathlib::Vector2;
using mathlib::Vector2Pair;

namespace {

/**
 * Java writes System.lineSeparator(), which is "\r\n" on Windows and "\n"
 * elsewhere. The port fixes it to "\n" so an SVG generated here is identical on
 * every platform; a test that compares against JVM output on Windows has to
 * normalise, which SvgTest does.
 */
constexpr const char *NL = "\n";

} // namespace

RendererSvg::RendererSvg(double width, double height, const Rgb &bg)
    : _decimal_format(mathlib::M::decimal_format()) {
    _2d_output_res = Vector2(width, height);
    _styles_color[static_cast<std::size_t>(Style::StyleBackground)] = bg;
    _styles_color[static_cast<std::size_t>(Style::StyleForeground)] = bg.negate();
    clear();
}

void RendererSvg::svg_begin_rect(double x1, double y1, double x2, double y2,
                                 bool terminate) {
    _out += "<rect ";
    _out += "x=\"" + format(x1) + "\" ";
    _out += "y=\"" + format(y1) + "\" ";
    _out += "width=\"" + format(x2 - x1) + "\" ";
    _out += "height=\"" + format(y2 - y1) + "\" ";
    if (terminate) {
        _out += " />";
        _out += NL;
    }
}

void RendererSvg::svg_begin_line(double x1, double y1, double x2, double y2,
                                 bool terminate) {
    _out += "<line ";
    _out += "x1=\"" + format(x1) + "\" ";
    _out += "y1=\"" + format(y1) + "\" ";
    _out += "x2=\"" + format(x2) + "\" ";
    _out += "y2=\"" + format(y2) + "\" ";
    if (terminate) {
        _out += " />";
        _out += NL;
    }
}

void RendererSvg::svg_begin_ellipse(double x, double y, double rx, double ry,
                                    bool terminate) {
    _out += "<ellipse ";
    _out += "cx=\"" + format(x) + "\" ";
    _out += "cy=\"" + format(y) + "\" ";
    _out += "rx=\"" + format(rx) + "\" ";
    _out += "ry=\"" + format(ry) + "\" ";
    if (terminate) {
        _out += " />";
        _out += NL;
    }
}

void RendererSvg::write_srgb(const Rgb &rgb) {
    char buf[32];
    // Java: String.format("#%02x%02x%02x", (int)(r*255.0), ...) -- the cast
    // truncates toward zero, so 1.0 gives ff and 0.999 gives fe.
    std::snprintf(buf, sizeof(buf), "#%02x%02x%02x", static_cast<int>(rgb.r * 255.0),
                  static_cast<int>(rgb.g * 255.0), static_cast<int>(rgb.b * 255.0));
    _out += buf;
}

void RendererSvg::svg_add_fill(const Rgb &rgb) {
    _out += " fill=\"";
    write_srgb(rgb);
    _out += "\"";
}

void RendererSvg::svg_end() {
    _out += " />";
    _out += NL;
}

void RendererSvg::clear() {
    _out.clear();
    svg_begin_rect(0.0, 0.0, _2d_output_res.x, _2d_output_res.y, false);
    svg_add_fill(get_style_color(Style::StyleBackground));
    svg_end();
    _out += "<defs>";
    _out += NL;
    _out += "<g id=\"dot\">";
    _out += NL;
    svg_begin_line(1, 1, 0, 0, true);
    _out += "</g>";
    _out += NL;
    _out += "<g id=\"cross\">";
    _out += NL;
    svg_begin_line(-3, 0, 3, 0, true);
    svg_begin_line(0, -3, 0, 3, true);
    _out += "</g>";
    _out += NL;
    _out += "<g id=\"square\">";
    _out += NL;
    svg_begin_line(-3, -3, -3, 3, true);
    svg_begin_line(-3, 3, 3, 3, true);
    svg_begin_line(3, 3, 3, -3, true);
    svg_begin_line(3, -3, -3, -3, true);
    _out += "</g>";
    _out += NL;
    _out += "<g id=\"round\">";
    _out += NL;
    svg_begin_ellipse(0, 0, 3, 3, false);
    _out += " fill=\"none\" />";
    _out += "</g>";
    _out += NL;
    _out += "<g id=\"triangle\">";
    _out += NL;
    svg_begin_line(0, -3, -3, 3, true);
    svg_begin_line(-3, 3, 3, 3, true);
    svg_begin_line(0, -3, +3, +3, true);
    _out += "</g>";
    _out += NL;
    _out += "</defs>";
    _out += NL;
}

void RendererSvg::group_begin(const std::string &name) {
    _out += "<g>";
    if (!name.empty())
        _out += "<title>" + name + "</title>";
    _out += NL;
}

void RendererSvg::group_end() {
    _out += "</g>";
    _out += NL;
}

void RendererSvg::svg_begin_use(const std::string &id, double x, double y,
                                bool terminate) {
    _out += "<use ";
    _out += "x=\"" + format(x) + "\" ";
    _out += "y=\"" + format(y) + "\" ";
    _out += "xlink:href=\"#" + id + "\" ";
    if (terminate) {
        _out += " />";
        _out += NL;
    }
}

void RendererSvg::svg_add_stroke(const Rgb &rgb) {
    _out += " stroke=\"";
    write_srgb(rgb);
    _out += "\"";
}

void RendererSvg::svg_add_dasharray() {
    if (_has_stroke_dasharray && !_stroke_dasharray.empty())
        _out += " stroke-dasharray=\"" + _stroke_dasharray + "\"";
}

void RendererSvg::svg_add_stroke_width() {
    if (_stroke_width != 1.0)
        _out += " stroke-width=\"" + format(_stroke_width) + "\"";
}

namespace {
const char *const kPointIds[] = {"dot", "cross", "round", "square", "triangle"};
constexpr std::size_t kPointIdCount = sizeof(kPointIds) / sizeof(kPointIds[0]);
} // namespace

void RendererSvg::draw_point(const Vector2 &p, const Rgb &rgb, PointStyle s) {
    if (static_cast<std::size_t>(s) >= kPointIdCount)
        s = PointStyle::PointStyleCross;
    Vector2 v2d = trans_pos(p);
    svg_begin_use(kPointIds[static_cast<std::size_t>(s)], v2d.x, v2d.y, false);
    svg_add_stroke(rgb);
    svg_end();
}

void RendererSvg::draw_segment(const Vector2Pair &l, const Rgb &rgb) {
    Vector2 v2da = trans_pos(l.v0);
    Vector2 v2db = trans_pos(l.v1);
    svg_begin_line(v2da.x, v2da.y, v2db.x, v2db.y, false);
    svg_add_stroke(rgb);
    svg_add_stroke_width();
    svg_add_dasharray();
    svg_end();
}

void RendererSvg::draw_circle(const Vector2 &c, double r, const Rgb &rgb, bool filled) {
    Vector2 v2d = trans_pos(c);
    svg_begin_ellipse(v2d.x, v2d.y, x_scale(r), y_scale(r), false);
    svg_add_stroke(rgb);
    if (filled)
        svg_add_fill(rgb);
    else
        _out += " fill=\"none\"";
    svg_end();
}

void RendererSvg::draw_text(const Vector2 &v, const Vector2 &dir, const std::string &str,
                            int a, int size, const Rgb &rgb) {
    int margin = size / 2;
    Vector2 v2d = trans_pos(v);
    double x = v2d.x;
    double y = v2d.y;
    double yo = y, xo = x;
    _out += "<text style=\"font-size:" + std::to_string(size) + ";";
    if (a & TextAlignLeft) {
        x += margin;
    } else if (a & TextAlignRight) {
        _out += "text-align:right;text-anchor:end;";
        x -= margin;
    } else {
        _out += "text-align:center;text-anchor:middle;";
    }
    if (a & TextAlignTop)
        y += size + margin;
    else if (a & TextAlignBottom)
        y -= margin;
    else
        y += size / 2.0;
    _out += "\" x=\"" + format(x) + "\" y=\"" + format(y) + "\"";
    double ra = mathlib::M::toDegrees(std::atan2(-dir.y, dir.x));
    if (ra != 0)
        _out += " transform=\"rotate(" + format(ra) + "," + format(xo) + "," +
                format(yo) + ")\"";
    svg_add_fill(rgb);
    _out += ">" + str + "</text>";
    _out += NL;
}

void RendererSvg::draw_polygon(const std::vector<Vector2> &array, const Rgb &rgb,
                               bool filled, bool closed) {
    if (array.size() < 3)
        return;
    closed = closed || filled;
    if (closed) {
        _out += "<polygon";
        if (filled) {
            svg_add_fill(rgb);
        } else {
            _out += " fill=\"none\"";
            svg_add_stroke(rgb);
        }
    } else {
        _out += "<polyline fill=\"none\"";
        svg_add_stroke(rgb);
    }
    if (!filled) {
        svg_add_stroke_width();
        svg_add_dasharray();
    }
    _out += " points=\"";
    for (std::size_t i = 0; i < array.size(); i++) {
        Vector2 v2d = trans_pos(array[i]);
        _out += format(v2d.x) + "," + format(v2d.y) + " ";
    }
    _out += "\" />";
    _out += NL;
}

Vector2 RendererSvg::trans_pos(const Vector2 &v) const {
    return Vector2(x_trans_pos(v.x), y_trans_pos(v.y));
}

double RendererSvg::y_trans_pos(double y) const {
    return (((y - _page.v1.y) / (_page.v0.y - _page.v1.y)) * _2d_output_res.y);
}

std::string &RendererSvg::write(std::string &s) const {
    s += "<?xml version=\"1.0\" standalone=\"no\"?>";
    s += NL;
    s += "<svg width=\"" + format(_2d_output_res.x) + "px\" height=\"" +
         format(_2d_output_res.y) + "px\" ";
    s += "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" ";
    s += "xmlns:xlink=\"http://www.w3.org/1999/xlink\">";
    s += NL;
    s += _out;
    s += "</svg>";
    s += NL;
    return s;
}

} // namespace redukti::render
