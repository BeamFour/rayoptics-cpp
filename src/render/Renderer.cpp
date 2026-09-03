// C++ port of org.redukti.render.rendering Renderer, RendererViewport and
// Renderer2d.
#include "redukti/render/Renderer.h"

#include "redukti/Exceptions.h"
#include "redukti/mathlib/M.h"

#include <cmath>

namespace redukti::render {

using mathlib::Transform3;
using mathlib::Vector2;
using mathlib::Vector2Pair;
using mathlib::Vector3;
using mathlib::Vector3Pair;

// ---------------------------------------------------------------------------
// Renderer
// ---------------------------------------------------------------------------

Renderer::PointStyle Renderer::point_style_of(int i) {
    switch (i) {
    case 0: return PointStyle::PointStyleDot;
    case 1: return PointStyle::PointStyleCross;
    case 2: return PointStyle::PointStyleRound;
    case 3: return PointStyle::PointStyleSquare;
    case 4: return PointStyle::PointStyleTriangle;
    default: throw IllegalArgumentException("PointStyle");
    }
}

Renderer::Renderer()
    : _feature_size(20.),
      _styles_color{Rgb(0.0, 0.0, 0.0, 1.0), Rgb(1.0, 1.0, 1.0, 1.0),
                    Rgb(1.0, 0.0, 0.0, 1.0), Rgb(0.5, 0.5, 1.0, 1.0),
                    Rgb(0.8, 0.8, 1.0, 1.0)},
      _ray_color_mode(RayColorMode::RayColorWavelen),
      _intensity_mode(IntensityMode::IntensityIgnore) {}

void Renderer::set_stroke_width(double width) {
    if (!(width > 0.0) || !std::isfinite(width))
        throw IllegalArgumentException("stroke width must be finite and positive");
    _stroke_width = width;
}

void Renderer::draw_polygon(const std::vector<Vector2> &array, const Rgb &rgb,
                            bool filled, bool closed) {
    (void)filled; // the base renderer draws the outline either way, as in Java
    std::size_t i;
    if (array.size() < 3)
        return;
    for (i = 0; i + 1 < array.size(); i++)
        draw_segment(Vector2Pair(array[i], array[i + 1]), rgb);
    if (closed)
        draw_segment(Vector2Pair(array[i], array[0]), rgb);
}

void Renderer::draw_circle(const Vector2 &v, double r, const Rgb &rgb, bool filled) {
    int count = std::min(100, std::max(6, static_cast<int>(2. * mathlib::M::PI * r /
                                                           _feature_size)));
    std::vector<Vector2> p(static_cast<std::size_t>(count), Vector2::vector2_0);
    double astep = 2. * mathlib::M::PI / count;
    double a = astep;
    // Java assigns p[0] = (r, 0) here and then overwrites it on the first
    // iteration below; kept so the geometry matches exactly.
    p[0] = Vector2(r, 0);
    for (int i = 0; i < count; i++, a += astep)
        p[static_cast<std::size_t>(i)] =
            v.plus(Vector2(r * std::cos(a), r * std::sin(a)));
    draw_polygon(p, rgb, filled, true);
}

void Renderer::draw_triangle(const mathlib::Triangle2 &t, bool filled, const Rgb &rgb) {
    auto arr = t.as_array();
    draw_polygon(std::vector<Vector2>(arr.begin(), arr.end()), rgb, filled, true);
}

void Renderer::draw_box(const Vector2Pair &c, const Rgb &rgb) {
    draw_segment(Vector2(c.v0.x, c.v0.y), Vector2(c.v1.x, c.v0.y), rgb);
    draw_segment(Vector2(c.v1.x, c.v1.y), Vector2(c.v1.x, c.v0.y), rgb);
    draw_segment(Vector2(c.v1.x, c.v1.y), Vector2(c.v0.x, c.v1.y), rgb);
    draw_segment(Vector2(c.v0.x, c.v0.y), Vector2(c.v0.x, c.v1.y), rgb);
}

// ---------------------------------------------------------------------------
// RendererViewport
// ---------------------------------------------------------------------------

RendererViewport::RendererViewport()
    : _window2d_fit(Vector2Pair::vector2_pair_00),
      _window2d(Vector2Pair::vector2_pair_00), _2d_output_res(Vector2::vector2_0),
      _margin_type(margin_type_e::MarginRatio),
      _margin(Vector2(0.13, 0.13), Vector2(0.13, 0.13)), _rows(1), _cols(1), _pageid(0),
      _page(Vector2Pair::vector2_pair_00), _fov(45.) {}

void RendererViewport::set_window(const Vector2 &center, const Vector2 &size,
                                  bool keep_aspect) {
    Vector2 s = size;
    if (keep_aspect) {
        double out_ratio = (_2d_output_res.x / _cols) / (_2d_output_res.y / _rows);
        if (std::abs(s.x / s.y) < out_ratio)
            s = Vector2(s.y * out_ratio, s.y);
        else
            s = Vector2(s.x, s.x / out_ratio);
    }
    Vector2 sby2 = s.divide(2.0);
    _window2d_fit = Vector2Pair(center.minus(sby2), center.plus(sby2));
    Vector2 ms0 = sby2;
    Vector2 ms1 = sby2;
    switch (_margin_type) {
    case margin_type_e::MarginLocal:
        ms0 = ms0.plus(_margin.v0);
        ms1 = ms1.plus(_margin.v1);
        break;
    case margin_type_e::MarginRatio:
        ms0 = ms0.plus(s.ebeTimes(_margin.v0));
        ms1 = ms1.plus(s.ebeTimes(_margin.v1));
        break;
    case margin_type_e::MarginOutput:
        ms0 = ms0.ebeDivide(
            Vector2::vector2_1.minus(_margin.v0.ebeDivide(_2d_output_res.times(2.0))));
        ms1 = ms1.ebeDivide(
            Vector2::vector2_1.minus(_margin.v1.ebeDivide(_2d_output_res.times(2.0))));
        break;
    }
    _window2d = Vector2Pair(center.minus(ms0), center.plus(ms1));
    update_2d_window();
    set_orthographic();
    set_page(_pageid);
}

void RendererViewport::set_window(const Vector2 &center, double radius,
                                  bool keep_aspect) {
    Vector2 size(radius, radius);
    set_window(center, size, keep_aspect);
}

void RendererViewport::set_window(const Vector2Pair &window, bool keep_aspect) {
    Vector2 center = window.v0.plus(window.v1).divide(2.0);
    Vector2 size(window.v1.x - window.v0.x, window.v1.y - window.v0.y);
    set_window(center, size, keep_aspect);
}

void RendererViewport::set_page(int page) {
    if (page >= _cols * _rows)
        throw IllegalArgumentException(
            "set_page: no such page number in current layout");
    _pageid = page;
    int row = page / _cols;
    int col = page % _cols;
    Vector2 size(_window2d.v1.x - _window2d.v0.x, _window2d.v1.y - _window2d.v0.y);
    Vector2 a(_window2d.v0.x - size.x * col,
              _window2d.v0.y - size.y * (_rows - 1 - row));
    Vector2 b(a.x + size.x * _cols, a.y + size.y * _rows);
    _page = Vector2Pair(a, b);
}

double RendererViewport::x_scale(double x) const {
    return ((x / (_page.v1.x - _page.v0.x)) * _2d_output_res.x);
}

double RendererViewport::y_scale(double y) const {
    return ((y / (_page.v1.y - _page.v0.y)) * _2d_output_res.y);
}

double RendererViewport::x_trans_pos(double x) const { return x_scale(x - _page.v0.x); }

double RendererViewport::y_trans_pos(double y) const { return y_scale(y - _page.v0.y); }

void RendererViewport::set_margin_output(double width, double height) {
    set_margin_output(width, height, width, height);
}

void RendererViewport::set_margin(double width, double height) {
    set_margin(width, height, width, height);
}

void RendererViewport::set_margin_ratio(double width, double height) {
    set_margin_ratio(width, height, width, height);
}

void RendererViewport::set_margin(double left, double bottom, double right, double top) {
    _margin_type = margin_type_e::MarginLocal;
    _margin = Vector2Pair(Vector2(left, bottom), Vector2(right, top));
    set_window(_window2d_fit, false);
}

void RendererViewport::set_margin_ratio(double left, double bottom, double right,
                                        double top) {
    _margin_type = margin_type_e::MarginRatio;
    _margin = Vector2Pair(Vector2(left, bottom), Vector2(right, top));
    set_window(_window2d_fit, false);
}

void RendererViewport::set_margin_output(double left, double bottom, double right,
                                         double top) {
    _margin_type = margin_type_e::MarginOutput;
    _margin = Vector2Pair(Vector2(left, bottom), Vector2(right, top));
    set_window(_window2d_fit, false);
}

void RendererViewport::draw_frame_2d() {
    std::vector<Vector2> fr(4, Vector2::vector2_0);
    fr[0] = _window2d_fit.v0;
    fr[1] = Vector2(_window2d_fit.v0.x, _window2d_fit.v1.y);
    fr[2] = _window2d_fit.v1;
    fr[3] = Vector2(_window2d_fit.v1.x, _window2d_fit.v0.y);
    draw_polygon(fr, get_style_color(Style::StyleForeground), false, true);
}

void RendererViewport::set_page_layout(int cols, int rows) {
    _cols = cols;
    _rows = rows;
    set_page(0);
}

void RendererViewport::set_camera_direction(const Vector3 &dir) {
    Transform3 t = get_camera_transform();
    t = t.set_direction(dir);
    set_camera_transform(t);
}

void RendererViewport::set_camera_position(const Vector3 &pos) {
    Transform3 t = get_camera_transform();
    t = t.set_translation(pos);
    set_camera_transform(t);
}

// ---------------------------------------------------------------------------
// Renderer2d
// ---------------------------------------------------------------------------

void Renderer2d::set_perspective() {
    double out_ratio = (_2d_output_res.y / _rows) / (_2d_output_res.x / _cols);
    if (out_ratio < 1.)
        _window2d = Vector2Pair(Vector2(-1. / out_ratio, -1.), Vector2(1. / out_ratio, 1.));
    else
        _window2d = Vector2Pair(Vector2(-1, -out_ratio), Vector2(1., out_ratio));
    _window2d_fit = _window2d;
    update_2d_window();
    set_page(_pageid);
    _projection_type = ProjectionType::Perspective;
    _eye_dist = 1. / std::tan(mathlib::M::toRadians(_fov) / 2.);
}

Vector2 Renderer2d::project(const Vector3 &v) const {
    switch (_projection_type) {
    case ProjectionType::Perspective:
        return projection_perspective(v);
    default:
        return projection_ortho(v);
    }
}

Vector2 Renderer2d::project_scale(const Vector3 &v) const {
    Vector2 v2d = project(v);
    return Vector2(x_trans_pos(v2d.x), y_trans_pos(v2d.y));
}

Vector2 Renderer2d::projection_ortho(const Vector3 &v) const {
    return _cam_transform.transform(v).project_xy();
}

Vector2 Renderer2d::projection_perspective(const Vector3 &v) const {
    Vector3 t = _cam_transform.transform(v);
    return Vector2(t.x * _eye_dist / -t.z, t.y * _eye_dist / -t.z);
}

void Renderer2d::draw_point(const Vector3 &p, const Rgb &rgb, PointStyle s) {
    draw_point(project(p), rgb, s);
}

void Renderer2d::draw_segment(const Vector3Pair &l, const Rgb &rgb) {
    draw_segment(Vector2Pair(project(l.point()), project(l.direction())), rgb);
}

void Renderer2d::draw_text(const Vector3 &pos, const Vector3 &dir, const std::string &str,
                           TextAlignMask a, int size, const Rgb &rgb) {
    draw_text(project(pos), project(dir), str, static_cast<int>(a), size, rgb);
}

} // namespace redukti::render
