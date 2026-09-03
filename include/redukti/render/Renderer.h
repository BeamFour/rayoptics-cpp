// C++ port of org.redukti.render.rendering: Rgb, Renderer, RendererViewport
// and Renderer2d.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RENDER_RENDERER_H
#define REDUKTI_RENDER_RENDERER_H

#include "redukti/mathlib/Transform3.h"
#include "redukti/mathlib/Triangle2.h"
#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector2Pair.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/mathlib/Vector3Pair.h"

#include <array>
#include <string>
#include <vector>

namespace redukti::render {

class Plot;

/** Java's `class Rgb`. */
class Rgb {
public:
    double r;
    double g;
    double b;
    double a;

    Rgb(double red, double green, double blue, double alpha)
        : r(red), g(green), b(blue), a(alpha) {}

    Rgb negate() const { return Rgb(1. - r, 1. - g, 1. - b, a); }

    static const Rgb rgb_black;
    static const Rgb rgb_red;
    static const Rgb rgb_green;
    static const Rgb rgb_blue;
    static const Rgb rgb_yellow;
    static const Rgb rgb_cyan;
    static const Rgb rgb_magenta;
    static const Rgb rgb_gray;
    static const Rgb rgb_white;
};

inline const Rgb Rgb::rgb_black{0.0f, 0.0f, 0.0f, 1.0f};
inline const Rgb Rgb::rgb_red{1.0f, 0.0f, 0.0f, 1.0f};
inline const Rgb Rgb::rgb_green{0.0f, 1.0f, 0.0f, 1.0f};
inline const Rgb Rgb::rgb_blue{0.0f, 0.0f, 1.0f, 1.0f};
inline const Rgb Rgb::rgb_yellow{1.0f, 1.0f, 0.0f, 1.0f};
inline const Rgb Rgb::rgb_cyan{0.0f, 1.0f, 1.0f, 1.0f};
inline const Rgb Rgb::rgb_magenta{1.0f, 0.0f, 1.0f, 1.0f};
inline const Rgb Rgb::rgb_gray{0.5f, 0.5f, 0.5f, 1.0f};
inline const Rgb Rgb::rgb_white{1.0f, 1.0f, 1.0f, 1.0f};

/** Java's abstract `Renderer`. */
class Renderer {
public:
    enum class IntensityMode {
        IntensityIgnore,
        IntensityShade,
        IntensityLogShade,
    };

    enum class RayColorMode {
        RayColorWavelen,
        RayColorFixed,
    };

    /** The int values are array indices into _styles_color, as in the Java. */
    enum class Style {
        StyleBackground = 0,
        StyleForeground = 1,
        StyleRay = 2,
        StyleSurface = 3,
        StyleGlass = 4,
        StyleLast = 5,
    };

    enum class PointStyle {
        PointStyleDot = 0,
        PointStyleCross = 1,
        PointStyleRound = 2,
        PointStyleSquare = 3,
        PointStyleTriangle = 4,
    };

    static PointStyle point_style_of(int i);

    /**
     * Java's `EnumSet<TextAlignMask>` becomes an int of these bits: the enum
     * values are already powers of two and the call sites only ever test
     * membership.
     */
    enum TextAlignMask {
        TextAlignCenter = 1,
        TextAlignLeft = 2,
        TextAlignRight = 4,
        TextAlignTop = 8,
        TextAlignBottom = 16,
        TextAlignMiddle = 32,
    };

    Renderer();
    virtual ~Renderer() = default;

    void set_stroke_dasharray(const std::string &dasharray) {
        _stroke_dasharray = dasharray;
        _has_stroke_dasharray = true;
    }

    /** Java passes null to clear it. */
    void clear_stroke_dasharray() {
        _stroke_dasharray.clear();
        _has_stroke_dasharray = false;
    }

    const std::string &get_stroke_dasharray() const { return _stroke_dasharray; }
    bool has_stroke_dasharray() const { return _has_stroke_dasharray; }

    void set_stroke_width(double width);
    double get_stroke_width() const { return _stroke_width; }

    Rgb get_style_color(Style s) const {
        return _styles_color[static_cast<std::size_t>(s)];
    }

    double get_feature_size() const { return _feature_size; }

    virtual void draw_point(const mathlib::Vector2 &p, const Rgb &rgb, PointStyle s) = 0;
    virtual void draw_text(const mathlib::Vector2 &pos, const mathlib::Vector2 &dir,
                           const std::string &str, int a, int size, const Rgb &rgb) = 0;
    virtual void draw_segment(const mathlib::Vector2Pair &s, const Rgb &rgb) = 0;
    virtual void draw_segment(const mathlib::Vector3Pair &s, const Rgb &rgb) = 0;

    void draw_point(const mathlib::Vector2 &p) {
        draw_point(p, Rgb::rgb_gray, PointStyle::PointStyleDot);
    }

    void draw_segment(const mathlib::Vector2Pair &s) { draw_segment(s, Rgb::rgb_gray); }

    void draw_segment(const mathlib::Vector2 &a, const mathlib::Vector2 &b,
                      const Rgb &rgb) {
        draw_segment(mathlib::Vector2Pair(a, b), rgb);
    }

    void draw_segment(const mathlib::Vector2 &a, const mathlib::Vector2 &b) {
        draw_segment(a, b, Rgb::rgb_gray);
    }

    void draw_segment(const mathlib::Vector3 &a, const mathlib::Vector3 &b,
                      const Rgb &rgb) {
        draw_segment(mathlib::Vector3Pair(a, b), rgb);
    }

    virtual void draw_polygon(const std::vector<mathlib::Vector2> &array, const Rgb &rgb,
                              bool filled, bool closed);

    virtual void draw_circle(const mathlib::Vector2 &v, double r, const Rgb &rgb,
                             bool filled);

    void draw_triangle(const mathlib::Triangle2 &t, bool filled, const Rgb &rgb);

    void draw_box(const mathlib::Vector2Pair &c, const Rgb &rgb);

    virtual void group_begin(const std::string &name) { (void)name; }
    virtual void group_end() {}

    virtual void draw_plot(const Plot &plot) { (void)plot; }

protected:
    double _feature_size;
    std::array<Rgb, static_cast<std::size_t>(Style::StyleLast)> _styles_color;
    RayColorMode _ray_color_mode;
    IntensityMode _intensity_mode;
    std::string _stroke_dasharray;
    bool _has_stroke_dasharray = false;
    double _stroke_width = 1.0;
};

/** Java's abstract `RendererViewport`. */
class RendererViewport : public Renderer {
public:
    enum class margin_type_e {
        MarginRatio,
        MarginLocal,
        MarginOutput,
    };

    RendererViewport();

    void set_2d_size(double width, double height) {
        _2d_output_res = mathlib::Vector2(width, height);
    }

    void set_window(const mathlib::Vector2 &center, const mathlib::Vector2 &size,
                    bool keep_aspect);
    void set_window(const mathlib::Vector2 &center, double radius, bool keep_aspect);
    void set_window(const mathlib::Vector2Pair &window, bool keep_aspect);

    virtual void update_2d_window() {}

    virtual void set_orthographic() = 0;
    virtual void set_perspective() = 0;

    void set_page(int page);

    double x_scale(double x) const;
    double y_scale(double y) const;
    double x_trans_pos(double x) const;
    virtual double y_trans_pos(double y) const;

    mathlib::Vector2Pair get_window2d_fit() const { return _window2d_fit; }
    mathlib::Vector2Pair get_window2d() const { return _window2d; }
    mathlib::Vector2 get_2d_output_res() const { return _2d_output_res; }

    void set_margin_output(double width, double height);
    void set_margin(double width, double height);
    void set_margin_ratio(double width, double height);
    void set_margin(double left, double bottom, double right, double top);
    void set_margin_ratio(double left, double bottom, double right, double top);
    void set_margin_output(double left, double bottom, double right, double top);

    void draw_frame_2d();

    void set_page_layout(int cols, int rows);

    void set_feature_size(double v) { _feature_size = v; }

    void set_camera_direction(const mathlib::Vector3 &dir);
    void set_camera_position(const mathlib::Vector3 &pos);

    virtual mathlib::Transform3 get_camera_transform() const = 0;
    virtual void set_camera_transform(const mathlib::Transform3 &t) = 0;

protected:
    mathlib::Vector2Pair _window2d_fit;
    mathlib::Vector2Pair _window2d;
    mathlib::Vector2 _2d_output_res;
    margin_type_e _margin_type;
    mathlib::Vector2Pair _margin;
    int _rows, _cols;
    int _pageid;
    mathlib::Vector2Pair _page;
    double _fov;
};

/** Java's abstract `Renderer2d`. */
class Renderer2d : public RendererViewport {
public:
    enum class ProjectionType {
        Ortho,
        Perspective,
    };

    void set_perspective() override;
    void set_orthographic() override { _projection_type = ProjectionType::Ortho; }

    mathlib::Vector2 project(const mathlib::Vector3 &v) const;
    mathlib::Vector2 project_scale(const mathlib::Vector3 &v) const;
    mathlib::Vector2 projection_ortho(const mathlib::Vector3 &v) const;
    mathlib::Vector2 projection_perspective(const mathlib::Vector3 &v) const;

    void draw_point(const mathlib::Vector3 &p, const Rgb &rgb, PointStyle s);
    void draw_segment(const mathlib::Vector3Pair &l, const Rgb &rgb) override;
    void draw_text(const mathlib::Vector3 &pos, const mathlib::Vector3 &dir,
                   const std::string &str, TextAlignMask a, int size, const Rgb &rgb);

    // The 2D overloads stay reachable alongside the 3D ones above.
    using Renderer::draw_point;
    using Renderer::draw_segment;
    using Renderer::draw_text;

    mathlib::Transform3 get_camera_transform() const override { return _cam_transform; }
    void set_camera_transform(const mathlib::Transform3 &t) override {
        _cam_transform = t;
    }

protected:
    ProjectionType _projection_type = ProjectionType::Ortho;
    mathlib::Transform3 _cam_transform;
    double _eye_dist = 0.0;
};

} // namespace redukti::render

#endif // REDUKTI_RENDER_RENDERER_H
