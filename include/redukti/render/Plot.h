// C++ port of org.redukti.render.plotting: PlotStyleMask, PlotData, PlotAxes,
// Plot and PlotRenderer.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RENDER_PLOT_H
#define REDUKTI_RENDER_PLOT_H

#include "redukti/Text.h"
#include "redukti/data/DataSet.h"
#include "redukti/render/Renderer.h"

#include <string>
#include <vector>

namespace redukti::render {

/** Java's `enum PlotStyleMask`; the values are combined as a bit mask. */
enum PlotStyleMask {
    LinePlot = 1,
    PointPlot = 2,
    InterpolatePlot = 4,
    ValuePlot = 8,
    Filled = 16,
};

class PlotData {
public:
    /** Java's nested `enum LineStyle`; Solid carries a null dash array. */
    enum class LineStyle {
        Solid,
        Dashed,
        Dotted,
        DashDot,
    };

    static const char *dasharray(LineStyle s);

    /** Borrowed; the caller owns the data set and outlives the plot. */
    explicit PlotData(data::DataSet *s);

    data::DataSet *get_set() const { return _set; }

    void set_label(const std::string &title) { _label = title; }
    const std::string &get_label() const { return _label; }

    void set_color(const Rgb &color) { _color = color; }
    Rgb get_color() const { return _color; }

    void enable_style(PlotStyleMask style) { _style |= static_cast<int>(style); }
    void disable_style(PlotStyleMask style) { _style &= ~static_cast<int>(style); }
    void set_style(int style) { _style = style; }
    int get_style() const { return _style; }

    void set_line_style(LineStyle line_style) { _line_style = line_style; }
    LineStyle get_line_style() const { return _line_style; }

private:
    data::DataSet *_set;
    Rgb _color;
    int _style;
    std::string _label;
    LineStyle _line_style = LineStyle::Solid;
};

class PlotAxes {
public:
    enum class AxisMask {
        X = 1,
        Y = 2,
        Z = 4,
        XY = 3,
        YZ = 6,
        XZ = 5,
        XYZ = 7,
    };

    enum class step_mode_e {
        step_interval,
        step_count,
        step_base,
    };

    class Axis {
    public:
        bool _axis = true;
        bool _tics = true;
        bool _values = true;
        step_mode_e _step_mode = step_mode_e::step_base;
        int _count = 5;
        double _step_base = 10.0;
        bool _si_prefix = false;
        bool _pow10_scale = true;
        int _pow10 = 0;
        std::string _unit;
        std::string _label;
        data::Range _range{0, 0};
    };

    std::array<Axis, 2> _axes;
    bool _grid = false;
    bool _frame = true;
    mathlib::Vector3 _pos = mathlib::Vector3::vector3_0;
    mathlib::Vector3 _origin = mathlib::Vector3::vector3_0;

    void set_tics_step(double step, AxisMask a);
    void set_tics_step(double step) { set_tics_step(step, AxisMask::XYZ); }

    void set_tics_count(int count, AxisMask a);
    void set_tics_count(int count) { set_tics_count(count, AxisMask::XYZ); }

    void set_tics_base(int min_count, double base, AxisMask a);
    void set_tics_base() { set_tics_base(5, 10.0, AxisMask::XYZ); }

    void set_origin(const mathlib::Vector3 &origin) { _origin = origin; }
    mathlib::Vector3 get_origin() const { return _origin; }

    void set_position(const mathlib::Vector3 &position) { _pos = position; }
    mathlib::Vector3 get_position() const { return _pos; }

    void set_show_grid(bool show) { _grid = show; }
    bool get_show_grid() const { return _grid; }

    void set_show_frame(bool show) { _frame = show; }
    bool get_show_frame() const { return _frame; }

    void set_show_axes(bool show, AxisMask a);
    bool get_show_axes(int axis) const { return _axes[static_cast<std::size_t>(axis)]._axis; }

    void set_show_tics(bool show, AxisMask a);
    bool get_show_tics(int axis) const { return _axes[static_cast<std::size_t>(axis)]._tics; }

    void set_show_values(bool show, AxisMask a);
    bool get_show_values(int axis) const {
        return _axes[static_cast<std::size_t>(axis)]._values;
    }

    void set_label(const std::string &label, AxisMask a);
    const std::string &get_label(int axis) const {
        return _axes[static_cast<std::size_t>(axis)]._label;
    }

    void set_unit(const std::string &unit, bool pow10_scale, bool si_prefix, int pow10,
                  AxisMask a);

    void set_range(const data::Range &r, AxisMask a);

    double get_tics_step(int index, const data::Range &r) const;
};

class Plot {
public:
    /** Borrowed; the caller owns the data set. */
    PlotData *add_plot_data(data::DataSet *data, const Rgb &color,
                            const std::string &label, int style);

    void add_plot_data(const PlotData &data) { _plots.push_back(data); }

    void erase_plot_data() { _plots.clear(); }

    int get_plot_count() const { return static_cast<int>(_plots.size()); }

    PlotData &get_plot_data(int index) { return _plots[static_cast<std::size_t>(index)]; }
    const PlotData &get_plot_data(int index) const {
        return _plots[static_cast<std::size_t>(index)];
    }

    void set_title(const std::string &title) { _title = title; }
    const std::string &get_title() const { return _title; }

    void set_color(const Rgb &color);
    void set_different_colors();
    void set_style(int style);

    void set_xy_swap(bool doswap) { _xy_swap = doswap; }
    bool get_xy_swap() const { return _xy_swap; }

    void fit_axes_range();

    PlotAxes &get_axes() { return _axes; }
    const PlotAxes &get_axes() const { return _axes; }

    int get_dimensions() const;

    data::Range get_x_data_range(int dimension) const;
    data::Range get_y_data_range() const;

private:
    std::string _title;
    std::vector<PlotData> _plots;
    PlotAxes _axes;
    bool _xy_swap = false;
};

class PlotRenderer {
public:
    PlotRenderer() : _decimal_format(mathlib::M::decimal_format(2)) {}

    void draw_plot(RendererViewport &r, Plot &plot);

    void draw_axes_2d(RendererViewport &renderer, PlotAxes &a);

private:
    void draw_plot_data_2d(RendererViewport &r, data::Set1d &data, PlotData &style);
    void draw_polyline(RendererViewport &r, const std::vector<mathlib::Vector2> &pts,
                       const Rgb &color);
    void draw_frame_2d(RendererViewport &r);
    void set_2d_plot_window(RendererViewport &r, Plot &plot);
    void draw_axes_tic2(RendererViewport &r, PlotAxes &a, int i, int pow10, bool oor,
                        double x);

    DecimalFormat _decimal_format;
};

} // namespace redukti::render

#endif // REDUKTI_RENDER_PLOT_H
