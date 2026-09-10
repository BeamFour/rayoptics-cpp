// The software is ported from Goptical, hence is licensed under the GPL.
// Copyright (c) 2021 Dibyendu Majumdar
// Goptical: Copyright (C) 2010-2011 Free Software Foundation, Inc; Author: Alexandre Becoulet
// Licensed under the GNU General Public License, version 3 or later; see LICENSE-GPL-3.0.txt
//
// C++ port of org.redukti.render.plotting: PlotStyleMask, PlotData, PlotAxes,
// Plot and PlotRenderer.
#ifndef REDUKTI_RENDER_PLOT_H
#define REDUKTI_RENDER_PLOT_H

#include "redukti/Text.h"
#include "redukti/data/DataSet.h"
#include "redukti/render/Renderer.h"

#include <string>
#include <vector>

namespace redukti::render {

/**
 * Specifies available styles for 2d and 3d plot data rendering
 */
/** Java's `enum PlotStyleMask`; the values are combined as a bit mask. */
enum PlotStyleMask {
    /**
     * Draw lines between knots
     */
    LinePlot = 1,
    /**
     * Draw points for each knot
     */
    PointPlot = 2,
    /**
     * Draw a smooth curve
     */
    InterpolatePlot = 4,
    /**
     * Print ploted values near knots
     */
    ValuePlot = 8,
    /**
     * Draw filled 3d surface
     */
    Filled = 16,
};

class PlotData {
public:
    /**
     * Line pattern used when drawing curves/lines for this data set.
     * The associated value is an SVG stroke-dasharray string (null = solid).
     */
    /** Java's nested `enum LineStyle`; Solid carries a null dash array. */
    enum class LineStyle {
        Solid,
        Dashed,
        Dotted,
        DashDot,
    };

    static const char *dasharray(LineStyle s);

    /**
     * Create a new data plot descriptor which describe the
     * specified dataset.
     */
    /**
     * Get the described data set
     */
    /** Borrowed; the caller owns the data set and outlives the plot. */
    explicit PlotData(data::DataSet *s);

    data::DataSet *get_set() const { return _set; }

    /**
     * Set data set plotting label
     */
    void set_label(const std::string &title) { _label = title; }
    /**
     * Get data set plotting label
     */
    const std::string &get_label() const { return _label; }

    /**
     * Set data set plotting color
     */
    void set_color(const Rgb &color) { _color = color; }
    /**
     * Set data set plotting color
     */
    Rgb get_color() const { return _color; }

    /**
     * Enable a plotting style
     */
    void enable_style(PlotStyleMask style) { _style |= static_cast<int>(style); }
    /**
     * Disable a plotting style
     */
    void disable_style(PlotStyleMask style) { _style &= ~static_cast<int>(style); }
    /**
     * Set the plotting style mask
     */
    void set_style(int style) { _style = style; }
    /**
     * Get the plotting style mask
     */
    int get_style() const { return _style; }

    /**
     * Set the line pattern (solid, dashed, ...) used for this data set.
     */
    void set_line_style(LineStyle line_style) { _line_style = line_style; }
    /**
     * Get the line pattern used for this data set.
     */
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
    /**
     * Specify axes
     */
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

    /**
     * This sets distance between axis tics to specified value.
     * see set_tics_count, set_tics_base
     */
    void set_tics_step(double step, AxisMask a);
    void set_tics_step(double step) { set_tics_step(step, AxisMask::XYZ); }

    /**
     * @This sets tics count. @see {set_tics_step, set_tics_base}
     */
    void set_tics_count(int count, AxisMask a);
    void set_tics_count(int count) { set_tics_count(count, AxisMask::XYZ); }

    /**
     * @This sets distance between axis tics to best fit power of
     * specified base divided by sufficient factor of 2 and 5 to
     * have at least @tt min_count tics. @see {set_tics_step,
     * set_tics_count}
     */
    void set_tics_base(int min_count, double base, AxisMask a);
    void set_tics_base() { set_tics_base(5, 10.0, AxisMask::XYZ); }

    /**
     * This sets axis tics values origin.
     */
    void set_origin(const mathlib::Vector3 &origin) { _origin = origin; }
    /**
     * This returns axes tics values origin.
     */
    mathlib::Vector3 get_origin() const { return _origin; }

    /**
     * This returns axis position
     */
    void set_position(const mathlib::Vector3 &position) { _pos = position; }
    /**
     * This returns axis position
     */
    mathlib::Vector3 get_position() const { return _pos; }

    /**
     * This sets grid visibility. Grid points use tic
     * step.
     */
    void set_show_grid(bool show) { _grid = show; }
    /**
     * see set_show_grid
     */
    bool get_show_grid() const { return _grid; }

    /**
     * @This sets frame visibility.
     */
    void set_show_frame(bool show) { _frame = show; }
    /**
     * see set_show_frame
     */
    bool get_show_frame() const { return _frame; }

    /**
     * @This sets axes visibility.
     */
    void set_show_axes(bool show, AxisMask a);
    /**
     * see set_show_axes
     */
    bool get_show_axes(int axis) const { return _axes[static_cast<std::size_t>(axis)]._axis; }

    /**
     * This sets tics visibility. Tics are located on axes and
     * frame. see {set_show_axes, set_show_frame}
     */
    void set_show_tics(bool show, AxisMask a);
    /**
     * see set_show_tics
     */
    bool get_show_tics(int axis) const { return _axes[static_cast<std::size_t>(axis)]._tics; }

    /**
     * @This sets tics value visibility. When frame is visible,
     * tics value is located on frame tics instead of axes tics.
     * @see {set_show_axes, set_show_frame}
     */
    void set_show_values(bool show, AxisMask a);
    /**
     * see set_show_values
     */
    bool get_show_values(int axis) const {
        return _axes[static_cast<std::size_t>(axis)]._values;
    }

    /**
     * This set axis label
     */
    void set_label(const std::string &label, AxisMask a);
    /**
     * Get axis label
     */
    const std::string &get_label(int axis) const {
        return _axes[static_cast<std::size_t>(axis)]._label;
    }

    /**
     * This sets axis unit.
     *
     * When @tt pow10_scale is set, value will be scaled to shorten
     * their length and appropriate power of 10 factor will be
     * displayed in axis label.
     *
     * If @tt si_prefix is set, SI letter decimal prefix is used
     * and the @tt pow10 parameter can be used to scale base unit
     * by power of ten (useful when input data use scaled SI base unit).
     */
    void set_unit(const std::string &unit, bool pow10_scale, bool si_prefix, int pow10,
                  AxisMask a);

    /**
     * Set value range for given axis. Default range is [0,0] which
     * means automatic range.
     */
    void set_range(const data::Range &r, AxisMask a);

    /**
     * get distance between axis tics
     */
    double get_tics_step(int index, const data::Range &r) const;
};

/**
 * data plots container
 *
 * This class is used to describe a data plot. It contains a list
 * of PlotData objects and describes some plot properties
 * (title, range, ...).
 *
 * Plots can be built from data sets or obtained directly from
 * various analysis functions.
 */
class Plot {
public:
    /**
     * Create and add plot data from specified data set.
     */
    /** Borrowed; the caller owns the data set. */
    PlotData *add_plot_data(data::DataSet *data, const Rgb &color,
                            const std::string &label, int style);

    /**
     * Add plot data
     */
    void add_plot_data(const PlotData &data) { _plots.push_back(data); }

    /**
     * Discard all plot data set
     */
    void erase_plot_data() { _plots.clear(); }

    /**
     * Get plot data set count
     */
    int get_plot_count() const { return static_cast<int>(_plots.size()); }

    PlotData &get_plot_data(int index) { return _plots[static_cast<std::size_t>(index)]; }
    /**
     * Get plot data set at given index
     */
    const PlotData &get_plot_data(int index) const {
        return _plots[static_cast<std::size_t>(index)];
    }

    /**
     * Set plot main title
     */
    void set_title(const std::string &title) { _title = title; }
    /**
     * Get plot main title
     */
    const std::string &get_title() const { return _title; }

    /**
     * Set color for all plots
     */
    void set_color(const Rgb &color);
    void set_different_colors();
    /**
     * Set plot style for all plot
     */
    void set_style(int style);

    /**
     * Swap x and y axis for 2d plots
     */
    void set_xy_swap(bool doswap) { _xy_swap = doswap; }
    /**
     * Get x and y axis swap state for 2d plots
     */
    bool get_xy_swap() const { return _xy_swap; }

    void fit_axes_range();

    /**
     * Get plot axes object
     */
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
    // Up to 2 fraction digits (trailing zeros suppressed) so fractional
    // axes such as 0, 0.2, 0.4 ... render correctly while integer axes
    // (0, 20, 40 ...) still print without a decimal point.
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
