// C++ port of org.redukti.plotter
#include "redukti/plotter/Plotter.h"

#include "redukti/data/DataSet.h"
#include "redukti/mathlib/M.h"
#include "redukti/rayoptics/util/Orientation.h"
#include "redukti/render/Plot.h"
#include "redukti/render/RendererSvg.h"

#include <cmath>
#include <memory>

namespace redukti::plotter {

using data::DiscreteSet;
using data::Interpolation;
using data::Range;
using mathlib::Vector2;
using mathlib::Vector2Pair;
using mathlib::Vector3;
using rayoptics::analysis::MTFResultByFreq;
using rayoptics::raytr::RayFanType;
using render::Plot;
using render::PlotAxes;
using render::PlotData;
using render::PlotRenderer;
using render::PlotStyleMask;
using render::Renderer;
using render::RendererSvg;
using render::Rgb;
namespace Orientation = rayoptics::util::Orientation;

// ---------------------------------------------------------------------------
// Colors
// ---------------------------------------------------------------------------

Rgb Colors::get_wavelen_color(double wl) {
    if (wl < 380.0 || wl > 780.0)
        return Rgb::rgb_black;
    double s = 1.0;
    if (wl < 420.0)
        s = 0.3 + 0.7 * (wl - 380.0f) / 40.0;
    else if (wl > 700.0)
        s = 0.3 + 0.7 * (780.0 - wl) / 80.0;
    if (wl < 510.0) {
        if (wl < 490.0) {
            if (wl < 440.0)
                return Rgb(s * -(wl - 440.0) / 60.0, 0.0, s, 1.0);
            else
                return Rgb(0.0, s * (wl - 440.0) / 50.0, s, 1.0);
        } else {
            return Rgb(0.0, s, s * -(wl - 510.0) / 20.0, 1.0);
        }
    } else {
        if (wl < 645.0) {
            if (wl < 580.0)
                return Rgb(s * (wl - 510.0) / 70.0, s, 0.0, 1.0);
            else
                return Rgb(s, s * -(wl - 645.0) / 65.0, 0.0, 1.0);
        } else {
            return Rgb(s, 0.0, 0.0, 1.0);
        }
    }
}

// ---------------------------------------------------------------------------
// SpotDiagram
// ---------------------------------------------------------------------------

std::string SpotDiagram::plot(std::optional<double> radius_in) const {
    RendererSvg r(640, 640, Rgb::rgb_black);
    double radius = radius_in.has_value() ? *radius_in : result->max_radius * 1000;
    r.set_window(Vector2Pair(Vector2(-radius, -radius), Vector2(radius, radius)), true);
    PlotAxes axes;
    axes.set_show_axes(false, PlotAxes::AxisMask::XY);
    axes.set_label("Sagittal distance", PlotAxes::AxisMask::X);
    axes.set_label("Tangential distance", PlotAxes::AxisMask::Y);
    axes.set_unit("m", true, true, -3, PlotAxes::AxisMask::XY);
    axes.set_tics_count(3, PlotAxes::AxisMask::XY);
    PlotRenderer plotRenderer;
    plotRenderer.draw_axes_2d(r, axes);
    for (const auto &intercepts : result->intercepts) {
        for (std::size_t i = 0; i < intercepts.x.size(); i++) {
            r.draw_point(Vector2(intercepts.x[i] * 1000, intercepts.y[i] * 1000),
                         Colors::get_wavelen_color(intercepts.wvl),
                         Renderer::PointStyle::PointStyleDot);
        }
    }
    std::string out;
    r.write(out);
    return out;
}

// ---------------------------------------------------------------------------
// GeoMTFPlot
// ---------------------------------------------------------------------------

std::string GeoMTFPlot::plot() const {
    int count = 0;
    for (std::size_t i = 0; i < geo_mtf->mtf.freq.size(); i++) {
        if (geo_mtf->mtf.freq[i] > 100.)
            break;
        count++;
    }
    Plot plot;
    plot.set_title("MTF for " + fld->toString() + " wvl " + doubleToString(geo_mtf->wvl));
    plot.get_axes().set_position(Vector3::vector3_0);
    plot.get_axes().set_range(Range(0, 100.0), PlotAxes::AxisMask::X);
    plot.get_axes().set_range(Range(0, 1.0), PlotAxes::AxisMask::Y);
    // Plot borrows its data sets, so they have to outlive the draw below.
    // deque-like stability is not needed, only that the addresses hold: the
    // vector is sized up front.
    std::vector<std::unique_ptr<DiscreteSet>> sets;
    for (int xy = 0; xy < Orientation::COUNT; xy++) {
        const auto &x_data = geo_mtf->mtf.freq;
        const auto &y_data =
            xy == Orientation::SAGITTAL ? geo_mtf->mtf.mag_x : geo_mtf->mtf.mag_y;
        auto set = std::make_unique<DiscreteSet>();
        set->set_interpolation(Interpolation::Linear);
        for (int i = 0; i < count; i++)
            set->add_data(x_data[static_cast<std::size_t>(i)],
                          y_data[static_cast<std::size_t>(i)]);
        plot.add_plot_data(set.get(),
                           xy == Orientation::SAGITTAL ? Rgb::rgb_black : Rgb::rgb_blue,
                           xy == Orientation::SAGITTAL ? "Sagittal" : "Tangential",
                           static_cast<int>(PlotStyleMask::InterpolatePlot));
        sets.push_back(std::move(set));
    }
    plot.get_axes().set_label("Spatial Frequency", PlotAxes::AxisMask::X);
    plot.get_axes().set_label("Modulation", PlotAxes::AxisMask::Y);
    plot.get_axes().set_unit("cycles/mm", false, false, 0, PlotAxes::AxisMask::X);
    RendererSvg r(640, 640);
    PlotRenderer plotRenderer;
    plotRenderer.draw_plot(r, plot);
    std::string out;
    r.write(out);
    return out;
}

// ---------------------------------------------------------------------------
// GeoMTFByFieldPlot
// ---------------------------------------------------------------------------

namespace {

const Rgb FREQ_COLORS[] = {Rgb::rgb_red,     Rgb::rgb_blue, Rgb::rgb_black,
                           Rgb::rgb_magenta, Rgb::rgb_cyan, Rgb::rgb_green};
constexpr std::size_t FREQ_COLOR_COUNT = sizeof(FREQ_COLORS) / sizeof(FREQ_COLORS[0]);

const char *const FREQ_COLOR_NAMES[] = {"red", "blue", "black", "magenta", "cyan",
                                        "green"};
constexpr std::size_t FREQ_COLOR_NAME_COUNT =
    sizeof(FREQ_COLOR_NAMES) / sizeof(FREQ_COLOR_NAMES[0]);

/** The Java holds one static DecimalFormat, M.decimal_format(). */
const DecimalFormat &df() {
    static const DecimalFormat instance = mathlib::M::decimal_format();
    return instance;
}

} // namespace

std::string GeoMTFByFieldPlot::freq_legend(const std::vector<int> &freqs) {
    std::string sb;
    for (std::size_t i = 0; i < freqs.size(); i++) {
        if (i > 0)
            sb += ",";
        sb += std::to_string(freqs[i]);
        sb += "=";
        sb += FREQ_COLOR_NAMES[i % FREQ_COLOR_NAME_COUNT];
    }
    return sb;
}

std::string GeoMTFByFieldPlot::plot() const {
    Plot plot;
    plot.set_title("MTF");
    plot.get_axes().set_position(Vector3::vector3_0);
    plot.get_axes().set_range(Range(0, 1.0), PlotAxes::AxisMask::X);
    plot.get_axes().set_range(Range(0, 100.0), PlotAxes::AxisMask::Y);
    std::vector<double> x_data = fields;
    std::vector<std::unique_ptr<DiscreteSet>> sets;
    for (std::size_t i = 0; i < mtfs_by_freq.size(); i++) {
        const auto &mtf = mtfs_by_freq[i];
        const Rgb &color = FREQ_COLORS[i % FREQ_COLOR_COUNT];
        for (int xy = 0; xy < Orientation::COUNT; xy++) {
            auto set = std::make_unique<DiscreteSet>();
            set->set_interpolation(Interpolation::Cubic);
            const auto &mtf_data = (xy == Orientation::SAGITTAL) ? mtf.sag_mtf_by_field
                                                                 : mtf.tan_mtf_by_field;
            for (std::size_t j = 0; j < mtf_data.size(); j++)
                set->add_data(x_data[j], mtf_data[j] * 100.0);
            std::string label = df().format(mtf.freq) +
                                (xy == Orientation::SAGITTAL ? " Sagittal" : " Tangential");
            PlotData *pd =
                plot.add_plot_data(set.get(), color, label,
                                   static_cast<int>(PlotStyleMask::InterpolatePlot));
            pd->set_line_style(xy == Orientation::SAGITTAL ? PlotData::LineStyle::Solid
                                                           : PlotData::LineStyle::Dashed);
            sets.push_back(std::move(set));
        }
    }
    plot.get_axes().set_label("Fields", PlotAxes::AxisMask::X);
    plot.get_axes().set_label("MTF", PlotAxes::AxisMask::Y);
    plot.get_axes().set_unit("", false, false, 0, PlotAxes::AxisMask::Y);
    plot.get_axes().set_unit("", false, false, 0, PlotAxes::AxisMask::X);
    RendererSvg r(1024, 640);
    PlotRenderer plotRenderer;
    plotRenderer.draw_plot(r, plot);
    std::string out;
    r.write(out);
    return out;
}

std::string GeoMTFByFieldPlot::toString() const {
    std::string sb;
    sb += ",";
    for (std::size_t i = 0; i < mtfs_by_freq[0].tan_mtf_by_field.size(); i++) {
        if (i > 0)
            sb += ",";
        sb += df().format(fields[i]);
    }
    sb += "\n";
    for (std::size_t i = 0; i < mtfs_by_freq.size(); i++) {
        const auto &mtf = mtfs_by_freq[i];
        for (int xy = 0; xy < Orientation::COUNT; xy++) {
            // Java builds a DiscreteSet here and never uses it; dropped.
            const auto &mtf_data = (xy == Orientation::SAGITTAL) ? mtf.sag_mtf_by_field
                                                                 : mtf.tan_mtf_by_field;
            sb += std::to_string(mtf.freq);
            sb += " ";
            sb += Orientation::name(xy);
            sb += ",";
            for (std::size_t j = 0; j < mtf_data.size(); j++) {
                if (j > 0)
                    sb += ",";
                sb += df().format(mtf_data[j]);
            }
            sb += "\n";
        }
    }
    return sb;
}

// ---------------------------------------------------------------------------
// RayAberrationPlot
// ---------------------------------------------------------------------------

double RayAberrationPlot::auto_y_scale() const {
    double yscale = 0.0;
    for (const auto &fan_result : ray_aberration_results->results) {
        if (yscale < fan_result.max_y_val)
            yscale = fan_result.max_y_val;
    }
    return yscale;
}

namespace {
/** Java prints the enum by name. */
const char *fan_type_name(RayFanType t) {
    switch (t) {
    case RayFanType::TransverseRayFan: return "TransverseRayFan";
    case RayFanType::OpticalPathDifference: return "OpticalPathDifference";
    }
    return "?";
}
} // namespace

std::string RayAberrationPlot::plot(const rayoptics::raytr::TraceFanResult &fan_result,
                                    double yscale) const {
    if (yscale == 0)
        yscale = auto_y_scale();
    Plot plot;
    // fan_result.type is set by eval_*_fan before this is ever called.
    plot.set_title(std::string(fan_type_name(*fan_result.type)) + " " +
                   fan_result.fld->toString());
    plot.get_axes().set_position(Vector3::vector3_0);
    plot.get_axes().set_range(Range(-1.0, 1.0), PlotAxes::AxisMask::X);
    plot.get_axes().set_tics_step(1.0, PlotAxes::AxisMask::X);
    plot.get_axes().set_range(Range(-yscale, yscale), PlotAxes::AxisMask::Y);
    std::vector<std::unique_ptr<DiscreteSet>> sets;
    for (const auto &fan : fan_result.fans) {
        const auto &x_data = fan.fan_x;
        const auto &y_data = fan.fan_y;
        auto set = std::make_unique<DiscreteSet>();
        set->set_interpolation(Interpolation::Cubic);
        for (std::size_t i = 0; i < x_data.size(); i++) {
            // Java unboxes the Double here; a null would NPE, and the fan
            // callbacks only answer null for a ray that failed to trace.
            set->add_data(x_data[i], y_data[i].value());
        }
        plot.add_plot_data(set.get(), Colors::get_wavelen_color(fan.wvl), "label",
                           static_cast<int>(PlotStyleMask::InterpolatePlot));
        sets.push_back(std::move(set));
    }
    std::string x_label;
    std::string y_label;
    if (fan_result.type == RayFanType::TransverseRayFan) {
        if (fan_result.xy == 1) {
            x_label = "Py";
            y_label = "eY";
        } else {
            x_label = "Px";
            y_label = "eX";
        }
    } else {
        y_label = "W";
        if (fan_result.xy == 1)
            x_label = "Py";
        else
            x_label = "Px";
    }
    plot.get_axes().set_label(x_label, PlotAxes::AxisMask::X);
    plot.get_axes().set_label(y_label, PlotAxes::AxisMask::Y);
    RendererSvg r(640, 640);
    PlotRenderer plotRenderer;
    plotRenderer.draw_plot(r, plot);
    std::string out;
    r.write(out);
    return out;
}

} // namespace redukti::plotter
