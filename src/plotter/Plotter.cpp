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
// Distinct colors cycled per frequency (10, 30, 50, ...)
namespace Orientation = rayoptics::util::Orientation;

// ---------------------------------------------------------------------------
// Colors
// ---------------------------------------------------------------------------

Rgb Colors::get_wavelen_color(double wl) {
    // based on algorithm from Dan Bruton
    // (www.physics.sfasu.edu/astro/color.html)
    // http://www.physics.sfasu.edu/astro/color/spectra.html
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
                // 380 to 440
                return Rgb(s * -(wl - 440.0) / 60.0, 0.0, s, 1.0);
            else
                // 440 to 490
                return Rgb(0.0, s * (wl - 440.0) / 50.0, s, 1.0);
        } else {
            // 490 to 510
            return Rgb(0.0, s, s * -(wl - 510.0) / 20.0, 1.0);
        }
    } else {
        if (wl < 645.0) {
            if (wl < 580.0)
                // 510 to 580
                return Rgb(s * (wl - 510.0) / 70.0, s, 0.0, 1.0);
            else
                // 580 to 645
                return Rgb(s, s * -(wl - 645.0) / 65.0, 0.0, 1.0);
        } else {
            // 645 to 780
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
    // MTF is plotted as a percentage on a 0-100 scale
    plot.get_axes().set_range(Range(0, 100.0), PlotAxes::AxisMask::Y);
    std::vector<double> x_data = fields;
    std::vector<std::unique_ptr<DiscreteSet>> sets;
    // for each freq
    for (std::size_t i = 0; i < mtfs_by_freq.size(); i++) {
        const auto &mtf = mtfs_by_freq[i];
        // color encodes the frequency
        const Rgb &color = FREQ_COLORS[i % FREQ_COLOR_COUNT];
        for (int xy = 0; xy < Orientation::COUNT; xy++) {
            auto set = std::make_unique<DiscreteSet>();
            set->set_interpolation(Interpolation::Cubic);
            const auto &mtf_data = (xy == Orientation::SAGITTAL) ? mtf.sag_mtf_by_field
                                                                 : mtf.tan_mtf_by_field;
            for (std::size_t j = 0; j < mtf_data.size(); j++)
                // scale 0..1 MTF to a 0..100 percentage
                set->add_data(x_data[j], mtf_data[j] * 100.0);
            std::string label = df().format(mtf.freq) +
                                (xy == Orientation::SAGITTAL ? " Sagittal" : " Tangential");
            PlotData *pd =
                plot.add_plot_data(set.get(), color, label,
                                   static_cast<int>(PlotStyleMask::InterpolatePlot));
            // line pattern encodes sagittal vs tangential
            pd->set_line_style(xy == Orientation::SAGITTAL ? PlotData::LineStyle::Solid
                                                           : PlotData::LineStyle::Dashed);
            sets.push_back(std::move(set));
        }
    }
    plot.get_axes().set_label("Fields", PlotAxes::AxisMask::X);
    plot.get_axes().set_label("MTF", PlotAxes::AxisMask::Y);
    // keep both axes on a plain scale (no x10^n factor)
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
    // for each freq
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
    //plot.get_axes().set_unit("",false,false,0, PlotAxes.AxisMask.X);
    //plot.get_axes().set_unit("",false,false,0, PlotAxes.AxisMask.Y);
    //plot.get_axes().set_unit("",true,false,0, PlotAxes.AxisMask.Y);
    RendererSvg r(640, 640);
    PlotRenderer plotRenderer;
    plotRenderer.draw_plot(r, plot);
    std::string out;
    r.write(out);
    return out;
}

// ---------------------------------------------------------------------------
// PupilMapPlot
// ---------------------------------------------------------------------------

namespace {

using PupilSample = rayoptics::analysis::PupilMapAnalysis::Sample;

const Rgb PUPIL_PASSED(0.30f, 0.69f, 0.31f, 1.0f);
const Rgb PUPIL_FAILED(0.93f, 0.93f, 0.93f, 1.0f);
const Rgb PUPIL_NOMINAL = Rgb::rgb_black;
const Rgb PUPIL_PIECEWISE(0.76f, 0.07f, 0.12f, 1.0f);
const Rgb PUPIL_ELLIPSE(0.12f, 0.47f, 0.71f, 1.0f);

/** A colour per blocking surface, assigned in the order the surfaces turn up. */
const std::vector<Rgb> &pupilPalette() {
    static const std::vector<Rgb> palette{
        Rgb(0.84f, 0.15f, 0.16f, 1.0f), Rgb(1.00f, 0.50f, 0.05f, 1.0f),
        Rgb(0.58f, 0.40f, 0.74f, 1.0f), Rgb(0.55f, 0.34f, 0.29f, 1.0f),
        Rgb(0.89f, 0.47f, 0.76f, 1.0f), Rgb(0.74f, 0.74f, 0.13f, 1.0f),
        Rgb(0.09f, 0.75f, 0.81f, 1.0f), Rgb(0.68f, 0.78f, 0.91f, 1.0f),
        Rgb(1.00f, 0.73f, 0.47f, 1.0f), Rgb(0.77f, 0.69f, 0.84f, 1.0f)};
    return palette;
}

/** Java's LinkedHashMap<Integer, Rgb>: surfaces in the order they first block a ray. */
using SurfaceColours = std::vector<std::pair<int, Rgb>>;

const Rgb *colourFor(const SurfaceColours &colours, int surface) {
    for (const auto &entry : colours)
        if (entry.first == surface)
            return &entry.second;
    return nullptr;
}

Rgb pupilColourOf(const PupilSample &s, const SurfaceColours &colours) {
    if (s.passed)
        return PUPIL_PASSED;
    const Rgb *colour = s.blocked_by >= 0 ? colourFor(colours, s.blocked_by) : nullptr;
    return colour != nullptr ? *colour : PUPIL_FAILED;
}

bool sameRgb(const Rgb &a, const Rgb &b) {
    return a.r == b.r && a.g == b.g && a.b == b.b;
}

} // namespace

std::string PupilMapPlot::plot(int size) const {
    using rayoptics::analysis::PupilMapAnalysis;
    RendererSvg r(size, size, Rgb::rgb_white);
    double reach = map->reach;
    r.set_window(Vector2Pair(Vector2(-reach, -reach), Vector2(reach, reach)), true);

    SurfaceColours colours;
    for (const auto &s : map->samples)
        if (!s.passed && s.blocked_by >= 0 && colourFor(colours, s.blocked_by) == nullptr)
            colours.emplace_back(
                s.blocked_by, pupilPalette()[colours.size() % pupilPalette().size()]);

    // The grid is regular and its regions are contiguous, so each row is drawn as a few
    // filled runs rather than as tens of thousands of points.
    double cell = 2 * reach / (map->num_samples - 1.0);
    for (int j = 0; j < map->num_samples; j++) {
        int start = 0;
        for (int i = 1; i <= map->num_samples; i++) {
            Rgb run = pupilColourOf(map->sample(start, j), colours);
            if (i < map->num_samples && sameRgb(run, pupilColourOf(map->sample(i, j), colours)))
                continue;
            double x0 = map->coordinate(start) - cell / 2;
            double x1 = map->coordinate(i - 1) + cell / 2;
            double y0 = map->coordinate(j) - cell / 2;
            double y1 = map->coordinate(j) + cell / 2;
            r.draw_polygon({Vector2(x0, y0), Vector2(x1, y0), Vector2(x1, y1), Vector2(x0, y1)},
                           run, true, true);
            start = i;
        }
    }

    r.draw_circle(Vector2(0.0, 0.0), 1.0, PUPIL_NOMINAL, false);
    draw_region(r, PUPIL_PIECEWISE, false);
    draw_region(r, PUPIL_ELLIPSE, true);

    const auto &fld = *map->fld;
    double top = reach * 0.93;
    Vector2 direction(1.0, 0.0);
    r.draw_text(Vector2(-reach * 0.97, top), direction, "field " + formatF(fld.yv(), 2),
                Renderer::TextAlignLeft, 16, Rgb::rgb_black);
    r.draw_text(Vector2(-reach * 0.97, top - reach * 0.09), direction,
                "vig scales  x " + formatF(PupilMapAnalysis::scale(fld.vlx), 3) + "/" +
                    formatF(PupilMapAnalysis::scale(fld.vux), 3) + "  y " +
                    formatF(PupilMapAnalysis::scale(fld.vly), 3) + "/" +
                    formatF(PupilMapAnalysis::scale(fld.vuy), 3),
                Renderer::TextAlignLeft, 12, Rgb::rgb_gray);
    r.draw_text(Vector2(-reach * 0.97, -top + reach * 0.09), direction,
                "piecewise sampled " + formatF(100 * map->piecewise_quality.sampled, 0) +
                    "% covered " + formatF(100 * map->piecewise_quality.covered, 0) + "%",
                Renderer::TextAlignLeft, 12, PUPIL_PIECEWISE);
    r.draw_text(Vector2(-reach * 0.97, -top), direction,
                "ellipse   sampled " + formatF(100 * map->ellipse_quality.sampled, 0) +
                    "% covered " + formatF(100 * map->ellipse_quality.covered, 0) + "%",
                Renderer::TextAlignLeft, 12, PUPIL_ELLIPSE);

    int row = 0;
    for (const auto &entry : colours) {
        r.draw_text(Vector2(reach * 0.97, top - row * reach * 0.075), direction,
                    "blocked by s" + intToString(entry.first), Renderer::TextAlignRight, 12,
                    entry.second);
        row++;
    }
    std::string out;
    r.write(out);
    return out;
}

void PupilMapPlot::draw_region(Renderer &r, const Rgb &rgb, bool ellipse) const {
    using rayoptics::analysis::PupilMapAnalysis;
    const auto &fld = *map->fld;
    const int steps = 180;
    std::optional<Vector2> previous;
    for (int i = 0; i <= steps; i++) {
        double t = 2 * mathlib::M::PI * i / steps;
        double x = std::cos(t), y = std::sin(t);
        Vector2 point = ellipse
            ? Vector2(PupilMapAnalysis::ellipse_offset(fld.vlx, fld.vux) +
                          x * PupilMapAnalysis::ellipse_scale(fld.vlx, fld.vux),
                      PupilMapAnalysis::ellipse_offset(fld.vly, fld.vuy) +
                          y * PupilMapAnalysis::ellipse_scale(fld.vly, fld.vuy))
            : Vector2(x * fld.vignetting_scale_x(x), y * fld.vignetting_scale_y(y));
        if (previous.has_value())
            r.draw_segment(*previous, point, rgb);
        previous = point;
    }
}

} // namespace redukti::plotter
