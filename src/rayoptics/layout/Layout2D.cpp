// C++ port of org.redukti.rayoptics.layout.Layout2D
#include "redukti/rayoptics/layout/Layout2D.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/parax/ParaxTypes.h"
#include "redukti/rayoptics/raytr/Trace.h"
#include "redukti/rayoptics/seq/SequentialModel.h"
#include "redukti/rayoptics/specs/OpticalSpecs.h"
#include "redukti/render/RendererSvg.h"

#include <cmath>
#include <limits>

namespace redukti::rayoptics::layout {

using elem::surface::Surface;
using mathlib::Vector2;
using mathlib::Vector2Pair;
using mathlib::Vector3;
using render::Rgb;

namespace {

const Rgb ELEMENT_COLOR = Rgb::rgb_black;
const Rgb STOP_COLOR = Rgb::rgb_black;
const Rgb AXIS_COLOR = Rgb::rgb_gray;
const Rgb FIELD_COLORS[] = {Rgb::rgb_red,     Rgb::rgb_blue, Rgb::rgb_green,
                            Rgb::rgb_magenta, Rgb::rgb_cyan, Rgb::rgb_yellow};
constexpr std::size_t FIELD_COLOR_COUNT = sizeof(FIELD_COLORS) / sizeof(FIELD_COLORS[0]);
constexpr double STOP_STROKE_WIDTH = 2.5;

} // namespace

// ---------------------------------------------------------------------------
// Bounds
// ---------------------------------------------------------------------------

Layout2D::Bounds::Bounds()
    : minX(std::numeric_limits<double>::infinity()),
      minY(std::numeric_limits<double>::infinity()),
      maxX(-std::numeric_limits<double>::infinity()),
      maxY(-std::numeric_limits<double>::infinity()) {}

void Layout2D::Bounds::add(const Vector2 &p) {
    if (!std::isfinite(p.x) || !std::isfinite(p.y))
        return;
    minX = std::fmin(minX, p.x);
    maxX = std::fmax(maxX, p.x);
    minY = std::fmin(minY, p.y);
    maxY = std::fmax(maxY, p.y);
}

void Layout2D::Bounds::pad(double margin) {
    double dx = maxX - minX, dy = maxY - minY;
    double pad_ = std::fmax(std::fmax(dx, dy) * margin, 1.0e-6);
    minX -= pad_;
    maxX += pad_;
    minY -= pad_;
    maxY += pad_;
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

std::string Layout2D::renderSvg(optical::OpticalModel *model, double width, double height,
                                const LayoutOptions *options) {
    render::RendererSvg renderer(width, height);
    render(renderer, model, options);
    std::string out;
    renderer.write(out);
    return out;
}

void Layout2D::render(render::RendererViewport &renderer, optical::OpticalModel *model,
                      const LayoutOptions *options_in) {
    LayoutOptions defaults;
    const LayoutOptions &options = options_in == nullptr ? defaults : *options_in;
    if (options.surfaceSamples < 2)
        throw IllegalArgumentException("surfaceSamples must be >= 2");
    if (options.margin < 0.0)
        throw IllegalArgumentException("margin must be >= 0");
    if (options.useTraceFan && options.fanRayCount == 1)
        throw IllegalArgumentException("Trace.trace_fan requires fanRayCount >= 2");
    ElementModel elementModel(model);
    std::vector<Polyline> geometry;
    if (options.drawElements)
        addElements(geometry, model, elementModel, options.surfaceSamples);
    if (options.drawReferenceRays || options.fanRayCount > 0)
        addRays(geometry, model, options);
    Bounds b = bounds(geometry);
    if (!b.valid())
        b = modelBounds(model);
    b.pad(options.margin);
    renderer.set_window(Vector2Pair(Vector2(b.minX, b.minY), Vector2(b.maxX, b.maxY)),
                        true);
    if (options.drawOpticalAxis)
        renderer.draw_segment(Vector2(b.minX, 0.0), Vector2(b.maxX, 0.0), AXIS_COLOR);
    for (const Polyline &line : geometry) {
        renderer.set_stroke_width(line.strokeWidth);
        if (line.points.size() == 2)
            renderer.draw_segment(line.points[0], line.points[1], line.color);
        else if (line.points.size() > 2)
            renderer.draw_polygon(line.points, line.color, false, false);
    }
    renderer.set_stroke_width(1.0);
}

// ---------------------------------------------------------------------------
// Elements
// ---------------------------------------------------------------------------

void Layout2D::addElements(std::vector<Polyline> &out, optical::OpticalModel *model,
                           ElementModel &elementModel, int samples) {
    seq::SequentialModel *sm = model->seq_model.get();
    for (const auto &element : elementModel.elements()) {
        switch (element->type()) {
        case ElementType::LENS:
            addLens(out, sm, static_cast<const LensElement &>(*element), samples);
            break;
        case ElementType::CEMENTED_LENS:
            addCementedElement(out, sm, static_cast<const CementedElement &>(*element),
                               samples);
            break;
        case ElementType::STOP:
            addStop(out, sm, static_cast<const Stop &>(*element));
            break;
        case ElementType::APERTURE:
            addAperture(out, sm, static_cast<const Aperture &>(*element));
            break;
        case ElementType::DUMMY_INTERFACE: {
            const auto &dummy = static_cast<const DummyInterface &>(*element);
            if (dummy.label() == "Image")
                addImagePlane(out, model, dummy);
            break;
        }
        case ElementType::MIRROR: {
            const auto &mirror = static_cast<const Mirror &>(*element);
            auto surface = std::dynamic_pointer_cast<Surface>(mirror.surface);
            if (surface != nullptr)
                addSurface(out, sm, mirror.surfaceIndex, *surface,
                           surfaceRadius(*surface), samples);
            break;
        }
        case ElementType::AIR_GAP:
            break;
        }
    }
}

Layout2D::LensDrawing Layout2D::lensDrawing(const LensElement &lens) {
    double od1 = semiDiameter(*lens.surface1);
    double od2 = semiDiameter(*lens.surface2);
    double mechanicalRadius =
        std::fmax(maxAperture(*lens.surface1), maxAperture(*lens.surface2));
    double cv1 = lens.surface1->profile->cv;
    double cv2 = lens.surface2->profile->cv;
    double drawRadius1;
    double drawRadius2;
    bool flat1 = false;
    bool flat2 = false;
    if (cv1 > 0.0 && cv2 < 0.0) {
        drawRadius1 = mechanicalRadius;
        drawRadius2 = mechanicalRadius;
    } else if ((cv1 > 0.0 && cv2 > 0.0) || (cv1 < 0.0 && cv2 < 0.0)) {
        if (cv1 - cv2 > 0.0) {
            drawRadius1 = mechanicalRadius;
            drawRadius2 = mechanicalRadius;
        } else if (od1 > od2) {
            drawRadius1 = mechanicalRadius;
            drawRadius2 = od2;
            flat2 = true;
        } else {
            drawRadius1 = od1;
            drawRadius2 = mechanicalRadius;
            flat1 = true;
        }
    } else {
        drawRadius1 = od1;
        drawRadius2 = od2;
        flat1 = true;
        flat2 = true;
    }
    return LensDrawing{drawRadius1, drawRadius2, flat1, flat2, mechanicalRadius};
}

void Layout2D::addLens(std::vector<Polyline> &out, seq::SequentialModel *sm,
                       const LensElement &lensElement, int samples) {
    // The Java repeats the lensDrawing logic inline here; it is the same
    // computation, so it is shared.
    LensDrawing dr = lensDrawing(lensElement);
    double od1 = semiDiameter(*lensElement.surface1);
    double od2 = semiDiameter(*lensElement.surface2);
    addSurface(out, sm, lensElement.firstSurfaceIndex, *lensElement.surface1,
               dr.drawRadius1, samples);
    addSurface(out, sm, lensElement.secondSurfaceIndex, *lensElement.surface2,
               dr.drawRadius2, samples);
    if (dr.flat1)
        addFlats(out, sm, lensElement.firstSurfaceIndex, *lensElement.surface1, od1,
                 dr.mechanicalRadius);
    if (dr.flat2)
        addFlats(out, sm, lensElement.secondSurfaceIndex, *lensElement.surface2, od2,
                 dr.mechanicalRadius);
    addCommonEdges(out, sm, lensElement, dr.drawRadius1, dr.drawRadius2,
                   dr.mechanicalRadius);
}

void Layout2D::mergeSurfaceDrawing(std::map<int, SurfaceDrawing> &drawings, int index,
                                   const std::shared_ptr<Surface> &surface,
                                   double profileRadius, bool flat,
                                   double mechanicalRadius) {
    auto it = drawings.find(index);
    if (it == drawings.end()) {
        drawings[index] = SurfaceDrawing{surface, profileRadius, flat, mechanicalRadius};
    } else {
        const SurfaceDrawing previous = it->second;
        drawings[index] = SurfaceDrawing{
            surface, std::fmax(previous.profileRadius, profileRadius),
            previous.flat || flat,
            std::fmax(previous.mechanicalRadius, mechanicalRadius)};
    }
}

void Layout2D::addCementedElement(std::vector<Polyline> &out, seq::SequentialModel *sm,
                                  const CementedElement &element, int samples) {
    std::vector<LensElement> lenses;
    std::vector<LensDrawing> drawings;
    // A std::map keyed by surface index: the Java uses a LinkedHashMap, and the
    // keys are inserted in ascending surface order, so iteration order matches.
    std::map<int, SurfaceDrawing> surfaces;
    for (std::size_t i = 0; i < element.gaps.size(); i++) {
        LensElement lens(element.surfaceIndices[i], element.surfaceIndices[i + 1],
                         element.surfaces[i], element.surfaces[i + 1], element.gaps[i]);
        LensDrawing drawing = lensDrawing(lens);
        lenses.push_back(lens);
        drawings.push_back(drawing);
        mergeSurfaceDrawing(surfaces, lens.firstSurfaceIndex, lens.surface1,
                            drawing.drawRadius1, drawing.flat1, drawing.mechanicalRadius);
        mergeSurfaceDrawing(surfaces, lens.secondSurfaceIndex, lens.surface2,
                            drawing.drawRadius2, drawing.flat2, drawing.mechanicalRadius);
    }
    for (const auto &entry : surfaces) {
        const SurfaceDrawing &drawing = entry.second;
        addSurface(out, sm, entry.first, *drawing.surface, drawing.profileRadius, samples);
        if (drawing.flat)
            addFlats(out, sm, entry.first, *drawing.surface, drawing.profileRadius,
                     drawing.mechanicalRadius);
    }
    for (std::size_t i = 0; i < lenses.size(); i++) {
        const LensDrawing &drawing = drawings[i];
        addCommonEdges(out, sm, lenses[i], drawing.drawRadius1, drawing.drawRadius2,
                       drawing.mechanicalRadius);
    }
}

void Layout2D::addSurface(std::vector<Polyline> &out, seq::SequentialModel *sm, int index,
                          const Surface &surface, double radius, int samples) {
    std::vector<Vector2> points;
    for (int i = 0; i < samples; i++) {
        double y = -radius + 2.0 * radius * i / (samples - 1.0);
        try {
            double sag = surface.profile->sag(0.0, y);
            if (std::isfinite(sag))
                points.push_back(toLayout(sm->gbl_tfrms[static_cast<std::size_t>(index)],
                                          Vector3(0.0, y, sag)));
        } catch (const RuntimeException &) {
            flush(out, points, ELEMENT_COLOR);
            points.clear();
        }
    }
    flush(out, points, ELEMENT_COLOR);
}

void Layout2D::addFlats(std::vector<Polyline> &out, seq::SequentialModel *sm,
                        int surfaceIndex, const Surface &surface, double profileRadius,
                        double mechanicalRadius) {
    if (mechanicalRadius <= profileRadius)
        return;
    addFlat(out, sm, surfaceIndex, surface, profileRadius, mechanicalRadius);
    addFlat(out, sm, surfaceIndex, surface, -profileRadius, -mechanicalRadius);
}

void Layout2D::addFlat(std::vector<Polyline> &out, seq::SequentialModel *sm,
                       int surfaceIndex, const Surface &surface, double fromY,
                       double toY) {
    try {
        double sag = surface.profile->sag(0.0, fromY);
        const math::Tfm3d &tfm = sm->gbl_tfrms[static_cast<std::size_t>(surfaceIndex)];
        out.push_back(Polyline({toLayout(tfm, Vector3(0.0, fromY, sag)),
                                toLayout(tfm, Vector3(0.0, toY, sag))},
                               ELEMENT_COLOR));
    } catch (const RuntimeException &) {
    }
}

void Layout2D::addCommonEdges(std::vector<Polyline> &out, seq::SequentialModel *sm,
                              const LensElement &lensElement, double profileRadius1,
                              double profileRadius2, double mechanicalRadius) {
    addCommonEdge(out, sm, lensElement, profileRadius1, profileRadius2, mechanicalRadius);
    addCommonEdge(out, sm, lensElement, -profileRadius1, -profileRadius2,
                  -mechanicalRadius);
}

void Layout2D::addCommonEdge(std::vector<Polyline> &out, seq::SequentialModel *sm,
                             const LensElement &lensElement, double profileY1,
                             double profileY2, double edgeY) {
    try {
        double sag1 = lensElement.surface1->profile->sag(0.0, profileY1);
        double sag2 = lensElement.surface2->profile->sag(0.0, profileY2);
        Vector2 p1 = toLayout(
            sm->gbl_tfrms[static_cast<std::size_t>(lensElement.firstSurfaceIndex)],
            Vector3(0.0, edgeY, sag1));
        Vector2 p2 = toLayout(
            sm->gbl_tfrms[static_cast<std::size_t>(lensElement.secondSurfaceIndex)],
            Vector3(0.0, edgeY, sag2));
        out.push_back(Polyline({p1, p2}, ELEMENT_COLOR));
    } catch (const RuntimeException &) {
    }
}

void Layout2D::addStop(std::vector<Polyline> &out, seq::SequentialModel *sm,
                       const Stop &stop) {
    addApertureMarker(out, sm->gbl_tfrms[static_cast<std::size_t>(stop.surfaceIndex)],
                      maxAperture(*stop.referenceSurface), STOP_STROKE_WIDTH);
}

void Layout2D::addAperture(std::vector<Polyline> &out, seq::SequentialModel *sm,
                           const Aperture &aperture) {
    addApertureMarker(out, sm->gbl_tfrms[static_cast<std::size_t>(aperture.surfaceIndex)],
                      maxAperture(*aperture.referenceSurface), 1.0);
}

void Layout2D::addApertureMarker(std::vector<Polyline> &out, const math::Tfm3d &tfm,
                                 double radius, double strokeWidth) {
    double outer = radius * 1.2;
    out.push_back(Polyline({toLayout(tfm, Vector3(0, radius, 0)),
                            toLayout(tfm, Vector3(0, outer, 0))},
                           STOP_COLOR, strokeWidth));
    out.push_back(Polyline({toLayout(tfm, Vector3(0, -radius, 0)),
                            toLayout(tfm, Vector3(0, -outer, 0))},
                           STOP_COLOR, strokeWidth));
}

void Layout2D::addImagePlane(std::vector<Polyline> &out, optical::OpticalModel *model,
                             const DummyInterface &image) {
    double radius = maxAperture(*image.surface);
    try {
        if (model->optical_spec->parax_data != nullptr)
            radius = std::fmax(radius,
                               std::abs(model->optical_spec->parax_data->fod.img_ht));
    } catch (const RuntimeException &) {
    }
    const math::Tfm3d &tfm =
        model->seq_model->gbl_tfrms[static_cast<std::size_t>(image.surfaceIndex)];
    out.push_back(Polyline({toLayout(tfm, Vector3(0, -radius, 0)),
                            toLayout(tfm, Vector3(0, radius, 0))},
                           ELEMENT_COLOR));
}

// ---------------------------------------------------------------------------
// Rays
// ---------------------------------------------------------------------------

void Layout2D::addRays(std::vector<Polyline> &out, optical::OpticalModel *model,
                       const LayoutOptions &options) {
    auto &fields = model->optical_spec->fov->fields;
    double wavelength = model->seq_model->central_wavelength();
    for (std::size_t fi = 0; fi < fields.size(); fi++) {
        specs::Field &field = *fields[fi];
        const Rgb &color = FIELD_COLORS[fi % FIELD_COLOR_COUNT];
        if (options.drawReferenceRays) {
            raytr::TraceOptions to = traceOptions(options);
            raytr::RayResult chief =
                raytr::Trace::trace_ray(model, Vector2(0, 0), field, wavelength, to);
            addRay(out, model->seq_model.get(), chief.pkg, color);
            raytr::TraceOptions to2 = traceOptions(options);
            auto boundary = raytr::Trace::trace_boundary_rays_at_field(model, field,
                                                                       wavelength, to2);
            auto named = raytr::Trace::boundary_ray_dict(model, boundary);
            addRay(out, model->seq_model.get(), named["+Y"], color);
            addRay(out, model->seq_model.get(), named["-Y"], color);
        }
        if (options.fanRayCount > 0) {
            if (options.useTraceFan)
                addTraceFan(out, model, field, wavelength, color, options);
            else
                addDirectFan(out, model, field, wavelength, color, options);
        }
    }
}

void Layout2D::addDirectFan(std::vector<Polyline> &out, optical::OpticalModel *model,
                            specs::Field &field, double wavelength, const Rgb &color,
                            const LayoutOptions &options) {
    int count = options.fanRayCount;
    for (int i = 0; i < count; i++) {
        double py = count == 1 ? 0.0 : -1.0 + 2.0 * i / (count - 1.0);
        raytr::TraceOptions to = traceOptions(options);
        raytr::RayResult result =
            raytr::Trace::trace_ray(model, Vector2(0, py), field, wavelength, to);
        addRay(out, model->seq_model.get(), result.pkg, color);
    }
}

void Layout2D::addTraceFan(std::vector<Polyline> &out, optical::OpticalModel *model,
                           specs::Field &field, double wavelength, const Rgb &color,
                           const LayoutOptions &options) {
    raytr::TraceFanDef fan(Vector2(0.0, -1.0), Vector2(0.0, 1.0), options.fanRayCount);
    double focus = model->optical_spec->defocus()->get_focus();
    auto to = traceOptions(options);
    auto traced = raytr::Trace::trace_fan(model, fan, field, wavelength, focus, false,
                                          nullptr, to);
    for (const auto &item : traced)
        addRay(out, model->seq_model.get(), item.ray_pkg, color);
}

raytr::TraceOptions Layout2D::traceOptions(const LayoutOptions &options) {
    raytr::TraceOptions result;
    result.check_apertures = options.clipRays;
    result.rayerr_filter = "full";
    return result;
}

void Layout2D::addRay(std::vector<Polyline> &out, seq::SequentialModel *sm,
                      const std::shared_ptr<const raytr::RayPkg> &ray, const Rgb &color) {
    if (ray == nullptr || ray->ray.size() < 2)
        return;
    std::size_t count = std::min(ray->ray.size(), sm->gbl_tfrms.size());
    std::vector<Vector2> points;
    for (std::size_t i = 1; i < count; i++)
        points.push_back(toLayout(sm->gbl_tfrms[i], ray->ray[i].p));
    flush(out, points, color);
}

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

double Layout2D::semiDiameter(const seq::Interface &surface) {
    return std::fmax(std::abs(surface.max_aperture), 1.0e-9);
}

double Layout2D::maxAperture(const seq::Interface &surface) {
    double radius = semiDiameter(surface);
    try {
        double surfaceOd = surface.surface_od();
        if (std::isfinite(surfaceOd))
            radius = std::fmax(radius, std::abs(surfaceOd));
    } catch (const RuntimeException &) {
    }
    return radius;
}

double Layout2D::surfaceRadius(const seq::Interface &surface) {
    return maxAperture(surface);
}

Vector2 Layout2D::toLayout(const math::Tfm3d &tfm, const Vector3 &local) {
    Vector3 global = tfm.rt->multiply(local).add(tfm.t);
    return Vector2(global.z, global.y);
}

void Layout2D::flush(std::vector<Polyline> &out, const std::vector<Vector2> &points,
                     const Rgb &color) {
    if (points.size() >= 2)
        out.push_back(Polyline(points, color));
}

Layout2D::Bounds Layout2D::bounds(const std::vector<Polyline> &lines) {
    Bounds b;
    for (const Polyline &line : lines)
        for (const Vector2 &point : line.points)
            b.add(point);
    return b;
}

Layout2D::Bounds Layout2D::modelBounds(optical::OpticalModel *model) {
    Bounds b;
    seq::SequentialModel *sm = model->seq_model.get();
    for (std::size_t i = 1; i < sm->ifcs.size(); i++) {
        double r = surfaceRadius(*sm->ifcs[i]);
        const math::Tfm3d &tfm = sm->gbl_tfrms[i];
        b.add(toLayout(tfm, Vector3(0, -r, 0)));
        b.add(toLayout(tfm, Vector3(0, r, 0)));
    }
    return b;
}

} // namespace redukti::rayoptics::layout
