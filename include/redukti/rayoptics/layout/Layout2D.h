// C++ port of org.redukti.rayoptics.layout.Layout2D
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_LAYOUT_LAYOUT2D_H
#define REDUKTI_RAYOPTICS_LAYOUT_LAYOUT2D_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/mathlib/Vector3.h"
#include "redukti/rayoptics/layout/ElementModel.h"
#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/render/Renderer.h"

#include <map>
#include <string>
#include <vector>

namespace redukti::rayoptics::layout {

/** Draws a 2D cross-section of the optical system. */
class Layout2D {
public:
    std::string renderSvg(optical::OpticalModel *model, double width, double height,
                          const LayoutOptions *options);

    void render(render::RendererViewport &renderer, optical::OpticalModel *model,
                const LayoutOptions *options);

private:
    class Polyline {
    public:
        std::vector<mathlib::Vector2> points;
        render::Rgb color;
        double strokeWidth;

        Polyline(std::vector<mathlib::Vector2> points_, const render::Rgb &color_,
                 double strokeWidth_ = 1.0)
            : points(std::move(points_)), color(color_), strokeWidth(strokeWidth_) {}
    };

    class LensDrawing {
    public:
        double drawRadius1, drawRadius2;
        bool flat1, flat2;
        double mechanicalRadius;
    };

    class SurfaceDrawing {
    public:
        std::shared_ptr<elem::surface::Surface> surface;
        double profileRadius;
        bool flat;
        double mechanicalRadius;
    };

    class Bounds {
    public:
        double minX, minY, maxX, maxY;

        Bounds();
        void add(const mathlib::Vector2 &p);
        bool valid() const { return minX <= maxX && minY <= maxY; }
        void pad(double margin);
    };

    void addElements(std::vector<Polyline> &out, optical::OpticalModel *model,
                     ElementModel &elementModel, int samples);
    void addLens(std::vector<Polyline> &out, seq::SequentialModel *sm,
                 const LensElement &lensElement, int samples);
    void addCementedElement(std::vector<Polyline> &out, seq::SequentialModel *sm,
                            const CementedElement &element, int samples);
    void addSurface(std::vector<Polyline> &out, seq::SequentialModel *sm, int index,
                    const elem::surface::Surface &surface, double radius, int samples);
    void addFlats(std::vector<Polyline> &out, seq::SequentialModel *sm, int surfaceIndex,
                  const elem::surface::Surface &surface, double profileRadius,
                  double mechanicalRadius);
    void addFlat(std::vector<Polyline> &out, seq::SequentialModel *sm, int surfaceIndex,
                 const elem::surface::Surface &surface, double fromY, double toY);
    void addCommonEdges(std::vector<Polyline> &out, seq::SequentialModel *sm,
                        const LensElement &lensElement, double profileRadius1,
                        double profileRadius2, double mechanicalRadius);
    void addCommonEdge(std::vector<Polyline> &out, seq::SequentialModel *sm,
                       const LensElement &lensElement, double profileY1, double profileY2,
                       double edgeY);
    void addStop(std::vector<Polyline> &out, seq::SequentialModel *sm, const Stop &stop);
    void addAperture(std::vector<Polyline> &out, seq::SequentialModel *sm,
                     const Aperture &aperture);
    void addApertureMarker(std::vector<Polyline> &out, const math::Tfm3d &tfm,
                           double radius, double strokeWidth);
    void addImagePlane(std::vector<Polyline> &out, optical::OpticalModel *model,
                       const DummyInterface &image);
    void addRays(std::vector<Polyline> &out, optical::OpticalModel *model,
                 const LayoutOptions &options);
    void addDirectFan(std::vector<Polyline> &out, optical::OpticalModel *model,
                      specs::Field &field, double wavelength, const render::Rgb &color,
                      const LayoutOptions &options);
    void addTraceFan(std::vector<Polyline> &out, optical::OpticalModel *model,
                     specs::Field &field, double wavelength, const render::Rgb &color,
                     const LayoutOptions &options);
    static raytr::TraceOptions traceOptions(const LayoutOptions &options);
    void addRay(std::vector<Polyline> &out, seq::SequentialModel *sm,
                const std::shared_ptr<const raytr::RayPkg> &ray, const render::Rgb &color);

    static LensDrawing lensDrawing(const LensElement &lens);
    static void mergeSurfaceDrawing(std::map<int, SurfaceDrawing> &drawings, int index,
                                    const std::shared_ptr<elem::surface::Surface> &surface,
                                    double profileRadius, bool flat,
                                    double mechanicalRadius);

    static double semiDiameter(const seq::Interface &surface);
    static double maxAperture(const seq::Interface &surface);
    static double surfaceRadius(const seq::Interface &surface);
    static mathlib::Vector2 toLayout(const math::Tfm3d &tfm, const mathlib::Vector3 &local);
    static void flush(std::vector<Polyline> &out, const std::vector<mathlib::Vector2> &points,
                      const render::Rgb &color);
    static Bounds bounds(const std::vector<Polyline> &lines);
    static Bounds modelBounds(optical::OpticalModel *model);
};

} // namespace redukti::rayoptics::layout

#endif // REDUKTI_RAYOPTICS_LAYOUT_LAYOUT2D_H
