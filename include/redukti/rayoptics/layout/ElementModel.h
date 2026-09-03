// C++ port of org.redukti.rayoptics.layout: ElementType, Element and its
// implementations, plus ElementModel and LayoutOptions.
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_LAYOUT_ELEMENTMODEL_H
#define REDUKTI_RAYOPTICS_LAYOUT_ELEMENTMODEL_H

#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Interface.h"

#include <memory>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::layout {

enum class ElementType {
    LENS,
    CEMENTED_LENS,
    STOP,
    APERTURE,
    DUMMY_INTERFACE,
    AIR_GAP,
    MIRROR,
};

/**
 * Java's `sealed`-ish set of records implementing Element.
 *
 * The Java uses `instanceof` pattern matching to dispatch; the port keeps one
 * base with a `type()` tag and downcasts on it, which is the same dispatch
 * written out.
 */
class Element {
public:
    virtual ~Element() = default;
    virtual std::string label() const = 0;
    virtual ElementType type() const = 0;
};

class DummyInterface : public Element {
public:
    int surfaceIndex;
    std::shared_ptr<seq::Interface> surface;
    std::string _label;

    DummyInterface(int surfaceIndex_, std::shared_ptr<seq::Interface> surface_,
                   std::string label_)
        : surfaceIndex(surfaceIndex_), surface(std::move(surface_)),
          _label(std::move(label_)) {}

    std::string label() const override { return _label; }
    ElementType type() const override { return ElementType::DUMMY_INTERFACE; }
};

class AirGap : public Element {
public:
    int gapIndex;
    std::shared_ptr<seq::Gap> gap;

    AirGap(int gapIndex_, std::shared_ptr<seq::Gap> gap_)
        : gapIndex(gapIndex_), gap(std::move(gap_)) {}

    std::string label() const override { return "Air " + std::to_string(gapIndex); }
    ElementType type() const override { return ElementType::AIR_GAP; }
};

class Aperture : public Element {
public:
    int surfaceIndex;
    std::shared_ptr<seq::Interface> referenceSurface;

    Aperture(int surfaceIndex_, std::shared_ptr<seq::Interface> referenceSurface_)
        : surfaceIndex(surfaceIndex_), referenceSurface(std::move(referenceSurface_)) {}

    std::string label() const override { return "Aperture " + std::to_string(surfaceIndex); }
    ElementType type() const override { return ElementType::APERTURE; }
};

class Mirror : public Element {
public:
    int surfaceIndex;
    std::shared_ptr<seq::Interface> surface;

    Mirror(int surfaceIndex_, std::shared_ptr<seq::Interface> surface_)
        : surfaceIndex(surfaceIndex_), surface(std::move(surface_)) {}

    std::string label() const override { return "Mirror " + std::to_string(surfaceIndex); }
    ElementType type() const override { return ElementType::MIRROR; }
};

class Stop : public Element {
public:
    int surfaceIndex;
    std::shared_ptr<seq::Interface> referenceSurface;

    Stop(int surfaceIndex_, std::shared_ptr<seq::Interface> referenceSurface_)
        : surfaceIndex(surfaceIndex_), referenceSurface(std::move(referenceSurface_)) {}

    std::string label() const override { return "Stop"; }
    ElementType type() const override { return ElementType::STOP; }
};

class LensElement : public Element {
public:
    int firstSurfaceIndex;
    int secondSurfaceIndex;
    std::shared_ptr<elem::surface::Surface> surface1;
    std::shared_ptr<elem::surface::Surface> surface2;
    std::shared_ptr<seq::Gap> gap;

    LensElement(int firstSurfaceIndex_, int secondSurfaceIndex_,
                std::shared_ptr<elem::surface::Surface> surface1_,
                std::shared_ptr<elem::surface::Surface> surface2_,
                std::shared_ptr<seq::Gap> gap_)
        : firstSurfaceIndex(firstSurfaceIndex_), secondSurfaceIndex(secondSurfaceIndex_),
          surface1(std::move(surface1_)), surface2(std::move(surface2_)),
          gap(std::move(gap_)) {}

    std::string label() const override { return "L" + std::to_string(firstSurfaceIndex); }
    ElementType type() const override { return ElementType::LENS; }
};

class CementedElement : public Element {
public:
    std::vector<int> surfaceIndices;
    std::vector<std::shared_ptr<elem::surface::Surface>> surfaces;
    std::vector<std::shared_ptr<seq::Gap>> gaps;

    CementedElement(std::vector<int> surfaceIndices_,
                    std::vector<std::shared_ptr<elem::surface::Surface>> surfaces_,
                    std::vector<std::shared_ptr<seq::Gap>> gaps_);

    std::string label() const override { return "CE" + std::to_string(surfaceIndices[0]); }
    ElementType type() const override { return ElementType::CEMENTED_LENS; }
};

class LayoutOptions {
public:
    bool drawOpticalAxis = true;
    bool drawElements = true;
    bool drawReferenceRays = true;
    int fanRayCount = 0;
    bool useTraceFan = false;
    bool clipRays = false;
    int surfaceSamples = 101;
    double margin = 0.05;

    LayoutOptions &fanRayCount_(int count) {
        fanRayCount = count;
        return *this;
    }
    LayoutOptions &useTraceFan_(bool value) {
        useTraceFan = value;
        return *this;
    }
    LayoutOptions &drawReferenceRays_(bool value) {
        drawReferenceRays = value;
        return *this;
    }
    LayoutOptions &clipRays_(bool value) {
        clipRays = value;
        return *this;
    }
};

class ElementModel {
public:
    explicit ElementModel(optical::OpticalModel *opticalModel);

    void updateModel();

    const std::vector<std::shared_ptr<Element>> &elements() const { return _elements; }

    static bool isAperture(seq::SequentialModel *sm, int surfaceIndex, double wavelength);
    static bool isAir(const std::shared_ptr<seq::Gap> &gap, double wavelength);

private:
    optical::OpticalModel *opticalModel;
    std::vector<std::shared_ptr<Element>> _elements;
};

} // namespace redukti::rayoptics::layout

#endif // REDUKTI_RAYOPTICS_LAYOUT_ELEMENTMODEL_H
