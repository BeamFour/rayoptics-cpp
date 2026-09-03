// C++ port of org.redukti.rayoptics.layout.ElementModel and CementedElement
#include "redukti/rayoptics/layout/ElementModel.h"

#include "redukti/Exceptions.h"
#include "redukti/rayoptics/optical/OpticalModel.h"
#include "redukti/rayoptics/seq/Medium.h"
#include "redukti/rayoptics/seq/SequentialModel.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace redukti::rayoptics::layout {

using elem::surface::Surface;
using seq::Air;
using seq::Gap;
using seq::InteractMode;
using seq::SequentialModel;

CementedElement::CementedElement(
    std::vector<int> surfaceIndices_,
    std::vector<std::shared_ptr<Surface>> surfaces_,
    std::vector<std::shared_ptr<Gap>> gaps_)
    : surfaceIndices(std::move(surfaceIndices_)), surfaces(std::move(surfaces_)),
      gaps(std::move(gaps_)) {
    if (gaps.size() < 2 || surfaces.size() != gaps.size() + 1 ||
        surfaceIndices.size() != surfaces.size())
        throw IllegalArgumentException(
            "a cemented element requires n gaps and n + 1 surfaces");
}

ElementModel::ElementModel(optical::OpticalModel *opticalModel_)
    : opticalModel(opticalModel_) {
    updateModel();
}

void ElementModel::updateModel() {
    SequentialModel *sm = opticalModel->seq_model.get();
    std::vector<std::shared_ptr<Element>> next;
    int lastSurface = static_cast<int>(sm->ifcs.size()) - 1;
    next.push_back(std::make_shared<DummyInterface>(0, sm->ifcs[0], "Object"));
    double wvl = sm->central_wavelength();
    for (int i = 1; i < lastSurface; i++) {
        auto ifc = sm->ifcs[static_cast<std::size_t>(i)];
        if (ifc->interact_mode == InteractMode::REFLECT)
            next.push_back(std::make_shared<Mirror>(i, ifc));
        if (sm->stop_surface.has_value() && *sm->stop_surface == i) {
            next.push_back(std::make_shared<Stop>(i, ifc));
        } else if (isAperture(sm, i, wvl)) {
            next.push_back(std::make_shared<Aperture>(i, ifc));
        }
    }
    for (int i = 0; i < static_cast<int>(sm->gaps.size());) {
        if (isAir(sm->gaps[static_cast<std::size_t>(i)], wvl)) {
            i++;
            continue;
        }
        int firstGap = i;
        while (i < static_cast<int>(sm->gaps.size()) &&
               !isAir(sm->gaps[static_cast<std::size_t>(i)], wvl))
            i++;
        int gapCount = i - firstGap;
        if (i >= static_cast<int>(sm->ifcs.size()))
            continue;
        if (gapCount == 1) {
            auto s1 = std::dynamic_pointer_cast<Surface>(
                sm->ifcs[static_cast<std::size_t>(firstGap)]);
            auto s2 = std::dynamic_pointer_cast<Surface>(
                sm->ifcs[static_cast<std::size_t>(firstGap + 1)]);
            if (s1 != nullptr && s2 != nullptr)
                next.push_back(std::make_shared<LensElement>(
                    firstGap, firstGap + 1, s1, s2,
                    sm->gaps[static_cast<std::size_t>(firstGap)]));
        } else if (gapCount > 1) {
            std::vector<int> indices;
            std::vector<std::shared_ptr<Surface>> surfaces;
            std::vector<std::shared_ptr<Gap>> gaps;
            bool allSurfaces = true;
            for (int j = firstGap; j <= i; j++) {
                indices.push_back(j);
                auto surface =
                    std::dynamic_pointer_cast<Surface>(sm->ifcs[static_cast<std::size_t>(j)]);
                if (surface != nullptr)
                    surfaces.push_back(surface);
                else
                    allSurfaces = false;
                if (j < i)
                    gaps.push_back(sm->gaps[static_cast<std::size_t>(j)]);
            }
            if (allSurfaces)
                next.push_back(std::make_shared<CementedElement>(indices, surfaces, gaps));
        }
    }
    next.push_back(std::make_shared<DummyInterface>(
        lastSurface, sm->ifcs[static_cast<std::size_t>(lastSurface)], "Image"));
    _elements = next;
}

bool ElementModel::isAperture(SequentialModel *sm, int surfaceIndex, double wavelength) {
    if (surfaceIndex < 1 || surfaceIndex >= static_cast<int>(sm->ifcs.size()) - 1 ||
        surfaceIndex >= static_cast<int>(sm->gaps.size()))
        return false;
    auto ifc = sm->ifcs[static_cast<std::size_t>(surfaceIndex)];
    // Java binds the preceding gap here and never reads it; kept out.
    auto followingGap = sm->gaps[static_cast<std::size_t>(surfaceIndex)];
    return ifc->interact_mode != InteractMode::REFLECT && ifc->profile != nullptr &&
           std::abs(ifc->profile->cv) < 1.0e-12 && std::abs(followingGap->thi) < 1.0e-12 &&
           isAir(followingGap, wavelength);
}

bool ElementModel::isAir(const std::shared_ptr<Gap> &gap, double wavelength) {
    if (gap == nullptr || gap->medium == nullptr)
        return true;
    if (std::dynamic_pointer_cast<Air>(gap->medium) != nullptr)
        return true;
    const auto &name = gap->medium->name();
    if (name.has_value()) {
        std::string lower;
        for (char c : *name)
            lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        if (lower == "air")
            return true;
    }
    return std::abs(gap->medium->rindex(wavelength) - 1.0) < 1.0e-12;
}

} // namespace redukti::rayoptics::layout
