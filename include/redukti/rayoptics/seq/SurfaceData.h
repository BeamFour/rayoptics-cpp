// C++ port of org.redukti.rayoptics.seq.{SurfaceData,PathSeg,NewSurfaceSpec}
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SEQ_SURFACEDATA_H
#define REDUKTI_RAYOPTICS_SEQ_SURFACEDATA_H

#include "redukti/rayoptics/elem/surface/Surface.h"
#include "redukti/rayoptics/math/Tfm3d.h"
#include "redukti/rayoptics/seq/Gap.h"
#include "redukti/rayoptics/seq/Interface.h"
#include "redukti/rayoptics/util/ZDir.h"

#include <memory>
#include <optional>
#include <string>

namespace redukti::rayoptics::seq {

/**
 * Builder for one surface of a prescription. Fields are package-private in the
 * Java and all the boxed ones are genuinely nullable -- mode() clears the glass
 * fields back to null.
 */
class SurfaceData {
public:
    std::optional<double> refractive_index;
    std::optional<double> v_number;
    double curvature;
    double thickness;
    std::optional<double> max_aperture_;
    std::optional<std::string> catalog_name;
    std::optional<std::string> glass_name;
    std::optional<InteractMode> interact_mode;

    SurfaceData(double curvature_, double thickness_)
        : curvature(curvature_), thickness(thickness_) {}

    SurfaceData *rindex(double index, double vd) {
        this->refractive_index = index;
        this->v_number = vd;
        this->glass_name = std::nullopt;
        this->catalog_name = std::nullopt;
        return this;
    }

    SurfaceData *rindex(double index, double vd, const std::string &glass_name_,
                        const std::string &catalog_name_) {
        this->refractive_index = index;
        this->v_number = vd;
        this->glass_name = glass_name_;
        this->catalog_name = catalog_name_;
        return this;
    }

    SurfaceData *mode(InteractMode mode_) {
        this->interact_mode = mode_;
        this->glass_name = std::nullopt;
        this->catalog_name = std::nullopt;
        this->refractive_index = std::nullopt;
        this->v_number = std::nullopt;
        return this;
    }

    SurfaceData *max_aperture(double map) {
        this->max_aperture_ = map;
        return this;
    }
};

/**
 * One step along the sequential path.
 *
 * SequentialModel holds its interfaces and gaps as shared_ptr (they are
 * polymorphic and get sliced into new lists by path()/reverse_path()), so a
 * PathSeg shares them rather than borrowing. Every field is nullable: the
 * five-argument zip_longest that builds these pads the short lists with null.
 */
class PathSeg {
public:
    std::shared_ptr<Interface> ifc;
    std::shared_ptr<Gap> gap;
    std::optional<math::Tfm3d> Tfrm;
    std::optional<double> Indx;
    std::optional<util::ZDir> Zdir;

    PathSeg(std::shared_ptr<Interface> ifc_, std::shared_ptr<Gap> gap_,
            std::optional<math::Tfm3d> Tfrm_, std::optional<double> Indx_,
            std::optional<util::ZDir> Zdir_)
        : ifc(std::move(ifc_)), gap(std::move(gap_)), Tfrm(std::move(Tfrm_)),
          Indx(Indx_), Zdir(Zdir_) {}
};

class NewSurfaceSpec {
public:
    std::shared_ptr<elem::surface::Surface> surface;
    std::shared_ptr<Gap> gap;
    double rndx;
    std::optional<math::Tfm3d> tfrm;
    util::ZDir z_dir;

    NewSurfaceSpec(std::shared_ptr<elem::surface::Surface> surface_,
                   std::shared_ptr<Gap> gap_, double rndx_,
                   std::optional<math::Tfm3d> tfrm_, util::ZDir z_dir_)
        : surface(std::move(surface_)), gap(std::move(gap_)), rndx(rndx_),
          tfrm(std::move(tfrm_)), z_dir(z_dir_) {}
};

} // namespace redukti::rayoptics::seq

#endif // REDUKTI_RAYOPTICS_SEQ_SURFACEDATA_H
