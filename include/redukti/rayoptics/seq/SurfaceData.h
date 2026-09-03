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
 * Every field is *borrowed*, and every one is nullable: the five-argument
 * zip_longest that builds these pads the short lists with null.
 *
 * The interface, gap and transform are raw pointers into the SequentialModel
 * that produced the path. That model owns them, they do not change once it is
 * built, and it outlives every trace run against it -- so sharing ownership
 * bought nothing and cost two atomic refcount operations per segment, while
 * holding the transform by value cost a 104-byte copy per segment. A path is
 * rebuilt for every single ray, so both were significant: see PERFORMANCE.md.
 *
 * The consequence is a lifetime rule. A PathSeg is only valid until the model
 * is rebuilt (update_model, add_surface, and anything else that touches ifcs,
 * gaps or lcl_tfrms). Do not hold one across such a call. reverse_path()
 * computes its transforms on the fly and parks them in the model to satisfy
 * this; see the note there.
 */
class PathSeg {
public:
    Interface *ifc;
    Gap *gap;
    const math::Tfm3d *Tfrm;
    std::optional<double> Indx;
    std::optional<util::ZDir> Zdir;

    PathSeg(Interface *ifc_, Gap *gap_, const math::Tfm3d *Tfrm_,
            std::optional<double> Indx_, std::optional<util::ZDir> Zdir_)
        : ifc(ifc_), gap(gap_), Tfrm(Tfrm_), Indx(Indx_), Zdir(Zdir_) {}
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
