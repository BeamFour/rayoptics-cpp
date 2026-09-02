// C++ port of org.redukti.rayoptics.specs.OpticalSpecs
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_SPECS_OPTICALSPECS_H
#define REDUKTI_RAYOPTICS_SPECS_OPTICALSPECS_H

#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/specs/FieldSpec.h"
#include "redukti/rayoptics/specs/PupilSpec.h"
#include "redukti/rayoptics/specs/SpecTypes.h"
#include "redukti/rayoptics/specs/WvlSpec.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <memory>
#include <optional>
#include <string>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::parax {
class ParaxData;
}

namespace redukti::rayoptics::specs {

/**
 * Aggregate the optical usage information: aperture, field of view, spectrum
 * and focus.
 *
 * `opt_model` is the back-reference to the owning OpticalModel. It is a raw
 * pointer: the OpticalModel owns this object, so it always outlives it, and
 * keeping it as a plain pointer means every method that reaches back through it
 * ports with its Java signature unchanged.
 */
class OpticalSpecs {
public:
    static bool do_aiming_default;

    std::unique_ptr<PupilSpec> pupil;
    std::unique_ptr<FieldSpec> fov;
    std::unique_ptr<WvlSpec> wvls;
    std::unique_ptr<FocusRange> focus;
    /** Null until update_optical_properties() has run. */
    std::shared_ptr<parax::ParaxData> parax_data;
    optical::OpticalModel *opt_model = nullptr;
    bool do_aiming = true;

    explicit OpticalSpecs(optical::OpticalModel *opt_model_);
    ~OpticalSpecs();

    void update_model();

    void update_optical_properties();

    void apply_scale_factor(double scale_factor);

    Coord obj_coords(Field &fld) { return fov->obj_coords(fld); }

    FocusRange *defocus() { return focus.get(); }

    /**
     * Returns the field, wavelength and defocus for the given indices.
     * `wl` selects a wavelength (null means the central one); `fr` is the focus
     * range parameter, -1.0 to 1.0.
     */
    util::Triple<Field *, double, double> lookup_fld_wvl_focus(
        int fi, std::optional<int> wl, std::optional<double> fr);

    util::Triple<Field *, double, double> lookup_fld_wvl_focus(int fi) {
        return lookup_fld_wvl_focus(fi, std::nullopt, 0.0);
    }

    ConjugateType conjugate_type(std::optional<ImageKey> space);

    /** Object and image space refractive indices, signed by z_dir. */
    util::Pair<double, double> obj_img_rindex();

    /** Start a ray from the object surface for the given pupil coordinate. */
    Coord ray_start_from_osp(const std::vector<double> &pupil, Field &fld,
                             raytr::PupilType pupil_type);

    WvlSpec *spectral_region() { return wvls.get(); }
    FieldSpec *field_of_view() { return fov.get(); }

    void list_str(std::string &sb) const;
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_OPTICALSPECS_H
