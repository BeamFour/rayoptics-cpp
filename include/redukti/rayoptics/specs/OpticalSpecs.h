// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.specs.OpticalSpecs
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
 * The OpticalSpecs class holds the optical usage definition of the model.
 * Aperture, field of view, wavelength, and focal position are all aspects of
 * the OpticalSpecs.
 * *
 * The first order properties are calculated and maintained by OpticalSpecs
 * in the parax_data variable. This is an instance of ParaxData that includes
 * the paraxial axial and chief rays, and the FirstOrderData that contains
 * first order properties.
 */
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

    /**
     * Aperture specification
     */
    std::unique_ptr<PupilSpec> pupil;
    /**
     * Field of view specification
     */
    std::unique_ptr<FieldSpec> fov;
    /**
     * Wavelengths
     */
    std::unique_ptr<WvlSpec> wvls;
    /**
     * Focal position
     */
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

    /* returns field, wavelength and defocus data
    Args:
        fi (int): index into the field_of_view list of Fields
        wl (int): index into the spectral_region list of wavelengths
        fr (float): focus range parameter, -1.0 to 1.0
    Returns:
        (**fld**, **wvl**, **foc**)

        - **fld** - :class:`Field` instance for field_of_view[fi]
        - **wvl** - wavelength in nm
        - **foc** - focus shift from image interface
    */
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

    /**
     * return the refractive indices in object and image space.
     */
    /** Object and image space refractive indices, signed by z_dir. */
    util::Pair<double, double> obj_img_rindex();

    /**
     * turn pupil and field specs into ray start specification.
     *
     *         Args:
     *             pupil: aperture coordinates of ray
     *             fld: instance of :class:`~.Field`
     *             pupil_type: controls how `pupil` data is interpreted
     *                 - 'rel pupil': relative pupil coordinates
     *                 - 'aim pt': aim point on pupil plane
     *                 - 'aim dir': aim direction in object space
     */
    /** Start a ray from the object surface for the given pupil coordinate. */
    Coord ray_start_from_osp(const std::vector<double> &pupil, Field &fld,
                             raytr::PupilType pupil_type);

    WvlSpec *spectral_region() { return wvls.get(); }
    FieldSpec *field_of_view() { return fov.get(); }

    void list_str(std::string &sb) const;
};

} // namespace redukti::rayoptics::specs

#endif // REDUKTI_RAYOPTICS_SPECS_OPTICALSPECS_H
