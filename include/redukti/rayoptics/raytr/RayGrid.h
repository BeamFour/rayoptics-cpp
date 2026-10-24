// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.RayGrid
#ifndef REDUKTI_RAYOPTICS_RAYTR_RAYGRID_H
#define REDUKTI_RAYOPTICS_RAYTR_RAYGRID_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/raytr/Analyses.h"
#include "redukti/rayoptics/raytr/RayTrace.h"

#include <array>
#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::specs {
class Field;
}

namespace redukti::rayoptics::raytr {

/**
 * Container for a square grid of rays.
 *
 *     Attributes:
 *         opt_model: :class:`~.OpticalModel` instance
 *         f: index into :class:`~.FieldSpec` or a :class:`~.Field` instance
 *         wl: wavelength (nm) to trace the fan, or central wavelength if None
 *         foc: focus shift to apply to the results
 *         image_pt_2d: base image point. if None, the chief ray is used
 *         image_delta: image offset to apply to image_pt_2d
 *         num_rays: number of samples along the side of the grid
 */
class RayGrid {
public:
    /** Borrowed. */
    optical::OpticalModel *opt_model;
    /** Borrowed; the trace sets its chief ray and reference sphere. */
    specs::Field *fld;
    double wvl;
    double foc;
    std::optional<mathlib::Vector2> image_pt_2d;
    std::optional<mathlib::Vector2> image_delta;
    int num_rays;
    double value_if_none;
    double oversize;
    /** output_filter, rayerr_filter and check_apertures (clip_rays) */
    TraceOptions rt_kwargs;

    /** Empty until the first update_data. */
    std::optional<Analyses::WavefrontGridPkg> grid_pkg;
    /** [3][num_rays][num_rays]: pupil x, pupil y and OPD in waves */
    std::vector<std::vector<std::vector<double>>> grid;

    RayGrid(optical::OpticalModel *opt_model, int f, std::optional<double> wl,
            std::optional<double> foc, std::optional<mathlib::Vector2> image_pt_2d,
            std::optional<mathlib::Vector2> image_delta, std::optional<std::string> output_filter,
            std::optional<std::string> rayerr_filter, int num_rays, bool clip_rays,
            double value_if_none, double oversize);

    RayGrid(optical::OpticalModel *opt_model, specs::Field &f, std::optional<double> wl,
            std::optional<double> foc, std::optional<mathlib::Vector2> image_pt_2d,
            std::optional<mathlib::Vector2> image_delta, std::optional<std::string> output_filter,
            std::optional<std::string> rayerr_filter, int num_rays, bool clip_rays,
            double value_if_none, double oversize);

    /** The upstream defaults: field 0, central wavelength, nominal focus, 21 rays, clipped. */
    explicit RayGrid(optical::OpticalModel *opt_model);

    /**
     * @param rebuild true for upstream build='rebuild': retrace the grid before refocusing
     */
    RayGrid &update_data(bool rebuild);

    /**
     * Calculate the input and output grid spacings.
     *
     * @param maxdim The total width of the sampling grid
     * @return delta_x: The linear grid spacing on the entrance pupil,
     *         delta_xp: The linear grid spacing on the image plane
     */
    std::array<double, 2> calc_psf_scaling(int maxdim) const;
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_RAYGRID_H
