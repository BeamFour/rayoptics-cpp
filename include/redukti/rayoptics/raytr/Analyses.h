// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.Analyses
#ifndef REDUKTI_RAYOPTICS_RAYTR_ANALYSES_H
#define REDUKTI_RAYOPTICS_RAYTR_ANALYSES_H

#include "redukti/mathlib/Vector2.h"
#include "redukti/rayoptics/raytr/PreOpdPkg.h"
#include "redukti/rayoptics/raytr/RayTrace.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::specs {
class Field;
}

namespace redukti::rayoptics::raytr {

/**
 * Aberration calculations for (fld, wvl, foc), including focus and image shift.
 *
 * Port of the rectilinear grid and PSF functions of upstream raytr/analyses.py.
 * The ray grid is traced once by #trace_wavefront, and #focus_wavefront
 * evaluates its OPD for a given focus and image point. #calc_psf computes the
 * monochromatic PSF of the resulting wavefront.
 *
 * Not ported: update_psf_data, which reads a maxdim attribute that upstream's RayGrid
 * does not have; the plotting code that supplies it passes maxdim to calc_psf directly.
 */
class Analyses {
public:
    /** Java's double[][][]: [num_rays][num_rays][3]. */
    using OpdGrid = std::vector<std::vector<std::array<double, 3>>>;

    /** Ray traced and focus-independent OPD data for a grid, as returned by #trace_wavefront. */
    /** Java's `record WavefrontGridPkg(List<List<GridItem>> grid, List<List<PreOpdPkg>> upd_grid)`. */
    class WavefrontGridPkg {
    public:
        std::vector<std::vector<GridItem>> grid;
        /** Null where the ray failed. */
        std::vector<std::vector<std::shared_ptr<const PreOpdPkg>>> upd_grid;

        WavefrontGridPkg(std::vector<std::vector<GridItem>> grid_,
                         std::vector<std::vector<std::shared_ptr<const PreOpdPkg>>> upd_grid_)
            : grid(std::move(grid_)), upd_grid(std::move(upd_grid_)) {}
    };

    /**
     * Trace a grid of rays at fld and wvl and return ray_pkgs in 2d list.
     *
     * Upstream defaults apply_vignetting to False unless the caller passes it; every
     * upstream caller relies on the default, and TraceOptions cannot express "unset", so it
     * is always off here. The grid itself is laid over the vignetted pupil.
     */
    static std::vector<std::vector<GridItem>> trace_ray_grid(
        optical::OpticalModel *opt_model, const TraceGridDef &grid_rng, specs::Field &fld,
        double wvl, double foc, bool append_if_none, const TraceOptions &trace_options);

    /**
     * Trace a grid of rays and evaluate the OPD across the wavefront.
     *
     * Upstream defaults check_apertures to True; here it is taken from
     * trace_options as given.
     *
     * @return [num_rays][num_rays][3] of pupil x, pupil y and OPD in waves
     */
    static OpdGrid eval_wavefront(optical::OpticalModel *opt_model, specs::Field &fld, double wvl,
                                  double foc, std::optional<mathlib::Vector2> image_pt_2d,
                                  std::optional<mathlib::Vector2> image_delta, int num_rays,
                                  double value_if_none, double oversize,
                                  const TraceOptions &trace_options);

    /**
     * Trace a grid of rays and pre-calculate data needed for rapid refocus.
     *
     * Upstream defaults check_apertures to True; here it is taken from
     * trace_options as given.
     */
    static WavefrontGridPkg trace_wavefront(optical::OpticalModel *opt_model, specs::Field &fld,
                                            double wvl, double foc,
                                            std::optional<mathlib::Vector2> image_pt_2d,
                                            std::optional<mathlib::Vector2> image_delta,
                                            int num_rays, double oversize,
                                            const TraceOptions &trace_options);

    /**
     * Given pre-traced rays and a ref. sphere, return the ray's OPD.
     *
     * @return [num_rays][num_rays][3] of pupil x, pupil y and OPD in waves
     */
    static OpdGrid focus_wavefront(optical::OpticalModel *opt_model,
                                   const WavefrontGridPkg &grid_pkg, specs::Field &fld,
                                   double wvl, double foc,
                                   std::optional<mathlib::Vector2> image_pt_2d,
                                   std::optional<mathlib::Vector2> image_delta,
                                   double value_if_none);

    // --- PSF calculation

    /**
     * Given 2 of 3 parameters, calculate the third.
     *
     * @param n       The total width of the sampling grid
     * @param n_pupil The sampling across the pupil
     * @param n_airy  The sampling across the central peak of the Airy disk
     * @return (n, n_pupil, n_airy)
     */
    static std::array<int, 3> psf_sampling(std::optional<int> n, std::optional<int> n_pupil,
                                           std::optional<int> n_airy);

    /**
     * Calculate the input and output grid spacings.
     *
     * @param ndim   The sampling across the wavefront
     * @param maxdim The total width of the sampling grid
     * @return delta_x: The linear grid spacing on the entrance pupil,
     *         delta_xp: The linear grid spacing on the image plane
     */
    static std::array<double, 2> calc_psf_scaling(optical::OpticalModel *opt_model,
                                                  const specs::Field &fld, double wvl, int ndim,
                                                  int maxdim);

    /**
     * Calculate the point spread function of wavefront W.
     * Unlike upstream, the pupil mask is kept separate from phase so zero OPD
     * remains illuminated. Non-finite samples and padding are excluded.
     *
     * @param wavefront ndim x ndim array of wavefront errors. No data
     *                  condition is indicated by nan
     * @param ndim      The sampling across the wavefront
     * @param maxdim    The total width of the sampling grid
     * @return AP, the PSF of the input wavefront
     * @throws IllegalArgumentException if the PSF has no finite positive peak
     */
    static std::vector<std::vector<double>> calc_psf(
        const std::vector<std::vector<double>> &wavefront, int ndim, int maxdim);

private:
    /** numpy.nan_to_num for a real value. */
    static double nan_to_num(double x);
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_ANALYSES_H
