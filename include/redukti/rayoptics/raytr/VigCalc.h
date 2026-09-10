// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
// Java version by Dibyendu Majumdar
// See LICENSE-ray-optics.txt
//
// C++ port of org.redukti.rayoptics.raytr.VigCalc
#ifndef REDUKTI_RAYOPTICS_RAYTR_VIGCALC_H
#define REDUKTI_RAYOPTICS_RAYTR_VIGCALC_H

#include "redukti/mathlib/ScalarObjectiveFunction.h"
#include "redukti/rayoptics/raytr/RayTypes.h"

#include <memory>
#include <optional>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}

namespace redukti::rayoptics::raytr {

// Vignetting and clear aperture setting operations
class VigCalc {
public:
    /** Null when any ray in the set is shorter than surface i. */
    static std::optional<double> max_aperture_at_surf(
        const std::vector<std::vector<std::shared_ptr<const RayPkg>>> &rayset, int i);

    /**
     * From existing fields and vignetting, calculate clear apertures.
     *
     *     Args:
     *         avoid_list: list of surfaces to skip when setting apertures.
     *         include_list: list of surfaces to include when setting apertures.
     *
     *     If specified, only one of either `avoid_list` or `include_list` should be specified. If neither is specified, all surfaces are set. If both are specified, the `avoid_list` is used.
     *
     *     If a surface is specified as the aperture stop, that surface's aperture is determined from the boundary rays of the first field.
     *
     *     The avoid_list idea and implementation was contributed by Quentin Bécar
     */
    /** Set the clear apertures of the model from traced marginal rays. */
    static void set_clear_apertures(optical::OpticalModel *opt_model,
                                    const std::vector<int> *avoid_list,
                                    const std::vector<int> *include_list);

    /**
     * From existing fields and vignetting, calculate clear apertures.
     *
     *     This function modifies the max_aperture maintained by the list of
     *     :class:`~.interface.Interface` in the
     *     :class:`~.sequential.SequentialModel`. For each interface, the smallest
     *     aperture that will pass all of the (vignetted) boundary rays, for each
     *     field, is chosen.
     *
     *     The change of the apertures is propagated to the
     *     :class:`~.elements.ElementModel` via
     *     :meth:`~.elements.ElementModel.sync_to_seq`.
     */
    static void set_ape(optical::OpticalModel *opm, const std::vector<int> *avoid_list,
                        const std::vector<int> *include_list);
    static void set_ape(optical::OpticalModel *opm) {
        set_ape(opm, nullptr, nullptr);
    }

    /**
     * From existing fields and clear apertures, calculate vignetting.
     */
    static void set_vig(optical::OpticalModel *opm, std::optional<bool> use_bisection);
    static void set_vig(optical::OpticalModel *opm) { set_vig(opm, std::nullopt); }

    static void set_stop_aperture(optical::OpticalModel *opm);

    static void set_pupil(optical::OpticalModel *opm, bool use_parax);
    static void set_pupil(optical::OpticalModel *opm) { set_pupil(opm, false); }

    static void calc_vignetting_for_field(optical::OpticalModel *opm, specs::Field &fld,
                                          double wvl, std::optional<bool> use_bisection,
                                          std::optional<int> max_iter_count);

    static VigResult calc_vignetted_ray(optical::OpticalModel *opm, int xy,
                                        const mathlib::Vector2 &start_dir,
                                        specs::Field &fld, double wvl,
                                        std::optional<int> max_iter_count);

    static VigResult calc_vignetted_ray_by_bisection(optical::OpticalModel *opm, int xy,
                                                     const mathlib::Vector2 &start_dir,
                                                     specs::Field &fld, double wvl,
                                                     std::optional<int> max_iter_count);

    /**
     * iterates a ray to r_target on interface indx, returns aim points on
     *     the paraxial entrance pupil plane
     *
     *     If indx is None, i.e. a floating stop surface, returns r_target.
     *
     *     If the iteration fails, a :class:`~.traceerror.TraceError` will be raised
     *
     *     Args:
     *         opm: :class:`~.OpticalModel` instance
     *         indx: index of interface whose edge is the iteration target
     *         xy: 0 or 1 depending on x or y axis as the pupil direction
     *         start_r0: iteration starting point
     *         r_target: clear aperture radius that is the iteration target.
     *         fld: :class:`~.Field` point for wave aberration calculation
     *         wvl: wavelength of ray (nm)
     *
     *     Returns:
     *         start_coords: pupil coordinates for ray thru r_target on ifc indx.
     */
    static mathlib::Vector2 iterate_pupil_ray(optical::OpticalModel *opt_model,
                                              std::optional<int> indx, int xy,
                                              double start_r0, double r_target,
                                              specs::Field &fld, double wvl);

    // FIXME this is same as R_pupil_coordnate except it returns null on error rather than throwing
    /**
     * Radial pupil coordinate error at surface `indx`.
     *
     * eval() returns null when the ray failed before that surface; find_edge
     * uses the null to steer its bracket.
     */
    class Fn_r_pupil_coordinate : public mathlib::ScalarObjectiveFunction {
    public:
        optical::OpticalModel *opt_model;
        int indx;
        int xy;
        specs::Field *fld;
        double wvl;
        double r_target;

        Fn_r_pupil_coordinate(optical::OpticalModel *opt_model_, int indx_, int xy_,
                              specs::Field *fld_, double wvl_, double r_target_)
            : opt_model(opt_model_), indx(indx_), xy(xy_), fld(fld_), wvl(wvl_),
              r_target(r_target_) {}

        std::optional<double> eval(double xy_coord) override;
    };
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_VIGCALC_H
