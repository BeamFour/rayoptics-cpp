// C++ port of org.redukti.rayoptics.raytr.VigCalc
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
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

class VigCalc {
public:
    /** Null when any ray in the set is shorter than surface i. */
    static std::optional<double> max_aperture_at_surf(
        const std::vector<std::vector<std::shared_ptr<const RayPkg>>> &rayset, int i);

    /** Set the clear apertures of the model from traced marginal rays. */
    static void set_clear_apertures(optical::OpticalModel *opt_model,
                                    const std::vector<int> *avoid_list,
                                    const std::vector<int> *include_list);

    static void set_ape(optical::OpticalModel *opm, const std::vector<int> *avoid_list,
                        const std::vector<int> *include_list);
    static void set_ape(optical::OpticalModel *opm) {
        set_ape(opm, nullptr, nullptr);
    }

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

    static mathlib::Vector2 iterate_pupil_ray(optical::OpticalModel *opt_model,
                                              std::optional<int> indx, int xy,
                                              double start_r0, double r_target,
                                              specs::Field &fld, double wvl);

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
