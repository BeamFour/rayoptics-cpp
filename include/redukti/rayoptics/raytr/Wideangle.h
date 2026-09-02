// C++ port of org.redukti.rayoptics.raytr.Wideangle
//
// Copyright 2017-2025 Michael J. Hayford
// Original software https://github.com/mjhoptics/ray-optics
#ifndef REDUKTI_RAYOPTICS_RAYTR_WIDEANGLE_H
#define REDUKTI_RAYOPTICS_RAYTR_WIDEANGLE_H

#include "redukti/mathlib/ScalarObjectiveFunction.h"
#include "redukti/rayoptics/raytr/RayTypes.h"
#include "redukti/rayoptics/util/Tuples.h"

#include <optional>
#include <string>
#include <vector>

namespace redukti::rayoptics::optical {
class OpticalModel;
}
namespace redukti::rayoptics::seq {
class SequentialModel;
}

namespace redukti::rayoptics::raytr {

class Wideangle {
public:
    /** An entrance pupil position and the resulting ray height at the stop. */
    class ZEnpStopHt {
    public:
        double z_enp;
        double ht_at_stop;

        ZEnpStopHt(double z_enp_, double ht_at_stop_)
            : z_enp(z_enp_), ht_at_stop(ht_at_stop_) {}

        std::string toString() const;
    };

    /** Trace from an assumed entrance pupil position to the stop surface. */
    static RayResultWithStopCoord enp_z_coordinate(double z_enp,
                                                   seq::SequentialModel *seq_model,
                                                   int stop_idx,
                                                   const mathlib::Vector3 &dir0,
                                                   double obj_dist, double wvl);

    static RayResultWithZEnp find_real_enp(optical::OpticalModel *opm,
                                           std::optional<int> stop_idx,
                                           specs::Field &fld, double wvl,
                                           const std::string &selector);

    /** Defaults to the "rev1" search. */
    static RayResultWithZEnp find_real_enp(optical::OpticalModel *opm,
                                           std::optional<int> stop_idx,
                                           specs::Field &fld, double wvl);

    static RayResultWithZEnp find_real_enp_rev1(optical::OpticalModel *opm,
                                                std::optional<int> stop_idx,
                                                specs::Field &fld, double wvl,
                                                std::optional<bool> check_direction);

    static RayResultWithZEnp find_real_enp_orig(optical::OpticalModel *opm,
                                                std::optional<int> stop_idx,
                                                specs::Field &fld, double wvl);

    /**
     * Bisect towards the edge of the region where the ray still traces.
     *
     * A null from f.eval means the ray failed at that position, and that is
     * what steers the bracket -- see the note on ScalarObjectiveFunction.
     */
    static ZEnpStopHt find_edge(mathlib::ScalarObjectiveFunction &f, double a, double b,
                                std::optional<int> max_iter);

    static util::Pair<mathlib::Vector3, RayResult> find_z_enp_on_interval(
        optical::OpticalModel *opt_model, std::optional<int> stop_idx, double start_z,
        double end_z, double z_estimate, specs::Field &fld, double wvl);

    static RayResultWithZEnp find_z_enp(optical::OpticalModel *opt_model,
                                        std::optional<int> stop_idx, double z_enp_0,
                                        specs::Field &fld, double wvl);

    static std::vector<double> linspace(double start, double end, int num);

    /** Find the object space ray that lands at the requested real image height. */
    static RayDataWithZ_Enp eval_real_image_ht(optical::OpticalModel *opt_model,
                                               specs::Field &fld, double wvl);

    /** Wraps enp_z_coordinate as a scalar function of the pupil position. */
    class Enp_z_coordinate_wrapper : public mathlib::ScalarObjectiveFunction {
    public:
        seq::SequentialModel *seq_model;
        int stop_idx;
        mathlib::Vector3 dir0;
        double obj_dist;
        double wvl;

        Enp_z_coordinate_wrapper(seq::SequentialModel *seq_model_, int stop_idx_,
                                 const mathlib::Vector3 &dir0_, double obj_dist_,
                                 double wvl_)
            : seq_model(seq_model_), stop_idx(stop_idx_), dir0(dir0_),
              obj_dist(obj_dist_), wvl(wvl_) {}

        /** Null when the ray failed to reach the stop. */
        std::optional<double> eval(double z_enp) override;
    };

    class Eval_Z_Enp_Function : public mathlib::ScalarObjectiveFunction {
    public:
        seq::SequentialModel *seq_model;
        int stop_idx;
        mathlib::Vector3 dir0;
        double obj_dist;
        double wvl;
        double y_target;
        RayResult rr;
        bool rr_set = false;

        Eval_Z_Enp_Function(seq::SequentialModel *seq_model_, int stop_idx_,
                            const mathlib::Vector3 &dir0_, double obj_dist_, double wvl_,
                            double y_target_)
            : seq_model(seq_model_), stop_idx(stop_idx_), dir0(dir0_),
              obj_dist(obj_dist_), wvl(wvl_), y_target(y_target_) {}

        std::optional<double> eval(double z_enp) override;
    };
};

} // namespace redukti::rayoptics::raytr

#endif // REDUKTI_RAYOPTICS_RAYTR_WIDEANGLE_H
