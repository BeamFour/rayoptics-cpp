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

    /**
     * Locate the z center of the real pupil for `fld`, wrt 1st ifc
     *
     *     This function implements a 2 step process to finding the chief ray
     *     for `fld` and `wvl` for wide angle systems. `fld` should be of type
     *     ('object', 'angle'), even for finite object distances.
     *
     *     The first phase searches for the window of pupil locations by sampling the
     *     z coordinate starting from the paraxial pupil location. The real pupil can move either inward or outward from the paraxial pupil location. As soon as 2 successful rays are traced, the search direction is updated if needed. The search continues until z_enp values are found giving rays that straddle the stop center. If no interval is found that contains the central ray, a finer sampled search is done to find the edges more accurately. If only a single successful trace is in hand, a second, more finely subdivided search is conducted around the successful point.
     *
     *     The outcome is a range, start_z -> end_z, an estimate of where the crossing point is (z_estimate), and a ray iteration (using :func:`~.raytr.wideangle.find_z_enp_on_interval`) to find the center of the stop surface.
     */
    static RayResultWithZEnp find_real_enp_rev1(optical::OpticalModel *opm,
                                                std::optional<int> stop_idx,
                                                specs::Field &fld, double wvl,
                                                std::optional<bool> check_direction);

    //logger.info(f"fld: {fld.yv:3.1f}:   {z_enp=:8.4f}  {ht_at_stop=:10.2e}")
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

    /**
     * iterates a ray to [0, 0] on interface stop_ifc, returning aim info
     *
     *     This function finds the entrance pupil location, z_enp, inside a range of pupil locations. The rays in the interval must be trace without throwing TraceError exceptions (ignoring aperture clipping).
     *
     *     Args:
     *
     *         opt_model:  input OpticalModel
     *         stop_idx:   index of the aperture stop interface
     *         start_z:    lower bound of the z_enp interval to be searched
     *         end_z:      upper bound of the z_enp interval to be searched
     *         z_estimate: estimate of pupil location. this estimate must support
     *                     a raytrace up to stop_ifc
     *         fld:        field point
     *         wvl:        wavelength of raytrace (nm)
     *
     *     Returns z distance from 1st interface to the entrance pupil.
     *
     *     If stop_ifc is None, i.e. a floating stop surface, returns paraxial
     *     entrance pupil.
     *
     *     If the iteration fails, a TraceError will be raised
     */
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

    /**
     * returns the function value or None, if fct failed to evalute.
     */
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
